/*
 * accel_stream — stream raw accelerometer events from the Apple SPU to stdout.
 *
 * Uses the private IOHIDEventSystemClient API (no root, no entitlement needed):
 * match services with PrimaryUsagePage 0xff00 / PrimaryUsage 3, request a
 * ReportInterval of 1000 us, and forward every kIOHIDEventTypeAccelerometer
 * (type 13) event as one little-endian binary record on stdout:
 *
 *     struct '<Qddd'  =  uint64 timestamp_ns (mach_absolute_time converted
 *                        to nanoseconds), double x, y, z  (units of g)
 *
 * Output is block-buffered and flushed every 50 events or 100 ms. Runs until
 * stdin reaches EOF, SIGTERM/SIGINT arrives, or the stdout pipe breaks.
 * One line of service identity goes to stderr at start-up.
 *
 * Build:  cc -O2 -o accel_stream accel_stream.c -framework CoreFoundation -framework IOKit
 */
#include <CoreFoundation/CoreFoundation.h>
#include <mach/mach_time.h>
#include <poll.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct __IOHIDEventSystemClient *IOHIDEventSystemClientRef;
typedef struct __IOHIDServiceClient *IOHIDServiceClientRef;
typedef struct __IOHIDEvent *IOHIDEventRef;
typedef void (*IOHIDEventCallback)(void *, void *, void *, IOHIDEventRef);

IOHIDEventSystemClientRef IOHIDEventSystemClientCreate(CFAllocatorRef);
int IOHIDEventSystemClientSetMatching(IOHIDEventSystemClientRef, CFDictionaryRef);
CFArrayRef IOHIDEventSystemClientCopyServices(IOHIDEventSystemClientRef);
Boolean IOHIDServiceClientSetProperty(IOHIDServiceClientRef, CFStringRef, CFTypeRef);
CFTypeRef IOHIDServiceClientCopyProperty(IOHIDServiceClientRef, CFStringRef);
void IOHIDEventSystemClientScheduleWithRunLoop(IOHIDEventSystemClientRef, CFRunLoopRef, CFStringRef);
void IOHIDEventSystemClientRegisterEventCallback(IOHIDEventSystemClientRef, IOHIDEventCallback, void *, void *);
int64_t IOHIDEventGetType(IOHIDEventRef);
double IOHIDEventGetFloatValue(IOHIDEventRef, int32_t);
uint64_t IOHIDEventGetTimeStamp(IOHIDEventRef);

#define USAGE_PAGE 0xff00
#define USAGE 3
#define EVENT_TYPE 13          /* kIOHIDEventTypeAccelerometer */
#define REPORT_INTERVAL_US 1000
#define FLUSH_EVERY 50
#define FLUSH_PERIOD_S 0.1

#pragma pack(push, 1)
typedef struct { uint64_t t_ns; double x, y, z; } record_t;
#pragma pack(pop)

static mach_timebase_info_data_t g_tb;
static volatile sig_atomic_t g_stop = 0;
static int g_pending = 0;

static void on_signal(int sig) { (void)sig; g_stop = 1; }

static void on_event(void *t, void *s, void *c, IOHIDEventRef e) {
    (void)t; (void)s; (void)c;
    if (IOHIDEventGetType(e) != EVENT_TYPE) return;
    record_t r;
    r.t_ns = IOHIDEventGetTimeStamp(e) * g_tb.numer / g_tb.denom;
    r.x = IOHIDEventGetFloatValue(e, (EVENT_TYPE << 16) | 0);
    r.y = IOHIDEventGetFloatValue(e, (EVENT_TYPE << 16) | 1);
    r.z = IOHIDEventGetFloatValue(e, (EVENT_TYPE << 16) | 2);
    if (fwrite(&r, sizeof r, 1, stdout) != 1) exit(0);   /* stdout gone: parent died */
    if (++g_pending >= FLUSH_EVERY) { fflush(stdout); g_pending = 0; }
}

static int stdin_closed(void) {
    struct pollfd p = { .fd = STDIN_FILENO, .events = POLLIN };
    if (poll(&p, 1, 0) <= 0) return 0;
    if (p.revents & (POLLHUP | POLLERR | POLLNVAL)) return 1;
    char buf[64];
    return read(STDIN_FILENO, buf, sizeof buf) == 0;      /* 0 = EOF; data is ignored */
}

static void on_tick(CFRunLoopTimerRef timer, void *info) {
    (void)timer; (void)info;
    if (g_pending) { fflush(stdout); g_pending = 0; }
    if (g_stop || stdin_closed() || ferror(stdout)) CFRunLoopStop(CFRunLoopGetCurrent());
}

static CFDictionaryRef matching(int page, int usage) {
    CFNumberRef p = CFNumberCreate(NULL, kCFNumberIntType, &page);
    CFNumberRef u = CFNumberCreate(NULL, kCFNumberIntType, &usage);
    const void *k[] = { CFSTR("PrimaryUsagePage"), CFSTR("PrimaryUsage") };
    const void *v[] = { p, u };
    CFDictionaryRef d = CFDictionaryCreate(NULL, k, v, 2, &kCFTypeDictionaryKeyCallBacks,
                                           &kCFTypeDictionaryValueCallBacks);
    CFRelease(p); CFRelease(u);
    return d;
}

static void describe(IOHIDServiceClientRef svc, char *out, size_t n) {
    const char *keys[] = { "Product", "Transport", "VendorID", "ProductID" };
    size_t used = 0;
    for (size_t i = 0; i < 4 && used < n; i++) {
        CFStringRef key = CFStringCreateWithCString(NULL, keys[i], kCFStringEncodingUTF8);
        CFTypeRef val = IOHIDServiceClientCopyProperty(svc, key);
        CFRelease(key);
        if (!val) continue;
        char buf[128] = "?";
        long long num;
        if (CFGetTypeID(val) == CFStringGetTypeID())
            CFStringGetCString((CFStringRef)val, buf, sizeof buf, kCFStringEncodingUTF8);
        else if (CFGetTypeID(val) == CFNumberGetTypeID() && CFNumberGetValue(val, kCFNumberLongLongType, &num))
            snprintf(buf, sizeof buf, "%lld", num);
        used += snprintf(out + used, n - used, "%s=%s ", keys[i], buf);
        CFRelease(val);
    }
}

int main(void) {
    mach_timebase_info(&g_tb);
    signal(SIGTERM, on_signal);
    signal(SIGINT, on_signal);
    signal(SIGPIPE, SIG_DFL);
    setvbuf(stdout, NULL, _IOFBF, 1 << 16);

    IOHIDEventSystemClientRef client = IOHIDEventSystemClientCreate(kCFAllocatorDefault);
    if (!client) { fprintf(stderr, "IOHIDEventSystemClientCreate failed\n"); return 2; }
    CFDictionaryRef m = matching(USAGE_PAGE, USAGE);
    IOHIDEventSystemClientSetMatching(client, m);
    CFRelease(m);
    CFArrayRef svcs = IOHIDEventSystemClientCopyServices(client);
    long n = svcs ? CFArrayGetCount(svcs) : 0;
    if (n == 0) { fprintf(stderr, "no IOHID service with usage page 0x%x usage %d\n", USAGE_PAGE, USAGE); return 3; }

    int ri = REPORT_INTERVAL_US;
    CFNumberRef riN = CFNumberCreate(NULL, kCFNumberIntType, &ri);
    char ident[512] = "";
    for (long i = 0; i < n; i++) {
        IOHIDServiceClientRef svc = (IOHIDServiceClientRef)CFArrayGetValueAtIndex(svcs, i);
        IOHIDServiceClientSetProperty(svc, CFSTR("ReportInterval"), riN);
        if (i == 0) describe(svc, ident, sizeof ident);
    }
    CFRelease(riN);
    fprintf(stderr, "accel_stream: services=%ld %sreport_interval_us=%d timebase=%u/%u\n",
            n, ident, ri, g_tb.numer, g_tb.denom);
    fflush(stderr);

    IOHIDEventSystemClientScheduleWithRunLoop(client, CFRunLoopGetCurrent(), kCFRunLoopDefaultMode);
    IOHIDEventSystemClientRegisterEventCallback(client, on_event, NULL, NULL);
    CFRunLoopTimerRef tick = CFRunLoopTimerCreate(NULL, CFAbsoluteTimeGetCurrent() + FLUSH_PERIOD_S,
                                                  FLUSH_PERIOD_S, 0, 0, on_tick, NULL);
    CFRunLoopAddTimer(CFRunLoopGetCurrent(), tick, kCFRunLoopDefaultMode);
    CFRunLoopRun();
    fflush(stdout);
    return 0;
}
