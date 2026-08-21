# REG calibration instrument

A GCP / PEAR-style random-event-generator built from this Mac's own physical
noise sources, with the statistical discipline that makes such an instrument
meaningful: per-device empirical calibration, a CSPRNG control run through
byte-identical code, continuous health tests, and covariate logging so every
excursion can be explained or stand on its own.

It replaces the upstream `coherence_monitor`, whose Health metric turned out to
be a step function of a matrix rank (see the project history for the review).

## Eggs (independent physical channels)

| Egg | Physics | Whitening | Rate |
|---|---|---|---|
| `mic` | Mic/ADC thermal noise + room acoustics, float32 at the device's native 44.1 kHz, requantised to 24 bit | bit planes 0–3 (measured fair to bit 6), each XOR-folded 8:1 | ~22 kbit/s |
| `accel` | Built-in accelerometer via IOKit HID (~795 Hz native, 2⁻¹⁶ g steps, ~50 LSB rms noise) | per-axis LSB parity, XOR-folded 4:1 | ~600 bit/s |
| `camera` | Image-sensor noise through the ISP (raw LSB plane is temporally smeared; a 16×16 XOR-fold restores fair bits) | LSB plane folded 16×16, seeded block permutation | ~36 kbit/s |
| `control` | `/dev/random` — the kernel's Fortuna CSPRNG; zero physical content by construction | none | 6 kbit/s |

Bit rate is irrelevant beyond a few hundred bits per second — the GCP uses 200
bits per trial, one trial per second, per device. What matters is residual
serial correlation below ~10⁻³, stationarity over hours, and independence
between eggs. CPU timing jitter was evaluated and excluded (deterministic
function of the 24 MHz clock phase and of machine load).

## What the daemon records

Every wall-clock second, per egg: the canonical 200-bit trial sum (bits
XOR-masked with a fixed alternating template, GCP-style), every further complete
200-bit block as *extra* trials (calibration only — reaches ~10⁶ trials in
hours), bias / lag-1..10 serial correlation / runs test / mask-phase bias on the
pre-mask bits, and the egg's health (levels, device identity, input volume,
frame statistics, helper liveness). Host covariates: load, wall-vs-monotonic
clock drift, loop lateness. SQLite, WAL mode, crash-consistent.

## Running it

From **Terminal.app** (macOS authorises the camera per launching application):

```bash
cd ~/coherence-visualizer/reg && ./run_calibration.sh        # 24 h, detached; safe to close Terminal
tail -f calibration.log                                       # progress
./report.sh                                                   # analysis, any time (writes data/reg_calibration_report.md)
./run_calibration.sh stop                                     # clean stop
```

The report gives per-egg normalisation constants (mean / variance vs the
theoretical 100 / 50), whiteness and the implied trial-variance inflation,
stationarity in 1-min / 5-min / 1-h blocks, cross-egg correlation and the
Stouffer network-variance statistic (physical eggs vs control), covariate
coupling (does the acoustic envelope predict the mic's z²?), and a Monte-Carlo
simultaneous band from the control — the band a live display must use instead
of the classic pointwise parabola, which is exited about half the time under
the null.

## Honest expectations

The GCP's published effect is ~0.3σ per pre-registered event pooled over ~60
devices. A one-to-three-egg home instrument cannot detect an effect that size.
What this instrument can do is show, with a defined null and logged covariates,
whether *your* field ever departs from chance — and exactly what the room was
doing when it did.
