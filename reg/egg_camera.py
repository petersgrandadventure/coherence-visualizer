"""
Camera egg: sensor noise of the built-in camera, taken from the least
significant bit of the luma plane and spatially XOR-folded.

Physics.  Per-pixel temporal noise on this sensor is ~1.6-2.4 gray levels
(median/mean at 640x480, ~30 fps, gray mean ~93), so the LSB of every pixel is
nominally driven by photon shot noise and read noise.  But the frames are
ISP-processed (demosaic, denoise, gamma, YUV->RGB->gray), which smears that
noise across time and neighbours: the raw LSB plane has lag-1 temporal
correlation 0.65, spatial neighbour correlation ~0.05, and ~49% of pixels show
a fixed pattern (ideal 5.8%).  The plane is not usable as-is.

Whitening.  XOR-folding each 16x16 block (256 pixels) into one bit multiplies
the per-bit bias and correlation terms (piling-up lemma), so 256 weakly
correlated, individually biased bits give one bit that measures clean:
p1 = 0.5013, temporal r = 0.005, spatial r < 0.002, fixed-pattern fraction
5.0% vs 5.8% ideal.  8x8 is marginal (temporal r 0.035, FPN 6.9%), so 16 is
the default: 1200 bits/frame, ~36 kbit/s at 30 fps, ample for the daemon.
Region sums across the sensor are NOT independent (quadrant correlations up to
-0.28: shared ISP gain/exposure), so the whole sensor is ONE egg.  Regions are
also not interchangeable (a ceiling light or reflection gives saturated or
fixed-pattern blocks), so the folded bits of each frame are emitted in a fixed
seeded permutation of the block grid rather than row-major: every 200-bit
trial the daemon forms then samples the whole sensor, the canonical trial is
representative of the stream, and block-column parity is not locked to the
daemon's alternating mask phase.

No XOR mask is applied here; the daemon does that.  The capture thread also
keeps per-second covariates: exposure level, saturated/black fractions, a
lens-covered flag, and the raw (unfolded) LSB bias and lag-1 temporal
correlation on a fixed random pixel subset -- the direct monitor of the ISP
behaviour the fold has to defeat.

macOS: camera access is granted per launching application (TCC).  Start the
daemon from Terminal.app; from an unauthorised context the capture does not
open and the egg reports itself unavailable.
"""
import math
import threading
import time

import cv2
import numpy as np

from eggs import Egg


class CameraEgg(Egg):
    name = "camera"

    def __init__(self, index: int = 0, width: int = 640, height: int = 480, fold: int = 16,
                 first_frame_timeout: float = 3.0, n_probe: int = 20000, seed: int = 12345):
        super().__init__()
        if fold < 1:
            raise ValueError("fold must be >= 1")
        self.index, self.width, self.height, self.fold = index, width, height, int(fold)
        self.first_frame_timeout, self.n_probe, self.seed = first_frame_timeout, n_probe, seed
        self._cap = None
        self._thread = None
        self._stop = threading.Event()
        self._first = threading.Event()
        self._hlock = threading.Lock()
        self._acc = self._new_acc()
        self._t_health = None
        self._backend = None
        self._frame_shape = None          # (H, W) of the gray plane
        self._hb = self._wb = 0           # folded block grid
        self._perm = None                 # fixed seeded emission order of the blocks
        self._probe_idx = None            # fixed random pixel subset (flat indices)
        self._prev_probe = None
        self._frames_total = 0
        self._failures_total = 0
        self._process_errors = 0
        self._last_error = None

    # ------------------------------------------------------------ lifecycle
    def _fail(self, reason: str):
        self.available = False
        self.unavailable_reason = reason

    def start(self):
        hint = "macOS camera authorisation is per launching app: run the daemon from Terminal.app"
        try:
            cap = cv2.VideoCapture(self.index, cv2.CAP_AVFOUNDATION)
            if not cap.isOpened():
                cap.release()
                self._fail(f"VideoCapture({self.index}) did not open (no device, busy, or not authorised); {hint}")
                return
            cap.set(cv2.CAP_PROP_FRAME_WIDTH, self.width)
            cap.set(cv2.CAP_PROP_FRAME_HEIGHT, self.height)
            self._backend = cap.getBackendName()
        except Exception as e:
            try:
                cap.release()
            except Exception:
                pass
            self._fail(f"VideoCapture open raised {e!r}; {hint}")
            return
        self._cap = cap
        self._stop.clear()
        self._first.clear()
        self._thread = threading.Thread(target=self._loop, name="camera-egg", daemon=True)
        self._thread.start()
        if not self._first.wait(self.first_frame_timeout):
            self.stop()
            self._fail(f"capture opened ({self._backend}) but no frame arrived within "
                       f"{self.first_frame_timeout:g} s"
                       + (f" ({self._process_errors} frame(s) failed: {self._last_error})"
                          if self._process_errors else "") + f"; {hint}")
            return
        self._t_health = time.monotonic()
        h, w = self._frame_shape
        print(f"camera: {self._backend} {w}x{h} LSB, fold {self.fold}x{self.fold} -> "
              f"{self._hb * self._wb} bits/frame", flush=True)

    def stop(self):
        self._stop.set()
        t = self._thread
        if t is not None and t is not threading.current_thread():
            t.join(timeout=6.0)     # AVFoundation grab times out at ~5 s, so the thread exits by then
            if t.is_alive():
                self._last_error = "stop(): capture thread did not exit; capture not released"
        self._thread = None

    # ------------------------------------------------------------ capture
    def _loop(self):
        cap = self._cap
        consecutive = 0
        try:
            while not self._stop.is_set():
                ok, frame = cap.read()
                if not ok or frame is None:
                    consecutive += 1
                    self._failures_total += 1
                    with self._hlock:
                        self._acc["fail"] += 1
                    time.sleep(0.005 if consecutive < 30 else 0.1)
                    continue
                consecutive = 0
                try:
                    self._process(frame)
                except Exception as e:
                    self._process_errors += 1
                    self._last_error = repr(e)
        finally:
            cap.release()       # the thread owns the capture: no cross-thread release race
            self._cap = None

    def _set_geometry(self, shape):
        h, w = shape
        self._frame_shape = (h, w)
        self._hb, self._wb = h // self.fold, w // self.fold
        self._perm = np.random.default_rng(self.seed + 1).permutation(self._hb * self._wb)
        rng = np.random.default_rng(self.seed)
        n = min(self.n_probe, h * w)
        self._probe_idx = np.sort(rng.choice(h * w, size=n, replace=False))
        self._prev_probe = None

    def _process(self, frame):
        if frame.ndim == 2:
            gray = frame
        else:
            gray = cv2.cvtColor(frame, cv2.COLOR_BGRA2GRAY if frame.shape[2] == 4 else cv2.COLOR_BGR2GRAY)
        if gray.shape != self._frame_shape:
            self._set_geometry(gray.shape)
        lsb = gray & 1                                   # uint8 0/1, unfolded plane
        k, hb, wb = self.fold, self._hb, self._wb
        blk = lsb[:hb * k, :wb * k].reshape(hb, k, wb, k)
        folded = np.bitwise_xor.reduce(blk, axis=(1, 3)).ravel()
        self._push(folded[self._perm])                   # whole-sensor order, fixed for the run
        self._frames_total += 1
        self._first.set()

        # covariates for health(): levels and the raw-LSB statistics the fold must defeat
        mean, std = cv2.meanStdDev(gray)
        probe = lsb.ravel()[self._probe_idx]
        prev, self._prev_probe = self._prev_probe, probe
        s1 = int(np.count_nonzero(probe))
        with self._hlock:
            a = self._acc
            a["frames"] += 1
            a["gray_sum"] += float(mean[0, 0])
            a["std_sum"] += float(std[0, 0])
            a["sat"] += int(np.count_nonzero(gray == 255))
            a["blk"] += int(np.count_nonzero(gray == 0))
            a["npix"] += gray.size
            a["lsb_n"] += probe.size
            a["lsb_1"] += s1
            if prev is not None:
                a["pairs"] += probe.size
                a["sx"] += s1
                a["sy"] += int(np.count_nonzero(prev))
                a["sxy"] += int(np.count_nonzero(probe & prev))

    # ------------------------------------------------------------ health
    @staticmethod
    def _new_acc():
        return dict(frames=0, fail=0, gray_sum=0.0, std_sum=0.0, sat=0, blk=0, npix=0,
                    lsb_n=0, lsb_1=0, pairs=0, sx=0, sy=0, sxy=0)

    @staticmethod
    def _pearson01(n, sx, sy, sxy):
        """Pearson r between two 0/1 sequences from their sums (x^2 = x for bits)."""
        den = (n * sx - sx * sx) * (n * sy - sy * sy)
        if n < 2 or den <= 0:
            return None
        return round((n * sxy - sx * sy) / math.sqrt(den), 5)

    def health(self) -> dict:
        now = time.monotonic()
        with self._hlock:
            a, self._acc = self._acc, self._new_acc()
        dt = now - (self._t_health or now)
        self._t_health = now
        n = a["frames"]
        h = {
            "frames_in_second": n,
            "fps": round(n / dt, 2) if dt > 0 else None,
            "gray_mean": None, "gray_std": None, "frac_saturated": None, "frac_black": None,
            "lens_covered": None, "raw_lsb_p1": None, "temporal_r": None,
            "read_failures": a["fail"],
            "read_failures_total": self._failures_total,
            "frames_total": self._frames_total,
            "frame_shape": list(self._frame_shape) if self._frame_shape else None,
            "backend": self._backend,
            "device_index": self.index,
            "fold": self.fold,
            "bits_per_frame": self._hb * self._wb,
            "block_order": "seeded permutation",
            "thread_alive": bool(self._thread and self._thread.is_alive()),
        }
        if self._process_errors:
            h["process_errors"] = self._process_errors
        if self._last_error:
            h["last_error"] = self._last_error
        if n:
            gm = a["gray_sum"] / n
            h.update(
                gray_mean=round(gm, 2),
                gray_std=round(a["std_sum"] / n, 2),
                frac_saturated=round(a["sat"] / a["npix"], 5),
                frac_black=round(a["blk"] / a["npix"], 5),
                # AE maxes gain under a cover, so a covered sensor reads a FLAT mid-gray,
                # not black: key on the flat field (spatial std) as well as darkness
                lens_covered=bool(gm < 8.0 or (a["std_sum"] / n) < 15.0),
                raw_lsb_p1=round(a["lsb_1"] / a["lsb_n"], 5),
                temporal_r=self._pearson01(a["pairs"], a["sx"], a["sy"], a["sxy"]),
            )
        return h
