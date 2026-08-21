"""
Egg interface for the REG calibration daemon.

An *egg* (GCP terminology) is one physically independent random bit source.
Every egg produces a time-ordered stream of whitened (folded / debiased) bits;
the daemon handles XOR masking, trial formation, statistics, and storage.

Contract — subclasses implement:

    name : str                      short identifier, e.g. "mic", "accel", "camera", "control"
    start()                         begin acquisition (own thread/process); must not block
    stop()                          release hardware; idempotent
    take_bits() -> np.ndarray       uint8 array of 0/1, all whitened bits produced since the
                                    previous call, in production order; empty array if none
    health() -> dict                JSON-serialisable per-second health/covariate readings
                                    (levels, rates, device identity, anything that could
                                    explain an excursion). Called once per second right after
                                    take_bits(). Return {} if nothing new.
    available : bool                False if start() could not open the hardware; the daemon
                                    then skips the egg and records why in the run metadata.
    unavailable_reason : str

Whitening is the egg's job because it depends on the physics of the source
(which bit planes are noise-backed, how many samples to fold). Eggs must NOT
apply the XOR mask — the daemon does that identically for every egg.

Helpers below (xor_fold, bits_from_ints) are shared by the hardware eggs.
"""
import os
import threading
import time

import numpy as np


def xor_fold(bits: np.ndarray, k: int) -> np.ndarray:
    """XOR k consecutive bits into one. Drops the incomplete tail. bits: uint8 0/1."""
    n = (len(bits) // k) * k
    if n == 0:
        return np.zeros(0, dtype=np.uint8)
    return np.bitwise_xor.reduce(bits[:n].reshape(-1, k), axis=1).astype(np.uint8)


def bits_from_ints(values: np.ndarray, plane: int) -> np.ndarray:
    """Extract bit-plane `plane` (0 = LSB) from an integer array as uint8 0/1."""
    return ((values >> plane) & 1).astype(np.uint8)


class Egg:
    name = "egg"

    def __init__(self):
        self.available = True
        self.unavailable_reason = ""
        self._lock = threading.Lock()
        self._buf = []          # list of uint8 arrays awaiting collection

    def _push(self, bits: np.ndarray):
        if len(bits):
            with self._lock:
                self._buf.append(np.asarray(bits, dtype=np.uint8))

    def take_bits(self) -> np.ndarray:
        with self._lock:
            chunks, self._buf = self._buf, []
        return np.concatenate(chunks) if chunks else np.zeros(0, dtype=np.uint8)

    def start(self):
        pass

    def stop(self):
        pass

    def health(self) -> dict:
        return {}


class ControlEgg(Egg):
    """
    Pseudo-egg: the kernel CSPRNG (/dev/random — Fortuna on macOS) at a fixed
    bit rate. Zero physical content by construction; it is the null against
    which every physical egg is judged, run through byte-identical code.
    """
    name = "control"

    def __init__(self, bits_per_second: int = 6000):
        super().__init__()
        self.rate = bits_per_second
        self._last = None

    def start(self):
        self._last = time.monotonic()

    def take_bits(self) -> np.ndarray:
        now = time.monotonic()
        n = int((now - self._last) * self.rate)
        if n > 2 * self.rate:                    # after a stall/sleep synthesise at most 2 s, never the backlog
            n = 2 * self.rate
            self._last = now - n / self.rate
        self._last += n / self.rate
        if n <= 0:
            return np.zeros(0, dtype=np.uint8)
        raw = np.frombuffer(os.urandom((n + 7) // 8), dtype=np.uint8)
        return np.unpackbits(raw)[:n]

    def health(self) -> dict:
        return {"source": "/dev/random (Fortuna CSPRNG)", "rate": self.rate}
