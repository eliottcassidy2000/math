#!/usr/bin/env python3
"""Primitive Fourier maxima of the 3-adic Syracuse law to level 17 (S20, 2026-09-29).

By consistency of the laws, mu_hat_n(3^(n-h) u) = mu_hat_h(u) for u a unit mod 3^h, so the conductor profile
M(h) = max_(u unit mod 3^h) |mu_hat_h(u)| is level-independent; Mazur's (2.3) (Tao's fine-scale mixing) forces
M(h) <= C_A (h-1)^-A for every A.  Here M(h) for h <= 17 from one FFT of the level-17 law, with the argmax and
the Fourier mass per conductor level, and the values at the powers of two t = 2^s.
Run: python 04-computation/experiments/collatz_five_mirrors_fourier_deep_20260929.py   (about 8 GB, a few minutes)
"""
from __future__ import annotations

import sys
import time

import numpy as np

sys.path.insert(0, "04-computation/experiments")
from mazur_harmonic_mass_deep_20260928 import level_up  # noqa: E402

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 17
    t0 = time.time()
    mu = np.array([0.0, 1 / 3, 2 / 3])
    for n in range(2, N + 1):
        mu = level_up(mu, n)
    print(f"level {N} law [{time.time() - t0:.0f}s]")
    mod = 3 ** N
    amp = np.abs(np.fft.fft(mu))
    del mu
    print(f"FFT done [{time.time() - t0:.0f}s]")
    t = np.arange(mod, dtype=np.int64)
    v3 = np.zeros(mod, dtype=np.int8)
    m = t.copy()
    for _ in range(N):
        div = (m % 3 == 0) & (m > 0)
        v3[div] += 1
        m[div] //= 3
    del m
    lev = np.where(t == 0, 0, N - v3.astype(np.int64))
    del v3
    print("h, M(h) = max |mu_hat| over conductor 3^h, argmax/3^(N-h), Fourier mass at level h, mean |mu_hat|^2 at level h * 3^h:")
    prev = None
    for h in range(1, N + 1):
        sel = np.flatnonzero(lev == h)
        a = amp[sel]
        k = int(np.argmax(a))
        mx = float(a[k])
        am = int(t[sel[k]]) // 3 ** (N - h)
        mass = float((a ** 2).sum())
        ratio = f"{mx / prev:.4f}" if prev else "-"
        print(f"   h={h:2d}: M(h) = {mx:.6f} (ratio {ratio}), argmax u = {am} (= 2^{np.log2(am):.2f} if a power of two), mass {mass:.5f}, typical |mu_hat|^2 3^h = {mass / len(sel) * 3 ** h:.4f}")
        prev = mx
        del sel, a
    print("   |mu_hat_N(2^s)| for s = 0..60: " + ", ".join(f"{amp[pow(2, s, mod)]:.4f}" for s in range(0, 61)))
    print("DONE")
