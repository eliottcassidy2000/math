#!/usr/bin/env python3
"""Full-period closed recursion, lean version: the maximal Fourier coefficient of the 3-adic Syracuse law over ALL units
to level 19 (S23, 2026-09-29).

m_n(k) = mu_hat_n(2^k mod 3^n) on the whole cycle Z/L_n, L_n = 2 3^(n-1), by the exact one-pole filter
m_n(k) = (g(k-1) + m_n(k-1))/2, g(k) = omega_n(k) m_(n-1)(k mod L_(n-1)) (no valuation truncation; the wrap is handled
by a first pass to obtain the periodic state and a second pass from it).  Since 2 generates the units, this is the
complete primitive Fourier profile: M(n) = max_k |m_n(k)| over all units.  Memory: two complex128 vectors of lengths
L_(n-1), L_n and one chunk buffer (level 19: 4.1 + 12.4 GB).  Prints per level: M(n), the argmax k and -k mod L_n,
the offset k - n log2 3 of whichever of +-2^k is the power of two side, the five largest |m| with their exponents,
the ell^2 mass 3^n sum_k |m_n(k)|^2 / L_n (Parseval: the Fourier mass per level), and the S20 FFT value.
Run: python 04-computation/experiments/collatz_three_mirrors_fullperiod_max_20260929.py [NMAX] [NPRINT_FROM]
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np
from scipy.signal import lfilter, lfilter_zi

M_S20 = {1: 0.5773503, 2: 0.3779236, 3: 0.2522368, 4: 0.1769989, 5: 0.1292736, 6: 0.0961064, 7: 0.0758700,
         8: 0.0608907, 9: 0.0480262, 10: 0.0382783, 11: 0.0319442, 12: 0.0264582, 13: 0.0220524, 14: 0.0191280,
         15: 0.0162845, 16: 0.0144095, 17: 0.0125107, 18: 0.0111873}
LOG23 = math.log2(3)
B = [0.0, 0.5]
A = [1.0, -0.5]


def omega_chunks(n: int, chunk: int):
    """Yield (lo, hi, omega_n[lo:hi]) for k in [0, L_n) using exact int64 arithmetic (valid for n <= 19)."""
    mod = 3 ** n
    L = 2 * 3 ** (n - 1)
    assert mod * mod < 2 ** 63
    base = np.empty(chunk, dtype=np.int64)
    x = 1
    for k in range(chunk):
        base[k] = x
        x = (2 * x) % mod
    step = pow(2, chunk, mod)
    mult = 1
    lo = 0
    while lo < L:
        hi = min(L, lo + chunk)
        r = (base[: hi - lo] * np.int64(mult)) % np.int64(mod)
        yield lo, hi, np.exp(2j * np.pi * (r.astype(np.float64) / mod))
        mult = (mult * step) % mod
        lo = hi


def level_up(prev: np.ndarray, n: int, chunk: int = 1 << 22) -> np.ndarray:
    L = 2 * 3 ** (n - 1)
    Lp = len(prev)
    out = np.empty(L, dtype=np.complex128)
    # pass 1: run the filter around the cycle from zero state until the transient 2^-(steps) is below 2^-80
    zi = np.zeros(1, dtype=np.complex128)
    for _ in range(max(1, -(-80 // L))):
        for lo, hi, om in omega_chunks(n, chunk):
            g = om * prev[np.arange(lo, hi) % Lp] if Lp < L else om * prev[lo:hi]
            _, zi = lfilter(B, A, g, zi=zi)
    # pass 2: from the periodic state (the transient of pass 1 is 2^-L, zero in float64 for L >= 60)
    for lo, hi, om in omega_chunks(n, chunk):
        g = om * prev[np.arange(lo, hi) % Lp] if Lp < L else om * prev[lo:hi]
        y, zi = lfilter(B, A, g, zi=zi)
        out[lo:hi] = y
    return out


if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    NFROM = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    t0 = time.time()
    # level 1 by hand: L = 2, m_1(k) = sum_a 2^-a e(2^(k-a) mod 3 / 3)
    m = np.array([sum(2.0 ** (-a) * np.exp(2j * np.pi * (pow(2, k - a, 3) / 3)) for a in range(1, 80)) for k in range(2)], dtype=np.complex128)
    for n in range(2, NMAX + 1):
        prev = m
        m = level_up(prev, n)
        del prev
        L = len(m)
        absm = np.abs(m)
        top = np.argpartition(-absm, 5)[:5]
        top = top[np.argsort(-absm[top])]
        kmax = int(top[0])
        Mn = float(absm[kmax])
        mass = float(np.sum(absm ** 2)) * 3 ** n / L
        if n >= NFROM:
            desc = []
            for k in top:
                k = int(k)
                kk = k if k <= L // 2 else k - L      # the exponent nearest zero: k or k - L (= -(L - k))
                mirror = (k - L // 2) % L            # -2^k = 2^(k + L/2)
                mm = mirror if mirror <= L // 2 else mirror - L
                desc.append(f"k={k} (|m|={absm[k]:.6f}; nearest-zero exponent {kk}, mirror {mm})")
            print(f"n={n:2d}: L={L}: M(n) = {Mn:.7f} (S20 FFT {M_S20.get(n, float('nan')):.7f}); Fourier mass per level {mass:.4f}; top: " + "; ".join(desc) + f"   [{time.time() - t0:.0f}s]")
            sys.stdout.flush()
    print("DONE")
