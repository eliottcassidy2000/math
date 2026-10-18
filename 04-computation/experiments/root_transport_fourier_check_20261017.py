#!/usr/bin/env python3
"""Check of the Fourier form of the ideal ROOT transport (THM-4667 statement 10), S41.

Recomputes K_d from K_(d-1) by clock_duality_oct15 R8 (step) and compares its additive Fourier transform over
Z/3^(d+1) with Khat_d(xi, b) = sum_(o <= b) (-1)^(o+1) e(-mu_o xi / 3^(d+1)) Khat_(d-1)(2^-(o+1) xi mod 3^d, b - o),
mu_o = 2^-(o+1) - (1 - 4^-j)/3, j = floor(o/2). Output: maximal discrepancy per depth (machine precision expected).
Usage: python root_transport_fourier_check_20261017.py
Session opus-2026-10-17-S41.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_transport_depth_ladder_20261017 import K0, step
K = K0.copy()
for d in range(1, 5):
    Kn = step(K, d)                     # [period, mod]
    q = 3 ** (d + 1); qo = 3 ** d; P = 2 * q
    Fold = np.fft.fft(K, axis=1)        # Khat(xi) = sum_r K(r) e(-r xi / qo)
    Fnew = np.fft.fft(Kn, axis=1)
    pred = np.zeros_like(Fnew)
    xi = np.arange(q)
    for b in range(P):
        for o in range(b + 1):
            j = o // 2
            inv2 = pow(2, -(o + 1), q)
            # (1 - 4^-j)/3 as an integer mod q: compute 4^-j mod 3q, subtract from 1, divide by 3
            t = (1 - pow(4, -j, 3 * q)) % (3 * q); assert t % 3 == 0
            mu = (inv2 - t // 3) % q
            src = (pow(2, -(o + 1), qo) * xi) % qo
            pred[b] += (-1) ** (o + 1) * np.exp(-2j * np.pi * mu * xi / q) * Fold[(b - o) % (2 * qo), src]
    print('d =', d, 'max |Khat_d - Fourier recursion| =', np.abs(pred - Fnew).max(), ' scale', np.abs(Fnew).max())
    K = Kn
