"""HYP-9370 at n = 23: zeros of the forest polynomial of the 13 doubly regular 23-tournament classes, from mac-mini's exact
forest numbers (../sierpinski_ocf_20261022/drt23_forest_exact.json, field N = N_L).  Session opus-2026-10-22-S54.
Usage: python phi_zeros_drt23.py"""
import json, os
import numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, '..', 'sierpinski_ocf_20261022', 'drt23_forest_exact.json')))
n = 23
for d in data:
    N = [int(x) for x in d['N']]
    assert all(N[L] == 0 for L in range(1, len(N), 2))
    coef = [float(Fraction(N[2 * k], n ** (2 * k))) for k in range((len(N) + 1) // 2)]
    r = np.roots(coef[::-1])
    w0 = r[np.argmin(np.abs(r))]
    print(f"{d['label']}: H = {d['H']}; smallest zero w = z^2 = {w0.real:+.6f}{w0.imag:+.6f}i, |z| = {np.sqrt(abs(w0)):.6f}; "
          f"zeros in |z| <= 1: {int(np.sum(np.abs(r) <= 1))}")
