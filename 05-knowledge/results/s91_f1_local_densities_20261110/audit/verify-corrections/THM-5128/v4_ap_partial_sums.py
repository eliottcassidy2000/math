#!/usr/bin/env python3
"""Verifier check 4 (THM-5128 (D2), (D4)): partial sums A_p^sigma(n) = sum_{h < 2^n} chi(h) w^sigma(h) for n = 16, 20, 22,
from the definitions in the theorem (chi = [p | gcd(x, y)] - 1/p^2, w^0 = (1/x) ln(2(x+y)/(x+2y)),
w^1 = (1/x) ln(1 + x/(2y))), against ap_constants.out / the audit's a10_ap_26.out; the non-absolute convergence;
the census value (F_small - 1/p^2) #T / (2X) is not recomputed here."""
import numpy as np

LOGH = 22
H = 1 << LOGH
psi = np.zeros(H, dtype=np.int64)
for n in range(1, LOGH + 1):
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    par = h >> 1
    psi[lo:hi] = np.where(h & 1 == 0, psi[par] + par, psi[par] + (par >> 1) + 1)
hh = np.arange(H, dtype=np.int64)
x = (2 * psi + 1).astype(np.float64)
y = (2 * hh + 1).astype(np.float64)
g = np.gcd(2 * psi + 1, 2 * hh + 1)
w0 = np.log(2 * (x + y) / (x + 2 * y)) / x
w1 = np.log1p(x / (2 * y)) / x
want = {  # (n = 16, 20, 22) from a10_ap_26.out / ap_constants.out
    3: ((0.05972, 0.05535, 0.05386), (0.08418, 0.07856, 0.07656)),
    5: ((0.02963, 0.02946, 0.02945), (0.04276, 0.04255, 0.04252)),
    7: ((0.00520, 0.00786, 0.00778), (0.00565, 0.00923, 0.00916)),
    13: ((-0.00362, -0.00381, -0.00379), (-0.00469, -0.00493, -0.00490)),
    31: ((-0.00184, -0.00192, -0.00194), (-0.00248, -0.00257, -0.00260)),
    43: ((0.00236, 0.00235, 0.00235), (0.00325, 0.00323, 0.00323)),
}
ok = True
for p, (a0, a1) in want.items():
    chi = (g % p == 0) - 1.0 / p ** 2
    c0 = np.cumsum(chi * w0)
    c1 = np.cumsum(chi * w1)
    got0 = [c0[(1 << n) - 1] for n in (16, 20, 22)]
    got1 = [c1[(1 << n) - 1] for n in (16, 20, 22)]
    good = all(abs(a - b) < 6e-6 for a, b in zip(got0, a0)) and all(abs(a - b) < 6e-6 for a, b in zip(got1, a1))
    ok &= good
    absum = np.cumsum(np.abs(chi) * w0)
    per_block = absum[(1 << 22) - 1] - absum[(1 << 21) - 1]
    lim = 2 * (p * p - 1) / p ** 4 * 0.1615389173 * np.log(2)
    print("%s p = %2d: A^0(16, 20, 22) = %+.5f %+.5f %+.5f; A^1 = %+.5f %+.5f %+.5f; sum |chi| w^0 over B_22 = %.5f (limit %.5f)"
          % ("PASS" if good else "FAIL", p, *got0, *got1, per_block, lim))
    ok &= abs(per_block / lim - 1) < 0.2
print("ALL PASS" if ok else "SOME FAIL")
