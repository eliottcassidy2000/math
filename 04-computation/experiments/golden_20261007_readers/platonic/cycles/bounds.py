#!/usr/bin/env python3
"""Rigorous bounds on the minimal |element| of an integer cycle of the Terras map
T(x) = x/2 (x even), (3x+1)/2 (x odd), as a function of (L, k) = (period, # odd steps).

Product identity around a cycle:  prod_{odd x_i} (3 + 1/x_i) = 2^L.
 positive cycle (x_i >= 1):  2^L/3^k = prod (1 + 1/(3x_i)) <= (1 + 1/(3 x_min))^k
     => x_min <= 1 / (3 ((2^L/3^k)^(1/k) - 1)),   needs 3^k < 2^L <= 4^k
 negative cycle (x_i = -y_i, y_i >= 1): 3^k/2^L = prod 1/(1 - 1/(3y_i)) <= (1 - 1/(3 y_min))^(-k)
     => y_min <= 1 / (3 (1 - (2^L/3^k)^(1/k))),   needs 3^k > 2^L
We compute the max over all (L,k) with L <= LMAX using mpmath at 60 digits and report
the worst (L,k) on each side.  Exact integer comparison decides the sign of 2^L - 3^k.
"""
import sys
from mpmath import mp, mpf, power, log
mp.dps = 60
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
bestP = (mpf(0), None); bestN = (mpf(0), None)
rowsP = []; rowsN = []
p3 = [1]
for k in range(1, LMAX + 2):
    p3.append(p3[-1] * 3)
for L in range(1, LMAX + 1):
    twoL = 1 << L
    # positive side: k with 3^k < 2^L <= 4^k ; bound increasing in k -> take largest k
    # (we still loop over a small window for safety)
    kmax_pos = None
    k = 0
    # largest k with 3^k < 2^L
    lo = int(L / 1.5849625007211562) - 2
    lo = max(lo, 1)
    cands = [kk for kk in range(lo, lo + 6) if kk >= 1 and p3[kk] < twoL and twoL <= (1 << (2 * kk))]
    for kk in cands:
        r = mpf(twoL) / mpf(p3[kk])
        b = 1 / (3 * (power(r, mpf(1) / kk) - 1))
        if b > bestP[0]:
            bestP = (b, (L, kk))
    candsN = [kk for kk in range(lo, lo + 6) if kk >= 1 and kk <= L and p3[kk] > twoL]
    for kk in candsN:
        r = mpf(twoL) / mpf(p3[kk])
        b = 1 / (3 * (1 - power(r, mpf(1) / kk)))
        if b > bestN[0]:
            bestN = (b, (L, kk))
print("LMAX", LMAX)
print("positive side: max x_min bound = %s at (L,k)=%s" % (mp.nstr(bestP[0], 15), bestP[1]))
print("negative side: max y_min bound = %s at (L,k)=%s" % (mp.nstr(bestN[0], 15), bestN[1]))
