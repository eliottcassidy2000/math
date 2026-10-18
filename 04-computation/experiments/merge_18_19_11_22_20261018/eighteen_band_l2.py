#!/usr/bin/env python3
"""The l2 eighteen-band law (mac-mini-2026-10-18-sixcolour): exact first-transport energy E_1/N and its clock-lag split.

clock_duality's R7 (root_braid_logistic_20261017f.md section 3b) computes the l1 first-transport mass
A_2/N = 155896/87381 from ideal band coefficients c_r(x), r = t mod 9, constant on 3/2^b < x <= 6/2^b with exact period
18 in b (B has order 9 mod 9 and factor 4^9 = 2^18). Here the same coefficients are squared (l2) and resolved by the
literal step clock tau_j = a_j + 2j of each preimage s_j of t:
  E_0/N = sum over sources of h^2,   E_1/N = sum_t C(t)^2 / N,   E_1 - E_0 = 2 sum_D M_D,
where M_D sums (-1)^(a_j + a_j') h_j h_j' over preimage pairs of one node at clock lag D (THM-4690's depth-one meets).
Checks: l1 reproduces 155896/87381; E_0 = 34/7; lag-1 part = -5/7 (siblings)."""
from fractions import Fraction as Fr
from collections import defaultdict

kappa = (0, -3, 3); lam = (-1, -1, 2)
def h(r, y):                       # ideal source coefficient of a seed s = r mod 3 at s/N = y
    tot, k = 0, 0
    while 4 ** k * y <= 2:
        tot += kappa[(r + k) % 3]; k += 1
    if 4 ** k * y <= 4:
        tot += lam[(r + k) % 3]
    return tot

def preimages(r9, x):              # (a_j, j, residue of s_j mod 3, y_j = s_j/N) for t = r9 mod 9 at t/N = x
    out = []; rho = r9 % 9
    for j in range(0, 60):
        y_scale = Fr(4) ** j * x / 3
        if 2 * y_scale > 4: break   # all further preimages exceed 4N (a_j >= 1)
        if rho % 3 == 1: a = 2
        elif rho % 3 == 2: a = 1
        else: a = None
        if a is not None:
            sres = ((2 ** a * rho - 1) % 9) // 3
            out.append((a, j, sres, 2 ** a * y_scale))
        rho = (4 * rho + 1) % 9
    return out

import sys
BMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 18
# The signed coefficients c_r have exact period 18 in b (R7), so l1 and E_1 use one period. The diagonal and the
# lag-resolved pair sums are NOT periodic (squares do not cancel), so they are summed directly over b < BMAX.
l1 = Fr(0); E1 = Fr(0); lagsum = defaultdict(Fr); diag = Fr(0); E1_direct = Fr(0)
for b in range(max(18, BMAX)):
    x = Fr(9, 2 ** (b + 1))                         # representative point of the band 3/2^b < x <= 6/2^b
    w = Fr(1, 8 * 2 ** b)                           # density 1/24 per class mod 9 times band length 3/2^b
    for r9 in range(9):
        pre = preimages(r9, x)
        vals = [((-1) ** a, a + 2 * j, h(sres, y)) for (a, j, sres, y) in pre]
        C = sum(sg * hv for sg, _, hv in vals)
        if b < 18:
            l1 += w * abs(C); E1 += w * C * C
        E1_direct += w * C * C
        diag += w * sum(hv * hv for _, _, hv in vals)
        for i in range(len(vals)):
            for k in range(i + 1, len(vals)):
                D = abs(vals[i][1] - vals[k][1])
                lagsum[D] += w * vals[i][0] * vals[k][0] * vals[i][2] * vals[k][2]
per = 1 / (1 - Fr(1, 2 ** 18))
l1 *= per; E1 *= per
print('direct sums over b < %d (tail < 2^-%d * poly): E_1 direct = %.15f, periodic formula %.15f' % (max(18, BMAX), max(18, BMAX) - 10, float(E1_direct), float(E1)))
print('l1 first-transport mass A_2/N =', l1, '=', float(l1), '; clock_duality R7: 155896/87381 =', float(Fr(155896, 87381)), '; match:', l1 == Fr(155896, 87381))
print('diagonal part (= E_0/N), direct sum: %.15f ; 34/7 = %.15f' % (float(diag), 34 / 7))
print('E_1/N exact =', E1, '=', float(E1), '; THM-4667 ideal E_1 = 3.388403')
print('E_1 - E_0 (34/7) = %.15f ; 2 * sum of lag parts = %.15f' % (float(E1 - Fr(34, 7)), float(2 * sum(lagsum.values()))))
for D in sorted(lagsum):
    v = lagsum[D]
    if abs(float(v)) > 1e-9:
        print('  lag %2d: M_D = %.12f   (best small fit: %s)' % (D, float(v), Fr(float(v)).limit_denominator(200000)))
print('denominator of E_1:', E1.denominator, '=', __import__('sympy').factorint(E1.denominator))
