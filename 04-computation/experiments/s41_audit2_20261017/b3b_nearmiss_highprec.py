#!/usr/bin/env python3
"""Referee check B3b: evaluate the two ve near-miss candidates of b3 (and a few nearest non-candidates) at 60 digits.

A genuine unit has |x|^2 - 1 = 0 exactly; a near miss shows a deviation far above 10^-55 at 60-digit precision.
Usage: python b3b_nearmiss_highprec.py
"""
import json
import mpmath as mp
import audit2_field as F

mp.mp.dps = 60
S = F.golden_sets()
DA = [F.K(0)] + F.diffs(S['A']); DC = [F.K(0)] + F.diffs(S['C'])
D3 = F.uniq(F.diffs(S['C']) + F.diffs(S['L']))
zeta = mp.exp(2j * mp.pi / 15); jn = mp.mpc(0, mp.sqrt(11))


def hp(x):
    tot = mp.mpc(0)
    for part, jp in ((x[0], 1), (x[1], jn)):
        for i, c in enumerate(part):
            if c:
                tot += mp.mpf(c.numerator) / c.denominator * zeta ** i * jp
    return tot


w = [mp.re(r) for r in mp.polyroots([1, 2, -4, -4], maxsteps=200, extraprec=200) if 1 < mp.re(r) < 2][0]
u = (w + 3 - jn) / (w + 3 + jn)
d = 1 + mp.exp(2j * mp.pi / 5)
Z = d + u; N = abs(Z) ** 2; Dv = 1 - 4 / N; R = 4 - N
s = mp.mpc(0, mp.sqrt(-Dv)); h = mp.sqrt(R)
v = Z * (-1 + s) / 2; e = Z * s / h
Y = v * e
res = json.load(open('b3_census_numeric.json'))
for iA, i1, i2, i3 in res['summary']['ve']['others']:
    x = hp(DA[iA]) + u * hp(DC[i1]) + v * hp(DC[i2]) + Y * hp(D3[i3])
    print('ve candidate (dA %d, d1 %d, d2 %d, d3 %d): |x|^2 - 1 = %s  -> %s' % (
        iA, i1, i2, i3, mp.nstr(abs(x) ** 2 - 1, 8), 'NOT a unit' if abs(abs(x) ** 2 - 1) > mp.mpf(10) ** -40 else 'unit?'))
