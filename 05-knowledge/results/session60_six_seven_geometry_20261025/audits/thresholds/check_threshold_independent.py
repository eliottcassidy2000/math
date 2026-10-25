#!/usr/bin/env python3
"""Independent exact audit: branched strata and eigensheaves versus CI adjunction."""
from math import comb
from fractions import Fraction
from pathlib import Path
import json

def require(condition, label):
    if not condition:
        raise ValueError(label)

def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0

def hilbert(d, n, r, m):
    return sum(((-1) ** j * comb(n - d - 1, j) * choose(n - 1 + m - r * j, n - 1) for j in range(n - d)))

def character_pg(d, n, r):
    counts = [1]
    for _ in range(n - 1):
        nxt = [0] * (len(counts) + r - 1)
        for s, v in enumerate(counts):
            for e in range(r):
                nxt[s + e] += v
        counts = nxt
    return sum((v * choose((s + r - 1) // r - 1, d) for s, v in enumerate(counts)))

def euler_strata(d, n, r):
    return sum((comb(n, k) * r ** (n - 1 - k) * sum(((-1) ** j * comb(n - k, j) * (d - k - j + 1) for j in range(d - k + 1))) for k in range(d + 1)))

def euler_chern(d, n, r):
    coeff = [comb(n, k) for k in range(d + 1)]
    for _ in range(n - d - 1):
        nxt = [0] * (d + 1)
        for i, v in enumerate(coeff):
            for j in range(d - i + 1):
                nxt[i + j] += v * (-r) ** j
        coeff = nxt
    return coeff[d] * r ** (n - d - 1)
rows = []
for d in range(1, 6):
    for r in range(2, 9):
        for n in range(d + 1, 2 * d + 8):
            a = (r - 1) * n - r * (d + 1)
            pg = character_pg(d, n, r)
            e = euler_strata(d, n, r)
            require(pg == hilbert(d, n, r, a), (d, n, r, pg, a))
            require(e == euler_chern(d, n, r), (d, n, r, e))
            if a >= 0:
                require(pg >= 1, 'failed exact gate at original line 49')
            rows.append(dict(d=d, n=n, r=r, a=a, pg=pg, e=e))
cones = []
for r in range(2, 21):
    for m in range(3, 21):
        deg = r ** (m - 2)
        g = 1 + Fraction(r ** (m - 1) * (-2 + Fraction(m * (r - 1), r)), 2)
        require(g.denominator == 1, 'failed exact gate at original line 58')
        b = -1 - Fraction(2 * g - 2, deg)
        require(b == m - 1 - r * (m - 2), 'failed exact gate at original line 60')
        require((g == 0) == (r == 2 and m == 3), 'failed exact gate at original line 61')
        require((b >= 0) == (r == 2 and m == 3), 'failed exact gate at original line 62')
        cones.append(dict(r=r, m=m, genus=int(g), discrepancy=int(b)))
surface0 = [(r, n) for r in range(2, 101) for n in range(3, 101) if (r - 1) * n - 3 * r == 0]
surface1 = [(r, n) for r in range(2, 101) for n in range(3, 101) if (r - 1) * n - 3 * r == 1]
require(surface0 == [(2, 6), (4, 4)], 'failed exact gate at original line 67')
require(surface1 == [(2, 7), (3, 5), (5, 4)], 'failed exact gate at original line 68')
out = dict(status='PASS', method='branched-stratum Euler and character-eigensheaf geometric genus independently versus CI calculations', smooth_cases=len(rows), cone_cases=len(cones), surface_a0=surface0, surface_a1=surface1, smooth_rows=rows, cone_rows=cones)
Path(__file__).with_name('threshold_independent_results.json').write_text(json.dumps(out, indent=2) + '\n')
print(f'PASS: {len(rows)} independent smooth cases; {len(cones)} local RH/discrepancy cases; exact surface threshold lists.')
