#!/usr/bin/env python3
"""reader_cde: edge-case checks of THM-5141 (C) as re-derived by the reader.

(1) three-block formula with the explicit constant
      K(P, Z) = 2 P 2^(w+s) + 2 Z + 1 - 2 q(P) - 2 q(Z) - P 2^w rev_s(Z)
    for all P < 2^5 (P >= 1), w = 1..5, s = 0..5 (so s = 0, w = 1, t' = 1 are included).
(2) all collisions of n-bit words, n <= 16: w >= 1, oddness, the equation, |gamma| <= 2^w, sign(X) = sign(Xi)
    (THM-4890 (C)1 in this normalisation), the carries bounds with the exact rho of the pair
    (t' <= s + rho + 2, s <= t' + rho - 1) and the sharper THM-4890 form t' - s in [dXi - dX, dXi - dX + 2].
(3) per (n, t', s) cell: number of colliding ordered pairs <= sum over odd pairs of
    4 (1 + 2^max(t'-1,s) g / max(|X|,|Xi|)) restricted to the odd pairs that occur (pair-by-pair lattice-line bound).
"""
import sys
from math import gcd
from collections import defaultdict

def q(h):
    r = 0; i = 0
    while h >> i:
        if (h >> i) & 1:
            r += ((h >> (i + 1)) + 1) // 2
        i += 1
    return r

def xf(h): return 2 * h + 1 - 2 * q(h)

def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1: r |= 1 << (w - 1 - i)
    return r

bad = 0; cases = 0
for P in range(1, 32):
    for w in range(1, 6):
        for s in range(0, 6):
            for Z in range(1 << s):
                c = 2 ** (s + 1) - rev(Z, s)
                K = 2 * P * 2 ** (w + s) + 2 * Z + 1 - 2 * q(P) - 2 * q(Z) - P * 2 ** w * rev(Z, s)
                for x in range(1 << w):
                    h = (P << (w + s)) | (x << s) | Z
                    rhs = K + c * x - P * rev(x, w) - 2 * q(x) - (P & 1) * ((x >> (w - 1)) & 1) - (c & 1) * (x & 1)
                    cases += 1
                    if xf(h) != rhs: bad += 1
print(f"(1) three-block formula with explicit K(P,Z): cases {cases}, violations {bad}")

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16
for n in range(2, NMAX + 1):
    byval = defaultdict(list)
    for h in range(1 << (n - 1), 1 << n):
        byval[xf(h)].append(h)
    npairs = 0; viol = defaultdict(int); maxdiff = -99; mindiff = 99
    cell = defaultdict(lambda: defaultdict(set))
    for v, lst in byval.items():
        if len(lst) < 2: continue
        for h in lst:
            for h2 in lst:
                if h == h2: continue
                npairs += 1
                d = h ^ h2
                hi = d.bit_length() - 1
                lo = (d & -d).bit_length() - 1
                tp = n - 1 - hi; s = lo; w = hi - lo + 1
                P = h >> (hi + 1); Z = h & ((1 << s) - 1)
                x = (h >> s) & ((1 << w) - 1); x2 = (h2 >> s) & ((1 << w) - 1)
                X = x - x2; Xi = rev(x, w) - rev(x2, w)
                c = 2 ** (s + 1) - rev(Z, s)
                gam = 2 * (q(x) - q(x2)) + (P & 1) * (((x >> (w-1)) & 1) - ((x2 >> (w-1)) & 1)) + (c & 1) * ((x & 1) - (x2 & 1))
                if w < 1 or tp < 1: viol['w/t'] += 1
                if X % 2 == 0 or Xi % 2 == 0: viol['odd'] += 1
                if c * X - P * Xi != gam: viol['eq'] += 1
                if abs(gam) > 2 ** w: viol['gamma'] += 1
                if X * Xi <= 0: viol['sign'] += 1
                dX = w - abs(X).bit_length(); dXi = w - abs(Xi).bit_length()
                rho = max(dX, dXi) + 1          # least rho with 2^(w-rho) <= |X|, |Xi|
                if not (tp <= s + rho + 2 and s <= tp + rho - 1): viol['carries'] += 1
                if not (dXi - dX <= tp - s <= dXi - dX + 2): viol['thm4890'] += 1
                maxdiff = max(maxdiff, tp - s - (dXi - dX)); mindiff = min(mindiff, tp - s - (dXi - dX))
                cell[(tp, s)][(x, x2)].add((P, c))
    lat = 0
    for (tp, s), dd in cell.items():
        for (x, x2), sols in dd.items():
            w = n - tp - s
            X = x - x2; Xi = rev(x, w) - rev(x2, w); g = gcd(abs(X), abs(Xi))
            if len(sols) > 4 * (1 + 2 ** max(tp - 1, s) * g / max(abs(X), abs(Xi))): lat += 1
    print(f"(2,3) n = {n:2d}: ordered colliding pairs {npairs}; violations {dict(viol)}; t'-s-(dXi-dX) in [{mindiff}, {maxdiff}]; lattice-line violations {lat}")
