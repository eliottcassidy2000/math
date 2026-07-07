#!/usr/bin/env python3
"""
kps-2026-07-07-S54 (part 4) -- the COMPLETE one-swap ladder structure of the near-tight
LRC(14) single-scale families.

The census maps M({1..13}\\{j} u {jk}) (remove AP element j, add multiple jk) for all j, k:

  * LARGE j (7 <= j <= 13): M = k/(j*k + b_j), a FAREY LADDER -> 1/j as k->infinity.
      b_j = 8,7,5,7,3,5,1 for j = 7,8,9,10,11,12,13.  These are the near-tight ladders.
  * SMALL j (2 <= j <= 6): M is CONSTANT in k (removing a small element leaves a fixed rung):
      j=2,3 -> 2/17;  j=4,5 -> 2/19;  j=6 -> 2/23.
  * The UNIQUE one-swap family reaching the tight floor 1/14 is j=12, k=2 = GW {1..11,13,24}
      (the residue-liar's k=2 boundary; k>=3 gives k/(12k+5)).  This is exactly why the
      M-tight locus is {AP, GW}.
  * The first EXCITED rung 2/27 is hit by TWO one-swaps: j=13,k=2 and j=10,k=2.

Formalized so far (kernel-pure, GREEN): j=10 (tenSwap_lonely, k/(10k+7)), j=12
(residueLiar_lonely, k/(12k+5)), GW (gw_lonely).  j=7,8,9,11,13 are analogous residue-table
ladders (not yet formalized).
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce

def Mw(v):
    v = [x for x in v if x]
    S = sum(abs(x) for x in v); Q = min(4 * S, 2 * max(abs(x) for x in v) + 2)
    va = np.array(v, dtype=np.int64); bn, bd = 0, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q); r = np.outer(va, a) % q
        d = np.minimum(r, q - r); col = d.min(axis=0); qb = int(col.max())
        if qb * bd > bn * q: bn, bd = qb, q
    return Fraction(bn, bd)

print("M({1..13}\\{j} u {jk}) for j=1..13, k=2..7 -- the one-swap ladder table")
print(" j |  k=2     k=3     k=4     k=5     k=6     k=7   | closed form (k>=3)")
print("-" * 78)
for j in range(1, 14):
    rest = [x for x in range(1, 14) if x != j]
    Ms = []
    for k in range(2, 8):
        v = tuple(sorted(rest + [j * k]))
        if len(set(v)) < 13 or reduce(gcd, v) != 1:
            Ms.append(None); continue
        Ms.append(Mw(list(v)))
    good = [(k, m) for k, m in zip(range(2, 8), Ms) if m is not None and k >= 3]
    cf = ""
    if len(good) >= 2:
        (k1, m1), (k2, m2) = good[:2]
        d1 = Fraction(k1) / m1; d2 = Fraction(k2) / m2
        a = (d2 - d1) / (k2 - k1); b = d1 - a * k1
        if a.denominator == 1 and b.denominator == 1:
            a, b = int(a), int(b)
            if all(m == Fraction(k, a * k + b) for k, m in good):
                cf = "k/(%dk%+d)" % (a, b)
        if not cf and len({m for _, m in good}) == 1:
            cf = "CONST %s" % good[0][1]
    row = "  ".join(("%-6s" % (str(m) if m is not None else "--")) for m in Ms)
    tight = "  <- GW tight 1/14 at k=2" if (j == 12) else ""
    print(" %2d|  %s | %s%s" % (j, row, cf, tight))

print()
print("STRUCTURE: near-tight LRC(14) one-swaps = AP with a LARGE element j->jk, ladder")
print("k/(jk+b_j); only j=12,k=2 (GW) is tight (1/14); everything else >= 2/27 > 1/14 LONELY.")
print("Formalized GREEN: j=10 (tenSwap), j=12 (residue-liar), GW; j=7,8,9,11,13 analogous.")
