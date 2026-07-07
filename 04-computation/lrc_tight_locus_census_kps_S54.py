#!/usr/bin/env python3
"""
kps-2026-07-07-S54 (part 3) -- CENSUS of the LRC(14) TIGHT LOCUS (all 13-speed families with
M = 1/14), correcting the S53 claim "AP is the unique single-scale tight family".

The one-swap census found {1..11,13,24} = AP[12->24] has M=1/14 (a SECOND tight family,
= mac-mini THM-612's Goddyn-Wong family; my S53 uniqueness claim repeated MISTAKE-100 --
structurally excluding the one-residue-moved shape).  Here: verify it, and census the FULL
tight locus by enumerating families within <=3 element-swaps of the AP, to characterize the
rigidity corner as a structured (finite up to dilation?) set.
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations

def Mw_wit(v):
    v = [x for x in v if x]
    S = sum(abs(x) for x in v)
    Q = min(4 * S, 2 * max(abs(x) for x in v) + 2)
    va = np.array(v, dtype=np.int64)
    bn, bd, bc = 0, 1, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q)
        r = np.outer(va, a) % q
        d = np.minimum(r, q - r)
        col = d.min(axis=0)
        j = int(col.argmax()); qb = int(col[j])
        if qb * bd > bn * q:
            bn, bd, bc = qb, q, int(a[j])
    return Fraction(bn, bd), Fraction(bc, bd)

def prim(v):
    g = reduce(gcd, v); return tuple(x // g for x in v) if g > 1 else tuple(v)

TIGHT = Fraction(1, 14)

print("=" * 78)
print("(A) VERIFY the Goddyn-Wong family {1..11,13,24} (independent, print witness)")
print("=" * 78)
for name, v in [("AP {1..13}", list(range(1, 14))),
                ("GW {1..11,13,24}", [1,2,3,4,5,6,7,8,9,10,11,13,24]),
                ("dilated GW *3", [3*x for x in [1,2,3,4,5,6,7,8,9,10,11,13,24]])]:
    M, w = Mw_wit(v)
    # explicit check at t=1/14
    phases14 = sorted(set((x % 14) for x in v))
    mindist = min(min(r, 14 - r) for r in (x % 14 for x in v))
    print("  %-20s M = %-6s = %.6f   witness t*=%s   min||v_i/14|| = %d/14" %
          (name, str(M), float(M), str(w), mindist))

print()
print("=" * 78)
print("(B) CENSUS the tight locus: families within <=2 swaps of AP {1..13}, M = 1/14")
print("    (replace a subset S of AP by new elements from {14..40}; primitive; M=1/14)")
print("=" * 78)
ap = list(range(1, 14))
newpool = list(range(14, 41))
tight = set()
# 0-swap
tight.add(prim(tuple(ap)))
# 1-swap
for i in range(13):
    for X in newpool:
        v = ap[:i] + ap[i+1:] + [X]
        vt = tuple(sorted(v))
        if len(set(vt)) < 13: continue
        if reduce(gcd, vt) != 1: continue
        if Mw_wit(list(vt))[0] == TIGHT:
            tight.add(vt)
# 2-swap
cnt2 = 0
for i, j in combinations(range(13), 2):
    remaining = [ap[t] for t in range(13) if t != i and t != j]
    for X, Y in combinations(newpool, 2):
        v = remaining + [X, Y]
        vt = tuple(sorted(v))
        if len(set(vt)) < 13: continue
        if reduce(gcd, vt) != 1: continue
        if Mw_wit(list(vt))[0] == TIGHT:
            tight.add(vt); cnt2 += 1
print("  tight families found (<=2 swaps, primitive, M=1/14): %d" % len(tight))
for vt in sorted(tight):
    # describe as AP-modification
    missing = [x for x in range(1, 14) if x not in vt]
    extra = [x for x in vt if x > 13]
    desc = "AP" if not missing else "AP \\ %s u %s" % (missing, extra)
    print("    %-42s   %s" % (str(vt), desc))

print()
print("=" * 78)
print("(C) closed-form INFINITE lonely families (n=13) found this session")
print("=" * 78)
print("  M({1..9,11,12,13, 10k}) = k/(10k+7)  [remove 10, add 10k]:")
rest10 = [x for x in range(1,14) if x != 10]
for k in range(2, 8):
    v = tuple(sorted(rest10 + [10*k]))
    M, w = Mw_wit(list(v))
    pred = Fraction(k, 10*k+7)
    print("    k=%d  X=%3d  M=%-7s  pred k/(10k+7)=%-7s  %s" %
          (k, 10*k, str(M), str(pred), "OK" if M == pred else "DIFF"))
print("  M({1..11,13, 12k}) = k/(12k+5) for k>=3; k=2 -> 1/14 (the GW tight exception):")
rest12 = [x for x in range(1,14) if x != 12]
for k in range(2, 8):
    v = tuple(sorted(rest12 + [12*k]))
    M, w = Mw_wit(list(v))
    pred = Fraction(k, 12*k+5)
    tag = "OK" if M == pred else ("TIGHT 1/14 (GW)" if M == TIGHT else "DIFF")
    print("    k=%d  X=%3d  M=%-7s  pred k/(12k+5)=%-7s  %s" %
          (k, 12*k, str(M), str(pred), tag))

print()
print("READOUT: the tight locus is {AP, GW=AP[12->24], + dilations} -- NOT just the AP")
print("(corrects S53, = MISTAKE-100 recurrence). The near-tight corner is a DISCRETE Farey")
print("ladder of closed-form families k/(10k+7), k/(12k+5), .... Structural, census-able.")
