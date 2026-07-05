#!/usr/bin/env python3
"""
klein-2026-07-04-S131 (HYP-4093) - THE COMPRESSED FLOOR: is every COMPRESSED covering family
lonely with M >= 1/13 (> 1/14)?  This is THE sole open leaf of the LRC(14) formalization.

State (LRCEndgameAssembly, kps-S6): LRC(14) <= LRC(13) + hcomp, where
  hcomp := every COMPRESSED covering family is lonely.
  COMPRESSED := forall i, exists j!=i, |v_i| <= 13|v_j|   <=>   v_max <= 13 * v_second.
The DOMINANT case (v_max > 13*v_second) is discharged by the peel/hdom; the deep well
{1..12,182} is DOMINANT (182 > 13*12=156), EXCLUDED from hcomp.  mac-mini-S45 flagged
"compressed => M >= 1/13" as the open direction (compressed leaf floors 7/89).

klein-S129: deep well UNIQUELY attains 14/183; every OTHER covering family is >= 7/89.
=> If the deep well is the only sub-7/89 covering family and it is dominant, every COMPRESSED
covering family is >= 7/89 > 1/13 > 1/14.  Test: enumerate compressed covering families (extremal
minimal-tightener shapes, S128 method + compressed filter) + all-large consecutive blocks; find
the MIN M, and whether ANY lands in the danger zone (1/14, 1/13) or below 1/14.

Exact (Fractions), minimal-tightener shapes (extremal, contain the min) -- fast.
"""
from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations

N = 14
ONE13 = F(1, 13); ONE14 = F(1, 14); SEV89 = F(7, 89)

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval(S, Qcap):
    best = F(0)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m
    return best

def is_covering(S, n=N):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

def is_compressed(S):
    s = sorted(S)
    return s[-1] <= 13 * s[-2]

def missing(A, n=N):
    return [q for q in range(2, n + 1) if not any(a % q == 0 for a in A)]

def qcap(S):
    return min(2 * max(S) + 2, 900)

def partitions_into(items, k):
    items = list(items)
    if k <= 0: return
    if k == 1: yield [items]; return
    if len(items) < k: return
    first, rest = items[0], items[1:]
    for p in partitions_into(rest, k - 1): yield [[first]] + p
    for p in partitions_into(rest, k):
        for i in range(len(p)): yield p[:i] + [[first] + p[i]] + p[i + 1:]

def smallest_tightener(block, used, scale=1):
    L = 1
    for q in block: L = lcm(L, q)
    x = ((14 + L - 1) // L) * L * scale
    # smallest multiple of L*scale that is >= 14 and unused
    step = L
    while x in used or x < 14: x += step
    return x

print(f"1/14={float(ONE14):.6f}  1/13={float(ONE13):.6f}  7/89={float(SEV89):.6f}")
print(f"DANGER ZONE for hcomp = (1/14, 1/13)")
print("=" * 88)

AP = list(range(1, 14))
danger = []; below14 = []
gmin = None; gat = None
fams = set()
# (1) minimal-tightener covering shapes (also scaled tighteners x1,x2,x3 to sample larger killers)
for d in range(0, 5):
    for drop in combinations(AP, d):
        A = [v for v in AP if v not in drop]
        Qm = missing(A)
        if len(Qm) < d: continue
        for part in partitions_into(Qm, d):
            for scale in (1, 2, 3, 4):
                used = set(A); T = []; ok = True
                for block in part:
                    t = smallest_tightener(block, used, scale); T.append(t); used.add(t)
                S = tuple(sorted(A + T))
                if len(S) == 13 and not missing(list(S)):
                    fams.add(S)
n1 = 0
for S in fams:
    S = list(S)
    if not is_compressed(S): continue
    n1 += 1
    M = Mval(S, qcap(S))
    if gmin is None or M < gmin: gmin, gat = M, S
    if ONE14 < M < ONE13: danger.append((M, S))
    if M < ONE14: below14.append((M, S))
print(f"(1) compressed covering families (minimal+scaled-tightener shapes): {n1} of {len(fams)}")

# (2) all-large consecutive blocks {V-12..V}
n2 = 0; bmin = None; bat = None
for V in range(14, 500):
    S = list(range(V - 12, V + 1))
    if S[0] <= 0 or not is_covering(S) or not is_compressed(S): continue
    n2 += 1
    M = Mval(S, qcap(S))
    if bmin is None or M < bmin: bmin, bat = M, S
    if gmin is None or M < gmin: gmin, gat = M, S
    if ONE14 < M < ONE13: danger.append((M, S))
    if M < ONE14: below14.append((M, S))
print(f"(2) compressed covering consecutive-blocks {{V-12..V}}: {n2}; block min {bmin} (~{float(bmin) if bmin else 0:.6f}) at {bat}")

print("=" * 88)
print(f"GLOBAL MIN over compressed covering families tested: {gmin} (~{float(gmin):.6f}) at {gat}")
print(f"  >= 1/13 ? {gmin >= ONE13}   >= 7/89 ? {gmin >= SEV89}   >= 1/14 ? {gmin >= ONE14}")
print(f"DANGER ZONE (1/14,1/13) hits: {len(danger)}")
for M, S in sorted(danger)[:12]: print(f"    M={M} (~{float(M):.6f})  {S}")
print(f"BELOW 1/14 (would refute LRC14): {len(below14)}")
for M, S in sorted(below14)[:12]: print(f"    M={M} (~{float(M):.6f})  {S}  <-- !!!")
print()
print("READING: no compressed covering family in (1/14,1/13) and min >= 7/89 => hcomp floors")
print(">= 1/13 > 1/14 (mac-mini 'compressed => M>=1/13' holds); deep well 14/183<1/13 is DOMINANT.")
print("DONE")
