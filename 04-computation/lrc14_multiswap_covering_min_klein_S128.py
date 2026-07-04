#!/usr/bin/env python3
"""
klein-2026-07-04-S128 (HYP-4083) - IS THE DEEP WELL THE GLOBAL COVERING-MIN?
Extending the ONE-SWAP map (S127) to MULTI-SWAP 13-element covering systems.

BACKGROUND. LRC(14) covering-min crux: min over 13-element covering systems S (for each
q in {2,...,14}, some v in S with q|v) of M(S)=max_t min_{v in S} ||v t||.  Claim: min = 14/183
= n/Phi_6(n), the deep well {1,...,12,182}.  S127 mapped the ONE-SWAP stratum; this session
extends to MULTI-SWAP.

STRUCTURE.  Every 13-element set = A u T, A subset {1,...,13} (AP core), T subset {14,...}
(tighteners), |A|+|T|=13.  A d-swap drops d AP elements.  A = {1,...,12} forces the single
tightener divisible by lcm(13,14)=182 (deep-well ladder).  Multi-swap = drop >= 2.

TWO monotonicities (both from S127, both used here):
  (i) SCALE a tightener up (same count): M INCREASES (more room to dodge a finer comb).
      => the covering-MIN uses the SMALLEST tighteners realizing each covering role.
  (ii) ADD a tightener (more count, same others): M DECREASES (extra constraint).
So the min-M d-swap for a given A uses d tighteners, each the SMALLEST integer >=14 realizing
its covering role.  The roles = a partition of the missing q's into d blocks; tightener(block)
= smallest multiple of lcm(block) that is >= 14 and not already present.  We enumerate all such
partitions (exhaustive over minimal-tightener families, which contain the min), exact M.

Exact (Fractions).  Save to 05-knowledge/results/.
"""
from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations

N = 14
DW = F(14, 183)

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval(S, Qcap):
    best = F(0)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1:
                continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best:
                best = m
    return best

def Mval_float(S, Qcap):
    best = 0.0
    for Q in range(2, Qcap + 1):
        invQ = 1.0 / Q
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1:
                continue
            m = min(cdist_q(v * a, Q) for v in S) * invQ
            if m > best:
                best = m
    return best

def is_covering(S, n=N):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

def missing(A, n=N):
    return [q for q in range(2, n + 1) if not any(a % q == 0 for a in A)]

def qcap(S):
    return min(2 * max(S) + 2, 700)

def partitions_into(items, k):
    """All set-partitions of `items` (a list) into exactly k nonempty blocks."""
    items = list(items)
    if k <= 0:
        return
    if k == 1:
        yield [items]; return
    if len(items) < k:
        return
    first, rest = items[0], items[1:]
    # first alone, partition rest into k-1
    for p in partitions_into(rest, k - 1):
        yield [[first]] + p
    # first joins one of the k blocks of rest
    for p in partitions_into(rest, k):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]

def smallest_tightener(block, used):
    """Smallest integer >= 14, divisible by lcm(block), not already in `used`."""
    L = 1
    for q in block:
        L = lcm(L, q)
    x = ((14 + L - 1) // L) * L      # first multiple of L that is >= 14
    while x in used:
        x += L
    return x

print(f"deep well target 14/183 = {float(DW):.6f}")
DWset = list(range(1, 13)) + [182]
print(f"sanity M(deep well)      = {Mval(DWset, qcap(DWset))}")
print("=" * 84)

AP = list(range(1, 14))
global_min, global_at = None, None
records = []      # (M_exact, S, d)

for d in (1, 2, 3, 4):
    dmin, dmin_at = None, None
    fams = []
    for drop in combinations(AP, d):
        A = [v for v in AP if v not in drop]
        Qm = missing(A)
        if len(Qm) < d:
            continue          # would need free slots (non-extremal); skip (see note)
        for part in partitions_into(Qm, d):
            used = set(A)
            T = []
            ok = True
            for block in part:
                t = smallest_tightener(block, used)
                if t is None:
                    ok = False; break
                T.append(t); used.add(t)
            if not ok:
                continue
            S = sorted(A + T)
            if len(S) != 13 or missing(S):
                continue
            fams.append(tuple(S))
    fams = set(fams)
    # exact M for all (these are minimal-tightener families; count is small)
    best = None; best_S = None; below = []
    for s in fams:
        s = list(s)
        me = Mval(s, qcap(s))
        if best is None or me < best:
            best, best_S = me, s
        if me < DW:
            below.append((me, s))
        records.append((me, s, d))
    dmin, dmin_at = best, best_S
    print(f"\n--- d={d}-swap : {len(fams)} minimal-tightener covering families ---")
    print(f"  MIN M = {dmin} (~{float(dmin) if dmin else 0:.6f})  at {dmin_at}")
    print(f"  >= 14/183 ? {dmin >= DW if dmin else 'n/a'}")
    if below:
        print(f"  *** {len(below)} FAMILIES BELOW 14/183: ***")
        for me, s in sorted(below)[:10]:
            print(f"      M={me} (~{float(me):.6f})  {s}")
    if global_min is None or (dmin is not None and dmin < global_min):
        global_min, global_at = dmin, dmin_at

print("\n" + "=" * 84)
print(f"GLOBAL MIN over d=1..4 minimal-tightener covering systems:")
print(f"  M = {global_min} (~{float(global_min):.6f})  at {global_at}")
print(f"  equals deep well 14/183 ? {global_min == DW}")
print("=" * 84)

# The 12 smallest covering families overall, with their swap-depth:
records.sort(key=lambda r: r[0])
print("\nThe 15 SMALLEST-M covering families found (M, d-swap, family):")
seen = set()
shown = 0
for me, s, d in records:
    key = tuple(s)
    if key in seen:
        continue
    seen.add(key)
    print(f"  M={str(me):>12} (~{float(me):.6f})  d={d}  {s}")
    shown += 1
    if shown >= 15:
        break

print("\nDONE")
