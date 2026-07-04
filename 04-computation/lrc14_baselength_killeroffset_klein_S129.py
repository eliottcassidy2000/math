#!/usr/bin/env python3
"""
klein-2026-07-04-S129 (HYP-4090) - THE BASE-LENGTH LAW: covering-min = base-optimum minus
killer-offset, and the deep well (longest base) is the global minimum.

THM-618 (mac-mini): the single-killer family {1..12,X} hides at t*=1/13-eps, runner-1 vs killer-X
equioscillation, M = 1/13 - 1/(13(X+1)) = base-optimum(1/13) MINUS a killer-offset; min at X=182.

CLAIM (this session): this generalizes. For ANY covering family, let r = length of the longest
contiguous base {1,...,r} contained in it.  Then:
  (i)  M(family) <= 1/(r+1)   [monotone: dropping runners raises M; {1,...,r} sub-family has M=1/(r+1)]
  (ii) M(family) = 1/(r+1) - offset  with the covering-min approached as r -> 12 (longest base).
So the covering-min DECREASES as the base LENGTHENS, and the deep well (r=12, the longest possible
base since covering forces missing q=13,14 => can't have {1,...,13}) is the GLOBAL minimum 14/183.
This EXPLAINS klein-S128's "monotone increasing in swap-depth": swap-depth d shortens the base,
d = 12 - r_effective, so more swaps => shorter base => larger 1/(r+1) => larger covering-min.

Verify: per family, r=longest contiguous base, base-opt=1/(r+1), M, offset; confirm M<=1/(r+1)
and that min-M families have the longest base.  Exact.
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
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m
    return best

def is_covering(S, n=N):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

def missing(A, n=N):
    return [q for q in range(2, n + 1) if not any(a % q == 0 for a in A)]

def qcap(S):
    return min(2 * max(S) + 2, 700)

def base_len(S):
    """Longest r with {1,...,r} all in S."""
    r = 0
    Sset = set(S)
    while (r + 1) in Sset:
        r += 1
    return r

def partitions_into(items, k):
    items = list(items)
    if k <= 0: return
    if k == 1: yield [items]; return
    if len(items) < k: return
    first, rest = items[0], items[1:]
    for p in partitions_into(rest, k - 1):
        yield [[first]] + p
    for p in partitions_into(rest, k):
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]

def smallest_tightener(block, used):
    L = 1
    for q in block: L = lcm(L, q)
    x = ((14 + L - 1) // L) * L
    while x in used: x += L
    return x

AP = list(range(1, 14))
fams = set()
for d in (1, 2, 3):
    for drop in combinations(AP, d):
        A = [v for v in AP if v not in drop]
        Qm = missing(A)
        if len(Qm) < d: continue
        for part in partitions_into(Qm, d):
            used = set(A); T = []; ok = True
            for block in part:
                t = smallest_tightener(block, used)
                if t is None: ok = False; break
                T.append(t); used.add(t)
            if not ok: continue
            S = sorted(A + T)
            if len(S) == 13 and not missing(S):
                fams.add(tuple(S))

print(f"deep well 14/183 = {float(DW):.6f}")
print(f"{len(fams)} covering families. Testing: M <= 1/(base_len+1), and min-M => longest base.")
print("=" * 90)

# group by base length
by_r = {}
viol = 0
for S in fams:
    S = list(S)
    r = base_len(S)
    M = Mval(S, qcap(S))
    if M > F(1, r + 1):
        viol += 1
        print(f"  VIOLATION: M={M} > 1/(r+1)={F(1,r+1)} for base r={r}, {S}")
    by_r.setdefault(r, []).append((M, S))

print(f"\nM <= 1/(base_len+1) holds for all {len(fams)} families? {viol == 0}")
print("=" * 90)
print(f"{'base r':>6} {'base-opt 1/(r+1)':>16} {'min M at this r':>16} {'offset':>12} {'argmin family':>20}")
for r in sorted(by_r):
    rows = sorted(by_r[r])
    mM, mS = rows[0]
    off = F(1, r + 1) - mM
    print(f"{r:>6} {str(F(1,r+1)):>16} {str(mM):>16} {str(off):>12}  {mS}")

# global picture
allrows = sorted((M, base_len(list(S)), list(S)) for S in fams)
gM, gr, gS = allrows[0]
print("=" * 90)
print(f"GLOBAL min M = {gM} at base length r={gr}: {gS}")
print(f"  is it the LONGEST base (r=12, deep well)? {gr == 12}")
print(f"  base-optimum 1/(r+1) = 1/13 = {float(F(1,13)):.6f}; killer-offset = 1/13 - 14/183 = {F(1,13)-DW} = 1/2379")
print()
print("MONOTONICITY: min-M by base length (should DECREASE as r increases toward 12):")
prev = None; mono = True
for r in sorted(by_r):
    mM = min(m for m, _ in by_r[r])
    arrow = ""
    if prev is not None and mM > prev: arrow = "  (up)"
    else: arrow = "  (down/eq)" if prev is not None else ""
    print(f"  r={r:>2}: min M = {mM} (~{float(mM):.6f}){arrow}")
    prev = mM
print()
print("READING: if M <= 1/(base_len+1) universally and min-M sits at the LONGEST base r=12")
print("(the deep well), then covering-min = base-optimum - killer-offset (THM-618 mechanism) with")
print("the base-length as the order parameter. This EXPLAINS S128's swap-depth monotonicity")
print("(swap = base-shortening) via the killer-offset geometry, unifying THM-618 + S128.")
print("DONE")
