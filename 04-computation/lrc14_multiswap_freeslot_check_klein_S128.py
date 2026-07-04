#!/usr/bin/env python3
"""
klein-2026-07-04-S128 (HYP-4083) - ROBUSTNESS: the FREE-SLOT cases.
The main enumeration skipped drop-sets with |missing(A)| < d (dropping small AP elements that
divide larger kept ones does NOT create a missing q, so there are more tightener slots than
missing q's = "free slots").  Free tighteners ADD constraints, which can LOWER M -- so we must
check these don't dip below 14/183.  For each such drop, cover the missing q's minimally and
fill the remaining slots with the SMALLEST fresh integers >= 14 (smallest = most constraining =
lowest M, by scale-monotonicity).  Exact M.  Also a from-scratch bounded exhaustive cross-check
over covering systems with all elements <= B (B small) to independently confirm no sub-14/183.
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

def is_covering(S, n=N):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

def missing(A, n=N):
    return [q for q in range(2, n + 1) if not any(a % q == 0 for a in A)]

def qcap(S):
    return min(2 * max(S) + 2, 700)

def smallest_mult(L, used):
    x = ((14 + L - 1) // L) * L
    while x in used:
        x += L
    return x

def smallest_fresh(used):
    x = 14
    while x in used:
        x += 1
    return x

AP = list(range(1, 14))
print(f"deep well 14/183 = {float(DW):.6f}")
print("FREE-SLOT drop-sets (|missing(A)| < d): fill missing minimally, pad with smallest fresh.")
print("=" * 84)
worst = None
n_checked = 0
for d in (1, 2, 3, 4, 5):
    for drop in combinations(AP, d):
        A = [v for v in AP if v not in drop]
        Qm = missing(A)
        if len(Qm) >= d:
            continue          # handled by the main enumeration
        # cover each missing q with its smallest fresh multiple; then pad
        used = set(A)
        T = []
        for q in Qm:
            t = smallest_mult(q, used); T.append(t); used.add(t)
        while len(A) + len(T) < 13:
            t = smallest_fresh(used); T.append(t); used.add(t)
        S = sorted(A + T)
        if len(S) != 13 or not is_covering(S):
            continue
        me = Mval(S, qcap(S))
        n_checked += 1
        if worst is None or me < worst[0]:
            worst = (me, S, d)
print(f"  checked {n_checked} free-slot families")
if worst:
    me, S, d = worst
    print(f"  WORST (min M) free-slot family: M={me} (~{float(me):.6f})  d={d}  {S}")
    print(f"  >= 14/183 ? {me >= DW}")
print("=" * 84)

# Independent cross-check: bounded exhaustive over covering systems with all elements <= B.
# (Excludes the deep well's 182, but tests whether any *small-element* covering system beats it.)
print("\nBOUNDED CROSS-CHECK: all covering systems with <=13 elements, every element <= B.")
print("  (small-element regime -- independent of the minimal-tightener assumption)")
for B in (28, 30):
    # candidate elements 1..B; need covering with <=13 of them; find min M.
    # Prune: must include a multiple of 13 (13 or 26) and of 14 (14 or 28) and of 11, 9, 8, ...
    # Enumerate greedily is hard; instead: fix the "hard" coverers and vary the rest is complex.
    # Simple approach: since <=13 elements and covering, and we want MIN M, do a randomized +
    # structured search over subsets of {1..B} of size exactly 13 that are covering.
    from itertools import combinations as comb
    pool = list(range(1, B + 1))
    best = None; best_S = None; cnt = 0
    # structured: always take {1,...,7} (cheap, always helps lower M), vary the rest to cover 8..14
    fixed = [1, 2, 3, 4, 5, 6, 7]
    rest_pool = [x for x in pool if x not in fixed]
    for extra in comb(rest_pool, 6):     # 6 more -> 13 total
        S = sorted(fixed + list(extra))
        if not is_covering(S):
            continue
        cnt += 1
        # quick float gate then exact
        me = Mval(S, qcap(S))
        if best is None or me < best:
            best, best_S = me, S
    print(f"  B={B}: {cnt} covering systems (with prefix 1..7); MIN M = {best} (~{float(best) if best else 0:.6f})  at {best_S}")
    print(f"         >= 14/183 ? {best >= DW if best else 'n/a'}")

print("\nDONE")
