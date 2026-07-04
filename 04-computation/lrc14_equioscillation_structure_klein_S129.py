#!/usr/bin/env python3
"""
klein-2026-07-04-S129 (HYP-4090) - IS THE COVERING-MIN A 2-POINT (small-runner vs killer)
EQUIOSCILLATION FOR *EVERY* COVERING FAMILY?

mac-mini THM-618 PROVED the single-killer ladder {1..12,X} via a runner-1 <-> killer-X 2-point
equioscillation at t* = 1/13 - 1/(13(X+1)): the covering-min = base-optimum(1/13) minus a
killer-offset.  QUESTION (this session): does that SAME 2-point structure govern the covering-min
of EVERY covering family (multi-swap too, not just the single-killer ladder)?  If yes, THM-618's
mechanism generalizes toward the full covering-min bound -- the "parametric/geometric" route opus-S70
pointed to after killing the Delsarte dual.

For each covering family S (from klein-S128's multi-swap enumeration), compute exact M and an
optimal witness t*=a/Q, then find the BINDING runners (those v with ||v a/Q|| = M).  Characterize:
  - how many runners bind (2-point equioscillation? or more);
  - is a "small" runner (<=12) and a "killer" (the coverer of q=13, i.e. 13|v) among the binders?
  - report the binding pair per family.

Exact (Fractions).
"""
from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations

N = 14
DW = F(14, 183)

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval_and_argmax(S, Qcap):
    """Exact M and ONE optimal (a,Q). Returns (M, a, Q)."""
    best = F(0); bestaQ = (0, 1)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1:
                continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best:
                best = m; bestaQ = (a, Q)
    return best, bestaQ[0], bestaQ[1]

def binders(S, a, Q, M):
    """Runners v with ||v a/Q|| == M (the equioscillating set)."""
    return [v for v in S if F(cdist_q(v * a, Q), Q) == M]

def is_covering(S, n=N):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

def missing(A, n=N):
    return [q for q in range(2, n + 1) if not any(a % q == 0 for a in A)]

def qcap(S):
    return min(2 * max(S) + 2, 700)

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

# Build the same minimal-tightener covering families as S128, d=1,2,3.
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
print(f"testing {len(fams)} minimal-tightener covering families (d=1,2,3)")
print("=" * 92)
print("Q: is the covering-min a 2-point (small runner v<=12 + killer 13|v) equioscillation?")
print("=" * 92)

n_two_point = 0
n_has_small_and_killer = 0
n_has_runner1 = 0
examples = []
size_hist = {}
for S in sorted(fams, key=lambda s: Mval_and_argmax(list(s), qcap(list(s)))[0])[:60]:
    S = list(S)
    M, a, Q = Mval_and_argmax(S, qcap(S))
    B = binders(S, a, Q, M)
    size_hist[len(B)] = size_hist.get(len(B), 0) + 1
    has_small = any(v <= 12 for v in B)
    has_killer = any(v % 13 == 0 for v in B)          # killer = covers q=13
    has_r1 = 1 in B
    if len(B) == 2: n_two_point += 1
    if has_small and has_killer: n_has_small_and_killer += 1
    if has_r1: n_has_runner1 += 1
    if len(examples) < 22:
        examples.append((M, S, a, Q, B))

print(f"{'M':>10} {'t*=a/Q':>12} {'#bind':>6}  binders            family")
for M, S, a, Q, B in examples:
    print(f"{str(M):>10} {str(a)+'/'+str(Q):>12} {len(B):>6}  {str(B):<18} {S}")

print()
print(f"binder-set-size histogram (over the 60 lowest-M families): {dict(sorted(size_hist.items()))}")
print(f"2-point equioscillations: {n_two_point}/60")
print(f"has (small runner <=12 AND killer 13|v) binding: {n_has_small_and_killer}/60")
print(f"has runner 1 binding: {n_has_runner1}/60")
print()
print("READING: if the covering-min is ~always a 2-point equioscillation between a small runner")
print("and a killer (13|v), THM-618's single-killer mechanism is the GENERAL structure -- the")
print("covering-min = base-optimum - killer-offset for every covering family, and the extremal")
print("pair is (runner-1, smallest-killer-182). That is the parametric/geometric route (opus-S70).")
print("DONE")
