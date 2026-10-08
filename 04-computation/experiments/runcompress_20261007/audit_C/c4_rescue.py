#!/usr/bin/env python3
"""audit_C item 4: 3-adic rescue (Terras backward tree) of Mersenne numbers, independent implementation.

Backward Terras tree of n: children of x are 2x (always) and (2x-1)/3 (when x = 2 mod 3; it is then odd).
n is 'rescued within depth DMAX' if some node m < n occurs at backward depth <= DMAX (then T^d(m) = n).
Pruning (exact): along any backward path, y + 1 >= (2/3)(x + 1) per step (equality for the odd step), so a node x at
depth d can have a descendant < n within r = DMAX - d further steps only if 2^r (x + 1) < 3^r (n + 1).
(The session's orphan_rescue.py prunes with x < n*3^r//2^r + 1, which can only differ on a window of width (3/2)^r.)

Runs: all 53 orphans in [1000, 6400] and ALL 2700 odd K in [1000, 6400] (the session used a 300-sample).
"""
import sys, time
from collections import defaultdict, Counter

DMAX = 36


def rescue_depth(n):
    frontier = [n]
    for d in range(1, DMAX + 1):
        r = DMAX - d
        A, Bc = 3 ** r * (n + 1), 2 ** r
        nxt = []
        for x in frontier:
            y = 2 * x
            if Bc * (y + 1) < A:
                nxt.append(y)
            if x % 3 == 2:
                m = (2 * x - 1) // 3
                if m < n:
                    return d
                if Bc * (m + 1) < A:
                    nxt.append(m)
        frontier = nxt
        if not frontier:
            return None
    return None


data = {}
for line in open('c4_mersenne_12800_audit.txt'):
    K, o, t = map(int, line.split())
    data[K] = (o, t)
by_odd = defaultdict(list)
orphan = {}
for K in range(2, 6401):
    o, t = data[K]
    ps = [Kp for Kp in by_odd[o] if t - data[Kp][1] == K - Kp]
    if K % 2 == 1 and K >= 3:
        orphan[K] = not ps
    by_odd[o].append(K)

t0 = time.time()
Ks = [K for K in range(1001, 6401, 2)]
res = {}
for K in Ks:
    res[K] = rescue_depth((1 << K) - 1)
orph = [K for K in Ks if orphan[K]]
ro = sum(res[K] is not None for K in orph)
ra = sum(res[K] is not None for K in Ks)
print(f"odd K in [1000, 6400]: {len(Ks)}; orphans {len(orph)}")
print(f"orphans rescued within depth {DMAX}: {ro}/{len(orph)} = {ro/len(orph):.3f}")
print(f"all odd K rescued within depth {DMAX}: {ra}/{len(Ks)} = {ra/len(Ks):.3f}")
non = [K for K in Ks if not orphan[K]]
rn = sum(res[K] is not None for K in non)
print(f"non-orphans rescued: {rn}/{len(non)} = {rn/len(non):.3f}")
print("rescue depth distribution (all odd K):", sorted(Counter(res[K] for K in Ks if res[K] is not None).items()))
print("depth-3 rescues are exactly K = 5 mod 6:", all((res[K] == 3) == (K % 6 == 5) for K in Ks))
dbl = [K for K in orph if res[K] is None]
print(f"doubly uncovered (orphan and not rescued): {len(dbl)} of {len(Ks)} = {len(dbl)/len(Ks):.4f}:", dbl)
# Fisher-type comparison: hypergeometric p-value for rescued orphans <= observed
from math import comb
N, Kr, n_ = len(Ks), ra, len(orph)
p_le = sum(comb(Kr, x) * comb(N - Kr, n_ - x) for x in range(0, ro + 1)) / comb(N, n_)
print(f"hypergeometric P(rescued orphans <= {ro}) = {p_le:.3f}")
print(f"time {time.time()-t0:.1f}s")
