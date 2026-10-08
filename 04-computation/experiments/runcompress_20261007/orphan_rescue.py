#!/usr/bin/env python3
"""Are Mersenne deletion orphans rescued by 3-adic predecessor certificates? For n = 2^K - 1 search the Terras backward tree
(x -> 2x always; x -> (2x-1)/3 when x = 2 mod 3) to depth DMAX for a node m < n (then m -> ... -> n, a certificate).
Reports the rescue rate for orphans and for all odd K in [1000, 6400]."""
import random
data = {}
for line in open('mersenne_sigma_6400.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
by_odd = {}; orphans = []; allK = []
for K in sorted(data):
    o, t = data[K]
    has = any(t - data[Kp][1] == K - Kp for Kp in by_odd.get(o, []))
    if K % 2 == 1 and 1000 <= K <= 6400:
        allK.append(K)
        if not has: orphans.append(K)
    by_odd.setdefault(o, []).append(K)
DMAX = 36
def min_rescue_depth(n):
    frontier = [n]
    for d in range(1, DMAX + 1):
        nxt = []
        for x in frontier:
            y = 2 * x
            nxt.append(y)
            if x % 3 == 2:
                m = (2 * x - 1) // 3
                if m < n: return d
                nxt.append(m)
        # prune nodes that can no longer come back below n: a node x needs factor (2/3)^j with j <= DMAX - d
        lim = n * 3 ** (DMAX - d) // 2 ** (DMAX - d) + 1
        frontier = [x for x in nxt if x < lim]
        if not frontier: return None
    return None
rnd = random.Random(3)
samp = rnd.sample(allK, 300)
res_o = [min_rescue_depth((1 << K) - 1) for K in orphans]
res_a = [min_rescue_depth((1 << K) - 1) for K in samp]
ro = sum(r is not None for r in res_o) / len(res_o); ra = sum(r is not None for r in res_a) / len(res_a)
print(f"orphans: {len(orphans)}, rescued by a smaller backward node within depth {DMAX}: {ro:.3f}; depths {sorted(r for r in res_o if r)[:12]}")
print(f"all odd K sample: rescue rate {ra:.3f}")
print("doubly uncovered orphan exponents:", [K for K, r in zip(orphans, res_o) if r is None])
