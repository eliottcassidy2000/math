#!/usr/bin/env python3
"""Do Mersenne deletion orphans carry long post-run source ones-runs (positive-drift re-excursions, THM-4603 (3))?
Compare, for odd K in [1000, 6400]: the longest run of consecutive odd Terras steps in the orbit of 2^K - 1 after its initial
ones-run, for orphans vs a random sample of non-orphans."""
import random, statistics
data = {}
for line in open('mersenne_sigma_6400.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
by_odd = {}
orph = set(); non = []
for K in sorted(data):
    o, t = data[K]
    has = any(t - data[Kp][1] == K - Kp for Kp in by_odd.get(o, []))
    if K % 2 == 1 and 1000 <= K <= 6400:
        (orph.add(K) if not has else non.append(K))
    by_odd.setdefault(o, []).append(K)
rnd = random.Random(1)
sample = rnd.sample(non, 3 * len(orph))
def max_post_run(K):
    n = (1 << K) - 1
    # skip the initial ones-run (K odd Terras steps)
    for _ in range(K): n = (3*n + 1) >> 1
    best = cur = 0; steps = 0; long_runs = 0
    while n != 1:
        if n & 1:
            cur += 1; n = (3*n + 1) >> 1
        else:
            if cur >= 12: long_runs += 1
            best = max(best, cur); cur = 0; n >>= 1
        steps += 1
    return max(best, cur), long_runs, steps
res_o = [max_post_run(K) for K in sorted(orph)]
res_n = [max_post_run(K) for K in sample]
def summ(res):
    return (statistics.mean(r[0] for r in res), statistics.median(r[0] for r in res),
            statistics.mean(r[1] for r in res), statistics.mean(r[2] for r in res))
so, sn = summ(res_o), summ(res_n)
print(f"orphans (n={len(res_o)}): mean max post-run ones-run {so[0]:.2f}, median {so[1]}, mean #runs>=12 {so[2]:.2f}, mean Terras length {so[3]:.0f}")
print(f"non-orphans (n={len(res_n)}): mean max post-run ones-run {sn[0]:.2f}, median {sn[1]}, mean #runs>=12 {sn[2]:.2f}, mean Terras length {sn[3]:.0f}")
# rank test (Mann-Whitney U, normal approx)
allv = sorted([(r[0], 'o') for r in res_o] + [(r[0], 'n') for r in res_n])
ranks = {}
i = 0
while i < len(allv):
    j = i
    while j < len(allv) and allv[j][0] == allv[i][0]: j += 1
    for k in range(i, j): ranks.setdefault(k, (i + j + 1) / 2)
    i = j
Ro = sum(ranks[k] for k in range(len(allv)) if allv[k][1] == 'o')
n1, n2 = len(res_o), len(res_n)
U = Ro - n1*(n1+1)/2; mu = n1*n2/2; sd = (n1*n2*(n1+n2+1)/12) ** 0.5
print(f"Mann-Whitney z (orphans larger) = {(U - mu)/sd:.2f}")
