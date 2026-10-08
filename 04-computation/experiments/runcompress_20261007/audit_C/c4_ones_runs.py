#!/usr/bin/env python3
"""audit_C item 4: longest post-run source ones-run (consecutive odd Terras steps after the initial K-run) of 2^K - 1,
orphans vs non-orphans, odd K in [1000, 6400]: all 53 orphans vs ALL 2647 non-orphans (the session used 159 sampled).
Mann-Whitney z with tie correction and a permutation p-value."""
import random, math
from collections import defaultdict
data = {}
for line in open('c4_mersenne_12800_audit.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
by = defaultdict(list); orphan = {}
for K in range(2, 6401):
    o, t = data[K]
    ps = [Kp for Kp in by[o] if t - data[Kp][1] == K - Kp]
    if K % 2 == 1 and K >= 3: orphan[K] = not ps
    by[o].append(K)
def longest(K):
    x = 3 ** K - 1          # after the initial K odd steps
    best = cur = 0
    while x != 1:
        if x & 1:
            cur += 1
            if cur > best: best = cur
            x = (3 * x + 1) >> 1
        else:
            cur = 0
            x >>= 1
    return best
Ks = list(range(1001, 6401, 2))
val = {K: longest(K) for K in Ks}
o = [val[K] for K in Ks if orphan[K]]; n = [val[K] for K in Ks if not orphan[K]]
allv = sorted(o + n)
# midranks
rank = {}
i = 0
while i < len(allv):
    j = i
    while j < len(allv) and allv[j] == allv[i]: j += 1
    rank[allv[i]] = (i + 1 + j) / 2
    i = j
R1 = sum(rank[x] for x in o); n1, n2 = len(o), len(n)
U = R1 - n1 * (n1 + 1) / 2
from collections import Counter
ties = Counter(allv); Nn = n1 + n2
sd = math.sqrt(n1 * n2 / 12 * ((Nn + 1) - sum(t ** 3 - t for t in ties.values()) / (Nn * (Nn - 1))))
z = (U - n1 * n2 / 2) / sd
rnd = random.Random(8)
obs = sum(o) / n1 - sum(n) / n2
pool = o + n; cnt = 0
for _ in range(4000):
    rnd.shuffle(pool)
    d = sum(pool[:n1]) / n1 - sum(pool[n1:]) / n2
    cnt += abs(d) >= abs(obs)
print(f"orphans n={n1}: mean longest post-run ones-run {sum(o)/n1:.2f}, median {sorted(o)[n1//2]}")
print(f"non-orphans n={n2}: mean {sum(n)/n2:.2f}, median {sorted(n)[n2//2]}")
print(f"Mann-Whitney z (orphans larger) = {z:.2f}; permutation p (two-sided, difference of means) = {cnt/4000:.3f}")
