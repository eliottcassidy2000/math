#!/usr/bin/env python3
"""Hypothesis battery 1 on the exact Mersenne orbit data (K <= 12800).
H1 tuning: jumps D = sigma(M_K) - sigma(M_{K-1}) != 0 satisfy ||D log2 3|| small (vs random integers of the same size).
   Identity: for n reaching 1, eps(n) := sigma_T(n) - sigma(n) log2 3 - log2 n = sum over odd orbit values x>1 of log2(1 + 1/(3x)) in (0, ...).
H6 braid width: number of distinct rivers among K in sliding windows.
H2 orphans vs two-run length v2(K-1)."""
import math, random
from collections import defaultdict, Counter
L3 = math.log2(3)
data = {}
for line in open('../runcompress_20261007/mersenne_sigma_12800.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
Ks = sorted(data)
def dist(x): return abs(x - round(x))
jumps = [data[K][0] - data[K-1][0] for K in Ks if K-1 in data and data[K][0] != data[K-1][0]]
dj = [dist(abs(j) * L3) for j in jumps]
rnd = random.Random(1)
dr = [dist(rnd.randint(1, 3000) * L3) for _ in range(20000)]
print(f"H1: {len(jumps)} nonzero jumps; median ||D log2 3|| = {sorted(dj)[len(dj)//2]:.4f} (random integers: {sorted(dr)[len(dr)//2]:.4f}); share < 0.05: {sum(d<0.05 for d in dj)/len(dj):.3f} vs {sum(d<0.05 for d in dr)/len(dr):.3f}")
# eps per K and per river
eps = {K: data[K][1] - data[K][0]*L3 - math.log2((1 << K) - 1) if K < 1000 else data[K][1] - data[K][0]*L3 - K for K in Ks}
cls = defaultdict(list)
for K in Ks: cls[(data[K][0], data[K][1] - K)].append(K)
spread = max(max(eps[K] for K in v) - min(eps[K] for K in v) for v in cls.values())
allv = sorted(eps.values())
print(f"    eps(M_K) range [{allv[0]:.4f}, {allv[-1]:.4f}], median {allv[len(allv)//2]:.4f}; max spread of eps inside a river {spread:.2e}")
# check identity eps = sum log2(1+1/(3x)) for a few K
for K in (101, 1001, 5001):
    n = (1 << K) - 1; s = 0.0; x = n; o = 0; t = 0
    while x != 1:
        if x & 1:
            s += math.log2(1 + 1/(3*x)); x = (3*x + 1) >> 1; o += 1; t += 1
        else:
            x >>= 1; t += 1
    print(f"    K={K}: eps from (sigma_T, sigma) = {t - o*L3 - (K if K >= 1000 else math.log2(n)):.6f}; reciprocal sum = {s:.6f}")
# H6 braid width
for W in (50, 200, 1000):
    widths = []
    for a in range(1000, 12800 - W, W):
        widths.append(len({(data[K][0], data[K][1]-K) for K in range(a, a+W)}))
    print(f"H6: distinct rivers per window of {W} consecutive K (K in [1000,12800]): mean {sum(widths)/len(widths):.2f}, max {max(widths)}")
# H2 orphans vs v2(K-1)
def v2(x): return (x & -x).bit_length() - 1
by_odd = {}; rows = []
for K in Ks:
    o, t = data[K]
    has = any(t - data[Kp][1] == K - Kp for Kp in by_odd.get(o, []))
    if K % 2 == 1 and K >= 200: rows.append((min(v2(K-1), 5), not has))
    by_odd.setdefault(o, []).append(K)
c = Counter(m for m, _ in rows); co = Counter(m for m, o in rows if o)
print("H2: orphan rate by m = v2(K-1) (K >= 200):", {m: f"{co[m]}/{c[m]}={co[m]/c[m]:.4f}" for m in sorted(c)})
