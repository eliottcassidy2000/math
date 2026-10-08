#!/usr/bin/env python3
"""Probe: the odd-step count sigma(M_K) of M_K = 2^K - 1 as a function of K (data K <= 12800).
Is it monotone? What are the jumps? Are rivers (partner classes) intervals?"""
from collections import defaultdict, Counter
data = {}
for line in open('../runcompress_20261007/mersenne_sigma_12800.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
Ks = sorted(data)
dec = [(K, data[K-1][0], data[K][0]) for K in Ks if K-1 in data and data[K][0] < data[K-1][0]]
print(f"K in [{Ks[0]},{Ks[-1]}]: sigma(M_K) < sigma(M_(K-1)) occurs {len(dec)} times; first: {dec[:8]}")
# running maximum structure
inc = [data[K][0] - data[K-1][0] for K in Ks if K-1 in data]
c = Counter(inc)
print("distribution of sigma(M_K) - sigma(M_(K-1)) (most common):", c.most_common(12))
# rivers = classes of equal (odd count, Terras length - K)  [partners share sigma_T - K]
cls = defaultdict(list)
for K in Ks: cls[(data[K][0], data[K][1] - K)].append(K)
sizes = Counter(len(v) for v in cls.values())
nonint = [v for v in cls.values() if len(v) > 1 and v[-1] - v[0] + 1 != len(v)]
print(f"river classes: {len(cls)}; size histogram (top): {sorted(sizes.items())[:10]} ... max {max(len(v) for v in cls.values())}")
print(f"rivers that are not intervals: {len(nonint)}; examples {[v[:8] for v in nonint[:4]]}")
# fraction of K whose river is an interval
interval_K = sum(len(v) for v in cls.values() if v[-1]-v[0]+1 == len(v))
print(f"K lying in interval rivers: {interval_K/len(Ks):.3f}")
# is sigma(M_K) monotone in K when restricted to river representatives? sigma along river minima
mins = sorted((v[0], k[0]) for k, v in cls.items())
nonmono = sum(1 for (a, s1), (b, s2) in zip(mins, mins[1:]) if s2 < s1)
print(f"river minima in K order: {len(mins)}; sigma decreases between consecutive river minima {nonmono} times")
