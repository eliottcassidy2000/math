#!/usr/bin/env python3
"""Growth of the polygonal order on pairs (1, m), n = 3: N_min(1,m) = least N with some j: ||j/N|| >= 1/3 and ||m j/N|| >= 1/3.
Records the running maxima (m, N_min) for m <= 30000 and compares with the least prime factor >= 3 of m-1 or m+1 (an upper bound)."""
from fractions import Fraction as F
def dist3(num, N):  # distance of num/N to Z, as Fraction
    r = num % N; return min(F(r, N), F(N - r, N))
def N_min_pair(m, Nmax=5000):
    for N in range(2, Nmax + 1):
        lo = -(-N // 3); hi = N - lo  # j in [ceil(N/3), N - ceil(N/3)] gives ||j/N|| >= 1/3
        for j in range(lo, hi + 1):
            if dist3(m * j, N) >= F(1, 3): return N
    return None
def lpf3(x):
    for p in range(3, x + 1):
        if x % p == 0: return p
    return None
best = 0; records = []
for m in range(2, 30001):
    N = N_min_pair(m)
    if N > best:
        best = N; records.append((m, N, min(v for v in (lpf3(m - 1), lpf3(m + 1)) if v is not None)))
        print(f"m = {m}: N_min(1,m) = {N}; bound min lpf>=3 of m-1, m+1 = {records[-1][2]}", flush=True)
print("records (m, N_min, bound):", records)
