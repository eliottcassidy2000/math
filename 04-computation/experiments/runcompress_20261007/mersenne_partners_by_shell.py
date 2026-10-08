#!/usr/bin/env python3
"""Mersenne line, full orbits: for odd K in [600, 6400], which deletion depths D give an equal-time partner
(same odd-step count, Terras difference D), split by the shell m = v2(K-1) (two-run length J = 1 + floor(m/2))."""
from collections import defaultdict
data = {}
for line in open('mersenne_sigma_6400.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
def v2(x): return (x & -x).bit_length() - 1
rows = defaultdict(lambda: [0, 0, 0, 0])   # m -> [count, D in {1,2}, D in {3,4}, any D <= K-1]
for K in range(601, 6401, 2):
    m = v2(K - 1); key = min(m, 6)
    o, t = data[K]
    p12 = any(data[K-D] == (o, t - D) for D in (1, 2))
    p34 = any(data[K-D] == (o, t - D) for D in (3, 4))
    anyp = any(data[K-D] == (o, t - D) for D in range(1, K - 1))
    r = rows[key]; r[0] += 1; r[1] += p12; r[2] += p34; r[3] += anyp
print("shell m (J = 1+floor(m/2))   #K   D in{1,2}   D in{3,4}   any D")
for m in sorted(rows):
    c, a, b, d = rows[m]
    lab = f"m={m}" if m < 6 else "m>=6"
    print(f"  {lab:6s} (J={1+m//2}{'+' if m >= 6 else ''})   {c:5d}    {a/c:.3f}       {b/c:.3f}     {d/c:.3f}")
