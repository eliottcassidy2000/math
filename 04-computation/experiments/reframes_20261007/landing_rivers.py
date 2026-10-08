#!/usr/bin/env python3
"""The +1-barrier children y*_D = 2/3^D - 1 satisfy y*_(D-1) = 3 y*_D + 2 (the Mersenne D = 1 relation), so their orbits
coalesce in rivers and land together.  Key: (N_D, s_0(D)) -- y*_(D-1) = 3 y*_D + 2 holds at EQUAL times, so an equal-time merge
forces equal landing value and equal (landing time - D).  Count rivers, runs, and the tail per river."""
import math
from collections import defaultdict
def landing(D):
    a, m, s = 2 - 3**D, D, 0
    while m > 0:
        if a & 1: a = (a + 3**(m-1)) // 2; m -= 1
        else: a //= 2
        s += 1
    return a, s
DMAX = 4500
key = {}
for D in range(3, DMAX):
    N, s = landing(D); key[D] = (N, s)
riv = defaultdict(list)
for D in range(3, DMAX): riv[key[D]].append(D)
same_next = sum(1 for D in range(3, DMAX - 1) if key[D] == key[D + 1])
print(f"D in [3,{DMAX}): {len(riv)} rivers (distinct (N_D, s0)); consecutive D with equal key: {same_next}/{DMAX-4}")
sizes = sorted((len(v) for v in riv.values()), reverse=True)
print("largest rivers:", sizes[:12])
Nr = [k[0] for k in riv]
C = 2.867
for x in (4, 16, 64, 256):
    print(f"   per-river P(N > {x:3d}) = {sum(1 for N in Nr if N > x)/len(Nr):.4f}   stationary C/x = {C/x:.4f}")
big = sorted(((k[0], min(v), max(v), len(v)) for k, v in riv.items()), reverse=True)[:8]
print("rivers with the largest N (N, minD, maxD, size):", big)
