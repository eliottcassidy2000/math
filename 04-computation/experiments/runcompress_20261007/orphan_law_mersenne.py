#!/usr/bin/env python3
"""Orphan law on the Mersenne line (mac-mini-2026-10-07 continuation).
M_K = 2^K - 1. An odd K is a deletion ORPHAN if no K' < K has the same odd-step count to 1 and Terras difference K - K'
(THM-4556: equal odd counts <=> equal-time merge before 1, verified <= 6000). Orphan fraction by window vs K^(-1/2)."""
import math
data = {}
for line in open('mersenne_sigma_12800.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
by_odd = {}
orph = []
for K in sorted(data):
    o, t = data[K]
    partners = [Kp for Kp in by_odd.get(o, []) if t - data[Kp][1] == K - Kp]
    if K % 2 == 1 and K >= 3: orph.append((K, len(partners) == 0, partners[-1] if partners else None))
    by_odd.setdefault(o, []).append(K)
# windows
wins = [(3, 50), (50, 100), (100, 200), (200, 400), (400, 800), (800, 1600), (1600, 3200), (3200, 6400), (6400, 12801)]
print("window      odd K   orphans  fraction   fraction*sqrt(K_mid)")
for a, b in wins:
    sel = [x for x in orph if a <= x[0] < b]
    n = len(sel); no = sum(1 for x in sel if x[1])
    mid = math.sqrt(a * b)
    print(f"[{a:5d},{b:5d})  {n:5d}   {no:5d}    {no/n:.4f}     {no/n*math.sqrt(mid):.3f}")
# also: partner distance (least shift) distribution for non-orphans
dists = [x[0] - x[2] for x in orph if not x[1] and x[0] >= 1000]
import collections
c = collections.Counter(dists)
print("least shift D (largest partner) for odd K >= 1000:", sorted(c.items())[:15], "...")
print("last orphans:", [x[0] for x in orph if x[1]][-15:])
