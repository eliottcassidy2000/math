#!/usr/bin/env python3
"""prime_path_mod12.py (klein-2026-11-01) -- the mod-12 obstruction to prime square-sum Hamiltonian paths.
Two primes = 1 (mod 3) are never adjacent (THM-4860 (b)), and the graph is bipartite by residue mod 4 (THM-4860 (a)).
Hence, writing c_r for the number of primes = r (mod 12) below N:
  class 1  (= 1 mod 12, side A) is adjacent only to class 11' = {p = 11 mod 12} u {3};
  class 7  (= 7 mod 12, side B) is adjacent only to class 5'  = {p = 5 mod 12} u {2}.
In a Hamiltonian path every vertex has two path edges except the two endpoints, so counting the path edges at class 1
against the capacity of class 11' gives 2 c_1 - e_1 <= 2 (c_11 + 1) - e_11', i.e. c_1 <= c_11 + 2, and likewise
c_7 <= c_5 + 2; together with | |A| - |B| | <= 1 (A = {2} u {1 mod 4}).  This script (i) checks the adjacency claims by
brute force below 10^5, (ii) prints the four class counts and all count conditions at N = 615,731, (iii) scans every
prime bound N below BOUND for the count conditions (sieve only), and (iv) below 700,000 combines them with the
isolated/leaf conditions of prime_path_scan.py.   usage: python3 prime_path_mod12.py [BOUND=5000000]"""
import math, sys
BOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 5_000_000
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
# (i) brute-force adjacency classes below 10^5
P = primes_below(100000); Ps = set(P); top = P[-1]
viol = 0; edges = 0
for p in P:
    s = math.isqrt(p + 2)
    while s * s <= p + top:
        q = s * s - p
        if q > p and q in Ps:
            edges += 1
            for a, b in ((p, q), (q, p)):
                if a % 12 == 1 and not (b % 12 == 11 or b == 3): viol += 1
                if a % 12 == 7 and not (b % 12 == 5 or b == 2): viol += 1
        s += 1
print('(i) edges below 10^5: %d; violations of the class rules: %d' % (edges, viol)); assert viol == 0
def counts(Pn):
    c = {1: 0, 5: 0, 7: 0, 11: 0}
    for p in Pn:
        if p > 3: c[p % 12] += 1
    return c
def conds(c):
    A = c[1] + c[5] + 1; B = c[7] + c[11] + 1          # 2 in A, 3 in B
    return {'|A|-|B|': A - B, 'c1-c11': c[1] - c[11], 'c7-c5': c[7] - c[5],
            'mod4 ok': abs(A - B) <= 1, 'class1 ok': c[1] <= c[11] + 2, 'class7 ok': c[7] <= c[5] + 2}
# (ii) N = 615,731
P = primes_below(BOUND)
import bisect
k = bisect.bisect_left(P, 615731)
c = counts(P[:k]); print('(ii) N = 615731: class counts', c, conds(c))
# (iii) scan all prime bounds below BOUND for the three count conditions
c = {1: 0, 5: 0, 7: 0, 11: 0}; passing = []; first = {}
for i, p in enumerate(P[:-1]):
    if p > 3: c[p % 12] += 1
    N = P[i + 1]
    A = c[1] + c[5] + 1; B = c[7] + c[11] + 1
    ok4 = abs(A - B) <= 1; ok1 = c[1] <= c[11] + 2; ok7 = c[7] <= c[5] + 2
    if ok4 and ok1 and ok7 and N > 100:
        passing.append(N)
        if 'all' not in first: first['all'] = (N, dict(c), A - B)
print('(iii) bounds N in (100, %d) passing all three count conditions: %d; first: %s' % (BOUND, len(passing), first.get('all')))
if passing:
    # group into windows
    win = []; start = passing[0]; prev = passing[0]
    for N in passing[1:]:
        if N - prev > 20000: win.append((start, prev)); start = N
        prev = N
    win.append((start, prev)); print('     windows:', win[:12])
# (iv) below 700,000: combine with isolated/leaf conditions (incremental degrees)
LIM = 700000; P7 = [p for p in P if p < LIM]; Ps7 = set(P7)
deg = {}; iso = 0; leaves = 0; c = {1: 0, 5: 0, 7: 0, 11: 0}; full = []
def bump(v):
    global iso, leaves
    d = deg[v]
    if d == 0: iso -= 1
    elif d == 1: leaves -= 1
    deg[v] = d + 1
    if d + 1 == 1: leaves += 1
for i, p in enumerate(P7):
    deg[p] = 0; iso += 1
    if p > 3: c[p % 12] += 1
    s = math.isqrt(p) + 1
    while s * s < 2 * p:
        q = s * s - p
        if q in Ps7: bump(q); bump(p)
        s += 1
    if i + 1 < len(P7):
        N = P7[i + 1]; A = c[1] + c[5] + 1; B = c[7] + c[11] + 1
        if iso == 0 and leaves <= 2 and abs(A - B) <= 1 and c[1] <= c[11] + 2 and c[7] <= c[5] + 2:
            full.append((N, i + 1, leaves, A - B, c[1] - c[11], c[7] - c[5]))
print('(iv) bounds N < 700000 passing ALL five conditions (isolated, leaves, mod 4, class 1, class 7): %d; first 8: %s' % (len(full), full[:8]))
