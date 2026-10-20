#!/usr/bin/env python3
"""tribonacci_partner_depth_20261020.py -- opus-2026-10-20-S46
Depth profile of the tribonacci node T*_n (transitive with reversed Hamiltonian path): for every
Hamiltonian path P of T*_n (there are tribonacci(n) of them, and T*_n is rigid for n >= 4, so
the tilings of the node are in bijection with the paths), relabel along P, read the tiling, and
compute its symmetric-core depth (least k with the restriction to the k-th concentric
sub-staircase sigma-symmetric).  Usage: python3 ... NMAX
"""
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12

def tstar(n):
    adj = [[0] * n for _ in range(n)]
    for i in range(n - 1): adj[i][i + 1] = 1
    for i in range(n):
        for j in range(i + 2, n): adj[j][i] = 1
    return adj

def ham_paths(adj, n):
    out = []
    def rec(path, used):
        if len(path) == n: out.append(tuple(path)); return
        v = path[-1]
        for w in range(n):
            if not (used >> w) & 1 and adj[v][w]: rec(path + [w], used | (1 << w))
    for s in range(n): rec([s], 1 << s)
    return out

def depth_of_path(adj, n, p):
    # tile (i,j) (positions along the path, j >= i+2): bit = 1 iff p[i] -> p[j]
    bit = {}
    for i in range(n):
        for j in range(i + 2, n): bit[(i, j)] = adj[p[i]][p[j]]
    k = 0
    while True:
        cells = [(i, j) for (i, j) in bit if i >= k and j <= n - 1 - k]
        if all(bit[(i, j)] == bit[(n - 1 - j, n - 1 - i)] for (i, j) in cells): return k
        k += 1

for n in range(4, NMAX + 1):
    adj = tstar(n); paths = ham_paths(adj, n)
    prof = {}
    for p in paths:
        d = depth_of_path(adj, n, p); prof[d] = prof.get(d, 0) + 1
    print(f"n={n:2d} H={len(paths):5d} depth profile {dict(sorted(prof.items()))}")
    sys.stdout.flush()
