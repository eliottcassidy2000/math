#!/usr/bin/env python3
"""Decode the minimisers printed by tbal (gentourng ascii) and report score sequences, minimum feedback arc sets and the
transitive order, to see which near-transitive families attain the strong minima."""
import itertools, sys

def decode(s, n):
    adj = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    return adj

def fas(adj, n):
    best = None
    for perm in itertools.permutations(range(n)):
        pos = {v: i for i, v in enumerate(perm)}
        back = [(u, v) for u in range(n) for v in range(n) if adj[u][v] and pos[u] > pos[v]]
        if best is None or len(back) < len(best[1]):
            best = (perm, back)
            if len(back) <= 1: break
    return best

for line in sys.stdin:
    line = line.strip()
    if not line.startswith('minimiser:'): continue
    s = line.split()[1]
    n = int((1 + (1 + 8 * len(s)) ** 0.5) / 2)
    adj = decode(s, n)
    scores = sorted(sum(r) for r in adj)
    perm, back = fas(adj, n)
    print(f"n={n} {s}: scores {scores}; min feedback arc set {len(back)}: order {perm}, backward arcs {back}")
