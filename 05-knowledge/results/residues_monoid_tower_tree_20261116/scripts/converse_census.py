#!/usr/bin/env python3
"""converse_census.py n: for every tournament class of order n (gentourng), the Hamiltonian-path count H and whether the
class is self-converse (isomorphic to its arc reversal); per count value: classes, self-converse classes, converse pairs."""
import sys, subprocess, itertools
from collections import defaultdict
n = int(sys.argv[1])
lines = subprocess.run(['gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
def parse(s):
    adj = [[0]*n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i+1, n):
            if s[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    return adj
def ham(adj):
    N = 1 << n; dp = [[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v] = 1
    for m in range(1, N):
        for l in range(n):
            a = dp[m][l]
            if not a: continue
            for v in range(n):
                if not (m >> v) & 1 and adj[l][v]: dp[m | (1<<v)][v] += a
    return sum(dp[N-1])
def selfconverse(adj):
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    for p in itertools.permutations(range(n)):
        if all(adj[i][j] == rev[p[i]][p[j]] for i in range(n) for j in range(n) if i != j): return True
    return False
tab = defaultdict(lambda: [0, 0])
for s in lines:
    adj = parse(s); H = ham(adj); sc = selfconverse(adj)
    tab[H][0] += 1; tab[H][1] += sc
print(f"order {n}: {len(lines)} classes, {sum(v[1] for v in tab.values())} self-converse")
print("H: classes, self-converse, converse pairs")
for H in sorted(tab):
    c, sc = tab[H]; print(f"  {H}: {c} {sc} {(c - sc)//2}")
only_pairs = [H for H in tab if tab[H][1] == 0]; only_self = [H for H in tab if tab[H][1] == tab[H][0]]
print("counts attained only by converse pairs:", only_pairs)
print("counts attained only by self-converse classes:", only_self)
