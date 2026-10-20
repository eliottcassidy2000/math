#!/usr/bin/env python3
"""Points of symmetry = tournaments all of whose Hamiltonian paths are PALINDROMIC (the reversal v_i -> v_{n+1-i} of the
path is an anti-automorphism, i.e. an isomorphism onto the converse). Conjecture (from the n <= 7 metagraph): the
irreducible points of symmetry are K_1, C_3 and the regular 5-tournament, and every point of symmetry is a transitive
composition of these in a palindromic pattern. Here: test the strong candidates of orders 7..15 that are circulant
(connection sets on Z_n up to multipliers), the Cayley tournaments of Z_3 x Z_3, and all strong tournaments on 7 vertices;
report H, the number of palindromic Hamiltonian paths, and whether all are palindromic."""
import itertools, subprocess, sys
from math import gcd
def ham_paths_list(n, adj):
    # all directed Hamiltonian paths (as tuples) by DFS; feasible while H <= ~2e5
    out = []
    def dfs(path, used):
        if len(path) == n: out.append(tuple(path)); return
        v = path[-1]
        for w in range(n):
            if not (used >> w) & 1 and adj[v][w]: dfs(path + [w], used | (1 << w))
    for v in range(n): dfs([v], 1 << v)
    return out
def is_antiauto(n, adj, sigma):
    return all((adj[u][v] == 0) or adj[sigma[v]][sigma[u]] for u in range(n) for v in range(n) if u != v)
def palindromic_count(n, adj, paths):
    c = 0
    for P in paths:
        sigma = [0] * n
        for i, v in enumerate(P): sigma[v] = P[n - 1 - i]
        if is_antiauto(n, adj, sigma): c += 1
    return c
def circulant_sets(n):
    units = [a for a in range(1, n) if gcd(a, n) == 1]
    classes = list(range(1, (n + 1) // 2)); seen = set(); reps = []
    for signs in itertools.product([0, 1], repeat=len(classes)):
        S = frozenset((x if s == 0 else n - x) for x, s in zip(classes, signs))
        if any(x in S and (n - x) in S for x in S): continue
        if len(S) != (n - 1) // 2: continue
        canon = min(tuple(sorted(a * x % n for x in S)) for a in units)
        if canon in seen: continue
        seen.add(canon); reps.append(S)
    return reps
print("=== circulant tournaments Z_n, n odd, 5 <= n <= 13 ===")
for n in [5, 7, 9, 11, 13]:
    for S in circulant_sets(n):
        adj = [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
        paths = ham_paths_list(n, adj); H = len(paths); pal = palindromic_count(n, adj, paths)
        print(f"n = {n}, S = {sorted(S)}: H = {H}, palindromic = {pal}, all palindromic: {pal == H}", flush=True)
print("=== Cayley tournaments of Z_3 x Z_3 (connection sets with S u -S = G - 0) ===")
G = [(a, b) for a in range(3) for b in range(3)]; idx = {g: i for i, g in enumerate(G)}
nonzero = [g for g in G if g != (0, 0)]
pairs = []; seenp = set()
for g in nonzero:
    ng = ((-g[0]) % 3, (-g[1]) % 3)
    key = frozenset([g, ng])
    if key not in seenp: seenp.add(key); pairs.append((g, ng))
seen_t = set()
for signs in itertools.product([0, 1], repeat=len(pairs)):
    S = [p[s] for p, s in zip(pairs, signs)]
    adj = [[0] * 9 for _ in range(9)]
    for i, g in enumerate(G):
        for s_ in S:
            h = ((g[0] + s_[0]) % 3, (g[1] + s_[1]) % 3); adj[i][idx[h]] = 1
    key = tuple(sorted(sum(r) for r in adj)), tuple(sorted(S))
    paths = ham_paths_list(9, adj); H = len(paths); pal = palindromic_count(9, adj, paths)
    print(f"S = {sorted(S)}: H = {H}, palindromic = {pal}, all: {pal == H}", flush=True)
print("=== all strong tournaments on 7 vertices (gentourng -c): points of symmetry among them ===")
out = subprocess.run(['gentourng', '-c', '-q', '7'], capture_output=True, text=True).stdout.split()
pos = 0; tot = 0
for line in out:
    n = 7; adj = [[0] * n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if line[k] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    paths = ham_paths_list(n, adj); H = len(paths); pal = palindromic_count(n, adj, paths); tot += 1
    if pal == H: pos += 1; print(f"  strong 7-vertex point of symmetry: H = {H}")
print(f"strong 7-vertex tournaments: {tot}; points of symmetry among them: {pos}")
