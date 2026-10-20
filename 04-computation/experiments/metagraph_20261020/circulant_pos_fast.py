#!/usr/bin/env python3
"""Strong points of symmetry of odd order are circulant (THM-4729: rotations by 2 along any Hamiltonian cycle are products of
two path-reversal anti-automorphisms). Test every circulant tournament of odd order n <= 17 (connection sets up to multipliers):
H by bitmask DP; palindromic paths = sum over involutory anti-automorphisms sigma of the number of sigma-palindromic paths,
each determined by its first half (one vertex per sigma-orbit, ending at the fixed point of sigma): DP over subsets of orbits.
Anti-automorphisms of a circulant: x -> -a x + c with a in Stab(S) (multiplier automorphisms) -- we enumerate all affine maps
x -> u x + c (gcd(u,n)=1) that are involutory anti-automorphisms and verify them directly."""
import itertools, sys
from math import gcd
def circulant_sets(n):
    units = [a for a in range(1, n) if gcd(a, n) == 1]
    classes = list(range(1, (n + 1) // 2)); seen = set(); reps = []
    for signs in itertools.product([0, 1], repeat=len(classes)):
        S = frozenset((x if s == 0 else n - x) for x, s in zip(classes, signs))
        canon = min(tuple(sorted(a * x % n for x in S)) for a in units)
        if canon in seen: continue
        seen.add(canon); reps.append(S)
    return reps
def ham_count(n, adj):
    full = (1 << n) - 1; dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            nb = adj[v]
            for w in nb:
                if not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[full])
def pal_count(n, adjm, S):
    # involutory anti-automorphisms among affine maps
    total = 0
    for u in range(1, n):
        if gcd(u, n) != 1: continue
        for c in range(n):
            sig = [(u * x + c) % n for x in range(n)]
            if any(sig[sig[x]] != x for x in range(n)): continue
            if not all(adjm[sig[y]][sig[x]] for x in range(n) for y in range(n) if adjm[x][y]): continue
            fixed = [x for x in range(n) if sig[x] == x]
            if len(fixed) != 1: continue   # odd n: exactly one fixed point for a palindromic reversal
            f = fixed[0]
            orbits = []; seen = set()
            for x in range(n):
                if x == f or x in seen: continue
                seen |= {x, sig[x]}; orbits.append((x, sig[x]))
            k = len(orbits); oid = {}
            for i, (a, b) in enumerate(orbits): oid[a] = i; oid[b] = i
            # DP: number of paths v1 -> ... -> v_k -> f using one vertex from each orbit; state (set of used orbits, last vertex)
            dp = {}
            for a, b in orbits:
                for v in (a, b): dp[(1 << oid[v], v)] = dp.get((1 << oid[v], v), 0) + 1
            full = (1 << k) - 1; cnt = 0
            for mask in range(1, 1 << k):
                for v in range(n):
                    c0 = dp.get((mask, v), 0)
                    if not c0: continue
                    if mask == full:
                        if adjm[v][f]: cnt += c0
                        continue
                    for w in range(n):
                        if w == f: continue
                        if adjm[v][w] and not (mask >> oid[w]) & 1:
                            dp[(mask | (1 << oid[w]), w)] = dp.get((mask | (1 << oid[w]), w), 0) + c0
            total += cnt
    return total
for n in [5, 7, 9, 11, 13, 15, 17]:
    for S in circulant_sets(n):
        adjm = [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
        adj = [[j for j in range(n) if adjm[i][j]] for i in range(n)]
        H = ham_count(n, adj); pal = pal_count(n, adjm, S)
        print(f"n = {n:2d} S = {sorted(S)}: H = {H}, palindromic = {pal}, ratio = {pal/H:.4f}{'  <-- POINT OF SYMMETRY' if pal == H else ''}", flush=True)
