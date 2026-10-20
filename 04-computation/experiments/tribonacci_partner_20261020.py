#!/usr/bin/env python3
"""tribonacci_partner_20261020.py -- opus-2026-10-20-S46

The line partner of the transitive node.  The transitive tournament has one tiling (all tiles
forward); its complement tiling (all tiles backward) is the tournament T*_n with arcs i -> i+1
and j -> i for j >= i+2, i.e. the transitive tournament with its Hamiltonian path reversed.
Claims checked here for n <= NMAX:
  (1) H(T*_n) = tribonacci A000213: h(n) = h(n-1) + h(n-2) + h(n-3), h(1) = h(2) = 1, h(3) = 3
      (proof in the note: every Hamiltonian path is A, 0, 1, ..., k with A a Hamiltonian path on
      {k+1..n-1}; the paths ending at 0 number e(n) = h(n-1) - e(n-1) = A001590);
  (2) T*_n is rigid (|Aut| = 1) for n >= 4 and self-converse (i -> n-1-i is an anti-automorphism);
  (3) the blue multiplicity bm(T*_n) = number of Hamiltonian paths whose positional reversal is an
      anti-automorphism, against the guess 2^(floor(n/2)-1) + 1 (data 3,3,5,5,9 at n = 4..8).
"""
import sys, itertools
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 16

def tstar(n):
    adj = [[0] * n for _ in range(n)]
    for i in range(n - 1): adj[i][i + 1] = 1
    for i in range(n):
        for j in range(i + 2, n): adj[j][i] = 1
    return adj

def ham_paths(adj, n):
    paths = []
    def rec(path, used):
        if len(path) == n: paths.append(tuple(path)); return
        v = path[-1]
        for w in range(n):
            if not (used >> w) & 1 and adj[v][w]: rec(path + [w], used | (1 << w))
    for s in range(n): rec([s], 1 << s)
    return paths

def is_anti_aut(adj, n, rho):
    return all(adj[rho[u]][rho[v]] == adj[v][u] for u in range(n) for v in range(n) if u != v)

def aut_order(adj, n):
    return sum(1 for p in itertools.permutations(range(n)) if all(adj[p[u]][p[v]] == adj[u][v] for u in range(n) for v in range(n) if u != v))

trib = {1: 1, 2: 1, 3: 3}
for n in range(4, NMAX + 1): trib[n] = trib[n - 1] + trib[n - 2] + trib[n - 3]
print("n  H(T*)  tribonacci  e(n)=paths ending at 0  bm(T*)  2^(floor(n/2)-1)+1  SC  |Aut|")
for n in range(3, NMAX + 1):
    adj = tstar(n); paths = ham_paths(adj, n); h = len(paths)
    e = sum(1 for p in paths if p[-1] == 0)
    bm = 0
    for p in paths:
        rho = [0] * n
        for i, v in enumerate(p): rho[v] = p[n - 1 - i]
        if is_anti_aut(adj, n, rho): bm += 1
    rho0 = [n - 1 - i for i in range(n)]
    sc = is_anti_aut(adj, n, rho0)
    aut = aut_order(adj, n) if n <= 9 else None
    guess = 2 ** (n // 2 - 1) + 1
    print(f"{n:2d} {h:7d} {trib[n]:10d} {e:8d} {bm:14d} {guess:18d}  {sc}  {aut}")
    sys.stdout.flush()
