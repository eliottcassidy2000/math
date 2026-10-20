"""Arity-3 relational information: 3-hypertournaments (each triple carries one of the two cyclic orientations) up to
isomorphism, by Burnside over cycle types with the sign of the induced permutation on each triple-orbit; and k-coloured
complete graphs (pairs with k symmetric states, no twist)."""
from fractions import Fraction as F
from math import factorial, gcd
from itertools import combinations
from sympy.utilities.iterables import partitions
from collections import Counter
import sys, json
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 11
def zlam(parts):
    z = 1
    for L, m in Counter(parts).items(): z *= L**m * factorial(m)
    return z
def perm_from_parts(parts):
    sigma = []; start = 0
    for L in parts:
        for i in range(L): sigma.append(start + (i+1) % L)
        start += L
    return sigma
def hyper3_fixed(sigma, n):
    """number of 3-hypertournaments fixed by sigma: 2^(#orbits of triples whose stabiliser acts evenly), 0 if some orbit flips"""
    seen = set(); count = 0
    for t in combinations(range(n), 3):
        if t in seen: continue
        # orbit of the set t; track the induced permutation sign on return
        cur = list(t); order = list(t); k = 0
        while True:
            cur = [sigma[v] for v in cur]; k += 1
            if set(cur) == set(t):
                break
            seen.add(tuple(sorted(cur)))
        seen.add(t)
        # induced permutation of t after k steps: position i of t -> cur[i]; parity of this permutation on the 3-set
        perm = [t.index(cur[i]) for i in range(3)]
        # sign: a 3-cycle or identity is even (orientation preserved), a transposition is odd (flipped)
        inv = sum(1 for a in range(3) for b in range(a+1, 3) if perm[a] > perm[b])
        if inv % 2 == 1: return 0
        count += 1
    return 2**count
def graphs_k_fixed(sigma, n, k):
    """k-coloured complete graphs fixed by sigma: k^(#orbits of unordered pairs)"""
    seen = set(); orbits = 0
    for p in combinations(range(n), 2):
        if p in seen: continue
        cur = p
        while True:
            cur = tuple(sorted((sigma[cur[0]], sigma[cur[1]]))); seen.add(cur)
            if cur == p: break
        orbits += 1
    return k**orbits
H = []; G3 = []; G4 = []
for n in range(1, NMAX+1):
    h = F(0); g3 = F(0); g4 = F(0)
    for p in partitions(n):
        parts = [L for L, m in p.items() for _ in range(m)]
        sigma = perm_from_parts(parts); z = zlam(parts)
        h += F(hyper3_fixed(sigma, n), z); g3 += F(graphs_k_fixed(sigma, n, 3), z); g4 += F(graphs_k_fixed(sigma, n, 4), z)
    H.append(int(h)); G3.append(int(g3)); G4.append(int(g4))
print("3-hypertournaments (cyclic orientation on triples) n=1..:", H)
print("3-coloured complete graphs n=1..:", G3)
print("4-coloured complete graphs n=1..:", G4)
json.dump(dict(hyper3=H, graphs3=G3, graphs4=G4), open("hyper_and_coloured.json", "w"))
