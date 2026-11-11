#!/usr/bin/env python3
"""audit_tiling.py -- relation between Theorem T (3) (uniform labelled tournaments) and the closed form recorded in
THM-555 / HYP-2707 for the 'tiling' ensemble (tournaments on 0..n-1 containing the path 0 -> 1 -> ... -> n-1, the other
C(n-1,2) arcs free).  E[(-1)^c3] in both ensembles, exhaustively; c3 from the classical formula C(n,3) - sum C(s_i,2)
(validated against direct counts in audit_small.c / audit_rand.c).  Also: the score-parity vector is uniform on its coset
in the tiling ensemble iff the graph of free pairs (complement of a Hamiltonian path) is connected, i.e. iff n >= 4."""
from fractions import Fraction as F
from math import comb
for n in range(3, 9):
    free = [(i, j) for i in range(n) for j in range(i + 2, n)]
    tot = 0; sg = 0; par = {}
    for t in range(1 << len(free)):
        s = [0] * n
        for i in range(n - 1): s[i] += 1
        for b, (i, j) in enumerate(free):
            if (t >> b) & 1: s[i] += 1
            else: s[j] += 1
        c3 = comb(n, 3) - sum(x * (x - 1) // 2 for x in s)
        sg += -1 if c3 & 1 else 1; tot += 1
        key = tuple(x & 1 for x in s); par[key] = par.get(key, 0) + 1
    uni = len(par) == 2 ** (n - 1) and len(set(par.values())) == 1
    print(f"n={n}: tiling ensemble E[(-1)^c3] = {F(sg, tot)} ; uniform ensemble (Theorem T) 2^-floor((n-1)/2) = {F(1, 2**((n-1)//2))} ; "
          f"score-parity vector uniform on its coset in the tiling ensemble: {uni} ({len(par)} of {2**(n-1)} vectors occur)")
