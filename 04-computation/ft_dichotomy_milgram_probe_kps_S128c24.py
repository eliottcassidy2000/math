#!/usr/bin/env python3
"""ft_dichotomy_milgram_probe_kps_S128c24.py -- kind-pasteur S128 cont.24.
(A) Milgram probe: Gauss sums over the E4 level lattice, G_n(w) = sum_T e(w*E4(T)/16), w=1,2 --
    testing whether the mod-8/16 level laws carry discriminant-form (8th-root-of-unity) structure.
(B) Aut-order spectrum n<=6 (odd orders only -- the Feit-Thompson dichotomy's data side)."""
import sys, cmath
from math import comb, pi
from itertools import permutations
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)
for n in [4, 5, 6]:
    m = comb(n, 2)
    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]
    G1 = 0j; G2 = 0j
    aut = Counter()
    perms = list(permutations(range(n)))
    for mask in range(1 << m):
        B = [[False]*n for _ in range(n)]
        for k, (u, v) in enumerate(pairs):
            if (mask >> k) & 1: B[u][v] = True
            else: B[v][u] = True
        d = [sum(r) for r in B]
        E4 = sum((2*x - (n-1))**2 for x in d)
        G1 += cmath.exp(2j*pi*E4/16)
        G2 += cmath.exp(2j*pi*E4/8)
        if n <= 5:
            a = 0
            for pm in perms:
                if all(B[u][v] == B[pm[u]][pm[v]] for u in range(n) for v in range(n) if u != v):
                    a += 1
            aut[a] += 1
    N = 1 << m
    print("n=%d: G(1/16)/2^m = %.4f + %.4fi  |G|/2^m = %.4f ; G(1/8)/2^m = %.4f%+.4fi |.|=%.4f"
          % (n, G1.real/N, G1.imag/N, abs(G1)/N, G2.real/N, G2.imag/N, abs(G2)/N))
    if n <= 5:
        print("   Aut-order spectrum (all odd = Feit-Thompson solvable): %s" % dict(aut))
print("DONE")
