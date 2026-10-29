#!/usr/bin/env python3
"""THM-4848 (ii): explicit monomial conjugation. For m in H, mu_m(r) = tau^(c(r)) pi(r); check numerically
B_1[r, r'] = zeta^(m^-1 c(r)) B_(m^-1)[pi r, pi r'] zeta^(-m^-1 c(r')) for P7 (k = 2, 3) and P11 (k = 2)."""
import os, sys, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a3_runobs import circulant, observer, observer2_fast, blocks_data
def check(S, q, k):
    adj = circulant(S, q)
    M, tup, idx = (observer2_fast(q, adj) if k == 2 else observer(q, adj, k))
    R, Ms = blocks_data(M, tup, idx, q)
    Ridx = {r: i for i, r in enumerate(R)}
    zeta = cmath.exp(2j * cmath.pi / q)
    def B(j):
        return [[sum(Ms[a][b][s] * zeta ** (j * s) for s in range(q)) for b in range(len(R))] for a in range(len(R))]
    H = [m for m in range(1, q) if {(m * x) % q for x in S} == set(S)]
    worst = 0.0
    B1 = B(1)
    for m in H:
        mi = pow(m, -1, q)
        Bmi = B(mi)
        pi, c = {}, {}
        for r in R:
            mr = tuple((m * x) % q for x in r)
            cc = mr[0]                       # mu_m(r) = tau^cc (rep with first coordinate 0)
            rep = tuple((x - cc) % q for x in mr)
            pi[r], c[r] = rep, cc
        for a, r in enumerate(R):
            for b, rp in enumerate(R):
                pred = zeta ** (mi * c[r]) * Bmi[Ridx[pi[r]]][Ridx[pi[rp]]] * zeta ** (-mi * c[rp])
                worst = max(worst, abs(pred - B1[a][b]))
    return H, worst
for S, q, k in (({1, 2, 4}, 7, 2), ({1, 2, 4}, 7, 3), ({1, 3, 4, 5, 9}, 11, 2), ({1, 2, 3, 5, 6, 9}, 13, 2)):
    H, w = check(S, q, k)
    print(f"q={q} S={sorted(S)} k={k}: H={H}; max |B_1 - Delta P B_(m^-1) P^-1 Delta^-1| = {w:.2e}")
