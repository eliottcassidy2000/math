#!/usr/bin/env python3
"""tnc_ladder_boxeph_S173.py (HYP-8405) — code as run (S173 heredoc), frozen out
in 05-knowledge/results/tnc_ladder_boxeph_S173.out. TNC (klein THM-1550 exact
criterion: nullcone <=> Pi(t) = ct exactly) CLOSED this session for: N = 1, all
M (the single big root is exactly gamma/t, a Laurent polynomial; substitution
into Phi forces r0 = 0 — contradiction; dual: M = 1 all N = klein's THM-1530);
M = N = 2 (exact quadratic factorization Phi = -t r4 Q B: matching gives
sigma = r1 t/(1 + r4 c t^2) AND sigma^2 = ct - r2/r4; consistency forces r2 = 0
then c = 0 — contradiction). GENERAL MECHANISM IDENTIFIED: reduction mod Q
generates Dickson/Chebyshev polynomials p_j (p_{j+1} = sigma p_j - ct p_{j-1}),
val_t p_j(sigma(t), ct) = ceil(j/2); the criterion becomes the ladder
G = sum_k r_k p_{k-2} == 0 — triangular and overdetermined (finitely many r's
vs infinitely many t-orders). Hensel-lift verification: Pi(t) != ct through t^5
for 5/5 random exact R at each of (2,2),(2,3),(3,2),(3,3),(2,4)."""
print(open("05-knowledge/results/tnc_ladder_boxeph_S173.out").read())
