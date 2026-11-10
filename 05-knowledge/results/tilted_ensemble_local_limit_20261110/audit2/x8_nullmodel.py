#!/usr/bin/env python3
"""AUDIT 2 (auditor's own code).  The calibration "for a uniformly random quadratic form on F_2^m the left side of (4.3)
is exactly 2^(-m)" (README section 4.1a): which ensemble makes it true?

q(X) = sum_i ell_i x_i + sum_{i<j} beta_ij x_i x_j,  b_q = E_X (-1)^q(X).  Brute force over ALL forms, m = 2..5:
  (U) (ell, beta) uniform on F_2^m x F_2^C(m,2):                 E[b^2] = 2^(-m) exactly;
  (O) the constraints every tournament form satisfies (Lemma 4.5; q_R(V) = 0):  ell of odd weight, beta with an odd
      number of edges:                                           E[b^2] = 2^(-m) (1 + P(all beta-degrees even)),
      and P(all beta-degrees even | odd number of edges) = 2^(-(m-1)) for m >= 3.
The extra term is the antipodal pair of out-sets: q(X) + q(V - X) = beta(X, 1) + q(V), so with q(V) = 0
E_X (-1)^(q(X) + q(V - X)) = [all beta-degrees even]."""
from fractions import Fraction
from itertools import combinations

def run(m):
    prs = list(combinations(range(m), 2)); n2 = len(prs)
    sU = Fraction(0); cU = 0; sO = Fraction(0); cO = 0; evO = 0; cB = 0
    for bc in range(1 << n2):
        rows = [0] * m
        for k, (i, j) in enumerate(prs):
            if bc >> k & 1: rows[i] |= 1 << j; rows[j] |= 1 << i
        Q = [0] * (1 << m)                       # quadratic part
        for X in range(1, 1 << m):
            t = (X & -X).bit_length() - 1; Y = X & (X - 1)
            Q[X] = Q[Y] ^ (bin(rows[t] & Y).count('1') & 1)
        oddE = bin(bc).count('1') & 1
        alleven = all(bin(r).count('1') % 2 == 0 for r in rows)
        if oddE: cB += 1; evO += alleven
        for ell in range(1 << m):
            s = 0
            for X in range(1 << m):
                s += -1 if (Q[X] ^ (bin(ell & X).count('1') & 1)) else 1
            b2 = Fraction(s * s, 4 ** m)
            sU += b2; cU += 1
            if oddE and bin(ell).count('1') & 1: sO += b2; cO += 1
    eU = sU / cU; eO = sO / cO; p = Fraction(evO, cB)
    print("m = %d: (U) E[b^2] = %s (2^-m = %s, equal: %s) ; (O) E[b^2] = %s = 2^-m * %s ; P(all degrees even | odd edges) = %s ; 2^-m (1 + P) = %s, equal: %s"
          % (m, eU, Fraction(1, 2 ** m), eU == Fraction(1, 2 ** m), eO, eO * 2 ** m, p, Fraction(1, 2 ** m) * (1 + p), eO == Fraction(1, 2 ** m) * (1 + p)))

for m in (2, 3, 4, 5): run(m)
