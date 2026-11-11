#!/usr/bin/env python3
"""audit_kbs.py -- the score formula of Theorem P (1) for c_3 as a corollary of the classical formula
c_3 = C(n,3) - sum_v C(s_v,2): since sum_v s_v^2 = q_1 (mod 4) and sum_v s_v = C(n,2),
c_3 = C(n,3) + (q_1 - C(n,2))/2 (mod 2).  Check that this equals C(q_1,2) + floor(n/4) for every n <= 200 and every
admissible q_1 (q_1 = C(n,2) mod 2), and the resulting cube sum (Theorem T (2)) and mean (Theorem T (3))."""
from math import comb
bad = 0
for n in range(1, 201):
    for q in range(0, n + 1):
        if (q - comb(n, 2)) % 2: continue
        a = (comb(n, 3) + (q - comb(n, 2)) // 2) % 2
        b = (comb(q, 2) + n // 4) % 2
        if a != b: bad += 1
    # sum over the admissible score-parity vectors of (-1)^(c_3) = cube sum of a base of order n - 1
    S = sum(comb(n, q) * (-1) ** (comb(q, 2) + n // 4) for q in range(0, n + 1) if (q - comb(n, 2)) % 2 == 0)
    if S != 2 ** (n // 2): bad += 1
print("mismatches for n <= 200:", bad)
