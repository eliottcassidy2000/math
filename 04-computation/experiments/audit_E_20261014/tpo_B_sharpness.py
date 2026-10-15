#!/usr/bin/env python3
"""Audit E (2026-10-14), HYP-9279(B) sharpness for EVERY n >= 5: at the THM-4655 equality word w = 1 0 1^(m-3) 0 (m >= 4),
u_2 has P(u_2<z) = 2/3, P(u_2<u_1) = 1/3, and against every other vertex P(u_j<u_2) = 2^(m-4)/(3(1+2^(m-3))) (3 <= j <= m-1)
or 2^(m-3)/(3(1+2^(m-3))) (j = m), both < 1/3; so the vertex balance of u_2 is exactly 1/3.  Exact check m = 4..16."""
from fractions import Fraction as F
for m in range(4, 17):
    w = (1, 0) + (1,) * (m - 3) + (0,)
    full = (1 << m) - 1; H = 0; c2z = c21 = 0; cj = [0] * m
    for A in range(full + 1):
        B = full ^ A
        if A and w[A.bit_length() - 1] != 0: continue
        if B and w[(B & -B).bit_length() - 1] != 1: continue
        H += 1
        if A >> 1 & 1:
            c2z += 1
            if not A & 1: c21 += 1
        else:
            for j in range(2, m):
                if A >> j & 1: cj[j] += 1          # u_(j+1) before u_2
    D = 3 * (1 + 2 ** (m - 3))
    ok = (H == D and F(c2z, H) == F(2, 3) and F(c21, H) == F(1, 3)
          and all(F(cj[j], H) == F(2 ** (m - 4), D) for j in range(2, m - 1)) and F(cj[m - 1], H) == F(2 ** (m - 3), D))
    best = max([min(F(c2z, H), 1 - F(c2z, H)), min(F(c21, H), 1 - F(c21, H))] + [min(F(cj[j], H), 1 - F(cj[j], H)) for j in range(2, m)])
    print(f"m={m} (n={m+1}): closed forms {'ok' if ok else 'MISMATCH'}; vertex balance of u_2 = {best}")
