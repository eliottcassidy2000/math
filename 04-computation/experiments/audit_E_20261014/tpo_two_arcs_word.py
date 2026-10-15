#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4655 1(e) "two is attained for every m >= 3": the word w = 1 0^(m-1) (reverse of the nearly
transitive tournament, H = 2^(m-1) + 1) has exactly two balanced arcs on its Hamiltonian cycle z u_1 ... u_m z, namely
(u_1,u_2) with P(u_2<u_1) = 2^(m-2)/(2^(m-1)+1) and (u_m,z) with P(u_m<z) = (2^(m-2)+1)/(2^(m-1)+1); the other arcs have
P(z<u_1) = 2^(m-1)/(2^(m-1)+1) and P(u_(i+1)<u_i) = 2^(m-3)/(2^(m-1)+1) (i >= 2).  Exact check via the bijection, m = 3..16."""
from fractions import Fraction as F
for m in range(3, 17):
    w = (1,) + (0,) * (m - 1)
    full = (1 << m) - 1; H = 0
    zu1 = u2u1 = umz = 0; back = [0] * m
    for A in range(full + 1):
        B = full ^ A
        if A and w[A.bit_length() - 1] != 0: continue
        if B and w[(B & -B).bit_length() - 1] != 1: continue
        H += 1
        if not A & 1: zu1 += 1                      # u_1 in B: z before u_1
        if (A >> 1 & 1) and not (A & 1): u2u1 += 1  # u_2 before u_1
        if A >> (m - 1) & 1: umz += 1               # u_m before z
        for i in range(1, m - 1):                   # u_(i+2) before u_(i+1)  (1-based i+1 >= 2)
            if (A >> (i + 1) & 1) and not (A >> i & 1): back[i] += 1
    D = 2 ** (m - 1) + 1
    ok = (H == D and F(zu1, H) == F(2 ** (m - 1), D) and F(u2u1, H) == F(2 ** (m - 2), D) and F(umz, H) == F(2 ** (m - 2) + 1, D)
          and all(F(back[i], H) == F(2 ** (m - 3), D) for i in range(1, m - 1)))
    bal = lambda p: F(1, 3) <= p <= F(2, 3)
    nb = bal(F(zu1, H)) + bal(F(u2u1, H)) + bal(F(umz, H)) + sum(bal(F(back[i], H)) for i in range(1, m - 1))
    print(f"m={m}: H={H}; formulas {'ok' if ok else 'MISMATCH'}; balanced cycle arcs = {nb}")
