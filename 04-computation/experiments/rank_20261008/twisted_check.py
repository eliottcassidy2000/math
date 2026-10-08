#!/usr/bin/env python3
"""Independent check of audit G's twisted obstruction: if l does not divide d, every m_i = d (mod l), and
r_i = d psi(m_i) (mod l) for a homomorphism psi: <m_i> -> Z/l, then e_n - psi(M_n) = e_0 (mod l) along the pair chain.
Z_3 (1,1,5), r = (0,2,5): l = 2, psi(5) = 1, psi(1) = 0 (psi(M) = exponent of 5 in M, mod 2).  Odd offsets never merge."""
import random
from fractions import Fraction as Fr
d, m, r = 3, [1, 1, 5], [0, 2, 5]
assert all((m[i] * i + r[i]) % d == 0 for i in range(d))
def modd(q, n): return (q.numerator * pow(q.denominator, -1, n)) % n
def psi(M):
    num, den = M.numerator, M.denominator; k = 0
    while num % 5 == 0: num //= 5; k += 1
    while den % 5 == 0: den //= 5; k -= 1
    return k % 2
rnd = random.Random(1); viol = 0; steps = 0; merges = {1: 0, 2: 0, 3: 0, 4: 0}
for e0 in (1, 2, 3, 4):
    for _ in range(300):
        M, e = Fr(1), Fr(e0)
        for t in range(3000):
            j = rnd.randrange(d); i = (modd(M, d) * j + modd(e, d)) % d
            e = (m[i] * e + r[i] - Fr(m[i], m[j]) * M * r[j]) / d; M = M * Fr(m[i], m[j]); steps += 1
            if (modd(e, 2) - psi(M) - e0) % 2 != 0: viol += 1
            if M == 1 and e == 0: merges[e0] += 1; break
print(f"Z_3 (1,1,5), r=(0,2,5): {steps} chain steps, invariant violations {viol}; merges out of 300 by T=3000: {merges}")
