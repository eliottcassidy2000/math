#!/usr/bin/env python3
"""Audit E (2026-10-14), THM-4656 part 4: golden Paley primes (p = 3 mod 4 split in Q(sqrt5), residue root of x^2-x-1
generating QR_p) versus Shanks' Fibonacci primitive roots (a primitive root g mod p with g^2 = g + 1; OEIS A003147).
Claim (audit): for p = 3 mod 4 the two notions coincide, because the non-residue root is -1/(residue root), and an element of
odd order (p-1)/2 times -1 has order p-1."""
from sympy import primerange, is_primitive_root
def order(x, p):
    k, y = 1, x % p
    while y != 1:
        y = y * x % p; k += 1
    return k
gpp, fpr3, fpr_all = [], [], []
for p in primerange(3, 2000):
    roots = [x for x in range(p) if (x * x - x - 1) % p == 0]
    if len(roots) != 2:
        continue
    if any(is_primitive_root(r, p) for r in roots):
        fpr_all.append(p)
        if p % 4 == 3:
            fpr3.append(p)
    if p % 4 == 3:
        qr = [r for r in roots if pow(r, (p - 1) // 2, p) == 1]
        if order(qr[0], p) == (p - 1) // 2:
            gpp.append(p)
print("golden Paley primes < 2000:", len(gpp))
print("Fibonacci-primitive-root primes = 3 mod 4 < 2000:", len(fpr3), "; identical lists:", gpp == fpr3)
print("all FPR primes < 600 (compare OEIS A003147):", [p for p in fpr_all if p < 600])
