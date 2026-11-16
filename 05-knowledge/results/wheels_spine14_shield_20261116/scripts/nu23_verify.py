#!/usr/bin/env python3
"""nu23_verify.py: verify the factordb.com factorisation of Phi_736(2) (status FF, queried 2026-11-16 01:25 MST) locally:
each factor is prime (sympy.isprime), the product is Phi_736(2), and each prime has multiplicative order exactly 736, so
nu(736) = 3, c(5, 23) = 11 + 3 = 14 < 23 and a(23) = 6 (least covering period 2^6 * 23 = 1472)."""
from sympy import cyclotomic_poly, isprime, n_order
v = int(cyclotomic_poly(736, 2))
F = [76392570609857, 1335570346574631363954390476479681,
     89915373937922777877614505795576574280561874130802566460033]
assert F[0] * F[1] * F[2] == v, "product mismatch"
for p in F:
    print(p, "prime:", isprime(p), "ord_p(2) =", n_order(2, p))
    assert isprime(p) and n_order(2, p) == 736
print("Phi_736(2) has", len(str(v)), "digits; nu(736) =", len(F))
c4 = 11
print("c(5,23) =", c4 + len(F), "< 23 =>", "a(23) = 6, least period 1472")
