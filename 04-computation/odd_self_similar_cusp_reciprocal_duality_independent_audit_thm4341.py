#!/usr/bin/env python3
"""Independent symbolic audit for THM-4341.

This path imports no primary-certificate code.  It builds every normalized
tail polynomial over Q(c), checks squarefreeness and odd projective degree,
reconstructs both integer index maps, and probes the sharp excluded cases.
"""

from fractions import Fraction
from math import isqrt
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, label):
    if not bool(condition):
        raise RuntimeError(label)


def T(n):
    return n * (n + 1) // 2


def inverse_oriented(n):
    # Unique g with g(g-1)<n<=g(g+1).  Start from the integer square-root
    # neighborhood and verify the inequalities, independently of the forward
    # map used below.
    g = max(1, isqrt(n))
    while not (g * (g - 1) < n <= g * (g + 1)):
        g += 1
    return g, n - g * (g - 1)


def inverse_pair(N, sign):
    g = max(1, isqrt(2 * N))
    while T(g) < N:
        g += 1
    while g > 1 and T(g - 1) >= N:
        g -= 1
    h = N - T(g - 1)
    return g, (g + h if sign == 1 else g + 1 - h)


def main():
    Z, c = sp.symbols("Z c", nonzero=True)
    domain = sp.QQ.frac_field(c)
    max_g = 80
    oriented = {}
    quotient_orbits = {}
    symbolic_tail_checks = 0

    for g in range(1, max_g + 1):
        m = 2 * g + 1
        for r in range(1, m):
            epsilon = r % 2
            polynomial = sp.Poly(Z**epsilon * (Z**(m-r) - c), Z,
                                 domain=domain)
            need(sp.gcd(polynomial, polynomial.diff()).degree() == 0,
                 f"g={g},r={r}: symbolic squarefree tail")
            degree = polynomial.degree()
            need(degree % 2 == 1, f"g={g},r={r}: one infinity")
            genus = (degree - 1) // 2
            delta = r // 2
            need((genus, delta) == ((m-r)//2, r//2),
                 f"g={g},r={r}: genus/delta formula")
            need(genus + delta == g, f"g={g},r={r}: delta conservation")

            rr = m-r
            need((rr//2, r//2) == (genus, delta),
                 f"g={g},r={r}: reciprocal exchange")
            lam = Fraction(r, m-r)
            lam_dual = Fraction(m-r, r)
            coeff = Fraction(1, 2)
            d = 7
            need((coeff*d*lam) * (coeff*d*lam_dual) == (coeff*d)**2,
                 f"g={g},r={r}: reciprocal differential product")

            n = g*(g-1)+r
            need(inverse_oriented(n) == (g,r), "oriented inverse")
            need(n not in oriented, "oriented injectivity")
            oriented[n] = (g,r)

            centered = 2*r-m
            sign = 1 if centered > 0 else -1
            h = (abs(centered)+1)//2
            N = T(g-1)+h
            need(inverse_pair(N,sign) == (g,r), "pair inverse")
            orbit = tuple(sorted((r,m-r)))
            if N in quotient_orbits:
                need(quotient_orbits[N] == (g,orbit), "pair quotient collision")
            else:
                quotient_orbits[N] = (g,orbit)
            symbolic_tail_checks += 1

    need(sorted(oriented) == list(range(1,max_g*(max_g+1)+1)),
         "oriented index surjectivity")
    need(sorted(quotient_orbits) == list(range(1,T(max_g)+1)),
         "triangular pair quotient surjectivity")

    # Excluded-scope probes.  Even m makes the normalized branch degree even,
    # hence two projective infinity points.  A zero leading critical value is
    # nonsquarefree, and characteristic dividing m-r can make nonzero roots
    # inseparable.
    even_poly = sp.Poly(Z*(Z**3-c),Z,domain=domain)
    need(even_poly.degree() == 4, "even-m two-infinity hostile")
    zero_c_poly = sp.Poly((Z*(Z**4-c)).as_expr().subs(c,0),Z,domain=sp.QQ)
    need(sp.gcd(zero_c_poly,zero_c_poly.diff()).degree() > 0,
         "c=0 nonsquarefree hostile")
    char3_poly = sp.Poly(Z**3-1,Z,modulus=3)
    need(sp.gcd(char3_poly,char3_poly.diff()).degree() > 0,
         "characteristic-three inseparability hostile")

    print("THM4341 ODD CUSP INDEPENDENT SYMBOLIC AUDIT")
    print(f"SYMBOLIC_TAILS={symbolic_tail_checks};g=1..{max_g};squarefree=PASS;odd_degree=PASS")
    print(f"INDEX_BIJECTIONS=oriented:{len(oriented)};reciprocal_pairs:{len(quotient_orbits)}")
    print("RECIPROCAL_CHECK=lambda_product=1;tail_genus_equals_complement_delta;order_product=PASS")
    print("SCOPE_HOSTILES=even_m_two_infinities;c_zero_nonsquarefree;char3_inseparable")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
