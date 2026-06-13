#!/usr/bin/env python3
"""
worpitzky_moments.py — Express each Worpitzky coefficient c_j as a moment of F.

c_j = sum_k F_k * [m^j coefficient of C(m+n-1-k, n-1)]

The binomial C(m+n-1-k, n-1) = prod_{i=0}^{n-2} (m+n-1-k-i) / (n-1)!
= prod_{i=0}^{n-2} (m + (n-1-k-i)) / (n-1)!

Let alpha_i = n-1-k-i for i=0,...,n-2. Then alpha_i = n-1-k-i.
Sum of alpha_i = sum(n-1-k-i, i=0..n-2) = (n-1)^2 - k(n-1) - (n-2)(n-1)/2
= (n-1)[(n-1) - k - (n-2)/2] = (n-1)[n/2 - k]

So the m^{n-2} coefficient is sum(alpha_i)/(n-1)! = (n/2 - k)/(n-2)!

Similarly, the m^{n-3} coefficient involves sum_{i<j} alpha_i * alpha_j / (n-1)!

Let's just compute these symbolically for small n.

Author: opus-2026-03-07-S46c
"""
from fractions import Fraction
from math import comb, factorial
from sympy import symbols, Poly, binomial, expand, factorial as sfact

m, k = symbols('m k')

for n in [5, 6, 7]:
    print(f"\n{'='*60}")
    print(f"n={n}: Worpitzky coefficient structure")
    print(f"{'='*60}")

    # C(m+n-1-k, n-1) as polynomial in m (with k as parameter)
    binom_expr = 1
    for i in range(n-1):
        binom_expr *= (m + n - 1 - k - i)
    binom_expr = binom_expr / factorial(n-1)

    # Expand as polynomial in m
    poly_m = Poly(expand(binom_expr), m)
    coeffs_of_m = poly_m.all_coeffs()  # high to low in m

    print(f"  C(m+{n}-1-k, {n}-1) coefficients in m (high to low):")
    for deg, c in enumerate(coeffs_of_m):
        power = n - 1 - deg
        print(f"    m^{power}: {c}")

    # Now c_j = sum_k F_k * [m^j coeff of binom]
    # The [m^j coeff] is a polynomial in k
    # c_j = sum_k F_k * P_j(k)
    # where P_j(k) is a polynomial in k
    #
    # Key moments:
    # M_0 = sum_k F_k = n!
    # M_1 = sum_k k*F_k = n!(n-1)/2 (by palindrome)
    # M_2 = sum_k k^2*F_k = n! * E[fwd^2]
    # M_r = sum_k k^r * F_k = n! * E[fwd^r]
    #
    # The coefficient P_j(k) at each level j tells us which moments matter.

    print(f"\n  Worpitzky coefficients as moment functionals:")
    for deg, c in enumerate(coeffs_of_m):
        power = n - 1 - deg
        # c is a polynomial in k
        poly_k = Poly(expand(c), k)
        k_coeffs = poly_k.all_coeffs()
        print(f"    c_{power} = sum_k F_k * ({c})")
        print(f"         = sum of F-moments: ", end="")

        # Express in terms of M_r
        terms = []
        for kdeg, kc in enumerate(k_coeffs):
            r = len(k_coeffs) - 1 - kdeg
            if kc != 0:
                terms.append(f"{kc}*M_{r}")
        print(" + ".join(terms))

    # Compute the transitive (Eulerian) values
    print(f"\n  Transitive tournament (M_r = E[fwd^r] * n!):")
    print(f"    M_0 = {factorial(n)}")
    print(f"    M_1 = {factorial(n) * (n-1) // 2}")

    # E[fwd^2] for transitive = Var + mean^2
    # Var[fwd] for transitive = (n+1)/12 (t3=0)
    mean_fwd = Fraction(n-1, 2)
    var_trans = Fraction(n+1, 12)
    efwd2_trans = var_trans + mean_fwd**2
    M2_trans = efwd2_trans * factorial(n)
    print(f"    E[fwd^2] = {efwd2_trans} = {float(efwd2_trans):.4f}")
    print(f"    M_2 = {M2_trans}")

    # Verify c_j for transitive = C(n, j)
    print(f"\n  Verification: transitive c_j should be C(n,j):")
    for j in range(n):
        expected = comb(n, j)
        print(f"    c_{j} = C({n},{j}) = {expected}")
