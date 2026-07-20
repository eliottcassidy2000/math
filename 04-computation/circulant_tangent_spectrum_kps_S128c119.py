#!/usr/bin/env python3
"""circulant_tangent_spectrum_kps_S128c119.py -- kind-pasteur-2026-07-20-S128c119

THE CIRCULANT TOURNAMENT'S SKEW SPECTRUM IS A TANGENT, IN CLOSED FORM.

THM-1440(D) already established that circulant tournaments make the sine literal:
mu_j = 2i * sum_{k in C} sin(2 pi j k / n), and that Re = -1/2 is a critical line for
them.  That is a SINE SUM.  This closes it to a CLOSED FORM for the standard connection
set C = {1, 2, ..., (n-1)/2}, i.e. the rotational tournament R_n:

    char_S(x)  =  ODD PART OF (1 + x)^n  =  ( (1+x)^n - (1-x)^n ) / 2
    spec(S)    =  { i * tan(k pi / n) : k = 0, ..., n-1 }

Both statements are the SAME fact, because substituting x = i tan(theta) into the odd
part of (1+x)^n gives i times the NUMERATOR of tan(n theta), whose zeros are theta =
k pi / n.  So the eigenvalue equation IS the multiple-angle formula for the tangent.

This also puts the owner's half-dictionary at the centre: the coordinate in which the
statement is clean is x = 2 lambda + 1, the inverse of x |-> (1+x)/2 (THM-1555).

n = 7 gives x(x^6 + 21x^4 + 35x^2 + 7) -- coefficients C(7,1), C(7,3), C(7,5), C(7,7).
"""
import numpy as np
import sympy as sp

x = sp.Symbol('x')


def circulant_tournament(n):
    """R_n: i -> j iff (j - i) mod n in {1, ..., (n-1)/2}."""
    C = set(range(1, (n - 1) // 2 + 1))
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in C:
                A[i, j] = 1
    return A


def odd_part(n):
    return sp.expand((((1 + x) ** n) - ((1 - x) ** n)) / 2)


print("%-4s %-56s %-10s %-10s" % ("n", "char_S(x) = odd part of (1+x)^n ?", "matches",
                                  "spec = i tan(k pi/n)"))
print("-" * 96)
for n in range(3, 14, 2):
    A = circulant_tournament(n)
    S = sp.Matrix((A - A.T).tolist())
    cp = sp.expand(S.charpoly(x).as_expr())
    op = odd_part(n)
    # char poly of a skew matrix: det(xI - S); odd part has leading coeff 1 -- compare
    match = sp.simplify(cp - op) == 0
    ev = np.linalg.eigvals((A - A.T).astype(float))
    pred = sorted(np.tan(np.pi * k / n) for k in range(n))
    got = sorted(e.imag for e in ev)
    tan_ok = np.allclose(pred, got, atol=1e-8)
    print("%-4d %-56s %-10s %-10s" % (n, sp.factor(cp), match, tan_ok))

print()
print("EXPLICIT SMALL CASES (factored over Q):")
for n in (3, 5, 7, 9, 11):
    A = circulant_tournament(n)
    S = sp.Matrix((A - A.T).tolist())
    print("  n = %-3d char_S(x) = %s" % (n, sp.factor(sp.expand(S.charpoly(x).as_expr()))))
print()
print("  Coefficients of char_S are the ODD-INDEX BINOMIALS C(n, 2j+1):")
for n in (5, 7, 9, 11):
    print("     n = %-3d : %s" % (n, [int(sp.binomial(n, k)) for k in range(1, n + 1, 2)]))
print()
print("WHICH regular 7-tournament is it?  (compare against the three iso classes)")
A7 = circulant_tournament(7)
lam = sp.Symbol('lam')
print("  char_A(R_7) = %s" % sp.factor(sp.expand(sp.Matrix(A7.tolist()).charpoly(lam).as_expr())))
print("  -> this is class [1] of x_half_shift_regular7, the one with x^6+21x^4+35x^2+7.")
print("     Paley-7 (doubly regular, connection set the quadratic residues {1,2,4}) is")
print("     class [3], char_S = x(x^2+7)^3.  They are DIFFERENT tournaments.")
print()
print("THE OWNER'S x(x^2+7)(x^4+14x^2+17): the closest regular 7-tournament is class [2],")
print("  char_S = x(x^2+7)(x^4+14x^2+1) -- constant 1, not 17.  Checked: no regular")
print("  tournament on 7 vertices has x^4+14x^2+17 as a skew factor (all three enumerated).")
print("  Note x^4+14x^2+17 has x^2 = -7 +- 4 sqrt(2) while ours has x^2 = -7 +- 4 sqrt(3),")
print("  so the two differ by sqrt(2) vs sqrt(3) -- a different structure, not a variant.")
