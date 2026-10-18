#!/usr/bin/env python3
"""Two-step reciprocals of the six-colour field (THM-4700(2)): 1/n = v + conj(v) with v a unit iff -(4n^2 - 1) is a
square in F, iff 4n^2 - 1 = D y^2 with sqrt(-D) in F, D in {3, 11, 15, 55} (the imaginary quadratic subfields), i.e.
(2n)^2 - D y^2 = 1. Equivalent square-class test used here: (4n^2 - 1)/3 has square class in {1, 5, 33, 165}
(the real quadratic subfields of F+ are Q(sqrt5), Q(sqrt33), Q(sqrt165)). Also checks the Pell parity of D = 55."""
from fractions import Fraction as Fr
from sympy import factorint
def sqf(fr):
    x = fr.numerator * fr.denominator; s = 1
    for p, e in factorint(x).items():
        if e % 2: s *= p
    return s
ok = [(n, sqf(Fr(4 * n * n - 1, 3))) for n in range(2, 12001) if sqf(Fr(4 * n * n - 1, 3)) in (1, 5, 33, 165)]
print('n <= 12000 with L(1/n) = 2 (n, square class of (4n^2-1)/3):', ok)
D_of = {1: 3, 5: 15, 33: 11, 165: 55}
print('as Pell (2n)^2 - D y^2 = 1:', [(n, D_of[k], int(((4 * n * n - 1) // D_of[k]) ** 0.5)) for n, k in ok])
x, y = 89, 12; par = []
for m in range(6):
    par.append(x % 2); x, y = 89 * x + 55 * 12 * y, 12 * x + 89 * y
print('D = 55: x mod 2 along (89 + 12 sqrt55)^m, m = 1..6:', par, '(always odd: no n)')
