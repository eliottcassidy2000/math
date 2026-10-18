#!/usr/bin/env python3
"""Audit F2, check (1): discriminant of the tower order O = Z[zeta15][j][w] (power basis zeta^a j^b w^c).

disc(O) = det(Tr_{F/Q}(e_i e_k)) over the 48 basis monomials, with Tr computed as the trace of the multiplication
matrix (no tensor-product assumption).  Expected from the tensor-product formula:
  disc(O) = disc(Phi_15)^6 * disc(x^2+11)^24 * disc(x^3+2x^2-4x-4)^16 = (3^4 5^6)^6 (-44)^24 592^16.
If 19 does not divide disc(O), then 19 does not divide the index [O_F : O] (index^2 | disc(O)) and 19 is unramified
in F, so O is 19-maximal and O (x) Z_(19) = O_F (x) Z_(19): an element is integral at every prime above 19 iff its
coordinates have no 19 in their denominators.
"""
import sys, os
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f2_field import E, N, EXP
import flint
from sympy import factorint, discriminant, symbols, cyclotomic_poly

x = symbols('x')
d15 = discriminant(cyclotomic_poly(15, x), x)
d2 = discriminant(x ** 2 + 11, x)
d3 = discriminant(x ** 3 + 2 * x ** 2 - 4 * x - 4, x)
print('disc(Phi_15) = %s = %s ; disc(x^2+11) = %s ; disc(cubic) = %s = %s' % (d15, factorint(d15), d2, d3, factorint(d3)))

basis = []
for k in range(N):
    v = [0] * N; v[k] = 1
    basis.append(E(v))
# trace of each basis element = trace of its multiplication matrix
tr = []
for k in range(N):
    t = Fraction(0)
    for l in range(N):
        prod = basis[k] * basis[l]
        t += Fraction(prod.n[l], prod.d)
    tr.append(t)
assert all(t.denominator == 1 for t in tr)
tr = [int(t) for t in tr]
print('trace of 1 =', tr[0])
G = []
for i in range(N):
    row = []
    for k in range(N):
        pr = basis[i] * basis[k]
        s = Fraction(0)
        for l in range(N):
            if pr.n[l]: s += Fraction(pr.n[l], pr.d) * tr[l]
        assert s.denominator == 1
        row.append(int(s))
    G.append(row)
D = int(flint.fmpz_mat(G).det())
expected = (3 ** 4 * 5 ** 6) ** 6 * (-44) ** 24 * 592 ** 16
print('det(trace form) = disc(O): sign %+d, |disc| has %d digits' % ((D > 0) - (D < 0), len(str(abs(D)))))
print('equals the tensor-product prediction (3^4 5^6)^6 (-44)^24 592^16:', D == expected)
rest = abs(D); fac = {}
for p in (2, 3, 5, 11, 37):
    while rest % p == 0:
        rest //= p; fac[p] = fac.get(p, 0) + 1
print('factorisation over {2,3,5,11,37}:', fac, '; cofactor:', rest)
print('19 divides disc(O):', D % 19 == 0)
