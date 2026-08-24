#!/usr/bin/env python3
"""Exact gauge, centralizer, and elliptic-filling audit for the S6 source.

This checks only the displayed integer monodromy and affine-action lattices.
It does not construct the analytic quotient, extend an isomorphism through the
cusp, compute the global Cech class, or verify a complex structure on S6.
"""

from itertools import combinations
from math import gcd

from sympy import Matrix, eye


def require(label, condition):
    if not condition:
        raise RuntimeError(f"FAIL  {label}")
    print(f"PASS  {label}")


T1 = Matrix(
    [[1, 0, -6, 2], [0, -1, 1, 1], [0, -1, 0, 1], [0, 0, 0, 1]]
)
T2 = Matrix(
    [[1, 6, 0, -3], [0, 0, -1, 1], [0, 1, 0, 0], [0, 0, 0, 1]]
)
A1 = T1.inv().T
A2 = T2.inv().T
I4 = eye(4)

# Solve the simultaneous rational centralizer from a 32 by 16 exact response.
basis_matrices = []
for i in range(4):
    for j in range(4):
        elementary = Matrix.zeros(4)
        elementary[i, j] = 1
        basis_matrices.append(elementary)
columns = []
for elementary in basis_matrices:
    response = (elementary * T1 - T1 * elementary).col_join(
        elementary * T2 - T2 * elementary
    )
    columns.append(Matrix(list(response)))
centralizer_response = Matrix.hstack(*columns)
centralizer_kernel = centralizer_response.nullspace()
E14 = Matrix.zeros(4)
E14[0, 3] = 1
require("simultaneous rational centralizer has dimension two", len(centralizer_kernel) == 2)
require("identity centralizes T1,T2", I4 * T1 == T1 * I4 and I4 * T2 == T2 * I4)
require("E14 centralizes T1,T2", E14 * T1 == T1 * E14 and E14 * T2 == T2 * E14)
expected_span = Matrix.hstack(Matrix(list(I4)), Matrix(list(E14)))
actual_span = Matrix.hstack(*centralizer_kernel)
require("centralizer equals span(I,E14)", expected_span.row_join(actual_span).rank() == 2)
print("RESULT C_GL4(Z)(T1,T2)={+/- (I+bE14): b in Z}")

# The dual action is Q_b=I+bE41 and centralizes A1,A2.
E41 = E14.T
for b in range(-5, 6):
    Q = I4 + b * E41
    if Q * A1 != A1 * Q or Q * A2 != A2 * Q:
        raise RuntimeError(("dual centralizer", b))
print("PASS  dual Q_b centralizes A1,A2 for -5<=b<=5")

# Exact quotient of the raw Seifert twist lattice by standard gauge.
k1 = Matrix([1, 3, 0])
k2 = Matrix([1, 0, 4])
q = Matrix([0, -1, 1])
gauge_basis = Matrix.hstack(k1, k2, q)
phi = Matrix([[12, -4, -3]])
require("gauge-complement basis is unimodular", abs(int(gauge_basis.det())) == 1)
require("both standard gauges lie in ker(phi)", phi * k1 == Matrix([0]) and phi * k2 == Matrix([0]))
require("complement maps to signed p=1", phi * q == Matrix([1]))
print("RESULT Z^3/<gauge moves> is Z via signed p")

residue_table = {}
for ell1 in (1, 2):
    for ell2 in (1, 3):
        p_residue = int((-4 * ell1 - 3 * ell2) % 12)
        residue_table[(ell1, ell2)] = p_residue
require(
    "admissible residue table is bijective onto units mod 12",
    residue_table == {(1, 1): 5, (1, 3): 11, (2, 1): 1, (2, 3): 7},
)
print(f"RESULT admissible p residue decoder={residue_table}")

# Saturated invariant lattices and their primitive invariant covectors.
delta = Matrix([0, 0, 0, 1])
epsilon1 = Matrix([1, 2, -4, 0])
epsilon2 = Matrix([1, 3, -3, 0])
gamma = Matrix([[1, 0, 0, 0]])
psi1 = Matrix([[0, 2, 1, 3]])
psi2 = Matrix([[0, 1, 1, 2]])
require("epsilon1,delta are A1-invariant", A1 * epsilon1 == epsilon1 and A1 * delta == delta)
require("epsilon2,delta are A2-invariant", A2 * epsilon2 == epsilon2 and A2 * delta == delta)
for vectors, action, label in (
    ((epsilon1, delta), A1, "A1"),
    ((epsilon2, delta), A2, "A2"),
):
    lattice = Matrix.hstack(*vectors)
    minors = [
        abs(int(lattice.extract(rows, (0, 1)).det()))
        for rows in combinations(range(4), 2)
    ]
    minor_gcd = 0
    for value in minors:
        minor_gcd = gcd(minor_gcd, value)
    require(f"displayed {label}-invariant lattice is saturated", minor_gcd == 1 and (action - I4).rank() == 2)
require("gamma is invariant for both elliptic actions", gamma * A1 == gamma and gamma * A2 == gamma)
require("psi1 is A1-invariant with psi1(delta)=3", psi1 * A1 == psi1 and int((psi1 * delta)[0]) == 3)
require("psi2 is A2-invariant with psi2(delta)=2", psi2 * A2 == psi2 and int((psi2 * delta)[0]) == 2)

# Q_b shifts the invariant vectors by b*ell*delta.  The invariant-covector
# criterion for translation conjugacy of x -> A_j x+v_j/m_j is automatic at
# order 3 and forces b even at order 4 when ell2 is odd.
filling_checks = 0
for b in range(-12, 13):
    Q = I4 + b * E41
    if Q * epsilon1 - epsilon1 != b * delta:
        raise RuntimeError(("epsilon1 shift", b))
    if Q * epsilon2 - epsilon2 != b * delta:
        raise RuntimeError(("epsilon2 shift", b))
    for ell1 in range(-8, 9):
        if gcd(ell1, 3) != 1:
            continue
        order3_numerator = int((psi1 * (b * ell1 * delta))[0])
        if order3_numerator % 3:
            raise RuntimeError((b, ell1, order3_numerator))
        negative_gamma_numerator = -2 * ell1
        if negative_gamma_numerator % 3 == 0:
            raise RuntimeError(("negative sign unexpectedly extends at order 3", b, ell1))
        filling_checks += 1
    for ell2 in range(-9, 10, 2):
        order4_numerator = int((psi2 * (b * ell2 * delta))[0])
        if (order4_numerator % 4 == 0) != (b % 2 == 0):
            raise RuntimeError((b, ell2, order4_numerator))
        negative_gamma_numerator = -2 * ell2
        if negative_gamma_numerator % 4 == 0:
            raise RuntimeError(("negative sign unexpectedly extends at order 4", b, ell2))
        filling_checks += 1
print("PASS  Q_b shifts both primitive elliptic invariants by b delta")
print(f"PASS  elliptic filling congruence controls ({filling_checks} cases)")

# Marked local hostile: adding one delta to v2 changes psi2 residue by 2 mod 4
# while leaving ell2 and the discrete p coordinate unchanged.
v2 = -epsilon2
v2_hostile = v2 + delta
require("marked order-four hostile keeps ell2", int((gamma * v2)[0]) == int((gamma * v2_hostile)[0]) == -1)
require("marked order-four hostile changes invariant residue by 2", int((psi2 * (v2_hostile - v2))[0]) % 4 == 2)

print("RESULT punctured unipotent orbit shifts c0 by Z")
print("RESULT order-four filling permits only even shifts")
print("RESULT first filled marked invariant is [c0-c2] in C/(2Z)")
print("OPEN global overlap/cusp Cech class for the b=2 shift")
print("NONCONSEQUENCE no analytic quotient, homology, diffeomorphism, or S6 claim is verified")
print("ALL EXACT CHECKS PASSED")
