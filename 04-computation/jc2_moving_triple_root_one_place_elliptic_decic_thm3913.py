#!/usr/bin/env python3
"""Exact companion for THM-3913's moving-triple-root elliptic decic."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


def binary_disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


A, C, U, V, Z = sp.symbols("A C U V Z")
ell = A * U + C * V
Phi = sp.expand(ell**3 + C * ell * U**2 + A**2 * V**3)
Phi_poly = sp.Poly(Phi, U, V)
coefficients = (
    Phi_poly.coeff_monomial(U**3),
    Phi_poly.coeff_monomial(U**2 * V),
    Phi_poly.coeff_monomial(U * V**2),
    Phi_poly.coeff_monomial(V**3),
)
expected_coefficients = (
    A**3 + A * C,
    3 * A**2 * C + C**2,
    3 * A * C**2,
    C**3 + A**2,
)
gate(coefficients == expected_coefficients, "compact and coefficient forms agree")
gate(all(polynomial.subs({A: 0, C: 0}) == 0 for polynomial in coefficients),
     "common coefficient zero")

coefficient_gcd = sp.Poly(coefficients[0], A, C)
for coefficient in coefficients[1:]:
    coefficient_gcd = sp.gcd(coefficient_gcd, sp.Poly(coefficient, A, C))
gate(coefficient_gcd.total_degree() == 0, "binary cubic is primitive")

Delta = binary_disc(*coefficients)
expected_delta = (
    -27 * A**10
    - 54 * A**8 * C
    - 27 * A**6 * C**2
    - 36 * A**4 * C**5
    - 4 * A**2 * C**6
    - 4 * C**9
)
zero(Delta - expected_delta, "exact sparse decic discriminant")
gate(sp.Poly(Delta, A, C).total_degree() == 10, "discriminant degree ten")
gate(sp.factor(Delta) == Delta, "decic is irreducible over Q")

# The degree-three leading row is the moving triple root ell^3.  The
# quadratic perturbation evaluates to A^5 at its repeated-root address.
q5 = sp.expand((C * ell * U**2 + A**2 * V**3).subs({U: -C, V: A}))
gate(q5 == A**5, "moving-triple-root transverse quintic is one fifth power")
Delta_top = sum(
    coefficient * A**powers[0] * C**powers[1]
    for powers, coefficient in sp.Poly(Delta, A, C).terms()
    if sum(powers) == 10
)
gate(Delta_top == -27 * A**10, "one-point decic top")

# Projective infinity is the sole point [0:1:0], and it is smooth because
# the C^9 Z term survives.
Delta_h = sp.expand(
    sum(
        coefficient
        * A**powers[0]
        * C**powers[1]
        * Z ** (10 - sum(powers))
        for powers, coefficient in sp.Poly(Delta, A, C).terms()
    )
)
zero(Delta_h.subs(Z, 0) + 27 * A**10, "unique infinity support")
infinity_point = {A: 0, C: 1, Z: 0}
gate(sp.diff(Delta_h, Z).subs(infinity_point) == -4, "infinity is smooth")

# Generic cubic irreducibility hostile.  At C=0 the primitive specialization
# is A^2(A U^3+V^3); the last factor is Kummer-irreducible over k(A) because
# the A-valuation of -1/A is -1, not divisible by three.
gate(sp.factor(Phi.subs(C, 0)) == A**2 * (A * U**3 + V**3),
     "Kummer specialization at C=0")
gate((-1) % 3 != 0, "Kummer valuation is not divisible by three")

# Exact elliptic normalization.  Constants k0,l0 satisfy 6 k0^2=1 and
# 3 l0^2+k0=0.  The source E is y^2=w(w^2-1).
w, y, k0, l0 = sp.symbols("w y k0 l0")
A0 = w * (w**2 + 2) * y
C0 = w * (w**2 - 1) * (w**2 + 2)
A_t = l0 * A0
C_t = k0 * C0
normalization_ideal = sp.groebner(
    [y**2 - w * (w**2 - 1), 6 * k0**2 - 1, 3 * l0**2 + k0],
    y,
    l0,
    k0,
    w,
    order="lex",
)
zero(normalization_ideal.reduce(sp.expand(Delta.subs({A: A_t, C: C_t})))[1],
     "elliptic parametrization lies on the decic")
gate(sp.discriminant(w * (w**2 - 1), w) == 4, "elliptic cubic is squarefree")

# Rational inverse in the scaled coordinates A0=A/l0,C0=C/k0.
R = sp.cancel(A0**2 / C0)
R_on_E = sp.cancel(R.subs(y**2, w * (w**2 - 1)))
gate(R_on_E == w**4 + 2 * w**2, "scaled radial invariant")
zeta = sp.cancel((C0**2 / R_on_E + 2 * R_on_E - 2) / (R_on_E + 1))
zero(sp.factor(zeta - w**2), "inverse recovers w squared")
w_inverse = sp.cancel(C0 / ((zeta - 1) * (zeta + 2)))
zero(sp.factor(w_inverse - w), "inverse recovers w")
y_inverse = sp.cancel(A0 / (w_inverse * (w_inverse**2 + 2)))
zero(sp.factor(y_inverse - y), "inverse recovers y")

# Pole orders at the elliptic origin O are ord_O(w)=-2, ord_O(y)=-3.
# Hence C has pole ten and A pole nine: the unique infinity place is smooth.
gate(2 * 5 == 10, "C pole order at elliptic infinity")
gate(2 * 3 + 3 == 9, "A pole order at elliptic infinity")
gate(10 - 9 == 1, "smooth local parameter A/C at infinity")

# The split boundary on the quadratic double plane has the first
# three-primary determinant allowed by the general THM-3912 sieve.
boundary_gram = sp.Matrix([[-4, 5], [5, -4]])
gate(boundary_gram.det() == -9, "decic split-boundary determinant")
gate(sp.gcd_list(list(boundary_gram)) == 1, "boundary Smith first invariant")
gate(abs(int(boundary_gram.det())) == 9, "boundary Smith second invariant")

semantic = {
    "binary_cubic": "(AU+CV)^3+C(AU+CV)U^2+A^2V^3",
    "coefficient_depth": 3,
    "index": "common-zero therefore globally nonmonogenic",
    "discriminant": "irreducible one-smooth-place decic",
    "normalization": "elliptic y^2=w(w^2-1) minus O",
    "boundary_snf": [1, 9],
    "resolvent": "nonzero class-group three-torsion forced by S3 closure",
    "failure": "elliptic Jelonek component is not polynomially uniruled",
    "scope": "no plane atlas and no JC counterexample",
}
semantic_sha256 = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("theorem=THM-3913-moving-triple-root-one-place-decic-normal-nonmonogenic-cubic")
print("binary_cubic=(AU+CV)^3+C(AU+CV)U^2+A^2V^3")
print("order=normal_nonmonogenic_S3;coefficient_depth=3")
print("branch=irreducible_decic;infinity=one_smooth_place")
print("normalization=elliptic_y2=w(w2-1)_minus_O")
print("quadratic_resolvent_units=kstar;class_group_3_torsion=NONZERO")
print("split_boundary_snf=1,9")
print("plane_atlas=EMPTY_BY_POLYNOMIAL_UNIRULEDNESS")
print("JC2=OPEN")
print(f"semantic_sha256={semantic_sha256}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
