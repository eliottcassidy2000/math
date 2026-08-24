#!/usr/bin/env python3
"""Exact companion for THM-3977's simultaneous cusp-arm obstruction."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


x, t, a, L, u = sp.symbols("x t a L u", nonzero=True)
L_relation = 6 * L**2 - 1


def zero_mod_L(expression: sp.Expr, message: str) -> None:
    """Check a rational identity modulo 6L^2-1."""
    numerator = sp.together(expression).as_numer_denom()[0]
    remainder = sp.Poly(sp.expand(numerator), L, domain=sp.EX).rem(
        sp.Poly(L_relation, L, domain=sp.EX)
    ).as_expr()
    zero(remainder, message)


z = 1 + x**2 * t
p = sp.expand(z * t)
y = sp.expand(x * z * t**2)

A = sp.expand(y**2 + 2 * L * x + a * p)
C = sp.expand(
    y**3 + 3 * L * x * y + sp.Rational(3, 2) * a * p * y * (1 - z)
    - x * z / a
    + x**2 * (-sp.Rational(3, 2) * a * L
              + sp.Rational(3, 8) * a**2 * y)
)


def jac(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(f, x) * sp.diff(g, t)
                     - sp.diff(f, t) * sp.diff(g, x))


# ---------------------------------------------------------------------------
# Exact completion rows and simultaneous boundary jets.
# ---------------------------------------------------------------------------

expected_A = (
    2 * L * x + a * t + a * x**2 * t**2
    + x**2 * t**4 + 2 * x**4 * t**5 + x**6 * t**6
)
zero(A - expected_A, "source first-coordinate expansion")

X, Z, P, Y = sp.symbols("X Z P Y")
A_B = Y**2 + 2 * L * X + a * P
C_B = (
    Y**3 + 3 * L * X * Y
    + sp.Rational(3, 2) * a * P * Y * (1 - Z)
    - X * Z / a
    + X**2 * (-sp.Rational(3, 2) * a * L
              + sp.Rational(3, 8) * a**2 * Y)
)
zero(A_B.subs({X: 0, Z: 0, P: 0}) - Y**2,
     "boundary D first cusp row")
zero(C_B.subs({X: 0, Z: 0, P: 0}) - Y**3,
     "boundary D second cusp row")
zero(A_B.subs({X: 0, Z: 1, Y: 0}) - a * P,
     "interior arm first row")
zero(C_B.subs({X: 0, Z: 1, Y: 0}),
     "interior arm second row")

J = sp.together(jac(A, C))
zero_mod_L(J.subs(x, 0) - 1, "source x=0 Jacobian residue one")
J_minus_one_num = sp.together(a * (J - 1)).as_numer_denom()[0]
J_minus_one_red = sp.Poly(sp.expand(J_minus_one_num), L, domain=sp.EX).rem(
    sp.Poly(L_relation, L, domain=sp.EX)
).as_expr()
gate(sp.Poly(J_minus_one_red, x).eval(0) == 0,
     "Jacobian minus one divisible by source x")

# The genuine D-chart computation, not an illicit t-finite substitution.
A_jet = Y**2 + X * (2 * L - a * Y)
C_jet = (
    Y**3 + X * (3 * L * Y - sp.Rational(3, 2) * a * Y**2)
    + X**2 * (-sp.Rational(3, 2) * a * L
              + sp.Rational(3, 8) * a**2 * Y)
)
jet_det = sp.expand(
    sp.diff(A_jet, X) * sp.diff(C_jet, Y)
    - sp.diff(A_jet, Y) * sp.diff(C_jet, X)
)
jet_linear = sp.Poly(jet_det, X).coeff_monomial(X)
zero_mod_L(jet_linear - 1, "boundary D normalized Jacobian residue one")


# ---------------------------------------------------------------------------
# Exact derivative resultant and its reduced factor.
# ---------------------------------------------------------------------------

Ax = sp.diff(A, x)
At = sp.diff(A, t)
raw_residual = (
    15552 * L**5 * x**12
    + 3456 * L**4 * a * x**9
    - 1024 * L**4 * x**5
    + 96 * L**3 * a**2 * x**6
    + 32 * L**2 * a**3 * x**3
    + 4 * L * a**5 * x**4
    - 4 * L * a**4
    + a**6 * x
)
resultant = sp.factor(sp.resultant(Ax, At, t))
zero(resultant - 96 * x**24 * raw_residual,
     "unreduced Sylvester resultant factorization")

H = (
    3888 * L * x**12
    + 864 * a * x**9
    + 144 * L * a**2 * x**6
    - 256 * x**5
    + 36 * L * a**5 * x**4
    + 48 * a**3 * x**3
    + 9 * a**6 * x
    - 36 * L * a**4
)
zero_mod_L(resultant - sp.Rational(32, 3) * x**24 * H,
           "reduced derivative resultant factorization")
gate(sp.Poly(H, x).degree() == 12, "residual degree twelve")
gate(sp.Poly(H, x).LC() == 3888 * L, "residual leading coefficient")
gate(H.subs(x, 0) == -36 * L * a**4, "residual constant coefficient")
gate(sp.Poly(Ax, t).LC() == 6 * x**5,
     "partial-x stable t-leading coefficient")
gate(sp.Poly(At, t).LC() == 6 * x**6,
     "partial-t stable t-leading coefficient")


# ---------------------------------------------------------------------------
# Intrinsic u=x^2t critical-address cover and reconstruction.
# ---------------------------------------------------------------------------

f = u**4 * (u + 1)**2
g = u * (u + 1)
zero(sp.together(A.subs(t, u / x**2)
                 - (x**-6 * f + a * x**-2 * g + 2 * L * x)),
     "u-normalized first coordinate")

E1 = 2 * u**3 * (u + 1) * (3 * u + 2) + a * x**4 * (2 * u + 1)
E2 = -6 * f - 2 * a * x**4 * g + 2 * L * x**7
zero(E1 - (sp.diff(f, u) + a * x**4 * sp.diff(g, u)),
     "normalized u derivative")

X4 = -2 * u**3 * (u + 1) * (3 * u + 2) / (a * (2 * u + 1))
X7 = -u**4 * (u + 1)**2 / (L * (2 * u + 1))
zero(sp.together(E1.subs(x**4, X4)), "first critical power row")
zero(sp.together(E2.subs({x**4: X4, x**7: X7})),
     "second critical power row")

Phi = 9 * a**7 * (u + 1) * (2 * u + 1)**3 + 32 * u**5 * (3 * u + 2)**7
compatibility = (
    X4**7 - X7**4
    + 4 * u**16 * (u + 1)**7 * Phi / (a**7 * (2 * u + 1)**7)
)
zero_mod_L(compatibility, "power compatibility equals critical cover")
gate(sp.Poly(Phi, u).degree() == 12, "critical cover degree twelve")
gate(sp.Poly(Phi, u).LC() == 32 * 3**7,
     "critical cover stable leading coefficient")
gate(Phi.subs(u, 0) == 9 * a**7, "critical cover nonzero constant")
gate(Phi.subs(u, -1) == 32, "critical cover avoids u=-1")
gate(Phi.subs(u, -sp.Rational(1, 2)) == -sp.Rational(1, 128),
     "critical cover avoids u=-1/2")
gate(Phi.subs(u, -sp.Rational(2, 3)) == -a**7 / 9,
     "critical cover avoids u=-2/3")

x_reconstructed = sp.factor(X4**2 / X7)
zero(x_reconstructed
     + 4 * L * u**2 * (3 * u + 2)**2 / (a**2 * (2 * u + 1)),
     "explicit critical x reconstruction")
zero(sp.together(
    x_reconstructed**4 - X4
    - X4 * (X4**7 - X7**4) / X7**4
), "Bezout reconstruction fourth-power identity")
zero(sp.together(
    x_reconstructed**7 - X7
    - (X4**7 - X7**4) * (X4**7 + X7**4) / X7**7
), "Bezout reconstruction seventh-power identity")


# ---------------------------------------------------------------------------
# The delayed and actual formal boundary lifts remain globally critical.
# ---------------------------------------------------------------------------

A_lift_5 = sp.expand(A + sp.Rational(5, 6) * L * x**5 * (1 - z))
A_lift_full = sp.expand(
    A + sp.Rational(1, 4) * x**2 * (1 - z)
    + sp.Rational(5, 6) * L * x**5 * (1 - z)
)


def reduced_resultant_residual(first_coordinate: sp.Expr, label: str) -> tuple[sp.Expr, str]:
    derivative_x = sp.diff(first_coordinate, x)
    derivative_t = sp.diff(first_coordinate, t)
    lifted_resultant = sp.resultant(derivative_x, derivative_t, t)
    reduced_resultant = sp.Poly(
        sp.expand(lifted_resultant), L, domain=sp.EX
    ).rem(sp.Poly(L_relation, L, domain=sp.EX)).as_expr()
    residual = sp.cancel(reduced_resultant / x**24)
    gate(sp.expand(reduced_resultant - x**24 * residual) == 0,
         f"{label}: exact x^24 resultant factor")
    gate(sp.Poly(residual, x).degree() == 43,
         f"{label}: residual degree forty-three")
    gate(sp.Poly(residual, x).LC() == -93750,
         f"{label}: stable residual leading coefficient")
    gate(sp.factor(residual.subs(x, 0)) == -384 * L * a**4,
         f"{label}: stable residual constant coefficient")
    gate(sp.Poly(derivative_x, t).LC() == 6 * x**5,
         f"{label}: stable partial-x t-leading coefficient")
    gate(sp.Poly(derivative_t, t).LC() == 6 * x**6,
         f"{label}: stable partial-t t-leading coefficient")
    residual_hash = hashlib.sha256(
        sp.srepr(sp.expand(residual)).encode()
    ).hexdigest()
    return residual, residual_hash


_, lift_5_hash = reduced_resultant_residual(A_lift_5, "delayed x5 lift")
_, lift_full_hash = reduced_resultant_residual(A_lift_full, "actual full lift")


# The a=0 endpoint of the first coordinate remains critical.
u_endpoint = -sp.Rational(2, 3)
x7_endpoint = sp.Rational(16, 243) / L
zero(sp.diff(f, u).subs(u, u_endpoint), "a=0 endpoint u derivative")
zero((-6 * f + 2 * L * x7_endpoint).subs(u, u_endpoint),
     "a=0 endpoint x derivative")
zero_mod_L(
    H.subs(a, 0) - x**5 * (3888 * L * x**7 - 256),
    "a=0 resultant residual endpoint"
)


summary = {
    "checks": CHECKS,
    "family": "A=y^2+2Lx+ap, 6L^2=1, a!=0",
    "boundary": "D maps to cusp (y^2,y^3); L1 maps to (ap,0); both jets J=1",
    "resultant": "Res_t(Ax,At)=(32/3)x^24 H_a, deg H_a=12",
    "address_cover": "9a^7(u+1)(2u+1)^3+32u^5(3u+2)^7=0",
    "conclusion": "every A_a has an affine critical point; no polynomial mate",
    "formal_lifts": {
        "delayed_x5_residual_sha256": lift_5_hash,
        "actual_full_residual_sha256": lift_full_hash,
        "degree": 43,
        "leading": -93750,
        "constant": "-384La^4",
    },
    "endpoint": "a=0 first coordinate also critical; companion undefined",
    "scope": "displayed family closed; arbitrary B2 Darboux and JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3977 simultaneous cusp-arm critical-resultant companion")
print(f"CHECKS={CHECKS}")
print("BOUNDARY_JETS=D_CUSP_AND_L1_LINE;NORMALIZED_JACOBIAN_ONE")
print("RESULTANT=X24_TIMES_DEGREE12_RESIDUAL")
print("CRITICAL_COVER=DEGREE12_WITH_NONZERO_CONSTANT")
print("RECONSTRUCTION=EXACT_X4_X7_BEZOUT")
print("CONCLUSION=AFFINE_CRITICAL_POINT_FOR_EVERY_NONZERO_A")
print("FORMAL_LIFTS=DELAYED_X5_AND_ACTUAL_FULL_STILL_CRITICAL;DEGREE43")
print(f"DELAYED_X5_RESIDUAL_SHA256={lift_5_hash}")
print(f"ACTUAL_FULL_RESIDUAL_SHA256={lift_full_hash}")
print("ENDPOINT=A_ZERO_FIRST_COORDINATE_ALSO_CRITICAL")
print("SCOPE=DISPLAYED_FAMILY_CLOSED;ARBITRARY_B2_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
