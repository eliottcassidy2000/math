#!/usr/bin/env python3
"""Exact companion for THM-3923's four-ray address obstruction.

Reproduction:
  python3 04-computation/jc2_binary_cubic_four_ray_tangent_address_thm3923.py
  python3 -O 04-computation/jc2_binary_cubic_four_ray_tangent_address_thm3923.py
"""

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
    gate(sp.expand(expression) == 0, message)


def homogeneous_piece(expression: sp.Expr, variables: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), *variables)
    return sp.expand(
        sum(
            coefficient * sp.prod(variable**power for variable, power in zip(variables, monomial))
            for monomial, coefficient in polynomial.terms()
            if sum(monomial) == degree
        )
    )


A, C, X, Y, t, lam = sp.symbols("A C X Y t lambda")


def discriminant(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(b**2 * c**2 - 4 * a * c**3 - 4 * b**3 * d - 27 * a**2 * d**2 + 18 * a * b * c * d)


# ---------------------------------------------------------------------------
# The fixed THM-3808 linear packet and its four distinct tangent rays.
# ---------------------------------------------------------------------------

a1, b1, c1, d1 = A, C, 7 * A, -3 * A
Delta0 = sp.factor(discriminant(a1, b1, c1, d1))
expected = A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
zero(Delta0 - expected, "THM-3808 discriminant factorization")

ray_vectors = ((1, 0), (5, 1), (19, 4), (-17, 3))
for i, left in enumerate(ray_vectors):
    for j, right in enumerate(ray_vectors[:i]):
        gate(left[0] * right[1] - left[1] * right[0] != 0, f"rays {j},{i} are distinct")

dehomogenized = sp.Poly(sp.expand(Delta0.subs({A: t, C: 1})), t)
gate(dehomogenized.degree() == 4, "four finite projective roots in the C-chart")
gate(sp.gcd(dehomogenized, dehomogenized.diff()) == 1, "tangent quartic is squarefree")


# ---------------------------------------------------------------------------
# Arbitrary quadratic coefficient jets cannot change the degree-four packet.
# ---------------------------------------------------------------------------

u = sp.symbols("u0:12")
quadratics = (
    u[0] * A**2 + u[1] * A * C + u[2] * C**2,
    u[3] * A**2 + u[4] * A * C + u[5] * C**2,
    u[6] * A**2 + u[7] * A * C + u[8] * C**2,
    u[9] * A**2 + u[10] * A * C + u[11] * C**2,
)
Delta_jet = discriminant(
    a1 + quadratics[0],
    b1 + quadratics[1],
    c1 + quadratics[2],
    d1 + quadratics[3],
)
zero(homogeneous_piece(Delta_jet, (A, C), 4) - Delta0, "fixed linear row fixes tangent discriminant")
gate(homogeneous_piece(Delta_jet - Delta0, (A, C), 0) == 0, "no constant discriminant correction")
gate(homogeneous_piece(Delta_jet - Delta0, (A, C), 1) == 0, "no linear discriminant correction")
gate(homogeneous_piece(Delta_jet - Delta0, (A, C), 2) == 0, "no quadratic discriminant correction")
gate(homogeneous_piece(Delta_jet - Delta0, (A, C), 3) == 0, "no cubic discriminant correction")


# A tangent-identity base change also fixes the four-ray tangent cone.
p20, p11, p02, q20, q11, q02 = sp.symbols("p20 p11 p02 q20 q11 q02")
P = A + p20 * A**2 + p11 * A * C + p02 * C**2
Q = C + q20 * A**2 + q11 * A * C + q02 * C**2
Delta_base = sp.expand(P * (Q + 5 * P) * (4 * Q + 19 * P) * (3 * Q - 17 * P))
zero(homogeneous_piece(Delta_base, (A, C), 4) - Delta0, "tangent-identity base change preserves four rays")


# ---------------------------------------------------------------------------
# THM-3853's two inverse-discriminant targets: exact four-address controls.
# ---------------------------------------------------------------------------

target_rows: list[dict[str, object]] = []
for label, L, M in (("C", C, A), ("A+C", A + C, A)):
    # Write Delta0=L^4 D(M/L), then use L=-D/lambda, M=tL.
    if label == "C":
        D = -t * (5 * t + 1) * (17 * t - 3) * (19 * t + 4)
        A_param = sp.expand(-t * D / lam)
        C_param = sp.expand(-D / lam)
    else:
        D = -t * (4 * t + 1) * (15 * t + 4) * (20 * t - 3)
        A_param = sp.expand(-t * D / lam)
        L_param = sp.expand(-D / lam)
        C_param = sp.expand(L_param - A_param)

    delta_target = sp.expand(Delta0 + lam * L**5)
    zero(sp.together(delta_target.subs({A: A_param, C: C_param})), f"{label} target parametrization")
    D_poly = sp.Poly(D, t)
    gate(D_poly.degree() == 4, f"{label} has four finite addresses")
    gate(sp.gcd(D_poly, D_poly.diff()) == 1, f"{label} addresses are distinct")
    gate(sp.Poly(A_param, t).degree() == 5, f"{label} first coordinate degree five")
    gate(max(sp.Poly(A_param, t).degree(), sp.Poly(C_param, t).degree()) == 5,
         f"{label} has one polynomial infinity place")
    target_rows.append(
        {
            "orientation": label,
            "address_polynomial": str(sp.factor(D)),
            "address_discriminant": str(sp.factor(sp.discriminant(D, t))),
        }
    )


# Arbitrary order-five target perturbations preserve the quartic tangent cone.
v = sp.symbols("v0:6")
Phi5 = sum(v[i] * A ** (5 - i) * C**i for i in range(6))
zero(homogeneous_piece(Delta0 + Phi5, (A, C), 4) - Delta0, "all order-five targets preserve four rays")


# Scope control: three distinct tangent directions are not rejected by r <= 3.
three_ray_control = A**2 * C * (A - C)
three_ray_factors = sp.factor_list(three_ray_control)[1]
gate(len(three_ray_factors) == 3, "three-ray boundary has exactly three factors")
gate(sorted(exponent for _, exponent in three_ray_factors) == [1, 1, 2],
     "three-ray boundary has partition 2+1+1")


summary = {
    "checks": CHECKS,
    "delta0": str(Delta0),
    "fixed_linear_tangent_degree": sp.Poly(Delta0, A, C).total_degree(),
    "tangent_ray_count": len(ray_vectors),
    "targets": target_rows,
    "scope_boundary": "three distinct tangent rays meet, rather than violate, the cubic address cap",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3923 exact tangent-address companion")
print(f"CHECKS={CHECKS}")
print(f"DELTA0={Delta0}")
for row in target_rows:
    print(
        "TARGET="
        f"{row['orientation']};ADDR={row['address_polynomial']};"
        f"DISC={row['address_discriminant']}"
    )
print("BOUNDARY=three rays are not excluded by the rank-three address cap")
print(f"SEMANTIC_SHA256={semantic}")
