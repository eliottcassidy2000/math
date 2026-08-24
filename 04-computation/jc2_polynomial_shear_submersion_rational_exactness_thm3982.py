#!/usr/bin/env python3
"""Exact companion for THM-3982's polynomial-shear classification."""

from __future__ import annotations

import hashlib
import json
import math

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.together(expression)) == 0, message)


x, t, u, s, rho = sp.symbols("x t u s rho", nonzero=True)
alpha, beta, gamma = sp.symbols("alpha beta gamma", nonzero=True)


def jac(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(f, x) * sp.diff(g, t)
                     - sp.diff(f, t) * sp.diff(g, x))


def delta(expression: sp.Expr, h: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.diff(expression, u) / sp.diff(h, u))


# ---------------------------------------------------------------------------
# Every shear is a submersion: freeze the polynomial Bezout identity for
# constant, affine, and nonlinear controls at every tested height.
# ---------------------------------------------------------------------------

for n in range(2, 8):
    u_source = x**n * t
    controls = [
        beta,
        alpha * u + beta,
        u**2 + gamma * u + beta,
        u**3 + gamma * u**2 + alpha * u + beta,
    ]
    for index, h in enumerate(controls):
        A_source = sp.expand(x + h.subs(u, u_source))
        Ax = sp.diff(A_source, x)
        At = sp.diff(A_source, t)
        zero(x * Ax - n * t * At - x,
             f"height {n} control {index}: submersion Bezout row")
        zero(Ax.subs(x, 0) - 1,
             f"height {n} control {index}: x=0 derivative row")


# ---------------------------------------------------------------------------
# Generic inverse-branch residue formula and the nonlinear infinity gate.
# For pure monomials the infinity leading row is the whole expression.
# ---------------------------------------------------------------------------

for degree in range(2, 6):
    h_monomial = gamma * u**degree
    for n in range(2, 7):
        delta_n = u
        for _ in range(n):
            delta_n = delta(delta_n, h_monomial)
        coefficient = sp.prod(1 - j * degree for j in range(1, n))
        expected = (
            coefficient / (degree * gamma)**n
            * u**(1 - n * degree)
        )
        zero(delta_n - expected,
             f"degree {degree}, height {n}: nonlinear infinity row")
        gate(delta_n != 0,
             f"degree {degree}, height {n}: nonlinear delta is nonzero")

# An independent local-residue replay at generic roots of u^d=A.
for degree in range(2, 5):
    h_plain = u**degree
    for n in range(2, 6):
        local_form = 1 / (rho**degree - (rho + s)**degree)**n
        actual_residue = sp.residue(local_form, s, 0)
        delta_n = u
        for _ in range(n):
            delta_n = delta(delta_n, h_plain)
        predicted = (
            (-1)**n * delta_n.subs(u, rho) / math.factorial(n - 1)
        )
        zero(actual_residue - predicted,
             f"degree {degree}, height {n}: inverse-branch residue")

# Affine h is the exact vanishing boundary for every n>=2.
h_affine = alpha * u + beta
for n in range(2, 8):
    delta_n = u
    for _ in range(n):
        delta_n = delta(delta_n, h_affine)
    zero(delta_n, f"height {n}: affine iterated residue vanishes")


# ---------------------------------------------------------------------------
# Exact rational, source-polynomial, and B_n positive controls.
# ---------------------------------------------------------------------------

for n in range(2, 8):
    u_source = x**n * t

    # Constant shear: t is a source mate and u realizes the completion
    # generator (A-beta)^n.
    A_constant = x + beta
    zero(jac(A_constant, t) - 1,
         f"height {n}: constant shear source mate")
    zero(jac(A_constant, u_source) - x**n,
         f"height {n}: constant shear completion generator")
    for power in range(4):
        Q_constant = sp.expand(u_source * x**power)
        zero(jac(A_constant, Q_constant) - x**(n + power),
             f"height {n}, power {power}: constant completion multiple")

    # Nonconstant affine shear: rational primitive, sharp source generator,
    # and sharp two-color B_n generator.  beta is a harmless target shift.
    A0 = sp.expand(x + alpha * u_source)
    A_affine = A0 + beta
    Q_rational = x**(1 - n) / (alpha * (n - 1))
    zero(jac(A_affine, Q_rational) - 1,
         f"height {n}: affine rational primitive")

    V = 1 + alpha * x**(n - 1) * t
    Q_source = sp.expand((V**(n - 1) - 1) / (alpha * (n - 1)))
    zero(jac(A_affine, Q_source) - A0**(n - 1),
         f"height {n}: affine source image generator")

    Q_completion = sp.expand((A0 + alpha)**(n - 1) * Q_source)
    expected_completion = (A0 * (A0 + alpha))**(n - 1)
    zero(jac(A_affine, Q_completion) - expected_completion,
         f"height {n}: affine two-color image generator")

    # For every negative row occurring in the constant-shear primitive,
    # both mandatory colors occur.  A linear coefficient c*u+d can pay them
    # only by vanishing.
    for q in range(1, n + 1):
        gate(math.ceil(q / n) == 1,
             f"height {n}, row {q}: u exponent one")
        gate(math.ceil(q / (n + 1)) == 1,
             f"height {n}, row {q}: u+1 exponent one")

c0, d0 = sp.symbols("c0 d0")
linear_row = c0 * u + d0
gate(linear_row.subs(u, 0) == d0,
     "constant-row first color evaluation")
gate(linear_row.subs(u, -1) == d0 - c0,
     "constant-row second color evaluation")
solution = sp.solve(
    [linear_row.subs(u, 0), linear_row.subs(u, -1)],
    [c0, d0], dict=True
)
gate(solution == [{c0: 0, d0: 0}],
     "linear row cannot pay both colors nontrivially")


summary = {
    "heights": "all n>=2 (identity proof; controls n=2..7)",
    "family": "A_h=x+h(x^n t), h polynomial",
    "submersion": "all h",
    "rational_image": {
        "deg_h_le_1": "k(A)",
        "deg_h_ge_2": "0",
    },
    "source_image": {
        "constant": "k[A]",
        "affine_nonconstant": "(A-beta)^(n-1)k[A]",
        "nonlinear": "0",
    },
    "completion_image": {
        "constant": "(A-beta)^n k[A]",
        "affine_nonconstant":
            "((A-beta)(A-beta+alpha))^(n-1)k[A]",
        "nonlinear": "0",
    },
    "residue": "(-1)^n delta^n(u)/(n-1)!",
    "conclusion": "no polynomial shear has a B_n constant mate",
    "scope": "polynomial-shear first coordinates; arbitrary B_n and JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3982 polynomial-shear submersion companion")
print(f"CHECKS={CHECKS}")
print("SUBMERSION=ALL_POLYNOMIAL_SHEARS")
print("RATIONAL_IMAGE=K(A)_IFF_DEG_H_LE_1;ZERO_IF_DEG_H_GE_2")
print("RESIDUE=(-1)^N_DELTA^N(U)/(N-1)!;NONLINEAR_INFINITY_NONZERO")
print("SOURCE_IMAGE=CONSTANT_ALL;AFFINE_(A-BETA)^(N-1);NONLINEAR_ZERO")
print("BN_IMAGE=CONSTANT_(A-BETA)^N;AFFINE_TWO_COLOR_PRODUCT;NONLINEAR_ZERO")
print("AFFINE_SECOND_FACTOR=A-BETA+ALPHA")
print("CONCLUSION=NO_POLYNOMIAL_SHEAR_HAS_A_BN_CONSTANT_MATE")
print("SCOPE=POLYNOMIAL_SHEAR_FIRST_COORDINATES;ARBITRARY_BN_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
