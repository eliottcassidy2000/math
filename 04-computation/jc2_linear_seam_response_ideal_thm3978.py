#!/usr/bin/env python3
"""Exact companion for THM-3978's linear-seam response ideals."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    """Record one exact gate which survives optimized execution."""
    global CHECKS
    if not bool(condition):
        raise RuntimeError(f"gate failed: {label}")
    CHECKS += 1


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(expression) == 0, label)


def jacobian(first: sp.Expr, second: sp.Expr, x: sp.Symbol,
             t: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, t)
        - sp.diff(first, t) * sp.diff(second, x)
    )


x, t, c = sp.symbols("x t c")

for n in range(2, 10):
    z = 1 + x**n * t
    p = z * t
    y = x ** (n - 1) * z * t**2
    A = x + c * (z - 1)
    Ax = sp.diff(A, x)
    At = sp.diff(A, t)

    bezout = (
        (1 - c * n * x ** (n - 1) * t) * Ax
        + c * n**2 * x ** (n - 2) * t**2 * At
    )
    zero(bezout - 1, f"n={n}: source-gradient Bezout identity")

    rational_mate = 1 / (c * (n - 1) * x ** (n - 1))
    zero(jacobian(A, rational_mate, x, t) - 1,
         f"n={n}: exact rational constant mate")

    plane_factor = A / x
    plane_primitive = plane_factor ** (n - 1) / (c * (n - 1))
    zero(jacobian(A, plane_primitive, x, t) - A ** (n - 1),
         f"n={n}: minimal plane response")

    S = x + c * (2 * z - 1) + c**2 * x ** (n - 1) * p
    zero(A * (A + c) - x * S,
         f"n={n}: split-boundary response factorization")

    completion_primitive = S ** (n - 1) / (c * (n - 1))
    completion_target = (A * (A + c)) ** (n - 1)
    zero(jacobian(A, completion_primitive, x, t) - completion_target,
         f"n={n}: minimal completion response")

    # D_A does not preserve B_n: this source-polynomial identity has a term
    # of boundary order -1 when rewritten in completion generators.
    pole_expression = (
        (z - 1)**2 / x + 2 * x ** (n - 1) * p
        + c * (n + 1) * x ** (n - 1) * y
    )
    zero(jacobian(A, y, x, t) - pole_expression,
         f"n={n}: completion non-endomorphism pole identity")

# Exact boundary values and endpoint hostiles.
X0, Z0 = sp.symbols("X0 Z0")
A_abstract = X0 + c * (Z0 - 1)
zero(A_abstract.subs({X0: 0, Z0: 0}) + c,
     "added boundary value A=-c")
zero(A_abstract.subs({X0: 0, Z0: 1}),
     "retained arm value A=0")
zero((A_abstract + c).subs({X0: 0, Z0: 1}) - c,
     "opposite factor is a unit on retained arm")

zero(jacobian(x, t, x, t) - 1,
     "c=0 hostile has the source-plane mate t")
A_height_one = x + c * x * t
height_one_address = {x: 0, t: -1 / c}
gate(
    sp.cancel(sp.diff(A_height_one, x).subs(height_one_address)) == 0
    and sp.cancel(sp.diff(A_height_one, t).subs(height_one_address)) == 0,
    "n=1 hostile has a critical point",
)

summary = {
    "checks": CHECKS,
    "family": "A=x+c(z-1) on B_n, n>=2 and c nonzero",
    "submersion": "source Bezout identity and dA|D=dx",
    "rational_solutions": "Q=R(A)/(c(n-1)x^(n-1))+H(A)",
    "plane_response": "J(A,k[x,t]) intersect k[A]=(A^(n-1))",
    "completion_response": "J(A,B_n) intersect k[A]=([A(A+c)]^(n-1))",
    "minimal_completion_primitive": "S=A(A+c)/x in B_n",
    "hostiles": "c=0 plane mate; n=1 critical; J(A,-) not an endomorphism of B_n",
    "scope": "no polynomial constant mate; exact nonconstant response ideal; JC2 open",
}
semantic = hashlib.sha256(
    json.dumps(summary, sort_keys=True).encode()
).hexdigest()

print("THM-3978 linear-seam response-ideal companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=A_X_PLUS_C_Z_MINUS1;N_GE_2;C_NONZERO")
print("SUBMERSION=SOURCE_BEZOUT;BOUNDARY_DA_DX")
print("RATIONAL_SOLUTIONS=Q_R_OVER_C_NMINUS1_X_NMINUS1_PLUS_H")
print("PLANE_RESPONSE_IDEAL=A_NMINUS1")
print("COMPLETION_RESPONSE_IDEAL=A_APLUSC_NMINUS1")
print("MINIMAL_COMPLETION_PRIMITIVE=S_A_APLUSC_OVER_X")
print("HOSTILES=C0_PLANE_MATE;N1_CRITICAL;DA_NOT_B_ENDOMORPHISM")
print("SCOPE=NO_CONSTANT_POLYNOMIAL_MATE;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
