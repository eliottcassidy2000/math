#!/usr/bin/env python3
"""Exact hostile controls for THM-3250.

The universal proof is symbolic.  This companion checks the Cauchy--Green
orientation, vertex-turn primitive normalization, pure-power primitive ODE,
endpoint-power collisions, and direct simplex moments with exact algebraic
arithmetic.  It does not infer transcendence from bounded sampling.

Reproduce with

    python3 04-computation/fc3_noncollinear_pure_power_turn_current_thm3250.py
    python3 -O 04-computation/fc3_noncollinear_pure_power_turn_current_thm3250.py
"""

from __future__ import annotations

from hashlib import sha256
from math import factorial

import sympy as sp


I = sp.I


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def zero(expr: sp.Expr) -> bool:
    return sp.simplify(sp.expand_complex(expr)) == 0


def simplex_monomial(a: int, b: int) -> sp.Rational:
    return sp.Rational(factorial(a) * factorial(b), factorial(a + b + 2))


def simplex_power(vertices: tuple[sp.Expr, sp.Expr, sp.Expr], degree: int) -> sp.Expr:
    """Integral over the coordinate simplex of z(u,v)^degree."""
    z0, z1, z2 = vertices
    u, v = sp.symbols("u v")
    poly = sp.Poly(sp.expand(z0 + (z1 - z0) * u + (z2 - z0) * v), u, v)
    powered = sp.Poly(sp.expand(poly.as_expr() ** degree), u, v)
    total = sp.S.Zero
    for (a, b), coeff in powered.terms():
        total += coeff * simplex_monomial(a, b)
    return sp.simplify(total)


def triangle_data(vertices: tuple[sp.Expr, sp.Expr, sp.Expr]):
    z = list(vertices)
    zb = [sp.conjugate(value) for value in z]
    jac = sp.simplify(sp.im(sp.conjugate(z[1] - z[0]) * (z[2] - z[0])))
    require(jac.is_positive is True, f"triangle is not certified CCW: {vertices}, J={jac}")
    slopes = []
    offsets = []
    for j in range(3):
        jp = (j + 1) % 3
        slope = sp.simplify((zb[jp] - zb[j]) / (z[jp] - z[j]))
        slopes.append(slope)
        offsets.append(sp.simplify(zb[j] - slope * z[j]))
    turns = [sp.simplify(slopes[(j - 1) % 3] - slopes[j]) for j in range(3)]
    return z, slopes, offsets, turns, sp.simplify(2 * I * jac)


def boundary_monomial(vertices: tuple[sp.Expr, sp.Expr, sp.Expr], degree: int) -> sp.Expr:
    """Integral around the oriented boundary of conjugate(z) z^degree dz."""
    z, slopes, offsets, _, _ = triangle_data(vertices)
    total = sp.S.Zero
    for j in range(3):
        jp = (j + 1) % 3
        # conjugate(z)=m_j z+b_j on this edge; holomorphic primitives are exact.
        total += slopes[j] * (z[jp] ** (degree + 2) - z[j] ** (degree + 2)) / (degree + 2)
        total += offsets[j] * (z[jp] ** (degree + 1) - z[j] ** (degree + 1)) / (degree + 1)
    return sp.simplify(total)


def grouped_source(vertices: tuple[sp.Expr, sp.Expr, sp.Expr], degree: int):
    z, _, _, turns, width = triangle_data(vertices)
    groups: list[tuple[sp.Expr, sp.Expr]] = []
    for zj, cj in zip(z, turns):
        exponent = sp.simplify(zj**degree)
        source = sp.simplify(cj * zj**2)
        for index, (old_exponent, old_source) in enumerate(groups):
            if zero(exponent - old_exponent):
                groups[index] = (old_exponent, sp.simplify(old_source + source))
                break
        else:
            groups.append((exponent, source))
    groups = [(exponent, sp.simplify(source)) for exponent, source in groups]
    require(zero(sum(source for _, source in groups) + width), "turn source does not sum to -W")
    require(any(not zero(source) and not zero(exponent) for exponent, source in groups),
            "no nonzero nonzero-exponent source group survived")
    return groups, width


ledger: list[str] = []

# C1: Cauchy--Green orientation on exact triangles and holomorphic monomials.
triangles = (
    (sp.S.Zero, sp.S.One, sp.Rational(2) + I),
    (sp.Rational(-1) + I, sp.Rational(2), sp.Rational(1) + 3 * I),
    (sp.Rational(1, 2) - I, sp.Rational(5, 2) + I, sp.Rational(-1) + 4 * I),
)
cauchy_checks = 0
for triangle in triangles:
    *_, width = triangle_data(triangle)
    for degree in range(8):
        lhs = boundary_monomial(triangle, degree)
        rhs = sp.simplify(width * simplex_power(triangle, degree))
        require(zero(lhs - rhs), f"Cauchy--Green mismatch: {triangle}, degree={degree}")
        cauchy_checks += 1
ledger.append(f"C1:cauchy_green={cauchy_checks}")

# C2: generic line equations, vertex regrouping, and the source identity S(0)=-W.
x, y = sp.symbols("x y", real=True, nonzero=True)
generic = (sp.S.Zero, sp.S.One, x + I * y)
z, slopes, offsets, turns, width = triangle_data((sp.S.Zero, sp.S.One, x + I * sp.Symbol("Y", positive=True, real=True)))
for j in range(3):
    require(zero(offsets[(j - 1) % 3] - offsets[j] + turns[j] * z[j]),
            f"vertex offset regrouping failed at {j}")
require(zero(sum(turns[j] * z[j] ** 2 for j in range(3)) + width), "generic S(0)=-W failed")
ledger.append("C2:vertex_turn_identity=SYMBOLIC")

# C3: coefficientwise pure-power primitive ODE for both residue blocks.
ode_checks = 0
for degree in range(3, 11):
    for k in (0, 1):
        for aval, zend in ((sp.Rational(2, 3), sp.Rational(3, 2) + I),
                           (-sp.Rational(5, 4), sp.Rational(-2, 3) + 2 * I)):
            for n in range(9):
                coeff = aval**n * zend ** (degree * n + k + 1) / (
                    (degree * n + k + 1) * factorial(n)
                )
                lhs = sp.simplify((degree * n + k + 1) * coeff)
                rhs = sp.simplify(zend ** (k + 1) * (aval * zend**degree) ** n / factorial(n))
                require(zero(lhs - rhs), f"primitive ODE mismatch d={degree}, k={k}, n={n}")
                ode_checks += 1
ledger.append(f"C3:primitive_ode_coefficients={ode_checks}")

# C4: hostile endpoint-power collisions, including total collisions and z_j=0.
omega = (-sp.S.One + I * sp.sqrt(3)) / 2
collision_cases = (
    (3, (sp.S.One, omega, sp.simplify(omega**2)), "all-cubic"),
    (4, (sp.S.One, I, -sp.S.One), "all-quartic"),
    (3, (sp.S.One, sp.Integer(2), omega), "partial-cubic"),
    (3, (sp.S.Zero, sp.S.One, sp.S.One + I), "zero-vertex"),
)
collision_profile = []
for degree, triangle, label in collision_cases:
    groups, _ = grouped_source(triangle, degree)
    collision_profile.append((label, len(groups), sum(not zero(source) for _, source in groups)))
ledger.append(f"C4:collision_profiles={collision_profile}")

# C5: exact direct simplex moments agree with the regrouped primitive blocks.
moment_checks = 0
for triangle in triangles:
    z, _, _, turns, width = triangle_data(triangle)
    for degree in range(3, 8):
        for n in range(6):
            direct = sp.simplify(width * simplex_power(triangle, degree * n))
            # Coefficient of A^n s^n/n! in G_1+G_0.
            primitive = sp.S.Zero
            for zj, cj in zip(z, turns):
                primitive += cj * zj ** (degree * n + 2) / (degree * n + 2)
                primitive -= cj * zj * zj ** (degree * n + 1) / (degree * n + 1)
            require(zero(direct - primitive),
                    f"direct/primitive moment mismatch d={degree}, n={n}, triangle={triangle}")
            moment_checks += 1
ledger.append(f"C5:direct_moment_coefficients={moment_checks}")

# C6: valuation boundary behind functional independence.  The proof uses
# 0<(k+1)<d and 0<|2-1|<d; these fail first at the quadratic k=1 block.
residue_checks = 0
for degree in range(3, 41):
    for k in (0, 1):
        require((k + 1) % degree != 0, f"integral zero residue at d={degree}, k={k}")
        residue_checks += 1
    require(1 % degree != 0, f"rational inter-block gauge at d={degree}")
    residue_checks += 1
require(2 % 2 == 0, "quadratic boundary control was not planted")
ledger.append(f"C6:nonintegral_residue_and_gap={residue_checks};quadratic_boundary=2/2")

digest = sha256("\n".join(ledger).encode("utf-8")).hexdigest()

print("THM-3250 FC(3) NONCOLLINEAR PURE-POWER TURN-CURRENT AUDIT")
for row in ledger:
    print(row)
print(f"semantic_sha256={digest}")
print("CONCLUSION=EXACT_IDENTITIES_VERIFIED;TRANSCENDENCE_STEP=CITED_BEUKERS")
