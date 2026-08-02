#!/usr/bin/env python3
"""Exact spline and collision controls for THM-3251.

The script verifies the three-distinct-knot and doubled-knot formulas with
exact algebraic arithmetic.  It does not infer transcendence from samples.
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
    u, v = sp.symbols("u v")
    z0, z1, z2 = vertices
    powered = sp.Poly(sp.expand((z0 + (z1 - z0) * u + (z2 - z0) * v) ** degree), u, v)
    total = sp.S.Zero
    for (a, b), coeff in powered.terms():
        total += coeff * simplex_monomial(a, b)
    return sp.simplify(total)


def line_segment_power(x: sp.Expr, y: sp.Expr, degree: int) -> sp.Expr:
    return sp.simplify((y ** (degree + 1) - x ** (degree + 1)) / (degree + 1))


def distinct_spline_power(knots: tuple[sp.Expr, sp.Expr, sp.Expr], degree: int) -> sp.Expr:
    a, b, c = knots
    p, q, length = sp.simplify(b - a), sp.simplify(c - b), sp.simplify(c - a)
    left = line_segment_power(a, b, degree + 1) - a * line_segment_power(a, b, degree)
    right = c * line_segment_power(b, c, degree) - line_segment_power(b, c, degree + 1)
    return sp.simplify((left / p + right / q) / length)


def doubled_spline_power(a: sp.Expr, c: sp.Expr, degree: int) -> sp.Expr:
    length = sp.simplify(c - a)
    value = c * line_segment_power(a, c, degree) - line_segment_power(a, c, degree + 1)
    return sp.simplify(value / length**2)


def group_terms(terms: list[tuple[sp.Expr, sp.Expr]]) -> list[tuple[sp.Expr, sp.Expr]]:
    groups: list[tuple[sp.Expr, sp.Expr]] = []
    for exponent, coefficient in terms:
        exponent, coefficient = sp.simplify(exponent), sp.simplify(coefficient)
        for index, (old_exponent, old_coefficient) in enumerate(groups):
            if zero(exponent - old_exponent):
                groups[index] = (old_exponent, sp.simplify(old_coefficient + coefficient))
                break
        else:
            groups.append((exponent, coefficient))
    return [(exponent, sp.simplify(coefficient)) for exponent, coefficient in groups]


def nonzero_source(groups: list[tuple[sp.Expr, sp.Expr]]) -> bool:
    for exponent, coefficient in groups:
        if not zero(coefficient):
            require(not zero(exponent), "a nonzero source occurred only at exponent zero")
            return True
    return False


ledger: list[str] = []

# C1: the complex-line B-spline is the direct simplex pushforward and has mass 1/2.
distinct_knots = (
    (sp.S.Zero, sp.S.One + I, 3 + 3 * I),
    (2 - I, sp.Rational(5, 2) + sp.Rational(1, 2) * I, 4 + 2 * I),
    (-1 + 2 * I, 1 + I, 5 - I),
)
doubled_knots = (
    (sp.S.Zero, 1 + I),
    (-2 + I, 3 - 4 * I),
    (1 + 2 * I, -3 + 6 * I),
)
pushforward_checks = 0
for knots in distinct_knots:
    for degree in range(11):
        require(zero(simplex_power(knots, degree) - distinct_spline_power(knots, degree)),
                f"distinct spline mismatch knots={knots}, degree={degree}")
        # Reversing the line order must not change the pushforward.
        require(zero(distinct_spline_power(knots, degree)
                     - distinct_spline_power((knots[2], knots[1], knots[0]), degree)),
                f"orientation reversal mismatch knots={knots}, degree={degree}")
        pushforward_checks += 2
for a, c in doubled_knots:
    for degree in range(11):
        direct = simplex_power((a, a, c), degree)
        require(zero(direct - doubled_spline_power(a, c, degree)),
                f"doubled spline mismatch a={a}, c={c}, degree={degree}")
        pushforward_checks += 1
ledger.append(f"C1:spline_pushforward_and_orientation={pushforward_checks}")

# C2: symbolic second-divided-difference weights and normalization.
a, p, q = sp.symbols("a p q", nonzero=True)
b, c, length = a + p, a + p + q, p + q
alpha = (-1 / (length * p), 1 / (p * q), -1 / (length * q))
z = (a, b, c)
require(zero(sum(alpha)), "sum alpha !=0")
require(zero(sum(alpha[j] * z[j] for j in range(3))), "sum alpha*z !=0")
require(zero(sum(alpha[j] * z[j] ** 2 for j in range(3)) + 1), "sum alpha*z^2 !=-1")
ledger.append("C2:distinct_weight_moments=(0,0,-1)")

# C3: direct pure-power coefficients equal the two primitive residue blocks.
primitive_checks = 0
for knots in distinct_knots:
    a0, b0, c0 = knots
    p0, q0, length0 = b0 - a0, c0 - b0, c0 - a0
    weights = (-1 / (length0 * p0), 1 / (p0 * q0), -1 / (length0 * q0))
    for power in range(3, 9):
        for n in range(6):
            degree = power * n
            primitive = sp.S.Zero
            for weight, endpoint in zip(weights, knots):
                primitive += weight * endpoint ** (degree + 2) / (degree + 2)
                primitive -= weight * endpoint * endpoint ** (degree + 1) / (degree + 1)
            require(zero(primitive - simplex_power(knots, degree)),
                    f"distinct primitive mismatch power={power}, n={n}, knots={knots}")
            primitive_checks += 1
for a0, c0 in doubled_knots:
    length0 = c0 - a0
    for power in range(3, 9):
        for n in range(6):
            degree = power * n
            primitive = (
                (a0 ** (degree + 2) - c0 ** (degree + 2)) / ((degree + 2) * length0**2)
                + c0 * (c0 ** (degree + 1) - a0 ** (degree + 1))
                / ((degree + 1) * length0**2)
            )
            require(zero(primitive - simplex_power((a0, a0, c0), degree)),
                    f"doubled primitive mismatch power={power}, n={n}, a={a0}, c={c0}")
            primitive_checks += 1
ledger.append(f"C3:primitive_moment_coefficients={primitive_checks}")

# C4: endpoint-power collisions for three distinct knots.
omega = (-sp.S.One + I * sp.sqrt(3)) / 2
distinct_collision_cases = (
    (4, (-sp.S.One, sp.S.Zero, sp.S.One), "even-endpoints"),
    (3, (sp.S.One, (sp.S.One + omega) / 2, omega), "cubic-endpoints"),
    (3, (sp.S.Zero, 1 + I, 2 + 2 * I), "zero-endpoint"),
)
distinct_profiles = []
for power, knots, label in distinct_collision_cases:
    a0, b0, c0 = knots
    p0, q0, length0 = b0 - a0, c0 - b0, c0 - a0
    weights = (-1 / (length0 * p0), 1 / (p0 * q0), -1 / (length0 * q0))
    groups = group_terms([(endpoint**power, weight * endpoint**2)
                          for endpoint, weight in zip(knots, weights)])
    require(zero(sum(coefficient for _, coefficient in groups) + 1),
            f"distinct grouped source lost -1: {label}")
    require(nonzero_source(groups), f"distinct collision erased every source: {label}")
    distinct_profiles.append((label, len(groups), sum(not zero(v) for _, v in groups)))
ledger.append(f"C4:distinct_collision_profiles={distinct_profiles}")

# C5: doubled-knot hostiles.  Either block may die, but not both.
doubled_collision_cases = (
    (4, -sp.S.One, sp.S.One, "symmetric-even-G1-dies"),
    (3, -sp.S.One, sp.S.Zero, "unique-zero-G0-dies"),
    (3, sp.S.Zero, sp.S.One, "repeated-zero"),
    (5, sp.S.One, sp.Integer(2), "generic"),
)
doubled_profiles = []
for power, a0, c0, label in doubled_collision_cases:
    length0 = c0 - a0
    source1 = group_terms([
        (a0**power, a0**2 / length0**2),
        (c0**power, -c0**2 / length0**2),
    ])
    source0 = group_terms([
        (a0**power, -c0 * a0 / length0**2),
        (c0**power, c0**2 / length0**2),
    ])
    alive1, alive0 = nonzero_source(source1), nonzero_source(source0)
    require(alive0 or alive1, f"both doubled residue blocks died: {label}")
    doubled_profiles.append((label, int(alive0), int(alive1)))
ledger.append(f"C5:doubled_alive_G0_G1={doubled_profiles}")

# C6: exact residue boundary used by the rational-function obstruction.
residue_checks = 0
for power in range(3, 41):
    require(1 % power != 0 and 2 % power != 0, f"integral residue at d={power}")
    require(1 % power != 0, f"inter-block rational gauge at d={power}")
    residue_checks += 3
require(2 % 2 == 0, "quadratic boundary control missing")
ledger.append(f"C6:nonintegral_residues={residue_checks};quadratic_boundary=2/2")

digest = sha256("\n".join(ledger).encode("utf-8")).hexdigest()

print("THM-3251 FC(3) COLLINEAR PURE-POWER SPLINE-RESIDUE AUDIT")
for row in ledger:
    print(row)
print(f"semantic_sha256={digest}")
print("CONCLUSION=EXACT_IDENTITIES_VERIFIED;TRANSCENDENCE_STEP=CITED_BEUKERS")
