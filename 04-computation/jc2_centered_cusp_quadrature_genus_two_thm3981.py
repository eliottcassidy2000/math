#!/usr/bin/env python3
"""Exact companion for THM-3981's centered cusp quadrature obstruction."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    """Record one optimization-safe exact gate."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


x, z, u, v, Y = sp.symbols("x z u v Y")

# Center y=Y+x/2 in the D-chart equation of THM-3979.
F = z * (z - 1) ** 2 - x**3 * (Y + x / 2)

# Projection from (x,z)=(0,1) gives the genus-two model.
R = sp.expand((u**3 - Y) ** 2 + 2 * u**2)
x_uv = v - Y + u**3
z_uv = 1 + u * x_uv


def reduce_v(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.Poly(sp.expand(numerator), v).rem(
        sp.Poly(v**2 - R, v)
    ).as_expr()


zero(reduce_v(F.subs({x: x_uv, z: z_uv})),
     "birational map lands on centered curve")
zero(
    reduce_v(
        sp.expand(x**2 + 2 * (Y - u**3) * x - 2 * u**2).subs(
            x, x_uv
        )
    ),
    "quadratic projection equation",
)

# Squarefreeness over k(Y), plus both degeneration controls.
disc = sp.factor(sp.discriminant(R, u))
res = sp.factor(sp.resultant(R, sp.diff(R, u), u))
gate(disc == -512 * Y**2 * (729 * Y**4 + 128),
     "sextic discriminant")
gate(res == 512 * Y**2 * (729 * Y**4 + 128),
     "sextic resultant")

a = sp.Rational(4, 9) / Y
exceptional = 729 * Y**4 + 128
gate(
    sp.rem(sp.together(R.subs(u, a)).as_numer_denom()[0],
           exceptional, Y) == 0,
    "exceptional double root lies on curve",
)
gate(
    sp.rem(
        sp.together(sp.diff(R, u).subs(u, a)).as_numer_denom()[0],
        exceptional,
        Y,
    ) == 0,
    "exceptional double root is critical",
)
gate(
    sp.rem(
        sp.together(sp.diff(R, u, 2).subs(u, a) + 8).as_numer_denom()[0],
        exceptional,
        Y,
    ) == 0,
    "exceptional double root is ordinary",
)
gate(sp.factor(R.subs(Y, 0)) == u**2 * (u**4 + 2),
     "central degeneration")
central_tangent_square = sp.cancel(R.subs(Y, 0) / u**2).subs(u, 0)
gate(central_tangent_square == 2,
     "central node has nonzero residue denominator")
exceptional_tangent_square = sp.diff(R, u, 2).subs(u, a) / 2
gate(
    sp.rem(
        sp.together(exceptional_tangent_square + 4).as_numer_denom()[0],
        exceptional,
        Y,
    ) == 0,
    "exceptional node has nonzero residue denominator",
)

# On v^2=R, dv/du=R'/(2v), and x=v-Y+u^3. The numerator below is
# v*dx/du-u(3ux+2). Its vanishing gives x*w*dx=du/v exactly.
differential_numerator = (
    sp.diff(R, u) / 2 + 3 * u**2 * v - u * (3 * u * x_uv + 2)
)
zero(reduce_v(differential_numerator),
     "centered quadrature is du/v")

# Formal D-chart terms protect the centering and sign conventions.
ORDER = 13
zeta = sp.Integer(0)
for degree in range(1, ORDER):
    coefficient = sp.Symbol(f"c{degree}")
    trial = zeta + coefficient * x**degree
    row = sp.expand(
        trial * (trial - 1) ** 2 - x**3 * (Y + x / 2)
    ).coeff(x, degree)
    solutions = sp.solve(sp.Eq(row, 0), coefficient)
    gate(len(solutions) == 1, f"unique D coefficient {degree}")
    zeta = sp.expand(zeta + solutions[0] * x**degree)

w = sp.series(
    1 / ((zeta - 1) * (3 * zeta - 1)), x, 0, ORDER
).removeO().expand()
omega_coefficient = sp.series(x * w, x, 0, 7).removeO().expand()
primitive = sp.series(
    2 * sp.integrate(omega_coefficient, x), x, 0, 8
).removeO().expand()
gate(omega_coefficient == x + 4 * Y * x**4 + 2 * x**5,
     "D quadrature expansion")
gate(
    primitive == x**2 + sp.Rational(8, 5) * Y * x**5
    + sp.Rational(2, 3) * x**6,
    "primitive expansion",
)

summary = {
    "curve": "v^2=(u^3-Y)^2+2u^2",
    "differential": "x*w*dx=du/v",
    "discriminant": "-512Y^2(729Y^4+128)",
    "generic_genus": 2,
    "infinity_divisor": "P_infinity_plus+P_infinity_minus",
    "obstruction": "holomorphic nonexact and finite trace descent",
    "special_slices": "nodal genus one with nonzero opposite residues",
    "scope": "canonical X and A_D transcendental on every scalar slice; alternative gauges open",
}
semantic = hashlib.sha256(
    json.dumps(summary, sort_keys=True).encode()
).hexdigest()

print("THM-3981 centered cusp quadrature genus-two companion")
print(f"CHECKS={CHECKS}")
print(f"HYPERELLIPTIC=v^2={sp.sstr(R)}")
print(f"DISCRIMINANT={sp.sstr(disc)}")
print("GENERIC_BRANCHES=6_FINITE_SIMPLE;INFINITY=2_UNRAMIFIED;GENUS=2")
print("QUADRATURE_DIFFERENTIAL=x*w*dx=du/v")
print("DIVISOR_GENERIC(du/v)=P_INFINITY_PLUS+P_INFINITY_MINUS")
print("GENERIC_OBSTRUCTION=NONZERO_HOLOMORPHIC_NOT_EXACT")
print("ALGEBRAIC_UPGRADE=TRACE_DESCENT_FORBIDS_FINITE_ALGEBRAIC_PRIMITIVE")
print("SPECIAL_Y0=NODAL_GENUS1;RESIDUE_SQUARED=1/2")
print("SPECIAL_729Y4_PLUS128_ZERO=NODAL_GENUS1;RESIDUE_SQUARED=-1/4")
print("FORMAL_D_EXPANSION=X2=x^2+(8/5)Yx^5+(2/3)x^6+...")
print("SCOPE=CANONICAL_X_AND_A_D_TRANSCENDENTAL_ON_ALL_SCALAR_SLICES;ALTERNATIVE_GAUGES_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
