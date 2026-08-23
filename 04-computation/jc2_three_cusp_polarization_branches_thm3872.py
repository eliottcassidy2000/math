#!/usr/bin/env python3
"""Exact companion for THM-3872's cusp-polarization branch packet."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


t, x, y = sp.symbols("t x y")
x_t = t**4 - 2 * t**2
y_t = 3 * t**5 - 5 * t**3
delta = sp.expand(
    81 * x**5
    + 90 * x**4
    + 25 * x**3
    + 30 * x**2 * y**2
    + 30 * x * y**2
    - y**4
    + 8 * y**2
)

# THM-3864's three defect representatives and their cusp-value selectors.
h1 = t**2 * (t**2 - 1) * (9 * t**2 - 14)
h2 = t * (t**2 - 1) * (9 * t**2 - 4) * (3 * t**2 - 5)
h3 = t**3 * (t**2 - 1) * (3 * t**2 - 5) * (
    9 * t**4 + 24 * t**2 - 38
)
e0 = x + 1
e_plus = (-2 * x - y) / 4
e_minus = (-2 * x + y) / 4
cusp_addresses = (0, 1, -1)

evaluation_matrix = sp.Matrix(
    [
        [selector.subs({x: x_t, y: y_t}).subs(t, address) for selector in
         (e0, e_plus, e_minus)]
        for address in cusp_addresses
    ]
)
gate(evaluation_matrix == sp.eye(3), "cusp-value selector matrix")
for index, selector in enumerate((e0, e_plus, e_minus)):
    selector_t = selector.subs({x: x_t, y: y_t})
    gate(
        [sp.diff(selector_t, t).subs(t, address) for address in cusp_addresses]
        == [0, 0, 0],
        f"selector e{index} lies in the branch ring",
    )

derivative_matrix = sp.Matrix(
    [
        [sp.diff(h, t).subs(t, address) for h in (h1, h2, h3)]
        for address in cusp_addresses
    ]
)
gate(
    derivative_matrix
    == sp.Matrix([[0, -20, 0], [-10, -20, 20], [10, -20, 20]]),
    "defect derivative matrix",
)
gate(derivative_matrix.det() == -8000, "defect derivative coordinates")

# Canonical square/cube descents for the entire defect plane.
P1 = 81 * x**3 + 49 * x**2 + 8 * y**2
P2 = -648 * x**3 - 720 * x**2 + 81 * x * y**2 - 200 * x + 49 * y**2
P3 = (
    81 * x**3 * y**2
    + 657 * x**2 * y**2
    + 2356 * x * y**2
    + 84 * y**4
    + 1444 * y**2
)
P12 = 81 * x**2 * y + 137 * x * y + 56 * y
P13 = 81 * x**3 * y - 369 * x**2 * y - 266 * x * y + 46 * y**3
P23 = (
    3078 * x**4
    + 3420 * x**3
    + 81 * x**2 * y**2
    + 950 * x**2
    + 555 * x * y**2
    + 322 * y**2
)

Q1 = -243 * x**4 - 143 * x**3 + 81 * x**2 * y**2 + 120 * x * y**2 + 64 * y**2
Q2 = (
    6561 * x**4 * y
    + 17010 * x**3 * y
    + 18009 * x**2 * y
    + 243 * x * y**3
    + 8760 * x * y
    + 143 * y**3
    + 1600 * y
)
Q3 = (
    8991 * x**4 * y**3
    + 21897 * x**3 * y**3
    + 81 * x**2 * y**5
    + 92226 * x**2 * y**3
    + 5364 * x * y**5
    + 134292 * x * y**3
    + 5308 * y**5
    + 54872 * y**3
)
Q112 = -891 * x**3 * y - 1183 * x**2 * y + 81 * x * y**3 - 392 * x * y + 56 * y**3
Q122 = (
    4536 * x**4
    + 5040 * x**3
    - 1539 * x**2 * y**2
    + 1400 * x**2
    - 1247 * x * y**2
    + 81 * y**4
    - 256 * y**2
)
Q113 = (
    2835 * x**4 * y
    + 4797 * x**3 * y
    + 81 * x**2 * y**3
    + 1862 * x**2 * y
    + 424 * x * y**3
    + 368 * y**3
)
Q133 = (
    5913 * x**4 * y**2
    - 6147 * x**3 * y**2
    + 81 * x**2 * y**4
    - 22268 * x**2 * y**2
    + 2172 * x * y**4
    - 10108 * x * y**2
    + 2116 * y**4
)
Q223 = (
    -20088 * x**4 * y
    - 46944 * x**3 * y
    + 1539 * x**2 * y**3
    - 33560 * x**2 * y
    + 3693 * x * y**3
    - 7600 * x * y
    + 81 * y**5
    + 1606 * y**3
)
Q233 = (
    116964 * x**4 * y
    + 5265 * x**3 * y**3
    + 129960 * x**3 * y
    - 12051 * x**2 * y**3
    + 36100 * x**2 * y
    + 81 * x * y**5
    - 3052 * x * y**3
    + 1500 * y**5
    + 2812 * y**3
)
Q123 = (
    6561 * x**6
    + 11826 * x**5
    + 7065 * x**4
    + 1400 * x**3
    + 4617 * x**3 * y**2
    + 11211 * x**2 * y**2
    + 9270 * x * y**2
    + 2576 * y**2
)


def canonical_descents(c1: sp.Expr, c2: sp.Expr, c3: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    square = sp.expand(
        c1**2 * P1
        + c2**2 * P2
        + c3**2 * P3
        + 2 * c1 * c2 * P12
        + 2 * c1 * c3 * P13
        + 2 * c2 * c3 * P23
    )
    cube = sp.expand(
        c1**3 * Q1
        + c2**3 * Q2
        + c3**3 * Q3
        + 3 * c1**2 * c2 * Q112
        + 3 * c1 * c2**2 * Q122
        + 3 * c1**2 * c3 * Q113
        + 3 * c1 * c3**2 * Q133
        + 3 * c2**2 * c3 * Q223
        + 3 * c2 * c3**2 * Q233
        + 6 * c1 * c2 * c3 * Q123
    )
    return square, cube


def shifted_residual(
    c1: sp.Expr,
    c2: sp.Expr,
    c3: sp.Expr,
    branch_addition: sp.Expr,
    mixed_descent: sp.Expr,
) -> sp.Expr:
    """Residual after H=h_C+r, where mixed_descent pulls back to r*h_C."""

    square0, cube0 = canonical_descents(c1, c2, c3)
    square = sp.expand(square0 + 2 * mixed_descent + branch_addition**2)
    cube = sp.expand(
        cube0
        + 3 * branch_addition * square0
        + 3 * branch_addition * mixed_descent
        + branch_addition**3
    )
    quotient, remainder = sp.div(sp.expand(square**3 - cube**2), delta, x, y)
    zero(remainder, "shifted residual divisibility")
    return sp.expand(quotient)


# Mixed products needed by the seven minimal affine polarization slices.
A10 = -11 * x**2 - 7 * x + y**2                         # e0*h1
A30 = 27 * x**2 * y + 69 * x * y + y**3 + 38 * y       # e0*h3
A_plus_A = (                                                # e+*(-2h1+h2)
    -sp.Rational(81, 4) * x**3
    - sp.Rational(53, 4) * x**2
    + 4 * x * y
    - sp.Rational(7, 4) * y**2
    + 4 * y
)
A_plus_B = (                                                # e+*(2h1+h3)
    -sp.Rational(81, 4) * x**4
    - sp.Rational(9, 4) * x**3
    - sp.Rational(27, 2) * x**2 * y
    + 7 * x**2
    - sp.Rational(57, 4) * x * y**2
    - sp.Rational(31, 2) * x * y
    - sp.Rational(1, 2) * y**3
    - sp.Rational(23, 2) * y**2
)
A_minus_pair = (                                            # e-*(h2+h3)
    sp.Rational(81, 4) * x**4
    + sp.Rational(45, 2) * x**3
    - sp.Rational(27, 2) * x**2 * y
    + sp.Rational(25, 4) * x**2
    + sp.Rational(57, 4) * x * y**2
    - sp.Rational(23, 2) * x * y
    - sp.Rational(1, 2) * y**3
    + sp.Rational(53, 4) * y**2
    + 4 * y
)

mixed_controls = (
    (A10, e0, h1),
    (A30, e0, h3),
    (A_plus_A, e_plus, -2 * h1 + h2),
    (A_plus_B, e_plus, 2 * h1 + h3),
    (A_minus_pair, e_minus, h2 + h3),
)
for index, (descent, selector, normalized_element) in enumerate(mixed_controls):
    zero(
        descent.subs({x: x_t, y: y_t})
        - selector.subs({x: x_t, y: y_t}) * normalized_element,
        f"mixed product descent {index}",
    )

# Branch Z={0}: C2=0 and the only minimal value change is a*e0.
a, c = sp.symbols("a c")
R_h1_value = shifted_residual(1, 0, 0, a * e0, a * A10)
expected_h1_y_zero = sp.expand(
    -3 * a**4 * x**2
    - 6 * a**4 * x
    - 3 * a**4
    + 72 * a**3 * x**2
    + 24 * a**3 * x
    - 16 * a**3
    - 486 * a**2 * x**3
    + 234 * a**2 * x**2
    + 336 * a**2 * x
    - 3888 * a * x**3
    - 2352 * a * x**2
    + 6561 * x**4
    + 3888 * x**3
)
zero(R_h1_value.subs(y, 0) - expected_h1_y_zero,
     "h1 plus e0 specialized quartic")

quad2, quad1, quad0 = sp.symbols("quad2 quad1 quad0")
h1_square_difference = sp.Poly(
    expected_h1_y_zero - (quad2 * x**2 + quad1 * x + quad0) ** 2,
    x,
)
h1_square_groebner = sp.groebner(
    [h1_square_difference.coeff_monomial(x**degree) for degree in range(5)],
    quad2,
    quad1,
    quad0,
    a,
    order="lex",
)
expected_h1_square_groebner = sp.groebner(
    [
        256 * quad2 - 243 * quad0 * a**2 - 1944 * quad0 * a - 5184 * quad0,
        32 * quad1 - 21 * quad0 * a**2 - 168 * quad0 * a - 480 * quad0,
        quad0**2 + 96 * a**2 + 768 * a + 1280,
        (a + 4) ** 3,
    ],
    quad2,
    quad1,
    quad0,
    a,
    order="lex",
)
gate(h1_square_groebner == expected_h1_square_groebner,
     "h1 plus e0 quartic square ideal")
zero(R_h1_value.subs(a, -4) - (9 * x + 4) ** 4,
     "unique h1 affine-slice square residual")

# Add the h3 class coordinate and scale it to one.  If c!=0, the y=0
# quartic reduces homogeneously to the preceding calculation and forces
# a=-4c.  Its remaining cofactor never squares.  The c=0 seam is checked
# separately because the y=0 specialization then vanishes identically.
R_zero_branch = shifted_residual(
    c,
    0,
    1,
    a * e0,
    a * (c * A10 + A30),
)
expected_zero_y_zero = sp.expand(
    c**2
    * (
        -3 * a**4 * x**2
        - 6 * a**4 * x
        - 3 * a**4
        + 72 * a**3 * c * x**2
        + 24 * a**3 * c * x
        - 16 * a**3 * c
        - 486 * a**2 * c**2 * x**3
        + 234 * a**2 * c**2 * x**2
        + 336 * a**2 * c**2 * x
        - 3888 * a * c**3 * x**3
        - 2352 * a * c**3 * x**2
        + 6561 * c**4 * x**4
        + 3888 * c**4 * x**3
    )
)
zero(R_zero_branch.subs(y, 0) - expected_zero_y_zero,
     "full d0-zero branch y-zero quartic")

zero_seam = sp.cancel(R_zero_branch.subs(a, -4 * c) / (c + y) ** 2)
gate(sp.denom(zero_seam) == 1, "d0-zero seam square-factor division")
expected_zero_seam_x0 = sp.expand(
    -16
    * (
        -16 * c**4
        - 64 * c**3 * y
        + 10488 * c**2 * y**2
        + 37044 * c * y**5
        - 127072 * c * y**3
        + 37044 * y**6
        + 445835 * y**4
    )
)
zero(zero_seam.subs(x, 0) - expected_zero_seam_x0,
     "d0-zero seam x-zero sextic")
cube3, cube2, cube1, cube0 = sp.symbols("cube3 cube2 cube1 cube0")
zero_seam_square_difference = sp.Poly(
    expected_zero_seam_x0
    - (cube3 * y**3 + cube2 * y**2 + cube1 * y + cube0) ** 2,
    y,
)
zero_seam_square_groebner = sp.groebner(
    [zero_seam_square_difference.coeff_monomial(y**degree) for degree in range(7)],
    cube3,
    cube2,
    cube1,
    cube0,
    c,
    order="lex",
)
gate(zero_seam_square_groebner == sp.groebner(
    [1], cube3, cube2, cube1, cube0, c, order="lex"
), "d0-zero seam cofactor is never square")

h3_value_core = sp.cancel(R_zero_branch.subs(c, 0) / y**2)
gate(sp.denom(h3_value_core) == 1, "pure h3 value seam y-square division")
expected_h3_value_x0 = sp.expand(
    -3 * a**4
    - 8 * a**3 * y**3
    - 472 * a**3 * y
    - 1008 * a**2 * y**4
    - 27816 * a**2 * y**2
    - 42336 * a * y**5
    - 727776 * a * y**3
    - 592704 * y**6
    - 7133360 * y**4
)
zero(h3_value_core.subs(x, 0) - expected_h3_value_x0,
     "pure h3 value seam x-zero sextic")
h3_value_square_difference = sp.Poly(
    expected_h3_value_x0
    - (cube3 * y**3 + cube2 * y**2 + cube1 * y + cube0) ** 2,
    y,
)
h3_value_square_groebner = sp.groebner(
    [h3_value_square_difference.coeff_monomial(y**degree) for degree in range(7)],
    cube3,
    cube2,
    cube1,
    cube0,
    a,
    order="lex",
)
gate(h3_value_square_groebner == sp.groebner(
    [1], cube3, cube2, cube1, cube0, a, order="lex"
), "pure h3 value seam is never square")

# Branch Z={+}: d_+=0 gives (C1,C2,C3)=(-2+2s,1,s) after
# scaling the nonzero C2.  The square recurrence for the y=0 specialization
# has an empty exact coefficient ideal.  The Z={-} branch follows by
# t|->-t, y|->-y, h2|->-h2.
s, alpha = sp.symbols("s alpha")
R_plus_branch = shifted_residual(
    -2 + 2 * s,
    1,
    s,
    alpha * e_plus,
    alpha * (A_plus_A + s * A_plus_B),
)
plus_y_zero = sp.Poly(R_plus_branch.subs(y, 0), x)
gate(plus_y_zero.degree() == 7, "d-plus-zero specialized degree")
plus_coefficients = [
    plus_y_zero.coeff_monomial(x**degree) for degree in range(8)
]
zero(plus_coefficients[0] + 320000, "d-plus-zero constant coefficient")

# If p(x)=(A3*x^3+...+A0)^2 and u_i=A0*A_i, then u1,u2,u3
# are forced successively by p1,p2,p3.  The remaining four rows are exact
# necessary and sufficient compatibility equations.
p0 = plus_coefficients[0]
u1 = sp.cancel(plus_coefficients[1] / 2)
u2 = sp.cancel((plus_coefficients[2] - u1**2 / p0) / 2)
u3 = sp.cancel((plus_coefficients[3] - 2 * u1 * u2 / p0) / 2)
plus_square_conditions = [
    plus_coefficients[7],
    plus_coefficients[4] - (u2**2 + 2 * u1 * u3) / p0,
    plus_coefficients[5] - 2 * u2 * u3 / p0,
    plus_coefficients[6] - u3**2 / p0,
]
plus_square_numerators = [
    sp.together(condition).as_numer_denom()[0]
    for condition in plus_square_conditions
]
plus_square_groebner = sp.groebner(
    plus_square_numerators, s, alpha, order="lex"
)
gate(plus_square_groebner == sp.groebner([1], s, alpha, order="lex"),
     "d-plus-zero minimal slice has no square residual")
zero(delta.subs(y, -y) - delta, "normalization involution preserves branch curve")
zero(e_plus.subs(y, -y) - e_minus, "normalization involution swaps selectors")

# Branch Z={0,+}: the nonzero defect class is proportional to 2h1+h3.
# A nonzero e+ value creates an unavoidable odd x-degree seven.  If that
# value vanishes, the packet reduces to the already closed Z={0} slice.
value0, value_plus = sp.symbols("value0 value_plus")
R_zero_plus_branch = shifted_residual(
    2,
    0,
    1,
    value0 * e0 + value_plus * e_plus,
    value0 * (2 * A10 + A30) + value_plus * A_plus_B,
)
zero_plus_y_zero = sp.Poly(R_zero_plus_branch.subs(y, 0), x)
gate(zero_plus_y_zero.degree() == 7, "zero-plus pair specialized degree")
zero(
    zero_plus_y_zero.coeff_monomial(x**7)
    + sp.Rational(6561, 8) * value_plus**3,
    "zero-plus pair odd leading coefficient",
)
zero(
    R_zero_plus_branch.subs(value_plus, 0)
    - R_zero_branch.subs({c: 2, a: value0}),
    "zero-plus pair value-zero boundary",
)

# Branch Z={+,-}: the class is proportional to h2+h3.  Two successive
# odd-degree tests leave a single exceptional value pair; its x=0 quartic
# in q=y^2 has an empty square coefficient ideal.
value_minus = sp.symbols("value_minus")
A_plus_pair = sp.expand(A_plus_A + A_plus_B)
R_plus_minus_branch = shifted_residual(
    0,
    1,
    1,
    value_plus * e_plus + value_minus * e_minus,
    value_plus * A_plus_pair + value_minus * A_minus_pair,
)
plus_minus_y_zero = sp.Poly(R_plus_minus_branch.subs(y, 0), x)
gate(plus_minus_y_zero.degree() == 7, "plus-minus pair specialized degree")
zero(
    plus_minus_y_zero.coeff_monomial(x**7)
    + sp.Rational(6561, 8) * (value_plus - value_minus - 152) ** 3,
    "plus-minus first odd-degree row",
)
zero(
    plus_minus_y_zero.coeff_monomial(x**6)
    + sp.Rational(729, 2)
    * (value_plus - value_minus - 152) ** 2
    * (5 * value_plus - 5 * value_minus - 652),
    "plus-minus second leading row",
)
zero(
    plus_minus_y_zero.coeff_monomial(x**5).subs(
        value_minus, value_plus - 152
    )
    + sp.Rational(243 * 277248, 64) * (value_plus - 76) ** 2,
    "plus-minus second odd-degree seam",
)
exceptional_plus_minus = sp.expand(
    R_plus_minus_branch.subs({value_plus: 76, value_minus: -76})
)
zero(
    exceptional_plus_minus.subs(y, 0) + 512 * (9 * x + 5) ** 4,
    "plus-minus exceptional y-zero square control",
)
q = sp.symbols("q")
exceptional_x_zero_q = sp.Poly(exceptional_plus_minus.subs(x, 0), y)
gate(
    all(exponent[0] % 2 == 0 for exponent, _ in exceptional_x_zero_q.terms()),
    "plus-minus exceptional x-zero specialization is y-even",
)
exceptional_q = sp.expand(
    sum(
        coefficient * q ** (exponent[0] // 2)
        for exponent, coefficient in exceptional_x_zero_q.terms()
    )
)
expected_exceptional_q = (
    -592704 * q**4
    + 4946089 * q**3
    - 11376906 * q**2
    + 210000 * q
    - 320000
)
zero(exceptional_q - expected_exceptional_q,
     "plus-minus exceptional quartic")
q2, q1, q0 = sp.symbols("q2 q1 q0")
exceptional_square_difference = sp.Poly(
    exceptional_q - (q2 * q**2 + q1 * q + q0) ** 2, q
)
exceptional_square_groebner = sp.groebner(
    [exceptional_square_difference.coeff_monomial(q**degree) for degree in range(5)],
    q2,
    q1,
    q0,
    order="lex",
)
gate(exceptional_square_groebner == sp.groebner(
    [1], q2, q1, q0, order="lex"
), "plus-minus exceptional residual is not square")

# The three-derivative-zero stratum has trivial defect class because the
# derivative matrix is invertible.  It lies in R and gives only P=H^2,
# Q=H^3, hence the zero residual rather than a cubic extension.
gate(derivative_matrix.nullspace() == [], "all-zero derivative class is trivial")

# First omitted J-adic layer.  These three polynomials generate the exact
# ideal of the three cusp image points.  Along each individual generator ray
# through the unique positive representative h*=h1-4e0, no new square occurs.
j1 = x * (x + 1)
j2 = y * (x + 1)
j3 = y**2 + 4 * x
j_groebner = sp.groebner([j1, j2, j3], y, x, order="lex")
expected_j_groebner = sp.groebner(
    [4 * x + y**2, x * y + y, x**2 + x], y, x, order="lex"
)
gate(j_groebner == expected_j_groebner, "exact cusp-value-zero ideal generators")
zero(j_groebner.reduce(delta)[1], "branch equation belongs to cusp ideal")

h_star = h1 - 4 * e0.subs({x: x_t, y: y_t})
P_star = (x + 1) * (9 * x + 4) ** 2
Q_star = (
    y**2 - sp.Rational(1, 2) * (30 * x**2 + 30 * x + 8)
) * (9 * x + 4) ** 2
zero(P_star.subs({x: x_t, y: y_t}) - h_star**2,
     "positive representative square descent")
zero(Q_star.subs({x: x_t, y: y_t}) - h_star**3,
     "positive representative cube descent")
zero(P_star**3 - Q_star**2 - delta * (9 * x + 4) ** 4,
     "positive representative square residual")

A_j1 = -15 * x**3 - 15 * x**2 + x * y**2 - 4 * x
A_j2 = -15 * x**2 * y - 15 * x * y + y**3 - 4 * y
A_j3 = (
    81 * x**4
    + 9 * x**3
    - 44 * x**2
    + 15 * x * y**2
    - 16 * x
    + 4 * y**2
)
for index, (generator, mixed) in enumerate(
    ((j1, A_j1), (j2, A_j2), (j3, A_j3)), start=1
):
    zero(
        mixed.subs({x: x_t, y: y_t})
        - generator.subs({x: x_t, y: y_t}) * h_star,
        f"J-generator mixed descent {index}",
    )


def j_ray_residual(generator: sp.Expr, mixed: sp.Expr, parameter: sp.Symbol) -> sp.Expr:
    square = sp.expand(P_star + 2 * parameter * mixed + parameter**2 * generator**2)
    cube = sp.expand(
        Q_star
        + 3 * parameter * generator * P_star
        + 3 * parameter**2 * generator * mixed
        + parameter**3 * generator**3
    )
    quotient, remainder = sp.div(sp.expand(square**3 - cube**2), delta, x, y)
    zero(remainder, "J-generator ray residual divisibility")
    return sp.expand(quotient)


lam = sp.symbols("lam")
R_j1 = j_ray_residual(j1, A_j1, lam)
R_j1_y = sp.Poly(R_j1, y)
gate(R_j1_y.degree() == 2, "j1 ray y degree two")
zero(R_j1_y.coeff_monomial(y), "j1 ray has no linear y term")
zero(R_j1_y.coeff_monomial(y**2) + 8 * lam**3 * x**3,
     "j1 ray unavoidable y-square coefficient")
zero(R_j1.subs({x: 0, y: 0}) - 256, "j1 ray nonzero constant part")

R_j2 = j_ray_residual(j2, A_j2, lam)
R_j2_x_zero = sp.Poly(R_j2.subs(x, 0), y)
gate(R_j2_x_zero.degree() == 5, "j2 ray x-zero degree five")
zero(R_j2_x_zero.coeff_monomial(y**5) + 8 * lam**3,
     "j2 ray odd leading coefficient")

R_j3 = j_ray_residual(j3, A_j3, lam)
R_j3_y_zero = sp.Poly(R_j3.subs(y, 0), x)
gate(R_j3_y_zero.degree() == 7, "j3 ray y-zero degree seven")
zero(R_j3_y_zero.coeff_monomial(x**7) - 52488 * lam**3,
     "j3 ray odd leading coefficient")

semantic = {
    "polarization": "H in S; H2,H3 in R iff H(a)H'(a)=0 at 0,+1,-1",
    "selectors": "e0=x+1;e+=(-2x-y)/4;e-=(-2x+y)/4;evaluation identity",
    "fibers": "fixed derivative class; values free exactly at zero derivative coordinates;mod J",
    "branches": "eight linear polarizations;seven noncanonical minimal affine slices",
    "square_search": "only nontrivial square is h1-4e0;residual=(9x+4)^4",
    "J_rays": "x(x+1),y(x+1),y2+4x generate J;each individual ray has no new square",
    "ramification": "THM3869 leaves field discriminant Delta*(9x+4)^2",
    "scope": "arbitrary combinations/higher-degree additions from J remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3872-three-cusp-polarization-branches-and-minimal-affine-square-residual-gate")
print("polarization_branches=8;criterion=H(a)*H'(a)=0")
print("defect_derivative_rank=3;cusp_value_selector_rank=3")
print("minimal_noncanonical_branches=7")
print("minimal_slice_nontrivial_square_residuals=1")
print("unique_square_representative=h1-4*(x+1);residual=(9x+4)^4")
print("cusp_value_zero_generator_rays=3;new_square_residuals=0")
print("extra_cardano_ramification_removed=NO")
print("cusp_value_zero_ideal_combinations_and_higher_representatives=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
