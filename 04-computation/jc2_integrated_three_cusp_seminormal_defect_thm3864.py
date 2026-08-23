#!/usr/bin/env python3
"""Exact companion for THM-3864's three-direction seminormal defect."""

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


t, u, x, y = sp.symbols("t u x y")
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

# The exact conductor divisor: order two at each cusp address and order one
# at every normalization address over the three nodes.
node_opposite = 3 * t**2 - 5
node_quartic = 9 * t**4 - 18 * t**2 + 4
conductor = sp.expand(t**2 * (t**2 - 1) ** 2 * node_opposite * node_quartic)
conductor_in_ring = 27 * x**3 + 33 * x**2 + 10 * x + y**2
zero(
    conductor_in_ring.subs({x: x_t, y: y_t}) - conductor,
    "conductor generator belongs to branch ring",
)
gate(sp.degree(conductor, t) == 12, "conductor degree")
gate(
    sp.gcd(t * (t**2 - 1), node_opposite * node_quartic) == 1,
    "cusp and node addresses are disjoint",
)
gate(sp.gcd(node_opposite, node_quartic) == 1, "node packets are disjoint")
gate(sp.discriminant(node_opposite, t) != 0, "opposite node addresses simple")
gate(sp.discriminant(node_quartic, t) != 0, "quartic node addresses simple")


def node_remainders(polynomial: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    """Polynomial conditions for equality on all three node pairs."""

    opposite = sp.rem(
        sp.Poly(sp.expand(polynomial - polynomial.subs(t, -t)), t),
        sp.Poly(node_opposite, t),
    ).as_expr()
    partner = -sp.Rational(2, 3) / t
    cleared = sp.together(polynomial - polynomial.subs(t, partner)).as_numer_denom()[0]
    quartic = sp.rem(sp.Poly(cleared, t), sp.Poly(node_quartic, t)).as_expr()
    return sp.expand(opposite), sp.expand(quartic)


# In k[t]/(conductor), the node equalities have rank three.  Adding the
# three cusp-derivative rows has rank six.  Thus the seminormal quotient has
# dimension three over the branch ring.
coefficients = sp.symbols("c0:12")
generic = sum(coefficients[index] * t**index for index in range(12))
node_conditions = []
for remainder in node_remainders(generic):
    polynomial = sp.Poly(remainder, t)
    node_conditions.extend(
        polynomial.coeff_monomial(t**degree)
        for degree in range(max(polynomial.degree(), 0) + 1)
    )
node_matrix, _ = sp.linear_eq_to_matrix(node_conditions, coefficients)
gate(node_matrix.rank() == 3, "three independent node-gluing conditions")
cusp_conditions = [sp.diff(generic, t).subs(t, address) for address in (0, 1, -1)]
full_matrix, _ = sp.linear_eq_to_matrix(
    [*node_conditions, *cusp_conditions], coefficients
)
gate(full_matrix.rank() == 6, "node plus cusp conditions have rank six")
gate(12 - node_matrix.rank() == 9, "seminormal conductor quotient dimension")
gate(12 - full_matrix.rank() == 6, "branch conductor quotient dimension")
gate((12 - node_matrix.rank()) - (12 - full_matrix.rank()) == 3,
     "seminormal defect dimension three")

# Three square/cube-descent representatives.  The third is the missing
# direction beyond THM-3854's even and odd controls.
h_even = t**2 * (t**2 - 1) * (9 * t**2 - 14)
h_odd = t * (t**2 - 1) * (9 * t**2 - 4) * (3 * t**2 - 5)
h_third = t**3 * (t**2 - 1) * (3 * t**2 - 5) * (
    9 * t**4 + 24 * t**2 - 38
)
h_rows = (h_even, h_odd, h_third)

for index, h in enumerate(h_rows, start=1):
    opposite_remainder, quartic_remainder = node_remainders(h)
    zero(opposite_remainder, f"h{index} opposite-node equality")
    zero(quartic_remainder, f"h{index} quartic-node equalities")
    gate([h.subs(t, address) for address in (0, 1, -1)] == [0, 0, 0],
         f"h{index} vanishes at every cusp")

derivative_matrix = sp.Matrix(
    [[sp.diff(h, t).subs(t, address) for h in h_rows] for address in (0, 1, -1)]
)
gate(
    derivative_matrix
    == sp.Matrix([[0, -20, 0], [-10, -20, 20], [10, -20, 20]]),
    "exact cusp-derivative coordinate matrix",
)
gate(derivative_matrix.det() == -8000, "derivative coordinates span defect")

# Explicit square and cube descents for the third direction.
P_third = sp.expand(
    81 * x**3 * y**2
    + 657 * x**2 * y**2
    + 2356 * x * y**2
    + 84 * y**4
    + 1444 * y**2
)
Q_third = sp.expand(
    8991 * x**4 * y**3
    + 21897 * x**3 * y**3
    + 81 * x**2 * y**5
    + 92226 * x**2 * y**3
    + 5364 * x * y**5
    + 134292 * x * y**3
    + 5308 * y**5
    + 54872 * y**3
)
zero(P_third.subs({x: x_t, y: y_t}) - h_third**2,
     "third direction square descent")
zero(Q_third.subs({x: x_t, y: y_t}) - h_third**3,
     "third direction cube descent")

residual_core = sp.expand(
    6561 * x**4
    - 845640 * x**3
    - 2056104 * x**2
    - 10708704 * x
    - 592704 * y**2
    - 7133360
)
zero(
    P_third**3 - Q_third**2 - delta * y**6 * residual_core,
    "third direction residual quotient",
)
gate(sp.Poly(residual_core, y).degree() == 2, "residual core y degree")
gate(sp.Poly(residual_core, y).coeff_monomial(y) == 0,
     "residual core has no odd y term")
gate(sp.Poly(residual_core, y).coeff_monomial(y**2) == -592704,
     "residual core has nonzero y-square coefficient")
gate(residual_core.subs(y, 0) != 0, "residual core has nonzero constant part")

# Complete minimum-total-degree, parity-preserving lift family for h_third.
# P may change by a*delta; Q may change by delta*y*(b0+b1*x).
a, b0, b1, q = sp.symbols("a b0 b1 q")
P_mixed = P_third + a * delta
Q_mixed = Q_third + delta * y * (b0 + b1 * x)
mixed_numerator = sp.expand(P_mixed**3 - Q_mixed**2)
mixed_residual, mixed_remainder = sp.div(mixed_numerator, delta, x, y)
zero(mixed_remainder, "mixed family exact divisibility")
mixed_y_polynomial = sp.Poly(mixed_residual, y)
gate(
    all(exponent[0] % 2 == 0 for exponent, _ in mixed_y_polynomial.terms()),
    "mixed residual is y-even",
)
mixed_q = sp.Poly(
    sp.expand(
        sum(
            coefficient * q ** (exponent[0] // 2)
            for exponent, coefficient in mixed_y_polynomial.terms()
        )
    ),
    q,
)
gate(mixed_q.degree() == 4, "generic mixed q degree four")
zero(mixed_q.coeff_monomial(q**4) - (a - 84) ** 3,
     "mixed leading q coefficient")
gate(
    sp.Poly(mixed_q.coeff_monomial(q**3), x).coeff_monomial(x**4) == 6561,
    "unavoidable x4 q3 coefficient",
)
zero(
    mixed_q.coeff_monomial(q**0) - a**3 * x**6 * (9 * x + 5) ** 4,
    "mixed q-zero boundary",
)

# The complete canonical projective line spanned by h_even and h_third.
# Its residual has a visible square factor, but the quotient never squares.
L, N = sp.symbols("L N")
P_even = 81 * x**3 + 49 * x**2 + 8 * y**2
P_even_third = 81 * x**3 * y - 369 * x**2 * y - 266 * x * y + 46 * y**3
Q_even = -243 * x**4 - 143 * x**3 + 81 * x**2 * y**2 + 120 * x * y**2 + 64 * y**2
Q_even_even_third = (
    2835 * x**4 * y
    + 4797 * x**3 * y
    + 81 * x**2 * y**3
    + 1862 * x**2 * y
    + 424 * x * y**3
    + 368 * y**3
)
Q_even_third_third = (
    5913 * x**4 * y**2
    - 6147 * x**3 * y**2
    + 81 * x**2 * y**4
    - 22268 * x**2 * y**2
    + 2172 * x * y**4
    - 10108 * x * y**2
    + 2116 * y**4
)
P_line = sp.expand(L**2 * P_even + 2 * L * N * P_even_third + N**2 * P_third)
Q_line = sp.expand(
    L**3 * Q_even
    + 3 * L**2 * N * Q_even_even_third
    + 3 * L * N**2 * Q_even_third_third
    + N**3 * Q_third
)
h_line = L * h_even + N * h_third
zero(P_line.subs({x: x_t, y: y_t}) - h_line**2,
     "projective-line square descent")
zero(Q_line.subs({x: x_t, y: y_t}) - h_line**3,
     "projective-line cube descent")
line_residual, line_remainder = sp.div(sp.expand(P_line**3 - Q_line**2), delta, x, y)
zero(line_remainder, "projective-line residual divisibility")
line_core = sp.cancel(line_residual / (L + N * y) ** 2)
gate(sp.denom(line_core) == 1, "projective-line exact square-factor division")
zero(
    line_residual - (L + N * y) ** 2 * line_core,
    "projective-line symbolic residual factorization",
)
zero(
    line_core.subs(y, 0) - 243 * L**4 * x**3 * (27 * x + 16),
    "projective-line nonzero-L nonsquare specialization",
)
zero(
    line_core.subs(L, 0) - N**4 * y**4 * residual_core,
    "projective-line zero-L third-direction boundary",
)

# The second canonical projective line, spanned by h_even and h_odd.
U, V, zeta = sp.symbols("U V zeta")
P_even_odd = 81 * x**2 * y + 137 * x * y + 56 * y
P_odd = -648 * x**3 - 720 * x**2 + 81 * x * y**2 - 200 * x + 49 * y**2
Q_even_even_odd = (
    -891 * x**3 * y
    - 1183 * x**2 * y
    + 81 * x * y**3
    - 392 * x * y
    + 56 * y**3
)
Q_even_odd_odd = (
    4536 * x**4
    + 5040 * x**3
    - 1539 * x**2 * y**2
    + 1400 * x**2
    - 1247 * x * y**2
    + 81 * y**4
    - 256 * y**2
)
Q_odd = (
    6561 * x**4 * y
    + 17010 * x**3 * y
    + 18009 * x**2 * y
    + 243 * x * y**3
    + 8760 * x * y
    + 143 * y**3
    + 1600 * y
)
P_even_odd_line = sp.expand(
    U**2 * P_even + 2 * U * V * P_even_odd + V**2 * P_odd
)
Q_even_odd_line = sp.expand(
    U**3 * Q_even
    + 3 * U**2 * V * Q_even_even_odd
    + 3 * U * V**2 * Q_even_odd_odd
    + V**3 * Q_odd
)
h_even_odd_line = U * h_even + V * h_odd
zero(P_even_odd_line.subs({x: x_t, y: y_t}) - h_even_odd_line**2,
     "even-odd line square descent")
zero(Q_even_odd_line.subs({x: x_t, y: y_t}) - h_even_odd_line**3,
     "even-odd line cube descent")
even_odd_residual, even_odd_remainder = sp.div(
    sp.expand(P_even_odd_line**3 - Q_even_odd_line**2), delta, x, y
)
zero(even_odd_remainder, "even-odd line residual divisibility")

# On V!=0, homogeneity lets us set V=1 and zeta=(U/V)^2.  The y=0
# specialization is a quartic G_zeta(x).
even_odd_y_zero = sp.Poly(even_odd_residual.subs({y: 0, V: 1}), U)
gate(all(exponent[0] % 2 == 0 for exponent, _ in even_odd_y_zero.terms()),
     "even-odd y-zero specialization uses U squared")
quartic_zeta = sp.expand(
    sum(
        coefficient * zeta ** (exponent[0] // 2)
        for exponent, coefficient in even_odd_y_zero.terms()
    )
)
expected_quartic_zeta = sp.expand(
    (6561 * zeta**3 - 157464 * zeta**2 + 1259712 * zeta - 3359232) * x**4
    + (3888 * zeta**3 - 108864 * zeta**2 - 124416 * zeta - 7464960) * x**3
    + (-9576 * zeta**2 - 1304640 * zeta - 6220800) * x**2
    + (-470400 * zeta - 2304000) * x
    - 320000
)
zero(quartic_zeta - expected_quartic_zeta,
     "even-odd exact specialized quartic")

square_a, square_b, square_c = sp.symbols("square_a square_b square_c")
square_difference = sp.Poly(
    quartic_zeta - (square_a * x**2 + square_b * x + square_c) ** 2,
    x,
)
square_equations = [
    square_difference.coeff_monomial(x**degree) for degree in range(5)
]
square_groebner = sp.groebner(
    square_equations, square_a, square_b, square_c, zeta, order="lex"
)
expected_square_groebner = sp.groebner(
    [
        400 * square_a + 243 * square_c * zeta - 1296 * square_c,
        200 * square_b - 147 * square_c * zeta - 720 * square_c,
        square_c**2 + 320000,
        zeta**2,
    ],
    square_a,
    square_b,
    square_c,
    zeta,
    order="lex",
)
gate(square_groebner == expected_square_groebner,
     "even-odd quartic square coefficient ideal")

odd_residual_core = (
    41472 * x**2
    + 6561 * x * y**2
    + 46080 * x
    + 3888 * y**2
    + 12800
)
zero(
    even_odd_residual.subs(U, 0)
    + V**6 * (9 * x + 5) ** 2 * odd_residual_core,
    "even-odd h2 endpoint residual",
)
zero(
    even_odd_residual.subs(V, 0)
    - U**6 * (6561 * x**4 + 3888 * x**3 - 512 * y**2),
    "even-odd h1 endpoint residual",
)

semantic = {
    "conductor": "t2(t2-1)2(3t2-5)(9t4-18t2+4);exact",
    "seminormalization": "three node equalities;branch adds three cusp derivatives",
    "defect": "dimension 3;derivative determinant -8000",
    "third_direction": "t3(t2-1)(3t2-5)(9t4+24t2-38);square/cube descent",
    "bounded_no_go": "minimum-degree parity-preserving delta lifts have no square residual",
    "projective_line": "all canonical L*h1+N*h3 residuals nonsquare",
    "second_line": "all canonical U*h1+V*h2 residuals nonsquare",
    "scope": "general mixed/higher lifts and h2-h3/full-plane combinations open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3864-integrated-three-cusp-conductor-seminormal-three-direction-gate")
print("conductor_degree=12;conductor=EXACT")
print("node_conditions=3;cusp_derivative_conditions=3;seminormal_defect_dimension=3")
print("derivative_coordinate_determinant=-8000;third_direction=EXPLICIT")
print("third_square_cube_descent=YES;canonical_residual_square=NO")
print("minimum_degree_parity_preserving_mixed_lift_square=NONE")
print("canonical_projective_line_h1_h3_square_residual=NONE")
print("canonical_projective_line_h1_h2_square_residual=NONE")
print("general_mixed_higher_and_h2_h3_full_plane_lifts=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
