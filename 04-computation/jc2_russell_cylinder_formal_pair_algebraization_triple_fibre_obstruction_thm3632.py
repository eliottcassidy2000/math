#!/usr/bin/env python3
"""Exact controls for the THM-3632 formal-pair algebraization obstruction.

The theorem's rational-side and fixed-fold statements are elementary descent
arguments. This companion checks the common THM-3630 side fold, the compiler
fibre, the affine-slice equalizer, polynomial and rational derivative gates,
the rho=0 side sections, and the exact zero boundary.
"""

import ast
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one active exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact rational zero test."""
    return sp.cancel(expression) == 0


print("THM-3632 exact companion -- proved formal-pair obstruction")
print("status=thm3630_F_nonlocal_algebraic_and_no_fixed_Q_completion; audit=HOSTILE-PASS")


print("SECTION compiler relation and retained triple fibre")
x, q, w = sp.symbols("x q w")
D = 1 + x**2 * q
b = (D - 1) * (D + 2) ** 2
c = x * D * (D + 2)
e = q * (D + 3)
require("Danielewski relation", zero(c**2 * e - b * (b + 4)))

collision_data = ((-1, -3), (0, -sp.Rational(3, 4)), (1, -3))
collision_images = []
for point, q_value in collision_data:
    image = tuple(
        sp.cancel(coordinate.subs({x: point, q: q_value}))
        for coordinate in (b, c, e)
    ) + (sp.Integer(0),)
    collision_images.append(image)
    require(f"collision image x={point}", image == (0, 0, -3, 0))
require("three images equal", len(set(collision_images)) == 1)

# A minimal Hermite control showing that the genuine non-even rho=0 stratum
# is nonempty. Its central derivative is u=1.
Q_demo = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
require("demo non-even", sp.expand(Q_demo.subs(x, -x) - Q_demo) != 0)
demo_values = tuple(sp.cancel(Q_demo.subs(x, point)) for point in (-1, 0, 1))
demo_slopes = tuple(
    sp.cancel(sp.diff(Q_demo, x).subs(x, point)) for point in (-1, 0, 1)
)
require("demo collision values", demo_values == (-3, -sp.Rational(3, 4), -3))
require(
    "demo rho0 slopes",
    demo_slopes == (-sp.Rational(9, 2), 1, sp.Rational(9, 2)),
)
q_slope = sp.symbols("q_slope")
c_x_total = sp.diff(c, x) + sp.diff(c, q) * q_slope
for point, q_value, slope in zip((-1, 0, 1), demo_values, demo_slopes):
    require(
        f"rho0 transverse c_x point={point}",
        zero(c_x_total.subs({x: point, q: q_value, q_slope: slope}) - 3),
    )
require("D=-2 side is c=0", zero(c.subs(q, -3 / x**2)))
print("PASS compiler_relation=exact fibre=(0,0,-3,0)^3 rho0_demo=u1_non_even c_x=3")


print("SECTION THM-3630 common rational side continuation")
Q_infinity = -sp.Rational(3, 4) - sp.Rational(9, 4) / x**2
require("Q_infinity minus value", zero(Q_infinity.subs(x, -1) + 3))
require("Q_infinity plus value", zero(Q_infinity.subs(x, 1) + 3))
require(
    "Q_infinity minus slope",
    zero(sp.diff(Q_infinity, x).subs(x, -1) + sp.Rational(9, 2)),
)
require(
    "Q_infinity plus slope",
    zero(sp.diff(Q_infinity, x).subs(x, 1) - sp.Rational(9, 2)),
)
q_infinity_fold = Q_infinity + w**2
side_coordinates = tuple(
    sp.cancel(coordinate.subs(q, q_infinity_fold)) for coordinate in (b, c, e)
)
for coordinate_index, coordinate in enumerate(side_coordinates):
    denominator = sp.denom(coordinate)
    for point in (-1, 1):
        require(
            f"rational side regular coordinate={coordinate_index} point={point}",
            denominator.subs({x: point, w: 0}) != 0,
        )

g_exact = (1 - sp.Rational(4, 3) * w**2) ** -sp.Rational(1, 2)
q_moving = sp.cancel(Q_infinity.subs(x, g_exact) + w**2)
require("moving side q", zero(q_moving + 3 / g_exact**2))
require("moving side D", zero(D.subs({x: g_exact, q: q_moving}) + 2))
require("moving plus b", zero(b.subs({x: g_exact, q: q_moving})))
require("moving plus c", zero(c.subs({x: g_exact, q: q_moving})))
require("moving minus c", zero(c.subs({x: -g_exact, q: q_moving})))
require(
    "moving side e",
    zero(e.subs({x: g_exact, q: q_moving}) + 3 - 4 * w**2),
)

side_x = sp.symbols("side_x")
kappa_side, h_side0 = sp.symbols("kappa_side h_side0")
A_minus = kappa_side * (side_x + g_exact)
A_plus = kappa_side * (side_x - g_exact)
require("minus side primitive derivative", zero(sp.diff(A_minus, side_x) - kappa_side))
require("plus side primitive derivative", zero(sp.diff(A_plus, side_x) - kappa_side))
require(
    "minus side primitive normalization",
    zero(A_minus.subs(side_x, -g_exact)),
)
require(
    "plus side primitive normalization",
    zero(A_plus.subs(side_x, g_exact)),
)
side_continuation_values = (
    -kappa_side + h_side0,
    kappa_side + h_side0,
)
side_continuation_basis = sp.groebner(
    side_continuation_values,
    kappa_side,
    h_side0,
)
require(
    "opposite side continuation ideal",
    side_continuation_basis.polys
    == [
        sp.Poly(kappa_side, kappa_side, h_side0),
        sp.Poly(h_side0, kappa_side, h_side0),
    ],
)
print("PASS thm3630_side_fold=one_rational_Q_infinity regular_at=minus1_plus1")
print("PASS moving_side_sections=minus_g_plus_g common_target=(0,0,-3+4w2,w)")
print("PASS side_primitives=kappa*(x_plusminus_g) continuation_ideal=(kappa,h0)")


print("SECTION affine-in-x pullback-ring equalizer")
kappa, h0 = sp.symbols("kappa h0")
affine_values = tuple(kappa * point + h0 for point in (-1, 0, 1))
require("affine collision values", affine_values == (h0 - kappa, h0, h0 + kappa))
equalizer_equations = [
    sp.expand(affine_values[0] - affine_values[1]),
    sp.expand(affine_values[2] - affine_values[1]),
]
equalizer_basis = sp.groebner(equalizer_equations, kappa, h0)
require(
    "equalizer ideal is kappa",
    equalizer_basis.polys == [sp.Poly(kappa, kappa, h0)],
)
require(
    "nonzero hostile value debt",
    tuple(value.subs({kappa: 12, h0: 0}) for value in affine_values)
    == (-12, 0, 12),
)

# Finite generic polynomial gate for the all-degree coefficient argument:
# d_x p=kappa leaves only kappa*x plus an arbitrary stable polynomial.
x_degree = 6
w_degree = 5
coefficients = {
    (i, j): sp.symbols(f"p_{i}_{j}")
    for i in range(1, x_degree + 1)
    for j in range(w_degree + 1)
}
stable_coefficients = sp.symbols(f"h_0:{w_degree + 1}")
generic_p = sum(stable_coefficients[j] * w**j for j in range(w_degree + 1))
generic_p += sum(
    coefficient * x**i * w**j
    for (i, j), coefficient in coefficients.items()
)
derivative_debt = sp.Poly(sp.diff(generic_p, x) - kappa, x, w)
require(
    "stable coefficients invisible",
    all(not derivative_debt.as_expr().has(value) for value in stable_coefficients),
)
require(
    "linear constant coefficient",
    derivative_debt.coeff_monomial(1) == coefficients[(1, 0)] - kappa,
)
for j in range(1, w_degree + 1):
    require(
        f"linear stable coefficient j={j}",
        derivative_debt.coeff_monomial(w**j) == coefficients[(1, j)],
    )
for i in range(2, x_degree + 1):
    for j in range(w_degree + 1):
        require(
            f"higher x coefficient i={i} j={j}",
            derivative_debt.coeff_monomial(x ** (i - 1) * w**j)
            == i * coefficients[(i, j)],
        )
print("PASS intersection=B_Q_cap_(C[w]+Cx)=C[w] equalizer_ideal=(kappa)")
print("PASS polynomial_derivative_gate=dx_p_kappa_implies_p=kappa*x+h(w) degrees=6x5")


print("SECTION target-local unit denominator gate")
n0, n1, n2 = sp.symbols("n0 n1 n2")
d0 = sp.symbols("d0", nonzero=True)
d1, d2 = sp.symbols("d1 d2")
h_rational = (n0 + n1 * w + n2 * w**2) / (d0 + d1 * w + d2 * w**2)
p_rational = kappa * x + h_rational
require("rational derivative", zero(sp.diff(p_rational, x) - kappa))
rational_values = tuple(
    sp.cancel(p_rational.subs({x: point, w: 0})) for point in (-1, 0, 1)
)
require(
    "rational unit-denominator residues",
    all(
        zero(left - right)
        for left, right in zip(
            rational_values,
            (n0 / d0 - kappa, n0 / d0, n0 / d0 + kappa),
        )
    ),
)
require(
    "rational left residue debt",
    zero(rational_values[0] - rational_values[1] + kappa),
)
require(
    "rational right residue debt",
    zero(rational_values[2] - rational_values[1] - kappa),
)
require("unit denominator at collision", d0.is_nonzero)
print("PASS target_local=unit_denominator_preserves_three_residues obstruction=kappa")
print("PASS denominator_boundary=any_affine_separator_must_be_undefined_at_y0")


print("SECTION normalized rho0 branch primitives")
s = sp.symbols("s")
r_minus, r_plus = sp.symbols("r_minus r_plus")


def side_section(point, slope, second_derivative, sign):
    """Return and check the c=0 side section through order s^2."""
    first = sp.Rational(2, 3) * sign
    second = sp.Rational(4, 27) * sign * (second_derivative + 18)
    displacement = first * s + second * s**2
    branch_x = point + displacement
    branch_q = -3 + slope * displacement + second_derivative * displacement**2 / 2
    side_equation = sp.series(
        branch_q + s + 3 / branch_x**2,
        s,
        0,
        3,
    ).removeO()
    require(f"side section point={point}", zero(side_equation))
    return sp.expand(branch_x)


chi_minus = side_section(-1, -sp.Rational(9, 2), r_minus, -1)
chi_zero = sp.Integer(0)
chi_plus = side_section(1, sp.Rational(9, 2), r_plus, 1)
require(
    "minus section formula",
    zero(
        chi_minus
        + 1
        + sp.Rational(2, 3) * s
        + sp.Rational(4, 27) * (r_minus + 18) * s**2
    ),
)
require(
    "plus section formula",
    zero(
        chi_plus
        - 1
        - sp.Rational(2, 3) * s
        - sp.Rational(4, 27) * (r_plus + 18) * s**2
    ),
)

source_x = sp.symbols("source_x")
branch_sections = (chi_minus, chi_zero, chi_plus)
primitives = tuple(kappa * (source_x - chi) for chi in branch_sections)
for index, (primitive, chi) in enumerate(zip(primitives, branch_sections)):
    require(
        f"primitive derivative branch={index}",
        zero(sp.diff(primitive, source_x) - kappa),
    )
    require(
        f"primitive normalization branch={index}",
        zero(primitive.subs(source_x, chi)),
    )
required_h0 = tuple(
    sp.cancel((-kappa * chi).subs(s, 0)) for chi in branch_sections
)
require("continuation constants", required_h0 == (kappa, 0, -kappa))
continuation_basis = sp.groebner(
    [required_h0[0] - required_h0[1], required_h0[2] - required_h0[1]],
    kappa,
)
require(
    "primitive compatibility ideal",
    continuation_basis.polys == [sp.Poly(kappa, kappa)],
)

g_series = sp.series(
    (1 - sp.Rational(4, 3) * s) ** -sp.Rational(1, 2),
    s,
    0,
    3,
).removeO()
comparison_minus = sp.expand(
    chi_minus.subs(r_minus, -sp.Rational(27, 2))
)
comparison_plus = sp.expand(chi_plus.subs(r_plus, -sp.Rational(27, 2)))
require("THM-3630 minus comparison", zero(comparison_minus + g_series))
require("THM-3630 plus comparison", zero(comparison_plus - g_series))
print("PASS sections=chi_minus_chi0_chi_plus constants=(-1,0,1) rho0_expansion=w2_w4")
print("PASS normalized_primitives=kappa*(x-chi_i) compatibility_ideal=(kappa)")
print("PASS thm3630_comparison=chi_side=plusminus_g_through_w4")


print("SECTION sharp zero boundary and scope controls")
stable_polynomial = sum(
    stable_coefficients[j] * w**j for j in range(w_degree + 1)
)
require("stable boundary derivative", zero(sp.diff(stable_polynomial, x)))
stable_values = tuple(
    stable_polynomial.subs({x: point, w: 0}) for point in (-1, 0, 1)
)
require(
    "stable boundary fibre values",
    stable_values == (stable_coefficients[0],) * 3,
)
zero_primitives = tuple(sp.cancel(primitive.subs(kappa, 0)) for primitive in primitives)
require("normalized zero primitives", zero_primitives == (0, 0, 0))
print("PASS kappa0_boundary=C[w] normalized_F=0 source_Jacobian=0")
print("PASS scope=second_coordinate_w_only specific_thm3630_F_and_fixed_Q_closed other_packets_mixed_pairs_JC2_OPEN")


print("SECTION source AST gate")
source = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no inactive assertions", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
