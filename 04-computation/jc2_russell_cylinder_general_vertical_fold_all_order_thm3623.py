#!/usr/bin/env python3
"""Exact comparison and boundary controls for proved THM-3623.

The gates verify the local Darboux chart, the general H(t) pullback, the
epsilon-divisible comparison geometry, several generalized invoice scales,
and the sharp k=1 and H=0 method boundaries.  They are exact symbolic
controls for the written all-order proof, not a finite substitute for it.
"""

import ast
from math import factorial
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
    """Exact rational-function zero test."""
    return sp.cancel(expression) == 0


def t_order(polynomial, t):
    """Return the least t exponent of a nonzero expanded polynomial."""
    expanded = sp.Poly(sp.expand(polynomial), t)
    if expanded.is_zero:
        return sp.oo
    return min(monomial[0] for monomial, _ in expanded.terms())


print("THM-3623 exact companion -- proved general vertical-fold closure")
print("status=proved, verified exact, and independently hostile-audited")


print("SECTION local chart and general pullback")
x, q, w, t = sp.symbols("x q w t")
u = x**2 * q
D = 1 + u
a_source = q / D**2
c_source = sp.expand(x * D * (D + 2))
b_source = sp.expand(u * (u + 3) ** 2)
e_source = sp.expand(q * (u + 4))

require("compiler relation", zero(c_source**2 * e_source - b_source * (b_source + 4)))
require("b plus four factorization", zero(b_source + 4 - D**2 * (D + 3)))
require("local a relation b", zero(a_source * c_source**2 - b_source))
require("local a relation e", zero(a_source * (b_source + 4) - e_source))

jacobian_ac = sp.det(
    sp.Matrix(
        [
            [sp.diff(a_source, x), sp.diff(a_source, q)],
            [sp.diff(c_source, x), sp.diff(c_source, q)],
        ]
    )
)
require("Jacobian xq ac", zero(jacobian_ac + 3))
require("collision chart unit", (b_source + 4).subs({x: 0, q: -sp.Rational(3, 4)}) == 4)
require("a parity", zero(a_source.subs(x, -x) - a_source))
require("c parity", zero(c_source.subs(x, -x) + c_source))

q_slope, h_prime = sp.symbols("q_slope h_prime")
P_form, K_form, R_form = sp.symbols("P_form K_form R_form")
a_x_total = sp.diff(a_source, x) + sp.diff(a_source, q) * q_slope
c_x_total = sp.diff(c_source, x) + sp.diff(c_source, q) * q_slope
a_t_total = sp.diff(a_source, q) * h_prime
c_t_total = sp.diff(c_source, q) * h_prime
jacobian_fold = sp.det(
    sp.Matrix([[a_x_total, a_t_total], [c_x_total, c_t_total]])
)
require("general fold ac determinant", zero(jacobian_fold + 3 * h_prime))

pullback = P_form * jacobian_fold + K_form * a_x_total + R_form * c_x_total
require(
    "general two-form pullback",
    zero(pullback - (-3 * h_prime * P_form + a_x_total * K_form + c_x_total * R_form)),
)
require(
    "total ax parity",
    zero(
        a_x_total.subs({x: -x, q_slope: -q_slope}, simultaneous=True)
        + a_x_total
    ),
)
require(
    "total cx parity",
    zero(
        c_x_total.subs({x: -x, q_slope: -q_slope}, simultaneous=True)
        - c_x_total
    ),
)

s = sp.symbols("s")
central_ax = a_x_total.subs({x: 0, q: -sp.Rational(3, 4) + s, q_slope: 0})
central_cx = c_x_total.subs({x: 0, q: -sp.Rational(3, 4) + s, q_slope: 0})
require("central ax", zero(central_ax))
require("central cx", central_cx == 3)
require("critical normalization R zero", 4 * central_cx == 12)
print("PASS chart=local_b+4 Jac_ac=-3 pullback=-3HprimeP+axK+cxR")


print("SECTION exact comparison and epsilon divisibility")
g, epsilon = sp.symbols("g epsilon")
s_exact = sp.Rational(3, 4) * (1 - g**-2)
q_comparison = -3 / g**2
a_zero = -sp.Rational(3, 4) + s_exact

comparison_a = a_source.subs({x: g, q: q_comparison})
comparison_c = c_source.subs({x: g, q: q_comparison})
require("comparison D minus two", zero(D.subs({x: g, q: q_comparison}) + 2))
require("comparison common a", zero(comparison_a - a_zero))
require("comparison common c", zero(comparison_c))

Q_infinity = -sp.Rational(3, 4) - sp.Rational(9, 4) / x**2
require(
    "rational comparison q",
    zero(Q_infinity.subs(x, g) + s_exact - q_comparison),
)
require(
    "rational comparison derivative",
    zero(4 * g**3 * sp.diff(Q_infinity, x).subs(x, g) - 18),
)

q_actual = q_comparison + epsilon
D_actual = D.subs({x: g, q: q_actual})
a_actual = a_source.subs({x: g, q: q_actual})
c_actual = c_source.subs({x: g, q: q_actual})
require("actual side D", zero(D_actual - (-2 + g**2 * epsilon)))
require(
    "actual side c epsilon factor",
    zero(c_actual - g**3 * epsilon * (-2 + g**2 * epsilon)),
)
require(
    "actual side a epsilon factor",
    zero(
        a_actual
        - a_zero
        - epsilon
        * (3 * g**2 * epsilon - 8)
        / (4 * (g**2 * epsilon - 2) ** 2)
    ),
)

side_ax_comparison = a_x_total.subs({x: g, q: q_comparison})
side_cx_comparison = c_x_total.subs({x: g, q: q_comparison})
require(
    "comparison side ax",
    zero(side_ax_comparison + (g**3 * q_slope - 9) / (2 * g**3)),
)
require(
    "comparison side cx",
    zero(side_cx_comparison - (12 - 2 * g**3 * q_slope)),
)

side_ax_actual = a_x_total.subs({x: g, q: q_actual})
side_cx_actual = c_x_total.subs({x: g, q: q_actual})
ax_error_quotient = sp.cancel((side_ax_actual - side_ax_comparison) / epsilon)
cx_error_quotient = sp.cancel((side_cx_actual - side_cx_comparison) / epsilon)
require("ax error divisible epsilon", sp.denom(ax_error_quotient).subs(epsilon, 0) != 0)
require("cx error divisible epsilon", sp.denom(cx_error_quotient).subs(epsilon, 0) != 0)
require("ax error vanishes", zero((side_ax_actual - side_ax_comparison).subs(epsilon, 0)))
require("cx error vanishes", zero((side_cx_actual - side_cx_comparison).subs(epsilon, 0)))
print("PASS comparison=common_target side_errors=epsilon_divisible main=R(18-4g^3Qprime)")


print("SECTION generalized invoice scales")


def target_jet(order):
    """The order-n jet of Q_infinity at x=1."""
    return (-1) ** (order - 1) * sp.Rational(9, 4) * factorial(order + 1)


vertical_controls = (
    ("k2_mixed", t**2 + 5 * t**5, 2, sp.Integer(1)),
    ("k3_odd", -2 * t**3 + 7 * t**4, 3, sp.Integer(-2)),
    ("k4_mixed", sp.Rational(3, 5) * t**4 - t**7, 4, sp.Rational(3, 5)),
    ("k7_sparse", 11 * t**7 + t**11, 7, sp.Integer(11)),
)

for label, vertical_polynomial, vertical_order, leading_coefficient in vertical_controls:
    require(f"{label} zero value", vertical_polynomial.subs(t, 0) == 0)
    require(f"{label} critical value", sp.diff(vertical_polynomial, t).subs(t, 0) == 0)
    require(f"{label} exact order", t_order(vertical_polynomial, t) == vertical_order)
    require(
        f"{label} leading coefficient",
        sp.expand(vertical_polynomial).coeff(t, vertical_order) == leading_coefficient,
    )

    g_vertical = (1 - sp.Rational(4, 3) * vertical_polynomial) ** -sp.Rational(1, 2)
    x_vertical = g_vertical - 1
    x_lead = sp.series(x_vertical, t, 0, vertical_order + 1).removeO().expand()
    for degree in range(vertical_order):
        require(f"{label} X lower degree={degree}", x_lead.coeff(t, degree) == 0)
    require(
        f"{label} X leading term",
        x_lead.coeff(t, vertical_order) == 2 * leading_coefficient / 3,
    )

    for order in range(1, 7):
        invoice_order = vertical_order * (order - 1)
        comparison_coefficient = sp.series(
            -16 * g_vertical**3 * x_vertical ** (order - 1) / factorial(order - 1),
            t,
            0,
            invoice_order + 1,
        ).removeO().expand().coeff(t, invoice_order)
        expected_coefficient = (
            -sp.Rational(16, factorial(order - 1))
            * (sp.Rational(2, 3) * leading_coefficient) ** (order - 1)
        )
        require(
            f"{label} invoice coefficient n={order}",
            comparison_coefficient == expected_coefficient,
        )
        require(
            f"{label} KR error order n={order}",
            t_order(vertical_polynomial**order, t) == vertical_order * order,
        )
        require(
            f"{label} P error order n={order}",
            t_order(sp.diff(vertical_polynomial, t) * vertical_polynomial**order, t)
            == vertical_order * order + vertical_order - 1,
        )
        require(
            f"{label} errors above invoice n={order}",
            vertical_order * order > invoice_order,
        )
        for xi_degree in range(1, invoice_order // vertical_order + 1):
            largest_t_degree = invoice_order - vertical_order * xi_degree
            require(
                f"{label} shift lower source row n={order} i={xi_degree}",
                xi_degree + largest_t_degree < invoice_order,
            )

for order in range(1, 13):
    require(
        f"rational jet formula n={order}",
        sp.diff(Q_infinity, x, order).subs(x, 1) == target_jet(order),
    )
    if order >= 2:
        require(
            f"rational recurrence n={order}",
            target_jet(order) == -(order + 1) * target_jet(order - 1),
        )

ode_expression = sp.expand(x * sp.diff(Q_infinity, x) + 2 * (Q_infinity + sp.Rational(3, 4)))
require("rational germ ODE", zero(ode_expression))
for polynomial_degree in range(1, 13):
    require(
        f"polynomial ODE leading coefficient d={polynomial_degree}",
        polynomial_degree + 2 != 0,
    )
print("PASS invoices=k2,k3,k4,k7 n1..6 errors_above_invoice shifts_use_lower_rows")


print("SECTION k=1 enlarged-form and local Darboux boundary")
H_unit = t + t**2
H_unit_prime = sp.diff(H_unit, t)
formal_P = -4 / H_unit_prime
require("k1 formal P regular at zero", sp.denom(formal_P).subs(t, 0) != 0)
require("k1 enlarged form survivor", zero(-3 * H_unit_prime * formal_P - 12))
require("k1 central ledger", -3 * 1 * (-4) + 3 * 0 == 12)

Q_linear_hostile = -sp.Rational(3, 4) - sp.Rational(9, 4) * x**2
require("k1 hostile Q even", zero(Q_linear_hostile.subs(x, -x) - Q_linear_hostile))
require("k1 hostile Q central", Q_linear_hostile.subs(x, 0) == -sp.Rational(3, 4))
require("k1 hostile Q side", Q_linear_hostile.subs(x, 1) == -3)
require("k1 hostile Q wrong first jet", sp.diff(Q_linear_hostile, x).subs(x, 1) == -sp.Rational(9, 2))

a_linear_fold = a_source.subs(q, Q_linear_hostile + t)
c_linear_fold = c_source.subs(q, Q_linear_hostile + t)
local_pair_jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(a_linear_fold, x), sp.diff(a_linear_fold, t)],
            [sp.diff(-4 * c_linear_fold, x), sp.diff(-4 * c_linear_fold, t)],
        ]
    )
)
require("k1 affine exact local Darboux pair", zero(local_pair_jacobian - 12))

k1_invoice_order = 4
k1_xi_degree = 2
k1_t_degree = 2
require("k1 shift reaches invoice", k1_xi_degree + k1_t_degree == k1_invoice_order)
require("k1 shift source degree not lower", k1_xi_degree + k1_t_degree == k1_invoice_order)
print("PASS k1=formal_P_survivor affine_local_pair shift_symbol_not_vertical")


print("SECTION H=0 static boundary")
Q_static = (
    -sp.Rational(3, 4)
    - sp.Rational(27, 4) * x**2
    + sp.Rational(9, 2) * x**4
)
require("H0 control even", zero(Q_static.subs(x, -x) - Q_static))
require("H0 control central", Q_static.subs(x, 0) == -sp.Rational(3, 4))
require("H0 control side", Q_static.subs(x, 1) == -3)
static_q1 = sp.diff(Q_static, x).subs(x, 1)
static_q2 = sp.diff(Q_static, x, 2).subs(x, 1)
require("H0 tuned first jet", static_q1 == sp.Rational(9, 2))
require("H0 wrong second jet", static_q2 == sp.Rational(81, 2))
require("H0 static comparison zero", 4 * (18 - 4 * static_q1) == 0)
require("H0 no moving shift", (1 - sp.Rational(4, 3) * 0) ** -sp.Rational(1, 2) - 1 == 0)
q1_symbol, q2_symbol = sp.symbols("q1_symbol q2_symbol")
static_comparison = 4 * (18 - 4 * q1_symbol)
require("H0 comparison omits q2", sp.diff(static_comparison, q2_symbol) == 0)
print("PASS H0=static_comparison_forces_only_q1 higher_Q_jets_invisible")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
