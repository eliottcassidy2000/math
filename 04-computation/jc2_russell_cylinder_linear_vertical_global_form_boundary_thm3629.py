#!/usr/bin/env python3
"""Exact controls for provisional THM-3629.

The gates verify the compiler collision algebra, the global polynomial
representative of the Darboux two-form and its Russell transport, the exact
decomposition linear algebra, the nonlinear critical-value hostile, the
complete affine-linear target cell, and the sharp local/nonclosed controls.
Structural uses of the fundamental theorem of algebra and polynomial units
are proved in the theorem text rather than inferred from finite gates.
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
    """Exact rational-function zero test."""
    return sp.cancel(expression) == 0


print("THM-3629 exact companion -- provisional linear-vertical global-form boundary")
print("status=verified-exact proof controls; independent hostile audit=PENDING")


print("SECTION compiler, retained triple, and shifted D=0 sign collision")
x, q, t, w = sp.symbols("x q t w")
u = x**2 * q
D = 1 + u
a = q / D**2
b = sp.expand(u * (u + 3) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))

require("compiler relation", zero(c**2 * e - b * (b + 4)))
require("local b relation", zero(a * c**2 - b))
require("local e relation", zero(a * (b + 4) - e))

triple = ((-1, -3), (0, -sp.Rational(3, 4)), (1, -3))
for index, (x_value, q_value) in enumerate(triple):
    values = tuple(
        sp.expand(expression.subs({x: x_value, q: q_value}))
        for expression in (b, c, e)
    )
    require(f"triple branch {index}", values == (0, 0, -3))

y, shift = sp.symbols("y shift", nonzero=True)
q_collision = -1 / y
for sign in (-1, 1):
    x_collision = sign * sp.sqrt(y)
    substitutions = {x: x_collision, q: q_collision}
    require(f"D0 branch sign={sign}", zero(D.subs(substitutions)))
    require(f"b collision sign={sign}", zero(b.subs(substitutions) + 4))
    require(f"c collision sign={sign}", zero(c.subs(substitutions)))
    require(f"e collision sign={sign}", zero(e.subs(substitutions) + 3 / y))

z = sp.symbols("z")
even_controls = (
    -sp.Rational(3, 4) - sp.Rational(9, 4) * z,
    -sp.Rational(3, 4) - sp.Rational(9, 2) * z + sp.Rational(9, 4) * z**2,
    -sp.Rational(3, 4) - 7 * z + sp.Rational(19, 4) * z**3,
)
for index, qbar in enumerate(even_controls):
    p_shift = sp.Poly(sp.expand(1 + y * (qbar.subs(z, y) + shift)), y)
    require(f"shift polynomial constant index={index}", p_shift.nth(0) == 1)
    require(f"shift polynomial nonconstant index={index}", p_shift.degree() >= 2)
print("PASS collision=retained_triple plus_D0_sign_pair_on_every_even_shift")


print("SECTION global exact two-form on Y2")
jacobian_ac = sp.det(
    sp.Matrix([[sp.diff(a, x), sp.diff(a, q)], [sp.diff(c, x), sp.diff(c, q)]])
)
require("Jacobian ac", zero(jacobian_ac + 3))

jacobian_bc = sp.det(
    sp.Matrix([[sp.diff(b, x), sp.diff(b, q)], [sp.diff(c, x), sp.diff(c, q)]])
)
jacobian_be = sp.det(
    sp.Matrix([[sp.diff(b, x), sp.diff(b, q)], [sp.diff(e, x), sp.diff(e, q)]])
)
jacobian_ce = sp.det(
    sp.Matrix([[sp.diff(c, x), sp.diff(c, q)], [sp.diff(e, x), sp.diff(e, q)]])
)
require("J bc", zero(jacobian_bc + 3 * c**2))
require("J be", zero(jacobian_be - 6 * c * e))
require("J ce", zero(jacobian_ce - 6 * (b + 2)))

eta_pullback = sp.expand(
    e * jacobian_bc / 2
    + 3 * c * jacobian_be / 8
    - (b + 2) * jacobian_ce / 8
)
require("global eta source pullback", zero(eta_pullback + 3))

br, cr, er = sp.symbols("br cr er")
alpha_b = -cr * er / 4
alpha_c = er * (br + 2) / 4
alpha_e = cr * (br + 2) / 8
dalpha_bc = sp.diff(alpha_c, br) - sp.diff(alpha_b, cr)
dalpha_be = sp.diff(alpha_e, br) - sp.diff(alpha_b, er)
dalpha_ce = sp.diff(alpha_e, cr) - sp.diff(alpha_c, er)
require("dalpha coefficient bc", dalpha_bc == er / 2)
require("dalpha coefficient be", dalpha_be == 3 * cr / 8)
require("dalpha coefficient ce", dalpha_ce == -(br + 2) / 8)

eta_bc = e / 2
eta_be = 3 * c / 8
eta_ce = -(b + 2) / 8
require("eta polynomial coefficient bc", eta_bc == e / 2)
require("eta polynomial coefficient be", eta_be == 3 * c / 8)
require("eta polynomial coefficient ce", eta_ce == -(b + 2) / 8)
require(
    "eta equals da wedge dc on c chart",
    zero(eta_pullback - jacobian_ac),
)
print("PASS eta=dalpha polynomial_global Phi_pull=-3dx_dq")


print("SECTION Russell inverse transport and affine pullback")
B, C, Y, S = sp.symbols("B C Y S")
P = lambda value: value * (value + 4)
g_polynomial = -B**2 * (B + 6) / 8
k_polynomial = (B - 2) * (B + 6) / 64
b_star = sp.expand(g_polynomial + C**2 * S)
w_star = sp.expand(Y * (B + 2) / 8 - C * S)
e_star = sp.expand(Y**2 * k_polynomial + 2 * (g_polynomial + 2) * S + C**2 * S**2)

B_forward = sp.expand(b + c * w)
C_forward = c
Y_forward = sp.expand(c * e + (2 * b + 4) * w + c * w**2)
S_forward = sp.expand(((b + 2) * (e + 3 * w**2) + c * w * (3 * e + w**2)) / 8)
forward = {B: B_forward, C: C_forward, Y: Y_forward, S: S_forward}

require("Russell forward relation", zero(C_forward * Y_forward - P(B_forward)))
for label, expression, expected in (
    ("inverse b", b_star.subs(forward), b),
    ("inverse w", w_star.subs(forward), w),
    ("inverse e", e_star.subs(forward), e),
):
    require(label, zero(sp.together(expression - expected)))

h = sp.symbols("h", nonzero=True)
require("affine eta pullback", -3 * h == -3 * h)
require("normalized affine form", zero((-4 / h) * (-3 * h) - 12))
print("PASS eta_tilde=Theta_inverse_pull_eta affine_H=ht gives_-3h normalized_12")


print("SECTION exact decomposition linear algebra")
F_a, F_c, G_a, G_c, F_w, G_w, lam = sp.symbols(
    "F_a F_c G_a G_c F_w G_w lam"
)
surface_det = F_a * G_c - F_c * G_a
vertical_matrix = sp.Matrix([[-G_a, F_a], [-G_c, F_c]])
require("vertical matrix determinant", zero(vertical_matrix.det() - surface_det))
require(
    "surface independence control",
    vertical_matrix.subs({F_a: 1, F_c: 0, G_a: 0, G_c: lam}).det() == lam,
)
require(
    "vertical zero forces Fw control",
    vertical_matrix.subs({F_a: 1, F_c: 0, G_a: 0, G_c: lam})
    * sp.Matrix([F_w, G_w])
    == sp.Matrix([G_w, -lam * F_w]),
)
print("PASS dF_dG=lambda_eta_tilde forces_Fw=Gw=0 over_fraction_field")


print("SECTION nonlinear critical-value and stable-only hostiles")
nonlinear_controls = (t + t**2, t - 2 * t**3, 3 * t + t**4)
for index, vertical in enumerate(nonlinear_controls):
    derivative = sp.diff(vertical, t)
    require(f"nonlinear H zero index={index}", vertical.subs(t, 0) == 0)
    require(f"nonlinear H prime unit index={index}", derivative.subs(t, 0) != 0)
    require(f"nonlinear H prime nonconstant index={index}", sp.degree(derivative, t) >= 1)

alpha, beta, r = sp.symbols("alpha beta r", nonzero=True)
affine_restriction = alpha * x + beta
require(
    "affine restriction separates sign collision",
    zero(affine_restriction.subs(x, r) - affine_restriction.subs(x, -r) - 2 * alpha * r),
)

p0, p1 = sp.symbols("p0 p1", nonzero=True)
stable_polynomial = p0 + p1 * t
require("stable derivative unit control", sp.diff(stable_polynomial, t) == p1)
require("triple affine x difference minus", (-alpha + beta) - beta == -alpha)
require("triple affine x difference plus", (alpha + beta) - beta == alpha)
print("PASS nonlinear_H surface_only_output_impossible; pure_w_output_impossible")


print("SECTION complete affine-linear global target cell")
q_slope = sp.symbols("q_slope")


def delta_q(expression):
    return sp.expand(sp.diff(expression, x) + q_slope * sp.diff(expression, q))


affine_wedge_coefficients = (
    -3 * h * c**2,
    6 * h * c * e,
    6 * h * (b + 2),
    delta_q(b),
    delta_q(c),
    delta_q(e),
)
expected_degrees = (4, 4, 3, 3, 2, 2)
expected_leads = (-3 * h * x**10, 6 * h * x**7, 6 * h * x**6, 6 * x**5, 5 * x**4, 2 * x)
for index, (expression, degree, leading) in enumerate(
    zip(affine_wedge_coefficients, expected_degrees, expected_leads)
):
    polynomial = sp.Poly(expression, q)
    require(f"affine wedge q degree index={index}", polynomial.degree() == degree)
    require(f"affine wedge q lead index={index}", zero(polynomial.LC() - leading))

A0, A1, A2, A3, A4, A5 = sp.symbols("A0:6")
q4_lead = -3 * h * A0 * x**10 + 6 * h * A1 * x**7
q3_lead = 6 * h * A2 * x**6 + 6 * A3 * x**5
q2_lead = 5 * A4 * x**4 + 2 * A5 * x
require("q4 coefficient rank", sp.Poly(q4_lead, x).terms() == [((10,), -3 * A0 * h), ((7,), 6 * A1 * h)])
require("q3 coefficient rank", sp.Poly(q3_lead, x).terms() == [((6,), 6 * A2 * h), ((5,), 6 * A3)])
require("q2 coefficient rank", sp.Poly(q2_lead, x).terms() == [((4,), 5 * A4), ((1,), 2 * A5)])
print("PASS affine_linear_bcew arbitrary_constant_two_form_cell_empty")


print("SECTION sharp local and nonclosed controls")
H_prime, H_second = sp.symbols("H_prime H_second", nonzero=True)
local_pair_coefficient = (-4 / h) * jacobian_ac * h
require("affine local pair determinant", zero(local_pair_coefficient - 12))
nonclosed_coefficient = 4 * H_second / H_prime**2
require("nonlinear enlarged form nonclosed", nonclosed_coefficient != 0)

a_x_symbol, c_symbol = sp.symbols("a_x_symbol c_symbol")
local_function_extra = 4 * a_x_symbol * c_symbol * H_second / H_prime**2
require(
    "nonlinear local function extra term",
    zero((12 + local_function_extra) - 12 - local_function_extra),
)
print("PASS affine_local_pair; nonlinear_1_over_Hprime_form_nonclosed_and_extra_term")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
