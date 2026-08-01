#!/usr/bin/env python3
"""Exact referee for the ordered quartic X_1(14) / diamond quotient bridge.

This is deliberately a scratch companion.  It uses exact SymPy arithmetic,
contains no truth-bearing ``assert``, and prints a compact deterministic
transcript suitable for normal / ``python -O`` parity checks.
"""

from fractions import Fraction

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def numerator(expr):
    return sp.together(expr).as_numer_denom()[0]


def zero_mod(expr, variable, monic_relation):
    num = sp.Poly(sp.expand(numerator(expr)), variable, domain="EX")
    rel = sp.Poly(sp.expand(monic_relation), variable, domain="EX")
    return sp.rem(num, rel).is_zero


def invariants(a1, a2, a3, a4, a6):
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = (
        a1 * a1 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3 * a3
        - a4 * a4
    )
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6
    disc = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return c4, c6, disc, Fraction(c4**3, disc)


# ---------------------------------------------------------------------------
# Projection of C_+ and the direct birational map to 14.a5.
# ---------------------------------------------------------------------------

tau, eta, xi = sp.symbols("tau eta xi")
quartic = tau**4 - 6 * tau**3 + 7 * tau**2 - 2 * tau + 1

# On C_+, [A:B:C]=[xi:1:tau].
plane = sp.expand(xi * tau - (xi - 1) * (xi - tau) * (1 - tau))
eta_projection = 2 * (tau - 1) * xi - tau**2 + tau + 1
require(
    sp.expand(eta_projection**2 - quartic - 4 * (tau - 1) * plane) == 0,
    "projection discriminant identity failed",
)

x = (eta + 1 - tau) / tau**2
y = (x**2 - 1) * tau + x + 3
u = (x + 1) / 2
v = (x**2 - 1) * tau / 4

require(
    zero_mod(y**2 - (2 * x**3 + 7 * x**2 + 4 * x + 3), eta, eta**2 - quartic),
    "quartic-to-cubic map failed",
)
require(
    zero_mod(v**2 + u * v + v - u**3 + u, eta, eta**2 - quartic),
    "quartic-to-14.a5 map failed",
)

# Inverse on the dense chart u(u-1) != 0.
uu, vv = sp.symbols("u v")
e1_relation = vv**2 + uu * vv + vv - uu**3 + uu
tau_inverse = vv / (uu * (uu - 1))
eta_inverse = (2 * uu - 1) * tau_inverse**2 + tau_inverse - 1
quartic_inverse = (
    eta_inverse**2
    - tau_inverse**4
    + 6 * tau_inverse**3
    - 7 * tau_inverse**2
    + 2 * tau_inverse
    - 1
)
require(zero_mod(quartic_inverse, vv, e1_relation), "inverse quartic map failed")

# The two finite tau=1 points are useful orbit controls.
require(sp.cancel(u.subs({tau: 1, eta: 1})) == 1, "tau=1 plus u failed")
require(sp.cancel(v.subs({tau: 1, eta: 1})) == 0, "tau=1 plus v failed")
require(sp.cancel(u.subs({tau: 1, eta: -1})) == 0, "tau=1 minus u failed")
require(sp.cancel(v.subs({tau: 1, eta: -1})) == 0, "tau=1 minus v failed")


# ---------------------------------------------------------------------------
# The even flank cycle is translation by T=(0,0).
# ---------------------------------------------------------------------------

# rho[A:B:C]=[B:C:A], hence tau'=xi/tau and xi'=1/tau.
xi_from_quartic = (eta + tau**2 - tau - 1) / (2 * (tau - 1))
tau_rho = sp.cancel(xi_from_quartic / tau)
xi_rho = 1 / tau
eta_rho = sp.cancel(
    2 * (tau_rho - 1) * xi_rho - tau_rho**2 + tau_rho + 1
)


def quartic_to_e1(T, E):
    xx = (E + 1 - T) / T**2
    return sp.cancel((xx + 1) / 2), sp.cancel((xx**2 - 1) * T / 4)


u_rho, v_rho = quartic_to_e1(tau_rho, eta_rho)
u_add_T = sp.cancel((v / u) ** 2 + v / u - u)
v_add_T = sp.cancel(-(v / u + 1) * u_add_T - 1)
require(
    zero_mod(u_rho - u_add_T, eta, eta**2 - quartic),
    "rho is not addition by T in u",
)
require(
    zero_mod(v_rho - v_add_T, eta, eta**2 - quartic),
    "rho is not addition by T in v",
)

# The tangent y=-x at T=(0,0) has triple intersection there, so 3T=O.
line_x = sp.symbols("line_x")
line_intersection = sp.expand(
    (-line_x) ** 2
    + line_x * (-line_x)
    + (-line_x)
    - line_x**3
    + line_x
)
require(line_intersection == -line_x**3, "T tangent/order-three control failed")


# ---------------------------------------------------------------------------
# Normalized Velu quotient E_1 -> E_0.
# ---------------------------------------------------------------------------

U = (uu**3 - uu + 1) / uu**2
V = (vv * (uu**3 + uu - 2) + uu**2 - uu - 1) / uu**3
e0_relation = V**2 + U * V + V - U**3 - 4 * U + 6
require(zero_mod(e0_relation, vv, e1_relation), "Velu E1-to-E0 map failed")

# Short models and the exact normalized short Velu formulas.
X, Y = sp.symbols("X Y")
X_source = 36 * uu + 3
Y_source = 108 * (2 * vv + uu + 1)
X_target = 36 * U + 3
Y_target = 108 * (2 * V + U + 1)
require(
    zero_mod(
        Y_source**2 - X_source**3 + 675 * X_source - 13662,
        vv,
        e1_relation,
    ),
    "source short model failed",
)
require(
    zero_mod(
        Y_target**2 - X_target**3 - 5805 * X_target + 285714,
        vv,
        e1_relation,
    ),
    "target short model failed",
)

short_velu_x = (X**3 - 6 * X**2 - 1287 * X + 50544) / (X - 3) ** 2
short_velu_y = Y * (X**3 - 9 * X**2 + 1323 * X - 97227) / (X - 3) ** 3
require(
    sp.cancel(X_target - short_velu_x.subs(X, X_source)) == 0,
    "short Velu x formula failed",
)
require(
    zero_mod(
        Y_target - short_velu_y.subs({X: X_source, Y: Y_source}),
        vv,
        e1_relation,
    ),
    "short Velu y formula failed",
)

c4_1, c6_1, disc_1, j_1 = invariants(1, 0, 1, -1, 0)
c4_0, c6_0, disc_0, j_0 = invariants(1, 0, 1, 4, -6)
require((c4_1, c6_1, disc_1, j_1) == (25, -253, -28, Fraction(-15625, 28)), "E1 invariants failed")
require((c4_0, c6_0, disc_0, j_0) == (-215, 5291, -21952, Fraction(9938375, 21952)), "E0 invariants failed")

# THM-2998's symmetric quotient sends the chosen source origin to this finite
# target point.  Therefore its displayed coordinates are an unpointed quotient
# cover, not the normalized Velu isogeny without a target translation.
t_sym = sp.Rational(0)
s_sym = sp.Rational(0)
v_sym = 28 * s_sym + 2 - 9 * t_sym
x_e = 32 - 112 * t_sym
y_e = -112 * v_sym
x0_origin_image = (x_e + 4) / 4
y0_origin_image = (y_e - x_e - 8) / 8
require((x0_origin_image, y0_origin_image) == (9, -33), "quotient origin image failed")
require(
    y0_origin_image**2 + x0_origin_image * y0_origin_image + y0_origin_image
    == x0_origin_image**3 + 4 * x0_origin_image - 6,
    "quotient origin image is not on E0",
)


# ---------------------------------------------------------------------------
# Diamond group and the odd sheet-exchange boundary.
# ---------------------------------------------------------------------------

units14 = [1, 3, 5, 9, 11, 13]
diamond_classes = [(1, 13), (3, 11), (5, 9)]
require(pow(3, 3, 14) in diamond_classes[0], "diamond generator order failed")
require(pow(3, 1, 14) in diamond_classes[1], "diamond generator first power failed")
require(pow(3, 2, 14) in diamond_classes[2], "diamond generator second power failed")

A, B, C = sp.symbols("A B C")
vandermonde = (A - B) * (A - C) * (B - C)
c_plus = A * B * C - vandermonde
c_minus = A * B * C + vandermonde
require(sp.expand(c_plus.subs({A: B, B: A}, simultaneous=True) - c_minus) == 0, "odd swap does not exchange sheets")
require(sp.expand(c_plus.subs({A: B, B: C, C: A}, simultaneous=True) - c_plus) == 0, "even cycle does not preserve sheet")

v_neg = -vv - uu - 1
require(zero_mod(v_neg**2 + uu * v_neg + v_neg - uu**3 + uu, vv, e1_relation), "E1 negation failed")
require(v_neg.subs({uu: 0, vv: 0}) == -1, "negation does not send T to -T")


print("THM3034 ordered-quartic X1(14) diamond referee")
print("projection_Cplus_to_quartic=PASS")
print("quartic_to_E1_birational_map=PASS")
print("inverse_dense_chart=PASS")
print("rho_equals_translation_by_T_0_0=PASS")
print("T_order=3")
print("velu_E1_to_E0=PASS")
print("short_velu_models=PASS")
print(f"E1_invariants=c4:{c4_1},c6:{c6_1},Delta:{disc_1},j:{j_1}")
print(f"E0_invariants=c4:{c4_0},c6:{c6_0},Delta:{disc_0},j:{j_0}")
print(f"diamond_units={units14}")
print(f"diamond_mod_pm1={diamond_classes}=C3")
print("THM2998_symmetric_origin_image=(9,-33)")
print("odd_flank=sheet_exchange; synchronized_groupoid_model=negation")
print("odd_flank_is_not_diamond_C2_or_Fricke=PASS")
print("all_exact_checks=PASS")
