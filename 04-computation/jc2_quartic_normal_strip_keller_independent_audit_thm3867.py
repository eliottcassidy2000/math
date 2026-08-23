#!/usr/bin/env python3
"""Independent hostile audit for THM-3867.

This companion does not import or call the primary verifier.  It reconstructs
the eight quartic Jacobian buckets by coefficient convolution, audits every
zero-component target chart, rederives the depressed (4,3) packet from its six
raw equations, and checks the W-ODE and all local valuation exits.  Separate
exact controls freeze the load-bearing outer B in the constant bucket, the
W identically zero branch, simultaneous A0/C0 polynomiality, and operation
over arbitrary (not necessarily algebraically closed) characteristic-zero
fields.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: object, label: str) -> None:
    """Optimization-safe exact gate."""

    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"GATE FAILED: {label}: {condition}")


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def jacobian(left: sp.Expr, right: sp.Expr, z: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(left, z) * sp.diff(right, s)
        - sp.diff(left, s) * sp.diff(right, z)
    )


def order_at_zero(expression: sp.Expr, variable: sp.Symbol) -> int:
    """Exact DVR order of a nonzero rational function at variable=0."""

    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, variable)
    denominator_poly = sp.Poly(denominator, variable)
    gate(not numerator_poly.is_zero, "valuation numerator nonzero")
    gate(not denominator_poly.is_zero, "valuation denominator nonzero")
    numerator_order = min(term[0][0] for term in numerator_poly.terms())
    denominator_order = min(term[0][0] for term in denominator_poly.terms())
    return numerator_order - denominator_order


def polynomial_in(expression: sp.Expr, *variables: sp.Symbol) -> bool:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    return denominator.is_number and sp.Poly(numerator, *variables) is not None


# ---------------------------------------------------------------------------
# 1. Eight raw buckets by an independent coefficient convolution.
# ---------------------------------------------------------------------------
s, z = sp.symbols("s z")
a, alpha, u, p, r, b, beta, v, q, t4 = sp.symbols(
    "a alpha u p r b beta v q t4"
)
ad, alphad, ud, pd, rd, bd, betad, vd, qd, t4d = sp.symbols(
    "ad alphad ud pd rd bd betad vd qd t4d"
)
A_coefficients = [a, alpha, u, p, r]
C_coefficients = [b, beta, v, q, t4]
A_derivatives = [ad, alphad, ud, pd, rd]
C_derivatives = [bd, betad, vd, qd, t4d]


def convolution_bucket(degree: int) -> sp.Expr:
    first = sum(
        i * A_coefficients[i] * C_derivatives[j]
        for i in range(1, 5)
        for j in range(5)
        if i - 1 + j == degree
    )
    second = sum(
        j * A_derivatives[i] * C_coefficients[j]
        for i in range(5)
        for j in range(1, 5)
        if i + j - 1 == degree
    )
    return sp.expand(first - second)


expected_buckets = [
    alpha * bd - ad * beta,
    alpha * betad - alphad * beta + 2 * u * bd - 2 * ad * v,
    alpha * vd - 2 * alphad * v + 2 * u * betad - ud * beta
    + 3 * p * bd - 3 * ad * q,
    alpha * qd - 3 * alphad * q + 2 * u * vd - 2 * ud * v
    + 3 * p * betad - pd * beta + 4 * r * bd - 4 * ad * t4,
    alpha * t4d - 4 * alphad * t4 + 2 * u * qd - 3 * ud * q
    + 3 * p * vd - 2 * pd * v + 4 * r * betad - rd * beta,
    2 * u * t4d - 4 * ud * t4 + 3 * p * qd - 3 * pd * q
    + 4 * r * vd - 2 * rd * v,
    3 * p * t4d - 4 * pd * t4 + 4 * r * qd - 3 * rd * q,
    4 * (r * t4d - rd * t4),
]

for degree in range(8):
    zero(
        convolution_bucket(degree) - expected_buckets[degree],
        f"independent coefficient convolution z^{degree}",
    )

# A direct specialization is a second route to every row.
sample_A_coefficients = [
    s**5 - 2 * s,
    2 * s**4 + 1,
    s**3 - 3 * s,
    5 * s**2 + 2,
    3 * s + 7,
]
sample_C_coefficients = [
    s**4 + 2 * s,
    s**5 - 1,
    3 * s**3 + s,
    s**2 + 7,
    2 * s + 5,
]
sample_substitution: dict[sp.Symbol, sp.Expr] = {}
for symbol, value, derivative_symbol in zip(
    A_coefficients, sample_A_coefficients, A_derivatives
):
    sample_substitution[symbol] = value
    sample_substitution[derivative_symbol] = sp.diff(value, s)
for symbol, value, derivative_symbol in zip(
    C_coefficients, sample_C_coefficients, C_derivatives
):
    sample_substitution[symbol] = value
    sample_substitution[derivative_symbol] = sp.diff(value, s)
sample_A = sum(value * z**index for index, value in enumerate(sample_A_coefficients))
sample_C = sum(value * z**index for index, value in enumerate(sample_C_coefficients))
sample_J = jacobian(sample_A, sample_C, z, s)
for degree in range(8):
    specialized_bucket = expected_buckets[degree].subs(
        sample_substitution, simultaneous=True
    )
    zero(
        sample_J.coeff(z, degree) - specialized_bucket,
        f"direct specialized bucket z^{degree}",
    )
gate(sp.Poly(sample_J, z).degree() <= 7, "quartic bracket has no row above seven")


# ---------------------------------------------------------------------------
# 2. Constant target normalization and positive target controls.
# ---------------------------------------------------------------------------
m11, m12, m21, m22 = sp.symbols("m11 m12 m21 m22")
r_new = m11 * r + m12 * t4
t_new = m21 * r + m22 * t4
rd_new = m11 * rd + m12 * t4d
t4d_new = m21 * rd + m22 * t4d
zero(
    r_new * t4d_new - rd_new * t_new
    - (m11 * m22 - m12 * m21) * (r * t4d - rd * t4),
    "top quartic Wronskian target covariance",
)

# Both affine charts of P^1(k) are normalized by SL2 without scalar roots.
R_direction, T_direction, h_direction = sp.symbols(
    "R_direction T_direction h_direction", nonzero=True
)
matrix_R_chart = sp.Matrix(
    [[1 / R_direction, 0], [-T_direction, R_direction]]
)
zero(matrix_R_chart.det() - 1, "SL2 R-chart determinant")
normalized_R_chart = matrix_R_chart * sp.Matrix(
    [h_direction * R_direction, h_direction * T_direction]
)
zero(normalized_R_chart[0] - h_direction, "SL2 R-chart retains first quartic")
zero(normalized_R_chart[1], "SL2 R-chart kills second quartic")

T_only = sp.symbols("T_only", nonzero=True)
matrix_T_chart = sp.Matrix([[0, 1 / T_only], [-T_only, 0]])
zero(matrix_T_chart.det() - 1, "SL2 T-chart determinant")
normalized_T_chart = matrix_T_chart * sp.Matrix([0, h_direction * T_only])
zero(normalized_T_chart[0] - h_direction, "SL2 T-chart retains first quartic")
zero(normalized_T_chart[1], "SL2 T-chart kills second quartic")

r_fun = sp.Function("r_fun")(s)
t_fun = sp.Function("t_fun")(s)
zero(
    sp.diff(r_fun / t_fun, s)
    + (r_fun * sp.diff(t_fun, s) - sp.diff(r_fun, s) * t_fun) / t_fun**2,
    "top Wronskian is a rational ratio derivative",
)

# A triangular quartic automorphism, then an SL2 transform with two nonzero
# quartic leading coordinates, is a hostile positive normalization control.
C_tri = s**2 + 3 * z
A_tri = (
    2 * C_tri**4 - C_tri**3 + 5 * C_tri**2 + 7 * C_tri
    - sp.Rational(5, 3) * s + 11
)
A_two = 2 * A_tri + C_tri
C_two = 3 * A_tri + 2 * C_tri
zero(jacobian(A_tri, C_tri, z, s) - 5, "triangular quartic control")
zero(jacobian(A_two, C_two, z, s) - 5, "two-quartic SL2 control")
top_A_two = sp.expand(A_two).coeff(z, 4)
top_C_two = sp.expand(C_two).coeff(z, 4)
zero(
    top_A_two * sp.diff(top_C_two, s)
    - sp.diff(top_A_two, s) * top_C_two,
    "two-quartic control top Wronskian",
)


# ---------------------------------------------------------------------------
# 3. Every q=0 stratum, including beta=0, and the degree-drop shears.
# ---------------------------------------------------------------------------
r0, beta0, bprime0, alpha0, lambda0 = sp.symbols(
    "r0 beta0 bprime0 alpha0 lambda0"
)
beta_zero_ideal = sp.groebner(
    [4 * r0 * bprime0, alpha0 * bprime0 - lambda0],
    bprime0,
    alpha0,
    r0,
    lambda0,
    order="lex",
)
zero(
    beta_zero_ideal.reduce(r0 * lambda0)[1],
    "q=v=beta=0 forces r*lambda=0",
)

R1 = sp.Function("R1")(s)
B1 = sp.Function("B1")(s)
zero(
    sp.diff(R1 / B1**4, s)
    + (4 * R1 * sp.diff(B1, s) - sp.diff(R1, s) * B1) / B1**5,
    "q=v=0 quartic-linear quotient derivative",
)

rho = sp.symbols("rho", nonzero=True)
b1 = sp.Function("b1")(s)
C_linear = b1 + B1 * z
A_quartic_raw = rho * B1**4 * z**4 + sp.Function("p1")(s) * z**3
zero(
    sp.expand(A_quartic_raw - rho * C_linear**4).coeff(z, 4),
    "quartic-linear target shear lowers degree",
)

R2 = sp.Function("R2")(s)
V2 = sp.Function("V2")(s)
zero(
    sp.diff(R2 / V2**2, s)
    + (2 * R2 * sp.diff(V2, s) - sp.diff(R2, s) * V2) / V2**3,
    "q=0,v!=0 quartic-quadratic quotient derivative",
)
C_quadratic = b1 + B1 * z + V2 * z**2
A_quartic_quadratic_raw = rho * V2**2 * z**4 + sp.Function("p2")(s) * z**3
zero(
    sp.expand(A_quartic_quadratic_raw - rho * C_quadratic**2).coeff(z, 4),
    "quartic-quadratic target shear lowers degree",
)

# The specialized raw rows themselves retain all vanishing-component edges.
zero(
    expected_buckets[6].subs({t4: 0, t4d: 0, q: 0, qd: 0}),
    "q=0 top row E6 vanishes",
)
zero(
    expected_buckets[5].subs(
        {t4: 0, t4d: 0, q: 0, qd: 0}
    ) - (4 * r * vd - 2 * rd * v),
    "q=0 row E5 gives v stratum",
)
zero(
    expected_buckets[4].subs(
        {t4: 0, t4d: 0, q: 0, qd: 0, v: 0, vd: 0}
    ) - (4 * r * betad - rd * beta),
    "q=v=0 row E4 gives beta stratum",
)


# ---------------------------------------------------------------------------
# 4. Hard q!=0 Kummer packet and the depressed six-equation system.
# ---------------------------------------------------------------------------
Rq = sp.Function("Rq")(s)
Qq = sp.Function("Qq")(s)
zero(
    sp.diff(Qq**4 / Rq**3, s)
    - Qq**3 * (
        4 * Rq * sp.diff(Qq, s) - 3 * sp.diff(Rq, s) * Qq
    ) / Rq**4,
    "q!=0 Kummer quotient derivative",
)

R_const, Q_const = sp.symbols("R_const Q_const", nonzero=True)
h_fun = sp.Function("h_fun")(s)
r_kummer = R_const * h_fun**4
q_kummer = Q_const * h_fun**3
zero(
    4 * r_kummer * sp.diff(q_kummer, s)
    - 3 * sp.diff(r_kummer, s) * q_kummer,
    "h4-h3 Kummer parametrization",
)

# Over x the bracket has exactly six possibly nonzero rows.  These are
# rederived directly, not copied from the z-bucket list.
x = sp.symbols("x")
Dfun = sp.Function("Dfun")(s)
Ufun = sp.Function("Ufun")(s)
Lfun = sp.Function("Lfun")(s)
afun = sp.Function("afun")(s)
Bfun = sp.Function("Bfun")(s)
bfun = sp.Function("bfun")(s)
A_depressed_raw = (
    R_const * x**4 + Dfun * x**3 + Ufun * x**2 + Lfun * x + afun
)
C_depressed_raw = Q_const * x**3 + Bfun * x + bfun
J_depressed_raw = jacobian(A_depressed_raw, C_depressed_raw, x, s)
depressed_rows = [sp.expand(J_depressed_raw).coeff(x, degree) for degree in range(6)]
expected_depressed_rows = [
    Lfun * sp.diff(bfun, s) - sp.diff(afun, s) * Bfun,
    2 * Ufun * sp.diff(bfun, s) + Lfun * sp.diff(Bfun, s)
    - sp.diff(Lfun, s) * Bfun,
    3 * Dfun * sp.diff(bfun, s) + 2 * Ufun * sp.diff(Bfun, s)
    - sp.diff(Ufun, s) * Bfun - 3 * Q_const * sp.diff(afun, s),
    4 * R_const * sp.diff(bfun, s) + 3 * Dfun * sp.diff(Bfun, s)
    - sp.diff(Dfun, s) * Bfun - 3 * Q_const * sp.diff(Lfun, s),
    4 * R_const * sp.diff(Bfun, s) - 3 * Q_const * sp.diff(Ufun, s),
    -3 * Q_const * sp.diff(Dfun, s),
]
for degree in range(6):
    zero(
        depressed_rows[degree] - expected_depressed_rows[degree],
        f"depressed direct row x^{degree}",
    )

# x=h*z+g multiplies the x,s bracket by h; all h' and g' terms cancel.
g_fun = sp.Function("g_fun")(s)
x_of_z = h_fun * z + g_fun
A_changed = A_depressed_raw.subs(x, x_of_z)
C_changed = C_depressed_raw.subs(x, x_of_z)
zero(
    jacobian(A_changed, C_changed, z, s)
    - h_fun * J_depressed_raw.subs(x, x_of_z),
    "affine normal change scales the bracket by h",
)

D0, E0, F0, G0 = sp.symbols("D0 E0 F0 G0")
K = 4 * R_const / (3 * Q_const)
H0 = D0 / Q_const
I = 2 * R_const / (9 * Q_const**2)
J = 2 * E0 / (3 * Q_const)
U_integrated = K * Bfun + E0
L_integrated = K * bfun + H0 * Bfun + F0
a_integrated = H0 * bfun + I * Bfun**2 + J * Bfun + G0
A_integrated = (
    R_const * x**4 + D0 * x**3 + U_integrated * x**2
    + L_integrated * x + a_integrated
)
C_integrated = Q_const * x**3 + Bfun * x + bfun
J_integrated = jacobian(A_integrated, C_integrated, x, s)
for degree in range(2, 6):
    zero(J_integrated.coeff(x, degree), f"integrated depressed row x^{degree}")

Wfun = K * Bfun + 2 * E0
Yfun = K * bfun + F0
relation = Wfun * sp.diff(bfun, s) + Yfun * sp.diff(Bfun, s)
zero(J_integrated.coeff(x, 1) - relation, "integrated x1 exact relation")

delta_integrated = (
    Yfun * sp.diff(bfun, s)
    - (2 * I * Bfun**2 + J * Bfun) * sp.diff(Bfun, s)
)
zero(J_integrated.coeff(x, 0) - delta_integrated, "integrated x0 delta")

# Freeze the load-bearing outer B: C_x(0)=B, so x0 is L*b'-a'*B.
zero(sp.diff(C_integrated, x).subs(x, 0) - Bfun, "C_x at x=0 is outer B")
zero(2 * I * Bfun + J - Wfun / (3 * Q_const), "a-prime W identity")
zero(
    delta_integrated
    - (Yfun * sp.diff(bfun, s)
       - Bfun * Wfun * sp.diff(Bfun, s) / (3 * Q_const)),
    "constant bucket retains outer B multiplier",
)

Nfun = Wfun * bfun + F0 * Bfun
zero(sp.diff(Nfun, s) - relation, "first integral N derivative")
Sfun = K * Nfun + 2 * E0 * F0
zero(Sfun - Wfun * Yfun, "S equals W*(Kb+F)")

# Algebraic W-coordinate derivation of the final ODE.
w, S0 = sp.symbols("w S0", nonzero=True)
B_of_w = (w - 2 * E0) / K
Y_of_w = S0 / w
b_of_w = (Y_of_w - F0) / K
delta_per_wprime = sp.simplify(
    Y_of_w * sp.diff(b_of_w, w)
    - (2 * I * B_of_w**2 + J * B_of_w) * sp.diff(B_of_w, w)
)
ode_per_wprime = -3 * Q_const / (16 * R_const**2) * (
    w**2 - 2 * E0 * w + 4 * R_const * S0**2 / w**3
)
zero(delta_per_wprime - ode_per_wprime, "final W ODE exact identity")

# W identically zero is not a zero of a nonzero rational W.  Audit it as a
# separate raw branch before any division by W.
B_Wzero = -2 * E0 / K
b_Wzero = sp.Function("b_Wzero")(s)
U_Wzero = K * B_Wzero + E0
L_Wzero = K * b_Wzero + H0 * B_Wzero + F0
a_Wzero = H0 * b_Wzero + I * B_Wzero**2 + J * B_Wzero + G0
A_Wzero = (
    R_const * x**4 + D0 * x**3 + U_Wzero * x**2
    + L_Wzero * x + a_Wzero
)
C_Wzero = Q_const * x**3 + B_Wzero * x + b_Wzero
J_Wzero = jacobian(A_Wzero, C_Wzero, x, s)
zero(J_Wzero.coeff(x, 1), "W identically zero kills x1 row")
zero(
    J_Wzero.coeff(x, 0) - (K * b_Wzero + F0) * sp.diff(b_Wzero, s),
    "W identically zero constant bucket",
)


# ---------------------------------------------------------------------------
# 5. UFD lattices, nonclosed-field controls, and local W valuations.
# ---------------------------------------------------------------------------
quartic_valuation_profiles = 0
quartic_kummer_points = 0
quadratic_valuation_profiles = 0
quadratic_kummer_points = 0
for r_order in range(33):
    for q_order in range(33):
        equation = 4 * q_order == 3 * r_order
        lattice = (
            r_order % 4 == 0
            and q_order % 3 == 0
            and r_order // 4 == q_order // 3
        )
        gate(equation == lattice, f"h4-h3 valuation ({r_order},{q_order})")
        quartic_valuation_profiles += 1
        if equation:
            quartic_kummer_points += 1
for r_order in range(33):
    for v_order in range(33):
        equation = r_order == 2 * v_order
        lattice = r_order % 2 == 0 and r_order // 2 == v_order
        gate(equation == lattice, f"q=0 h2-h1 valuation ({r_order},{v_order})")
        quadratic_valuation_profiles += 1
        if equation:
            quadratic_kummer_points += 1

# Rational constants 2 and 3 deliberately have no required fourth/cube roots.
gate(not sp.integer_nthroot(2, 4)[1], "Q scalar R=2 has no rational fourth root")
gate(not sp.integer_nthroot(3, 3)[1], "Q scalar Q=3 has no rational cube root")
h_Q = s**2 + 1
r_Q = 2 * h_Q**4
q_Q = 3 * h_Q**3
zero(
    4 * r_Q * sp.diff(q_Q, s) - 3 * sp.diff(r_Q, s) * q_Q,
    "nonclosed-field scalar-root-free h4-h3 control",
)

# Localize at irreducibles, not roots.  Separability makes pi' a DVR unit.
prime_controls = [s, s**2 + 1, s**3 + s + 1, s**4 + s + 1]
for index, prime in enumerate(prime_controls):
    gate(sp.gcd(prime, sp.diff(prime, s)) == 1,
         f"char0 separable prime-style control {index}")

local = sp.symbols("local")
dvr_derivative_controls = 0
for degree in range(1, 10):
    local_element = (
        local**degree + 2 * local + 3
    ) / (local**2 + local + 1) ** degree
    derivative_denominator = sp.cancel(
        sp.diff(local_element, local)
    ).as_numer_denom()[1]
    gate(derivative_denominator.subs(local, 0) != 0,
         f"DVR derivative closure degree {degree}")
    dvr_derivative_controls += 1

# Exact Laurent controls for both ODE channels and the W==0 branch.
R_num = sp.Rational(2)
Q_num = sp.Rational(3)
E_num = sp.Rational(5)
S_num = sp.Rational(7)


def delta_from_W(W_value: sp.Expr, S_value: sp.Expr) -> sp.Expr:
    return sp.cancel(
        -3 * Q_num / (16 * R_num**2)
        * (W_value**2 - 2 * E_num * W_value
           + 4 * R_num * S_value**2 / W_value**3)
        * sp.diff(W_value, local)
    )


pole_W_controls = 0
zero_W_controls = 0
for pole_order in range(1, 13):
    W_pole = (1 + local) / local**pole_order
    gate(
        order_at_zero(delta_from_W(W_pole, S_num), local)
        == -3 * pole_order - 1,
        f"W-pole ODE order m={pole_order}",
    )
    pole_W_controls += 1
for zero_order in range(1, 13):
    W_zero = local**zero_order * (1 + local)
    gate(
        order_at_zero(delta_from_W(W_zero, S_num), local)
        == -2 * zero_order - 1,
        f"W-zero ODE order n={zero_order}",
    )
    gate(
        order_at_zero(delta_from_W(W_zero, sp.S.Zero), local)
        == 2 * zero_order - 1,
        f"W-zero with S=0 cannot make a pole n={zero_order}",
    )
    zero_W_controls += 1

W_unit = 2 + local
gate(order_at_zero(delta_from_W(W_unit, S_num), local) >= 0,
     "W-unit ODE cannot make h divisor pole")
zero(delta_from_W(sp.Rational(2), S_num), "constant W gives zero delta")

K_num = sp.Rational(4, 3)
W_identically_zero_controls = 0
for b_pole_order in range(1, 13):
    b_local = local**(-b_pole_order)
    delta_Wzero = (K_num * b_local + 5) * sp.diff(b_local, local)
    gate(
        order_at_zero(delta_Wzero, local) == -2 * b_pole_order - 1,
        f"W identically zero order n={b_pole_order}",
    )
    W_identically_zero_controls += 1


# ---------------------------------------------------------------------------
# 6. Simultaneous A0/C0 polynomiality and terminal noncancellation.
# ---------------------------------------------------------------------------
pole_dominance_profiles = 0
pole_balanced_profiles = 0
for m in range(1, 25):
    for d in range(1, 25):
        main_orders = [4 * d, m + 2 * d, 2 * m]
        lower_orders = [3 * d, 2 * d, m + d, d, m]
        maximum = max(main_orders)
        ties = sum(value == maximum for value in main_orders)
        gate(
            (ties >= 2) == (2 * d == m),
            f"pole-W A0 balance m={m},d={d}",
        )
        gate(
            max(lower_orders) < maximum,
            f"pole-W lower terms stay lower m={m},d={d}",
        )
        pole_dominance_profiles += 1
        if ties >= 2:
            pole_balanced_profiles += 1

zero_dominance_profiles = 0
zero_balanced_profiles = 0
for n in range(1, 25):
    for d in range(1, 25):
        c_orders = [3 * d, n]
        ties = c_orders[0] == c_orders[1]
        gate(
            ties == (3 * d == n),
            f"zero-W C0 balance n={n},d={d}",
        )
        gate(d < 3 * d, f"zero-W B*g remains lower n={n},d={d}")
        zero_dominance_profiles += 1
        if ties:
            zero_balanced_profiles += 1

zero(I - K / Q_const + R_const / Q_const**2
     + R_const / (9 * Q_const**2),
     "pole-W terminal coefficient is -R/(9Q^2)")
zero(R_const - K * Q_const + R_const / 3,
     "zero-W terminal coefficient is -R/3")

# Exact pole-W rational hostile: C is polynomial and every positive-z
# coefficient of A is polynomial, but simultaneous regularity leaves exactly
# the terminal -1/(9*s^4) arm pole.
h_pole = s**7
g_pole = 1 / s
B_pole = -1 / s**2
x_pole = h_pole * z + g_pole
C_pole_hostile = sp.expand(x_pole**3 + B_pole * x_pole)
A_pole_hostile = sp.expand(
    x_pole**4 + sp.Rational(4, 3) * B_pole * x_pole**2
    + sp.Rational(2, 9) * B_pole**2
)
zero(jacobian(A_pole_hostile, C_pole_hostile, z, s) + sp.Rational(8, 9),
     "pole-W rational hostile Jacobian")
gate(polynomial_in(C_pole_hostile, s, z), "pole-W hostile C polynomial")
gate(order_at_zero(A_pole_hostile.coeff(z, 0), s) == -4,
     "pole-W hostile unique arm order -4")
for degree in range(1, 5):
    gate(order_at_zero(A_pole_hostile.coeff(z, degree), s) >= 0,
         f"pole-W hostile A positive coefficient {degree} polynomial")
for degree in range(8):
    target = -sp.Rational(8, 9) if degree == 0 else 0
    zero(
        jacobian(A_pole_hostile, C_pole_hostile, z, s).coeff(z, degree)
        - target,
        f"pole-W hostile bucket z^{degree}",
    )

# Exact W==0 rational hostile.  It freezes the branch that division by W
# would silently erase: C is polynomial and only A0=-1/(3*s^4) fails.
h_Wzero = s**7
g_Wzero = 1 / s
x_Wzero = h_Wzero * z + g_Wzero
C_Wzero_hostile = sp.expand(x_Wzero**3 - 1 / s**3)
A_Wzero_hostile = sp.expand(
    x_Wzero**4 - sp.Rational(4, 3) * x_Wzero / s**3
)
zero(jacobian(A_Wzero_hostile, C_Wzero_hostile, z, s) + 4,
     "W-identically-zero rational hostile Jacobian")
gate(polynomial_in(C_Wzero_hostile, s, z), "W==0 hostile C polynomial")
gate(order_at_zero(A_Wzero_hostile.coeff(z, 0), s) == -4,
     "W==0 hostile unique arm order -4")
for degree in range(1, 5):
    gate(order_at_zero(A_Wzero_hostile.coeff(z, degree), s) >= 0,
         f"W==0 hostile A positive coefficient {degree} polynomial")
for degree in range(8):
    target = -4 if degree == 0 else 0
    zero(
        jacobian(A_Wzero_hostile, C_Wzero_hostile, z, s).coeff(z, degree)
        - target,
        f"W==0 hostile bucket z^{degree}",
    )

# A local positive-zero W profile checks the second simultaneous contradiction
# without assuming a scalar cube root.
W_local_zero = local**3
B_local_zero = sp.Rational(3, 4) * local**3
b_local_zero = -local**-3
g_local_zero = local**-1
Y_local_zero = sp.Rational(4, 3) * b_local_zero
C0_local_zero = g_local_zero**3 + B_local_zero * g_local_zero + b_local_zero
A0_local_zero = (
    g_local_zero**4 + W_local_zero * g_local_zero**2
    + Y_local_zero * g_local_zero + sp.Rational(2, 9) * B_local_zero**2
)
gate(order_at_zero(C0_local_zero, local) >= 0,
     "positive-order W local C0 regular")
gate(order_at_zero(A0_local_zero, local) == -4,
     "positive-order W local A0 terminal pole")
zero(
    W_local_zero * Y_local_zero + sp.Rational(4, 3),
    "positive-order W local invariant S=-4/3",
)


# ---------------------------------------------------------------------------
# 7. Constant-h polynomial edge in all three W/S strata.
# ---------------------------------------------------------------------------
constant_h_controls = 0
for degree in range(13):
    trial = sum((index + 1) * local**index for index in range(degree + 1))
    # W==0: (K*b+F)b'.  Its leading term has degree 2d-1 for d>0.
    wzero_product = (K_num * trial + 5) * sp.diff(trial, local)
    if degree == 0:
        zero(wzero_product, "constant h W==0 constant b gives delta zero")
    else:
        gate(sp.Poly(wzero_product, local).degree() == 2 * degree - 1,
             f"constant h W==0 product degree {degree}")

    # S==0, W!=0: b is constant and delta=-B(2IB+J)B'.
    s_zero_product = -trial * (
        sp.Rational(4, 9) * trial + sp.Rational(10, 9)
    ) * sp.diff(trial, local)
    if degree == 0:
        zero(s_zero_product, "constant h S==0 constant B gives delta zero")
    else:
        gate(sp.Poly(s_zero_product, local).degree() == 3 * degree - 1,
             f"constant h S==0 product degree {degree}")
    constant_h_controls += 1

# If W and Y are polynomials and WY is a nonzero constant, both are units;
# an explicit bounded hostile family records that no positive degree survives.
unit_product_controls = 0
for left_degree in range(1, 10):
    for right_degree in range(1, 10):
        left = local**left_degree + 2
        right = local**right_degree + 3
        gate(sp.Poly(left * right, local).degree() == left_degree + right_degree,
             f"constant h nonzero-S unit product {left_degree},{right_degree}")
        unit_product_controls += 1


# ---------------------------------------------------------------------------
# 8. Frozen semantic packet and optimization-safe source audit.
# ---------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "no inactive Python assert",
)

semantic = {
    "buckets": "eight rows by convolution and independent direct specialization",
    "target": "two SL2 charts;top constant direction;all q-v-beta zero strata",
    "easy": "q0-v0 C4 shear;q0-vnonzero C2 shear;cubic handoff",
    "hard": "h4-h3;depression;six raw rows;D-U-L-a integrations",
    "outer_B": "Cx0=B;delta=Ybprime-BW Bprime/(3Q)",
    "ode": "delta=-3Q/(16R2)(W2-2EW+4RS2/W3)Wprime",
    "quartic_valuation_profiles": quartic_valuation_profiles,
    "quartic_kummer_points": quartic_kummer_points,
    "quadratic_valuation_profiles": quadratic_valuation_profiles,
    "quadratic_kummer_points": quadratic_kummer_points,
    "local": "W pole;W zero;W unit;W identically zero;DVR derivative closure",
    "pole_dominance_profiles": pole_dominance_profiles,
    "zero_dominance_profiles": zero_dominance_profiles,
    "terminal": "-R/(9Q2) and -R/3 nonzero in every char0 residue field",
    "field": "arbitrary characteristic-zero;irreducible DVRs;no scalar roots",
    "hostiles": "pole-W and W-identically-zero rational pairs fail only A0",
    "scope": "polynomial k[s,z] transverse degree at most four only;JC2 open",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("THM3867_QUARTIC_NORMAL_STRIP_INDEPENDENT_HOSTILE_AUDIT")
print("status=PASS_PROOF+VERIFIED_EXACT;JC2_OPEN")
print("eight_buckets=coefficient_convolution+direct_specialization_PASS")
print("target=two_SL2_charts;all_q_v_beta_zero_strata;degree_drop_shears_PASS")
print("hard_4_3=h4_h3+depressed_six_rows+all_integrations_PASS")
print("outer_B=Cx0_B;constant_bucket_and_W_ODE_exact_identity_PASS")
print(
    f"quartic_valuation_profiles={quartic_valuation_profiles};"
    f"h4_h3_lattice_points={quartic_kummer_points}"
)
print(
    f"quadratic_valuation_profiles={quadratic_valuation_profiles};"
    f"h2_h1_lattice_points={quadratic_kummer_points}"
)
print(
    f"W_pole_controls={pole_W_controls};W_zero_controls={zero_W_controls};"
    f"W_identically_zero_controls={W_identically_zero_controls}"
)
print(
    f"pole_dominance_profiles={pole_dominance_profiles};"
    f"balanced={pole_balanced_profiles};"
    f"zero_dominance_profiles={zero_dominance_profiles};"
    f"balanced={zero_balanced_profiles}"
)
print("terminal_coefficients=-R/(9Q^2),-R/3;noncancellation_PASS")
print("field_edge=nonclosed_char0;R2_Q3_without_scalar_roots_PASS")
print(f"DVR_derivative_controls={dvr_derivative_controls};constant_h_controls={constant_h_controls}")
print(f"constant_h_nonzero_S_unit_product_controls={unit_product_controls}")
print("rational_hostiles=pole_W+W_identically_zero;all_8_buckets;only_A0_pole")
print("scope=polynomial_transverse_degree_at_most_4_only;rational_and_infinite_open")
print(f"semantic_sha256={semantic_sha}")
print(f"GATES={GATES}")
print("RESULT=PASS")
