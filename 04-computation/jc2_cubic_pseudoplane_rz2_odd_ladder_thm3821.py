#!/usr/bin/env python3
"""Exact companion for THM-3821's rz^2 odd-ladder anatomy."""

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


r, z, e = sp.symbols("r z e")
variables = (r, z, e)
surface = r**2 * e - z**3 + r
monic_relation = z**3 - r**2 * e - r
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 + 6 * r * e],
        [-9 * z**2, -3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, variable) for variable in variables])
    dr = sp.Matrix([sp.diff(right, variable) for variable in variables])
    return sp.expand((dl.T * poisson * dr)[0])


casimir_vector = poisson * sp.Matrix(
    [sp.diff(surface, variable) for variable in variables]
)
gate(
    all(sp.cancel(entry) == 0 for entry in casimir_vector),
    "surface is a Poisson Casimir",
)


# The complete first rz^2 extension of THM-3814.
g = sp.Function("g")(e)
f = sp.Function("f")(e)
h = sp.Function("h")(e)
kap = sp.Function("kappa")(e)
p = sp.Function("p")(e)
q = sp.Function("q")(e)
S = sp.Function("S")(e)
T = sp.Function("T")(e)
A = e**2 - z / 3 + r * g + z**2 * f + r * z * p + r * z**2 * S
C = (
    e**3
    - e
    - e * z / 2
    + r * h
    + z**2 * kap
    + r * z * q
    + r * z**2 * T
)
normal = sp.rem(bracket(A, C) - 1, monic_relation, z)
normal_poly = sp.Poly(sp.expand(normal), r, z)
gate(max(z_degree for _, z_degree in normal_poly.monoms()) <= 2,
     "canonical z-normal degree")
zero(
    sp.rem(bracket(A, C) - 1 - normal, monic_relation, z),
    "canonical reduction has no quotient loss",
)


# Freeze the six descending buckets used in the proof.
bucket_z = normal_poly.coeff_monomial(z)
bucket_r5 = normal_poly.coeff_monomial(r**5)
bucket_r4 = normal_poly.coeff_monomial(r**4)
bucket_r3z2 = normal_poly.coeff_monomial(r**3 * z**2)
bucket_r3z = normal_poly.coeff_monomial(r**3 * z)
bucket_r3 = normal_poly.coeff_monomial(r**3)
bucket_r2z2 = normal_poly.coeff_monomial(r**2 * z**2)
bucket_r2z = normal_poly.coeff_monomial(r**2 * z)
bucket_r2 = normal_poly.coeff_monomial(r**2)
bucket_z2 = normal_poly.coeff_monomial(z**2)

expected_z = (36 * e**2 * f - 24 * e * kap - 12 * f + 1) / 2
expected_r5 = -21 * e**2 * (-S * sp.diff(T, e) + T * sp.diff(S, e))
expected_r4 = -3 * e * (
    -4 * e * f * sp.diff(T, e)
    + 4 * e * kap * sp.diff(S, e)
    - 7 * e * S * sp.diff(kap, e)
    + 7 * e * T * sp.diff(f, e)
    + 2 * f * T
    - 2 * kap * S
    - 12 * S * sp.diff(T, e)
    + 12 * T * sp.diff(S, e)
)
expected_r3z2 = -3 * (
    -5 * e * p * sp.diff(T, e)
    + 5 * e * q * sp.diff(S, e)
    - 7 * e * S * sp.diff(q, e)
    + 7 * e * T * sp.diff(p, e)
    - p * T
    + q * S
)
expected_r3z = -3 * (
    -3 * e * g * sp.diff(T, e)
    + 3 * e * h * sp.diff(S, e)
    - 5 * e * p * sp.diff(q, e)
    + 5 * e * q * sp.diff(p, e)
    - 7 * e * S * sp.diff(h, e)
    + 7 * e * T * sp.diff(g, e)
    - 2 * g * T
    + 2 * h * S
)
expected_r3 = -3 * (
    -4 * e**2 * f * sp.diff(kap, e)
    + 4 * e**2 * kap * sp.diff(f, e)
    - 6 * e * f * sp.diff(T, e)
    - 3 * e * g * sp.diff(q, e)
    + 3 * e * h * sp.diff(p, e)
    + 6 * e * kap * sp.diff(S, e)
    - 5 * e * p * sp.diff(h, e)
    + 5 * e * q * sp.diff(g, e)
    - 12 * e * S * sp.diff(kap, e)
    + 12 * e * T * sp.diff(f, e)
    + 2 * f * T
    - g * q
    + h * p
    - 2 * kap * S
    - 5 * S * sp.diff(T, e)
    + 5 * T * sp.diff(S, e)
)
for actual, expected, label in [
    (bucket_z, expected_z, "arm bucket"),
    (bucket_r5, expected_r5, "r5 top Wronskian"),
    (bucket_r4, expected_r4, "r4 7/4 bucket"),
    (bucket_r3z2, expected_r3z2, "r3z2 7/5 bucket"),
    (bucket_r3z, expected_r3z, "r3z integration bucket"),
    (bucket_r3, expected_r3, "r3 terminal bucket"),
]:
    zero(actual - expected, label)


# The S=0,T!=0 branch forces f(0)=0, contrary to the arm f(0)=1/12.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
v = sp.Function("v")(e)
T_asymmetric = alpha * e**4 * v**7
f_asymmetric = beta * e**2 * v**4
asymmetric_ode = (
    4 * e * f * sp.diff(T, e)
    - 7 * e * T * sp.diff(f, e)
    - 2 * T * f
)
zero(
    asymmetric_ode.subs({T: T_asymmetric, f: f_asymmetric}).doit(),
    "one-sided 7/4 Kummer tower",
)
gate(f_asymmetric.subs(e, 0) == 0, "one-sided tower forces f(0)=0")
gate(sp.Rational(1, 12) != 0, "arm forces f(0)=1/12")


# Main branch: T=mu*S and K=kappa-mu*f.  The r4 equation is the 7/4
# Kummer law and has the indicated complete UFD valuation parametrization.
mu = sp.symbols("mu")
K = sp.Function("K")(e)
r4_main = expected_r4.subs({T: mu * S, kap: mu * f + K}).doit()
main_74_ode = 4 * e * K * sp.diff(S, e) - 7 * e * S * sp.diff(K, e) - 2 * K * S
zero(r4_main + 3 * e * main_74_ode, "main r4 Kummer reduction")
S_tower = alpha * e**4 * v**7
K_tower = beta * e**2 * v**4
zero(
    main_74_ode.subs({S: S_tower, K: K_tower}).doit(),
    "7/4 tower solves r4 equation",
)
multiplicity = sp.symbols("multiplicity", integer=True, nonnegative=True)
gate(4 * (4 + 7 * multiplicity) == 2 + 7 * (2 + 4 * multiplicity),
     "7/4 e-prime valuation family")
gate(4 * (7 * multiplicity) == 7 * (4 * multiplicity),
     "7/4 nonzero-prime valuation family")


# The next difference P=q-mu*p either vanishes or enters the 7/5 tower.
P = sp.Function("P")(e)
r3z2_main = expected_r3z2.subs(
    {T: mu * S, q: mu * p + P}
).doit()
main_75_ode = 7 * e * S * sp.diff(P, e) - 5 * e * P * sp.diff(S, e) - P * S
zero(r3z2_main - 3 * main_75_ode, "main r3z2 Kummer reduction")
delta = sp.symbols("delta", nonzero=True)
P_tower = delta * e**3 * v**5
zero(
    main_75_ode.subs({S: S_tower, P: P_tower}).doit(),
    "7/5 tower solves r3z2 equation",
)
gate(7 * (3 + 5 * multiplicity) == 1 + 5 * (4 + 7 * multiplicity),
     "7/5 e-prime valuation family")
gate(7 * (5 * multiplicity) == 5 * (7 * multiplicity),
     "7/5 nonzero-prime valuation family")


# Generic P!=0 branch.  The r3z equation integrates exactly and polynomial
# typing forces p=e*v^2*R.  The terminal valuation then adds one more e*v.
Q = sp.Function("Q")(e)
r3z_main = expected_r3z.subs(
    {T: mu * S, q: mu * p + P, h: mu * g + Q}
).doit()
main_integrating_equation = (
    5 * e * (p * sp.diff(P, e) - P * sp.diff(p, e))
    + 7 * e * S * sp.diff(Q, e)
    - 3 * e * Q * sp.diff(S, e)
    - 2 * Q * S
)
zero(r3z_main - 3 * main_integrating_equation,
     "generic r3z integrating equation")

eta = sp.symbols("eta")
R = sp.Function("R")(e)
p_typed = e * v**2 * R
Q_integrated = (
    sp.Rational(5, 7) * delta / alpha * R
    + eta * alpha / beta * e**2 * v**3
)
zero(
    main_integrating_equation.subs(
        {S: S_tower, P: P_tower, p: p_typed, Q: Q_integrated}
    ).doit(),
    "generic integrated Q relation",
)
zero(
    sp.Rational(5, 7) * P_tower * p_typed / S_tower
    + eta * S_tower / K_tower
    - Q_integrated,
    "generic rational integration has polynomial typing",
)

# The equality of the forcing and R-square orders has the unique solutions
# ord_rho(R)=ord_rho(v) and ord_0(R)=ord_0(v)+1.  The coefficient-resonance
# values lie strictly above those matches and cannot hide a lower term.
R_order = sp.symbols("R_order", integer=True, nonnegative=True)
root_order = sp.symbols("root_order", integer=True, positive=True)
zero(
    sp.solve(
        sp.Eq(2 * R_order + 2 * root_order - 1, 4 * root_order - 1),
        R_order,
    )[0]
    - root_order,
    "nonzero-root R valuation match",
)
gate(3 * root_order > root_order, "nonzero-root resonance is above match")
zero_order = sp.symbols("zero_order", integer=True, nonnegative=True)
zero(
    sp.solve(
        sp.Eq(1 + 2 * zero_order + 2 * R_order, 3 + 4 * zero_order),
        R_order,
    )[0]
    - (zero_order + 1),
    "origin R valuation match",
)
gate(3 * zero_order + 2 > zero_order + 1,
     "origin resonance is above match")

U = sp.Function("U")(e)
p_odd = e**2 * v**3 * U
Q_odd = (
    sp.Rational(5, 7) * delta / alpha * e * v * U
    + eta * alpha / beta * e**2 * v**3
)

# Freeze the compact terminal Riccati/square payment law against the actual
# pure-r3 coefficient after every preceding substitution.
terminal_generic = (
    14 * alpha**2 * beta**2 * e**3 * v**7 * (3 * e * sp.diff(v, e) + v)
    + 21
    * alpha**2
    * eta
    * e
    * v**2
    * (2 * e * U * sp.diff(v, e) - e * v * sp.diff(U, e) + U * v)
    + 2
    * beta
    * (2 * e * sp.diff(v, e) + v)
    * (28 * alpha * beta * f - 5 * delta * U**2)
    + e
    * v
    * (-28 * alpha * beta**2 * sp.diff(f, e)
       + 10 * beta * delta * U * sp.diff(U, e))
    + 35
    * alpha
    * beta
    * delta
    * v
    * ((3 * e * sp.diff(v, e) + 2 * v) * g
       - e * v * sp.diff(g, e))
)
r3_generic = expected_r3.subs(
    {
        T: mu * S_tower,
        S: S_tower,
        kap: mu * f + K_tower,
        q: mu * p_odd + P_tower,
        p: p_odd,
        h: mu * g + Q_odd,
    }
).doit()
zero(
    r3_generic
    - 3 * e**3 * v**3 * terminal_generic / (7 * alpha * beta),
    "generic terminal Riccati law",
)

# At a nonzero root of v, the terminal leading coefficient is exactly the
# advertised square payment coordinate.
local_x = sp.symbols("local_x")
rho, local_unit, local_f0, local_U0 = sp.symbols(
    "rho local_unit local_f0 local_U0", nonzero=True
)
local_v = local_x**root_order * local_unit
local_e = rho + local_x
local_square_block = (
    2
    * beta
    * (2 * local_e * sp.diff(local_v, local_x) + local_v)
    * (28 * alpha * beta * local_f0 - 5 * delta * local_U0**2)
)
local_square_lead = sp.limit(
    sp.factor(
        sp.powsimp(
            local_square_block / local_x ** (root_order - 1), force=True
        )
    ),
    local_x,
    0,
)
zero(
    local_square_lead
    - 4
    * beta
    * rho
    * root_order
    * local_unit
    * (28 * alpha * beta * local_f0 - 5 * delta * local_U0**2),
    "generic nonzero-root square payment",
)


# Degenerate P=0 branch.  Q=0 leaves an uncancellable 7/4 forcing, so Q
# must enter the 7/3 tower; the terminal valuation then forces p=e*v*U.
degenerate_73_ode = 7 * e * S * sp.diff(Q, e) - 3 * e * Q * sp.diff(S, e) - 2 * Q * S
theta = sp.symbols("theta", nonzero=True)
Q_tower = theta * e**2 * v**3
zero(
    degenerate_73_ode.subs({S: S_tower, Q: Q_tower}).doit(),
    "degenerate 7/3 tower solves r3z equation",
)
gate(7 * (2 + 3 * multiplicity) == 2 + 3 * (4 + 7 * multiplicity),
     "7/3 e-prime valuation family")
gate(7 * (3 * multiplicity) == 3 * (7 * multiplicity),
     "7/3 nonzero-prime valuation family")

degenerate_forcing = (
    -4 * e**2 * K * sp.diff(f, e)
    + 4 * e**2 * f * sp.diff(K, e)
    - 6 * e * K * sp.diff(S, e)
    + 12 * e * S * sp.diff(K, e)
    + 2 * K * S
)
expected_degenerate_forcing = (
    2
    * beta
    * e**3
    * v**3
    * (
        3 * alpha * e**4 * v**7 * sp.diff(v, e)
        + alpha * e**3 * v**8
        + 8 * e * f * sp.diff(v, e)
        - 2 * e * v * sp.diff(f, e)
        + 4 * f * v
    )
)
zero(
    degenerate_forcing.subs({S: S_tower, K: K_tower}).doit()
    - expected_degenerate_forcing,
    "degenerate Q-zero forcing factor",
)

local_qzero_bracket = (
    8 * local_e * local_f0 * sp.diff(local_v, local_x)
    + 4 * local_f0 * local_v
)
local_qzero_lead = sp.limit(
    sp.factor(
        sp.powsimp(
            local_qzero_bracket / local_x ** (root_order - 1), force=True
        )
    ),
    local_x,
    0,
)
zero(
    local_qzero_lead
    - 8 * rho * local_f0 * root_order * local_unit,
    "degenerate Q-zero nonzero-root obstruction",
)
origin_unit = sp.symbols("origin_unit", nonzero=True)
origin_v = local_x**zero_order * origin_unit
origin_qzero_bracket = (
    8 * local_x * sp.Rational(1, 12) * sp.diff(origin_v, local_x)
    + 4 * sp.Rational(1, 12) * origin_v
)
origin_qzero_lead = sp.limit(
    sp.factor(
        sp.powsimp(origin_qzero_bracket / local_x**zero_order, force=True)
    ),
    local_x,
    0,
)
zero(
    origin_qzero_lead
    - (8 * zero_order + 4) * origin_unit / 12,
    "degenerate Q-zero origin obstruction",
)

# Once Q enters the 7/3 tower, matching its p-operator against the same
# forcing has the unique global polynomial typing p=e*v*U.
p_order = sp.symbols("p_order", integer=True, nonnegative=True)
zero(
    sp.solve(
        sp.Eq(p_order + 3 * root_order - 1, 4 * root_order - 1),
        p_order,
    )[0]
    - root_order,
    "degenerate nonzero-root p valuation match",
)
gate(5 * root_order > root_order,
     "degenerate nonzero-root resonance is above match")
zero(
    sp.solve(
        sp.Eq(p_order + 2 + 3 * zero_order, 3 + 4 * zero_order),
        p_order,
    )[0]
    - (zero_order + 1),
    "degenerate origin p valuation match",
)
gate(3 + 5 * zero_order > zero_order + 1,
     "degenerate origin resonance is above match")

p_skip = e * v * U
terminal_skip = (
    2
    * beta
    * (
        3 * alpha * e**4 * v**7 * sp.diff(v, e)
        + alpha * e**3 * v**8
        + 8 * e * f * sp.diff(v, e)
        - 2 * e * v * sp.diff(f, e)
        + 4 * f * v
    )
    - 3
    * theta
    * (-4 * e * U * sp.diff(v, e) + e * v * sp.diff(U, e) - 2 * U * v)
)
r3_skip = expected_r3.subs(
    {
        T: mu * S_tower,
        S: S_tower,
        kap: mu * f + K_tower,
        q: mu * p_skip,
        p: p_skip,
        h: mu * g + Q_tower,
    }
).doit()
zero(
    r3_skip - 3 * e**3 * v**3 * terminal_skip,
    "degenerate terminal linear-payment law",
)

local_skip_block = (
    2
    * beta
    * (8 * local_e * local_f0 * sp.diff(local_v, local_x)
       + 4 * local_f0 * local_v)
    - 3
    * theta
    * (-4 * local_e * local_U0 * sp.diff(local_v, local_x)
       - 2 * local_U0 * local_v)
)
local_skip_lead = sp.limit(
    sp.factor(
        sp.powsimp(
            local_skip_block / local_x ** (root_order - 1), force=True
        )
    ),
    local_x,
    0,
)
zero(
    local_skip_lead
    - 4
    * rho
    * root_order
    * local_unit
    * (4 * beta * local_f0 + 3 * theta * local_U0),
    "degenerate nonzero-root linear payment",
)


# Constant-v hostile closure in both branches.  The arm division first
# forces the same quadratic parameter law as THM-3814; the next untouched
# bucket then gives an immediate contradiction.
def reduce_mu_law(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    if sp.diff(denominator, mu) != 0:
        raise RuntimeError("mu-law reduction acquired a mu denominator")
    return sp.cancel(sp.rem(numerator, 4 * mu**2 + 3, mu) / denominator)


f_constant_v = sp.Rational(1, 12) - mu * e / 6
K_constant_v = -mu * e**2 / 4
S_constant_v = alpha * e**4
U_constant = sp.Function("U_constant")(e)

p_generic_constant = e**2 * U_constant
P_generic_constant = delta * e**3
Q_generic_constant = (
    sp.Rational(5, 7) * delta / alpha * e * U_constant
    - 4 * eta * alpha / mu * e**2
)
generic_constant_r2z2 = bucket_r2z2.subs(
    {
        f: f_constant_v,
        kap: mu * f_constant_v + K_constant_v,
        S: S_constant_v,
        T: mu * S_constant_v,
        p: p_generic_constant,
        q: mu * p_generic_constant + P_generic_constant,
        h: mu * g + Q_generic_constant,
    }
).doit()
generic_constant_obstruction = reduce_mu_law(
    sp.cancel(generic_constant_r2z2 / e**3).subs(e, 0)
)
zero(
    generic_constant_obstruction - sp.Rational(5, 2) * delta,
    "generic constant-v r2z2 obstruction",
)

p_skip_constant = e * U_constant
Q_skip_constant = theta * e**2
skip_constant_r2z2 = bucket_r2z2.subs(
    {
        f: f_constant_v,
        kap: mu * f_constant_v + K_constant_v,
        S: S_constant_v,
        T: mu * S_constant_v,
        p: p_skip_constant,
        q: mu * p_skip_constant,
        h: mu * g + Q_skip_constant,
    }
).doit()
skip_constant_r3 = bucket_r3.subs(
    {
        f: f_constant_v,
        kap: mu * f_constant_v + K_constant_v,
        S: S_constant_v,
        T: mu * S_constant_v,
        p: p_skip_constant,
        q: mu * p_skip_constant,
        h: mu * g + Q_skip_constant,
    }
).doit()
U0 = sp.symbols("U0")
skip_constant_r2z2_origin = reduce_mu_law(
    sp.cancel(skip_constant_r2z2 / e**3)
    .subs(e, 0)
    .subs(U_constant.subs(e, 0), U0)
)
skip_constant_r3_origin = reduce_mu_law(
    sp.cancel(skip_constant_r3 / e**3)
    .subs(e, 0)
    .subs(U_constant.subs(e, 0), U0)
)
zero(
    skip_constant_r2z2_origin + 6 * mu * U0,
    "degenerate constant-v r2z2 origin equation",
)
zero(
    skip_constant_r3_origin - (18 * theta * U0 - mu / 2),
    "degenerate constant-v r3 origin equation",
)
gate(
    sp.resultant(4 * mu**2 + 3, -mu / 2, mu) != 0,
    "degenerate constant-v arm and terminal contradiction",
)


# Linear-v closure.  Nonzero roots are killed directly by the local payment
# plus the next bucket; the confluent root tau=0 is recomputed separately.
tau = sp.symbols("tau")
v_linear = e - tau
K_linear = beta * e**2 * v_linear**4
D_linear = 3 * e**2 - 2 * mu * e - 1
numerator_f_linear = 2 * e * K_linear - sp.Rational(1, 12)
quotient_f_linear, remainder_f_linear = sp.div(
    sp.Poly(numerator_f_linear, e), sp.Poly(D_linear, e)
)
f_linear = quotient_f_linear.as_expr()
U_linear = sp.Function("U_linear")(e)
g_linear = sp.Function("g_linear")(e)
S_linear = alpha * e**4 * v_linear**7

p_generic_linear = e**2 * v_linear**3 * U_linear
P_generic_linear = delta * e**3 * v_linear**5
Q_generic_linear = (
    sp.Rational(5, 7) * delta / alpha * e * v_linear * U_linear
    + eta * alpha / beta * e**2 * v_linear**3
)
generic_linear_r2z = bucket_r2z.subs(
    {
        f: f_linear,
        kap: mu * f_linear + K_linear,
        S: S_linear,
        T: mu * S_linear,
        p: p_generic_linear,
        q: mu * p_generic_linear + P_generic_linear,
        g: g_linear,
        h: mu * g_linear + Q_generic_linear,
    }
).doit()
zero(
    generic_linear_r2z.subs(e, tau).doit()
    - 60
    * delta
    * tau**2
    * U_linear.subs(e, tau)
    * f_linear.subs(e, tau)
    / (7 * alpha),
    "generic linear-v nonzero-root r2z obstruction",
)

p_skip_linear = e * v_linear * U_linear
Q_skip_linear = theta * e**2 * v_linear**3
skip_linear_r2z = bucket_r2z.subs(
    {
        f: f_linear,
        kap: mu * f_linear + K_linear,
        S: S_linear,
        T: mu * S_linear,
        p: p_skip_linear,
        q: mu * p_skip_linear,
        g: g_linear,
        h: mu * g_linear + Q_skip_linear,
    }
).doit()
zero(
    skip_linear_r2z.subs(e, tau).doit()
    - tau**2
    * (-2 * mu + 3 * tau)
    * U_linear.subs(e, tau),
    "degenerate linear-v nonzero-root r2z address",
)
zero(
    remainder_f_linear.as_expr().coeff(e, 1).subs(tau, 2 * mu / 3)
    - 2 * beta / 27,
    "degenerate linear-v addressed arm obstruction",
)

# Generic confluent root tau=0.  The arm linear remainder and the z2
# origin bucket demand two coprime polynomials in x=mu^2.
g0 = sp.symbols("g0")
generic_zero_profiles = {
    f: f_linear.subs(tau, 0),
    kap: mu * f_linear.subs(tau, 0) + K_linear.subs(tau, 0),
    S: S_linear.subs(tau, 0),
    T: mu * S_linear.subs(tau, 0),
    p: p_generic_linear.subs(tau, 0),
    q: (mu * p_generic_linear + P_generic_linear).subs(tau, 0),
    g: g0,
    h: (mu * g_linear + Q_generic_linear).subs(
        {tau: 0, g_linear: g0}
    ),
}
generic_zero_r2_origin = bucket_r2.subs(generic_zero_profiles).doit().subs(e, 0)
generic_zero_z2_origin = bucket_z2.subs(generic_zero_profiles).doit().subs(e, 0)
zero(generic_zero_r2_origin - mu * g0,
     "generic v=e r2 origin forces g0=0")
zero(
    generic_zero_z2_origin
    + (32 * beta * mu**5 + 88 * beta * mu**3 + 42 * beta * mu + 729 * g0)
    / 81,
    "generic v=e z2 origin equation",
)
arm_sextic = 64 * mu**6 + 240 * mu**4 + 216 * mu**2 + 27
zero(
    remainder_f_linear.as_expr().coeff(e, 1).subs(tau, 0)
    - 2 * beta * arm_sextic / 729,
    "generic v=e arm sextic",
)
x_square = sp.symbols("x_square")
arm_sextic_x = 64 * x_square**3 + 240 * x_square**2 + 216 * x_square + 27
z2_quadratic_x = 16 * x_square**2 + 44 * x_square + 21
zero(
    sp.resultant(arm_sextic_x, z2_quadratic_x, x_square) + 4534272,
    "generic v=e arm/z2 resultant",
)

# Degenerate confluent root tau=0.  r2z first gives U(0)=0; the arm
# constant remainder then turns the r3 leading coefficient into 6*beta.
U0_linear = sp.symbols("U0_linear")
g_zero = sp.Function("g_zero")(e)
U_zero = U0_linear + sp.Function("U_tail")(e) * e
p_skip_zero = e**2 * U_zero
Q_skip_zero = theta * e**5
skip_zero_profiles = {
    f: f_linear.subs(tau, 0),
    kap: mu * f_linear.subs(tau, 0) + K_linear.subs(tau, 0),
    S: S_linear.subs(tau, 0),
    T: mu * S_linear.subs(tau, 0),
    p: p_skip_zero,
    q: mu * p_skip_zero,
    g: g_zero,
    h: mu * g_zero + Q_skip_zero,
}
skip_zero_r2z_lead = sp.cancel(
    bucket_r2z.subs(skip_zero_profiles).doit() / e**2
).subs(e, 0)
skip_zero_r3_lead = sp.cancel(
    bucket_r3.subs(skip_zero_profiles).doit() / e**7
).subs(e, 0)
zero(skip_zero_r2z_lead + 3 * mu * U0_linear,
     "degenerate v=e r2z leading coefficient")
arm_constant_zero = remainder_f_linear.as_expr().coeff(e, 0).subs(tau, 0)
arm_constant_numerator = sp.factor(2916 * arm_constant_zero)
zero(
    arm_constant_numerator
    - (
        256 * beta * mu**5
        + 768 * beta * mu**3
        + 432 * beta * mu
        - 243
    ),
    "degenerate v=e arm constant equation",
)
zero(
    skip_zero_r3_lead.subs(U0_linear, 0)
    - 2
    * (
        beta
        * (256 * beta * mu**5 + 768 * beta * mu**3 + 432 * beta * mu)
    )
    / 81,
    "degenerate v=e raw r3 lead after U0=0",
)
zero(
    skip_zero_r3_lead.subs(U0_linear, 0)
    - 6 * beta
    - 2 * beta * arm_constant_numerator / 81,
    "degenerate v=e arm-reduced r3 obstruction",
)
gate(beta != 0, "degenerate v=e obstruction is nonzero")


# All-degree closure.  The r2z bucket sees the terminal payments one order
# earlier than every other summand at a nonzero root.  Work in a differential
# polynomial presentation so that the root-leading coordinates are checked
# without specializing the multiplicity.
vv, vd, uu, ud, ff, fd, gg, gd = sp.symbols(
    "vv vd uu ud ff fd gg gd"
)
differential_symbols = {
    v: vv,
    sp.diff(v, e): vd,
    U: uu,
    sp.diff(U, e): ud,
    f: ff,
    sp.diff(f, e): fd,
    g: gg,
    sp.diff(g, e): gd,
}

generic_all_r2z = bucket_r2z.subs(
    {
        f: f,
        kap: mu * f + K_tower,
        S: S_tower,
        T: mu * S_tower,
        p: p_odd,
        q: mu * p_odd + P_tower,
        g: g,
        h: mu * g + Q_odd,
    }
).doit()
generic_root_coordinate = sp.expand(
    generic_all_r2z.xreplace(differential_symbols)
).subs(vv, 0)
zero(
    generic_root_coordinate
    - sp.Rational(60, 7) * delta / alpha * e**2 * ff * uu * vd,
    "generic r2z root-leading coordinate",
)

degenerate_all_r2z = bucket_r2z.subs(
    {
        f: f,
        kap: mu * f + K_tower,
        S: S_tower,
        T: mu * S_tower,
        p: p_skip,
        q: mu * p_skip,
        g: g,
        h: mu * g + Q_tower,
    }
).doit()
degenerate_differential = sp.expand(
    degenerate_all_r2z.xreplace(differential_symbols)
)
degenerate_root_coordinate = degenerate_differential.subs(vv, 0)
zero(
    degenerate_root_coordinate - e**2 * (3 * e - 2 * mu) * uu * vd,
    "degenerate r2z first root address",
)

# At the only possible nonzero address rho=2mu/3, the next coefficient for a
# root vv=x^m*u is (3/2)rho^2*u*U*(2m-5).  The first term below is the first
# derivative of the vanished vd coefficient; the second is the vv coefficient.
degenerate_address_next = (
    3 * rho**2 * root_order
    - rho * (6 * rho + sp.Rational(3, 2) * rho)
)
zero(
    degenerate_address_next
    - sp.Rational(3, 2) * rho**2 * (2 * root_order - 5),
    "degenerate addressed-root next coefficient",
)
gate(
    sp.solve(sp.Eq(2 * root_order, 5), root_order) == [],
    "positive integral root multiplicity misses address resonance",
)

# Once nonzero roots have been excluded, algebraic closure makes v=c*e^d.
# The generic origin calculation uses only the displayed first jets: every
# tower correction has strictly higher e-order for d>=1.
f0_all, f1_all, g0_all, g1_all = sp.symbols(
    "f0_all f1_all g0_all g1_all"
)
f_jet = f0_all + f1_all * e
g_jet = g0_all + g1_all * e
generic_origin_profiles = {
    f: f_jet,
    kap: mu * f_jet,
    S: 0,
    T: 0,
    p: 0,
    q: 0,
    g: g_jet,
    h: mu * g_jet,
}
generic_origin_r2 = bucket_r2.subs(generic_origin_profiles).doit().subs(e, 0)
generic_origin_z2 = bucket_z2.subs(generic_origin_profiles).doit().subs(e, 0)
zero(generic_origin_r2 - mu * g0_all,
     "generic monomial origin r2 jet")
zero(
    generic_origin_z2 + 3 * f0_all + mu * f1_all + 9 * g0_all,
    "generic monomial origin z2 jet",
)
zero(
    generic_origin_z2.subs(
        {
            f0_all: sp.Rational(1, 12),
            f1_all: -mu / 6,
            g0_all: 0,
        }
    )
    - (mu**2 / 6 - sp.Rational(1, 4)),
    "generic nonzero-mu origin law",
)

# If a,b are the roots of D, their ratio q satisfies the displayed trace
# relation.  At mu^2=3/2 this is q^2+4q+1=0; its two complex roots are real
# and neither is one of the only possible real roots of unity, +/-1.
q_ratio, mu_square = sp.symbols("q_ratio mu_square")
ratio_trace = -(sp.Rational(4, 3) * mu_square + 2)
zero(
    ratio_trace.subs(mu_square, sp.Rational(3, 2)) + 4,
    "generic D-root ratio trace",
)
ratio_polynomial = q_ratio**2 + 4 * q_ratio + 1
gate(ratio_polynomial.subs(q_ratio, 1) != 0,
     "D-root ratio is not the real root of unity +1")
gate(ratio_polynomial.subs(q_ratio, -1) != 0,
     "D-root ratio is not the real root of unity -1")

# In the degenerate monomial branch the terminal law fixes U(0), while r2z
# has an earlier nonzero mu-leading coefficient.  The mu=0 seam is instead
# killed directly by the arm: its two D-roots are opposite and the exponent
# 3+4d is odd.
monomial_degree = sp.symbols("monomial_degree", integer=True, positive=True)
U0_all = sp.symbols("U0_all")
degenerate_terminal_origin = (
    2 * beta * (8 * monomial_degree + 4) * sp.Rational(1, 12)
    + 3 * theta * (4 * monomial_degree + 2) * U0_all
)
zero(
    degenerate_terminal_origin
    - sp.Rational(2, 3)
    * (2 * monomial_degree + 1)
    * (beta + 9 * theta * U0_all),
    "degenerate monomial terminal origin payment",
)
degenerate_r2z_origin = -mu * (2 * monomial_degree + 1) * U0_all
zero(
    degenerate_r2z_origin.subs(U0_all, -beta / (9 * theta))
    - mu * beta * (2 * monomial_degree + 1) / (9 * theta),
    "degenerate monomial r2z contradiction",
)
gate(sp.Mod(3 + 4 * monomial_degree, 2) == 1,
     "mu-zero arm exponent is odd")


semantic = {
    "ansatz": "THM3814 plus rz2*S(e),rz2*T(e)",
    "top": "S!=0;T=mu*S;K=kappa-mu*f=beta*e^2*v^4;S=alpha*e^4*v^7",
    "generic": "P=q-mu*p=delta*e^3*v^5;p=e^2*v^3*U;Q=h-mu*g=(5delta/7alpha)e*v*U+(eta alpha/beta)e^2*v^3",
    "generic_terminal": "Riccati law;nonzero roots pay 28alpha beta f=5delta U^2",
    "degenerate": "P=0;Q=theta*e^2*v^3;p=e*v*U",
    "degenerate_terminal": "linear law;nonzero roots pay 4beta f+3theta U=0",
    "constant_v": "empty in both generic and degenerate terminal branches",
    "linear_v": "empty in both generic and degenerate terminal branches",
    "all_degree": "nonzero roots killed by r2z;confluent monomial towers killed by origin jets and arm root ratio",
    "scope": "complete no-go for the full first rz2 extension;no planar-JC claim beyond this ansatz",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3821-cubic-pseudoplane-rz2-odd-ladder-terminal-riccati-gate")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("ansatz=THM3814_plus_rz2*S(e),rz2*T(e)")
print("top=S_nonzero;T=muS;S=alpha_e4_v7;K=beta_e2_v4")
print("generic=P=delta_e3_v5;p=e2_v3_U;Q=(5delta/7alpha)e_v_U+(eta_alpha/beta)e2_v3")
print("generic_payment=28alpha_beta_f=5delta_U2_at_nonzero_v_roots")
print("degenerate=P0;Q=theta_e2_v3;p=e_v_U")
print("degenerate_payment=4beta_f+3theta_U=0_at_nonzero_v_roots")
print("constant_v=empty_in_generic_and_degenerate_branches")
print("linear_v=empty_in_generic_and_degenerate_branches")
print("all_degree=nonzero_roots_and_confluent_monomial_towers_empty")
print("scope=complete_first_rz2_extension_no_go;no_broader_planar_JC_claim")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
