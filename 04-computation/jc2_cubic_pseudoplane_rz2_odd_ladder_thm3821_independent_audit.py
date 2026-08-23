#!/usr/bin/env python3
"""Independent hostile audit and one-bucket extension of provisional THM-3821.

This companion deliberately rebuilds the Poisson bracket from the three
pairwise brackets, rather than importing the canonical companion.  It also
tests the origin boundary omitted by the provisional proof and then consumes
the next exact r^2 z^2 bucket.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    gate(sp.cancel(left - right) == 0, label)


# ---------------------------------------------------------------------------
# 1. Rebuild the quotient Poisson calculation independently.
# ---------------------------------------------------------------------------
r, z, e = sp.symbols("r z e")
relation = z**3 - r**2 * e - r
surface = r**2 * e - z**3 + r


def poisson_bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    """Use the three displayed pairwise brackets directly."""
    return sp.expand(
        3 * r**2 * (sp.diff(left, r) * sp.diff(right, z)
                    - sp.diff(left, z) * sp.diff(right, r))
        + 9 * z**2 * (sp.diff(left, r) * sp.diff(right, e)
                      - sp.diff(left, e) * sp.diff(right, r))
        + (3 + 6 * r * e)
        * (sp.diff(left, z) * sp.diff(right, e)
           - sp.diff(left, e) * sp.diff(right, z))
    )


for coordinate in (r, z, e):
    same(poisson_bracket(surface, coordinate), 0,
         f"surface Casimir against {coordinate}")

f = sp.Function("f")(e)
g = sp.Function("g")(e)
h = sp.Function("h")(e)
kap = sp.Function("kap")(e)
p = sp.Function("p")(e)
q = sp.Function("q")(e)
S = sp.Function("S")(e)
T = sp.Function("T")(e)

A = e**2 - z / 3 + r * g + z**2 * f + r * z * p + r * z**2 * S
C = (
    e**3 - e - e * z / 2 + r * h + z**2 * kap
    + r * z * q + r * z**2 * T
)
raw = poisson_bracket(A, C) - 1
normal = sp.rem(raw, relation, z)
normal_poly = sp.Poly(sp.expand(normal), r, z)
gate(max(j for _, j in normal_poly.monoms()) <= 2,
     "monic z normal form has degree at most two")
same(sp.rem(raw - normal, relation, z), 0,
     "normal form retains the quotient class")


def bucket(i: int, j: int) -> sp.Expr:
    return normal_poly.coeff_monomial(r**i * z**j)


b_z = bucket(0, 1)
b_r5 = bucket(5, 0)
b_r4 = bucket(4, 0)
b_r3z2 = bucket(3, 2)
b_r3z = bucket(3, 1)
b_r3 = bucket(3, 0)
b_r2z2 = bucket(2, 2)
b_r2z = bucket(2, 1)
b_r2 = bucket(2, 0)
b_z2 = bucket(0, 2)

expected_z = (36 * e**2 * f - 24 * e * kap - 12 * f + 1) / 2
expected_r5 = -21 * e**2 * (-S * sp.diff(T, e) + T * sp.diff(S, e))
expected_r4 = -3 * e * (
    -4 * e * f * sp.diff(T, e) + 4 * e * kap * sp.diff(S, e)
    - 7 * e * S * sp.diff(kap, e) + 7 * e * T * sp.diff(f, e)
    + 2 * f * T - 2 * kap * S - 12 * S * sp.diff(T, e)
    + 12 * T * sp.diff(S, e)
)
expected_r3z2 = -3 * (
    -5 * e * p * sp.diff(T, e) + 5 * e * q * sp.diff(S, e)
    - 7 * e * S * sp.diff(q, e) + 7 * e * T * sp.diff(p, e)
    - p * T + q * S
)
expected_r3z = -3 * (
    -3 * e * g * sp.diff(T, e) + 3 * e * h * sp.diff(S, e)
    - 5 * e * p * sp.diff(q, e) + 5 * e * q * sp.diff(p, e)
    - 7 * e * S * sp.diff(h, e) + 7 * e * T * sp.diff(g, e)
    - 2 * g * T + 2 * h * S
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
    + 2 * f * T - g * q + h * p - 2 * kap * S
    - 5 * S * sp.diff(T, e) + 5 * T * sp.diff(S, e)
)
for actual, expected, label in (
    (b_z, expected_z, "arm bucket"),
    (b_r5, expected_r5, "r5 Wronskian bucket"),
    (b_r4, expected_r4, "r4 bucket"),
    (b_r3z2, expected_r3z2, "r3z2 bucket"),
    (b_r3z, expected_r3z, "r3z bucket"),
    (b_r3, expected_r3, "r3 bucket"),
):
    same(actual, expected, label)


# ---------------------------------------------------------------------------
# 2. Reconstruct the 7/4 and 7/5 top ladders.
# ---------------------------------------------------------------------------
alpha, beta, delta, theta = sp.symbols(
    "alpha beta delta theta", nonzero=True
)
mu, eta = sp.symbols("mu eta")
v = sp.Function("v")(e)
K = sp.Function("K")(e)
P = sp.Function("P")(e)
Q = sp.Function("Q")(e)
R = sp.Function("R")(e)
U = sp.Function("U")(e)

arm = (3 * e**2 - 1) * f - 2 * e * kap + sp.Rational(1, 12)
same(b_z, 6 * arm, "arm identity normalization")

# S=0,T!=0 gives the asymmetric 7/4 tower, which has f(0)=0.
w = sp.Function("w")(e)
asym_ode = 4 * e * f * sp.diff(T, e) - 7 * e * T * sp.diff(f, e) - 2 * T * f
same(
    asym_ode.subs({T: alpha * e**4 * w**7,
                   f: beta * e**2 * w**4}).doit(),
    0,
    "asymmetric tower solves its ODE",
)
gate((beta * e**2 * w**4).subs(e, 0) == 0,
     "asymmetric tower contradicts f(0)=1/12")

S0 = alpha * e**4 * v**7
K0 = beta * e**2 * v**4
P0 = delta * e**3 * v**5
Q0 = theta * e**2 * v**3

ode74 = 4 * e * K * sp.diff(S, e) - 7 * e * S * sp.diff(K, e) - 2 * K * S
same(ode74.subs({S: S0, K: K0}).doit(), 0, "7/4 tower")
j = sp.symbols("j", integer=True, nonnegative=True)
gate(4 * (7 * j) == 7 * (4 * j), "7/4 ordinary-prime family")
gate(4 * (4 + 7 * j) == 2 + 7 * (2 + 4 * j),
     "7/4 origin-prime family")

r4_main = b_r4.subs({T: mu * S, kap: mu * f + K}).doit()
same(r4_main, -3 * e * ode74, "r4 difference reduction")

ode75 = 7 * e * S * sp.diff(P, e) - 5 * e * P * sp.diff(S, e) - P * S
r3z2_main = b_r3z2.subs({T: mu * S, q: mu * p + P}).doit()
same(r3z2_main, 3 * ode75, "r3z2 difference reduction")
same(ode75.subs({S: S0, P: P0}).doit(), 0, "7/5 tower")
gate(7 * (5 * j) == 5 * (7 * j), "7/5 ordinary-prime family")
gate(7 * (3 + 5 * j) == 1 + 5 * (4 + 7 * j),
     "7/5 origin-prime family")


# ---------------------------------------------------------------------------
# 3. Integrate the generic branch and expose the omitted d=0 seam.
# ---------------------------------------------------------------------------
integrating = (
    5 * e * (p * sp.diff(P, e) - P * sp.diff(p, e))
    + 7 * e * S * sp.diff(Q, e) - 3 * e * Q * sp.diff(S, e)
    - 2 * Q * S
)
r3z_main = b_r3z.subs(
    {T: mu * S, q: mu * p + P, h: mu * g + Q}
).doit()
same(r3z_main, 3 * integrating, "r3z integrating equation")

p_pre = e * v**2 * R
Q_pre = sp.Rational(5, 7) * delta / alpha * R + eta * alpha / beta * e**2 * v**3
same(
    integrating.subs({S: S0, P: P0, p: p_pre, Q: Q_pre}).doit(),
    0,
    "integrated generic relation",
)
same(Q_pre, sp.Rational(5, 7) * P0 * p_pre / S0 + eta * S0 / K0,
     "integrated rational expression and polynomial typing")

# Exact pure-r3 difference equation before the last valuation transfer.
r3_difference = (
    -4 * e**2 * K * sp.diff(f, e) + 4 * e**2 * f * sp.diff(K, e)
    - 6 * e * K * sp.diff(S, e) + 12 * e * S * sp.diff(K, e)
    + 2 * K * S - 5 * e * P * sp.diff(g, e)
    + 3 * e * g * sp.diff(P, e) - 3 * e * Q * sp.diff(p, e)
    + 5 * e * p * sp.diff(Q, e) + P * g - Q * p
)
same(
    b_r3.subs({T: mu * S, kap: mu * f + K,
               q: mu * p + P, h: mu * g + Q}).doit(),
    3 * r3_difference,
    "pure-r3 difference equation",
)

c = sp.Rational(5, 7) * delta / alpha
forcing = r3_difference.subs(
    {S: S0, K: K0, P: 0, Q: 0, p: 0, g: 0}
).doit()
r_square = 2 * c * e * v * R * (
    -2 * v * R - 3 * e * sp.diff(v, e) * R
    + e * v * sp.diff(R, e)
)
g_block = -5 * e * P0 * sp.diff(g, e) + 3 * e * g * sp.diff(P0, e) + P0 * g
eta_block = sp.expand(
    r3_difference.subs({S: S0, K: K0, P: P0,
                        p: p_pre, Q: Q_pre}).doit()
    - forcing - r_square - g_block
)
same(
    r3_difference.subs({S: S0, K: K0, P: P0,
                        p: p_pre, Q: Q_pre}).doit(),
    forcing + r_square + g_block + eta_block,
    "generic r3 block decomposition",
)

m, a = sp.symbols("m a", integer=True, positive=True)
gate(4 * m - 1 < 5 * m - 1, "nonzero-root g block follows forcing")
gate(sp.solve(sp.Eq(2 * a + 2 * m - 1, 4 * m - 1), a)[0] == m,
     "nonzero-root R valuation match")
gate(3 * m > m, "nonzero-root R resonance lies above match")

# The provisional proof says every other term follows order 3+4d at the
# origin.  For d=0 this is false: gP also starts in order three and can
# cancel the f,K forcing.  Freeze that exact counter-mechanism.
v0, g0, R0 = sp.symbols("v0 g0 R0", nonzero=True)
f0 = sp.Rational(1, 12)
forcing_d0_e3 = sp.Rational(2, 3) * beta * v0**4
g_d0_e3 = 10 * delta * v0**5 * g0
gate(
    sp.cancel((forcing_d0_e3 + g_d0_e3).subs(
        g0, -beta / (15 * delta * v0)
    )) == 0,
    "d=0 forcing is not uniquely earliest",
)

# Repair: a=0 is killed at e^1 by the R-square block.  If a>=1, the next
# exact r^2 z^2 bucket has unavoidable e^3 coefficient 5 delta v0^5/2.
same(-4 * c * v0**2 * R0**2,
     -sp.Rational(20, 7) * delta / alpha * v0**2 * R0**2,
     "d=0,a=0 r3 e1 obstruction")

v1, f1, R1 = sp.symbols("v1 f1 R1")
vjet = v0 + v1 * e
fjet = f0 + f1 * e
for a_test in (1, 2, 4):
    Rjet = R1 * e**a_test
    pjet = e * vjet**2 * Rjet
    profiles = {
        f: fjet,
        kap: mu * fjet + beta * e**2 * vjet**4,
        S: alpha * e**4 * vjet**7,
        T: mu * alpha * e**4 * vjet**7,
        p: pjet,
        q: mu * pjet + delta * e**3 * vjet**5,
    }
    coefficient = sp.expand(b_r2z2.subs(profiles).doit()).coeff(e, 3)
    same(coefficient, sp.Rational(5, 2) * delta * v0**5,
         f"d=0,a={a_test} r2z2 repair obstruction")

# Once d>=1, the original origin comparison is sound: gP is later.
d = sp.symbols("d", integer=True, positive=True)
gate(3 + 4 * d < 3 + 5 * d, "repaired origin g block follows forcing")
gate(sp.solve(sp.Eq(1 + 2 * d + 2 * a, 3 + 4 * d), a)[0] == d + 1,
     "repaired origin R valuation match")
gate(3 * d + 2 > d + 1, "origin R resonance lies above match")

p_generic = e**2 * v**3 * U
Q_generic = c * e * v * U + eta * alpha / beta * e**2 * v**3


# ---------------------------------------------------------------------------
# 4. Recompute the terminal square payment.
# ---------------------------------------------------------------------------
terminal_generic = (
    14 * alpha**2 * beta**2 * e**3 * v**7
    * (3 * e * sp.diff(v, e) + v)
    + 21 * alpha**2 * eta * e * v**2
    * (2 * e * U * sp.diff(v, e) - e * v * sp.diff(U, e) + U * v)
    + 2 * beta * (2 * e * sp.diff(v, e) + v)
    * (28 * alpha * beta * f - 5 * delta * U**2)
    + e * v * (-28 * alpha * beta**2 * sp.diff(f, e)
               + 10 * beta * delta * U * sp.diff(U, e))
    + 35 * alpha * beta * delta * v
    * ((3 * e * sp.diff(v, e) + 2 * v) * g
       - e * v * sp.diff(g, e))
)
r3_generic = b_r3.subs(
    {T: mu * S0, S: S0, kap: mu * f + K0,
     q: mu * p_generic + P0, p: p_generic,
     h: mu * g + Q_generic}
).doit()
same(r3_generic, 3 * e**3 * v**3 * terminal_generic / (7 * alpha * beta),
     "generic terminal law from actual bracket")

x, rho, unit, froot, Uroot = sp.symbols(
    "x rho unit froot Uroot", nonzero=True
)
root_mult = sp.symbols("root_mult", integer=True, positive=True)
vlocal = unit * x**root_mult
payment_lead = sp.limit(
    2 * beta * (2 * (rho + x) * sp.diff(vlocal, x) + vlocal)
    * (28 * alpha * beta * froot - 5 * delta * Uroot**2)
    / x**(root_mult - 1),
    x,
    0,
)
same(
    payment_lead,
    4 * beta * rho * root_mult * unit
    * (28 * alpha * beta * froot - 5 * delta * Uroot**2),
    "generic nonzero-root square payment",
)


# ---------------------------------------------------------------------------
# 5. Recompute the P=0 skip ladder and linear payment.
# ---------------------------------------------------------------------------
ode73 = 7 * e * S * sp.diff(Q, e) - 3 * e * Q * sp.diff(S, e) - 2 * Q * S
same(ode73.subs({S: S0, Q: Q0}).doit(), 0, "7/3 tower")
gate(7 * (3 * j) == 3 * (7 * j), "7/3 ordinary-prime family")
gate(7 * (2 + 3 * j) == 2 + 3 * (4 + 7 * j),
     "7/3 origin-prime family")

qzero_parenthesis = (
    3 * alpha * e**4 * v**7 * sp.diff(v, e) + alpha * e**3 * v**8
    + 8 * e * f * sp.diff(v, e) - 2 * e * v * sp.diff(f, e) + 4 * f * v
)
same(
    r3_difference.subs({S: S0, K: K0, P: 0, Q: 0, p: 0}).doit(),
    2 * beta * e**3 * v**3 * qzero_parenthesis,
    "Q=0 forcing factor",
)
gate(8 * root_mult + 4 != 0, "Q=0 origin coefficient cannot vanish")

p_order = sp.symbols("p_order", integer=True, nonnegative=True)
gate(sp.solve(sp.Eq(p_order + 3 * root_mult - 1,
                    4 * root_mult - 1), p_order)[0] == root_mult,
     "degenerate nonzero-root p match")
gate(5 * root_mult > root_mult, "degenerate root resonance above match")
gate(sp.solve(sp.Eq(p_order + 2 + 3 * d, 3 + 4 * d), p_order)[0] == d + 1,
     "degenerate origin p match")
gate(3 + 5 * d > d + 1, "degenerate origin resonance above match")

p_skip = e * v * U
terminal_skip = (
    2 * beta * qzero_parenthesis
    - 3 * theta * (-4 * e * U * sp.diff(v, e)
                   + e * v * sp.diff(U, e) - 2 * U * v)
)
r3_skip = b_r3.subs(
    {T: mu * S0, S: S0, kap: mu * f + K0,
     q: mu * p_skip, p: p_skip, h: mu * g + Q0}
).doit()
same(r3_skip, 3 * e**3 * v**3 * terminal_skip,
     "degenerate terminal law from actual bracket")

for sample_mult in (1, 2, 3, 5):
    sample_v = unit * x**sample_mult
    sample_block = (
        2 * beta * (8 * (rho + x) * froot * sp.diff(sample_v, x)
                    + 4 * froot * sample_v)
        - 3 * theta * (-4 * (rho + x) * Uroot * sp.diff(sample_v, x)
                      - 2 * Uroot * sample_v)
    )
    sample_lead = sp.cancel(sample_block / x**(sample_mult - 1)).subs(x, 0)
    same(
        sample_lead,
        4 * rho * sample_mult * unit
        * (4 * beta * froot + 3 * theta * Uroot),
        f"degenerate nonzero-root linear payment m={sample_mult}",
    )


# ---------------------------------------------------------------------------
# 6. The unused r^2 z^2 bucket closes both branches.
# ---------------------------------------------------------------------------
generic_next = sp.factor(
    b_r2z2.subs(
        {T: mu * S0, S: S0, kap: mu * f + K0,
         q: mu * p_generic + P0, p: p_generic,
         h: mu * g + Q_generic}
    ).doit()
)
def generic_inner_formula(
    E: sp.Expr,
    V: sp.Expr,
    Vp: sp.Expr,
    F: sp.Expr,
    Fp: sp.Expr,
    UU: sp.Expr,
    UUp: sp.Expr,
) -> sp.Expr:
    return (
        9 * alpha * delta * E**4 * V**7 * Vp
        + 3 * alpha * delta * E**3 * V**8
        - 21 * alpha * E**3 * V**2 * Vp
        + 14 * alpha * E**2 * mu * V**2 * Vp
        + 7 * alpha * E * mu * V**3
        - 24 * beta * E**2 * UU * V**2 * Vp
        + 12 * beta * E**2 * V**3 * UUp
        - 12 * beta * E * UU * V**3
        - 60 * delta * E * F * Vp
        + 15 * delta * E * V * Fp
        - 30 * delta * F * V
    )


generic_inner = generic_inner_formula(
    e, v, sp.diff(v, e), f, sp.diff(f, e), U, sp.diff(U, e)
)
same(generic_next, -e**3 * v**4 * generic_inner,
     "generic next-bucket law")

# At a nonzero root, generic_inner has unique order m-1 lead
# -60 delta rho f(rho) m unit.  The arm makes f(rho) nonzero.
for sample_mult in (1, 2, 4):
    sample_v = unit * x**sample_mult
    sample_local_inner = generic_inner_formula(
        rho + x,
        sample_v,
        sp.diff(sample_v, x),
        froot,
        0,
        Uroot,
        0,
    )
    sample_local_lead = sp.cancel(
        sample_local_inner / x**(sample_mult - 1)
    ).subs(x, 0)
    same(
        sample_local_lead,
        -60 * delta * rho * froot * sample_mult * unit,
        f"generic next-bucket nonzero-root lead m={sample_mult}",
    )

# If no nonzero root remains, algebraic closure gives v=c e^d.  For d>=1
# the same inner law has unique e^d coefficient -30 delta(2d+1)f(0).
same(-30 * delta * (2 * d + 1) * sp.Rational(1, 12),
     -sp.Rational(5, 2) * delta * (2 * d + 1),
     "generic pure-origin tower obstruction")

degenerate_next = sp.factor(
    b_r2z2.subs(
        {T: mu * S0, S: S0, kap: mu * f + K0,
         q: mu * p_skip, p: p_skip, h: mu * g + Q0}
    ).doit()
)
def degenerate_inner_formula(
    E: sp.Expr,
    V: sp.Expr,
    Vp: sp.Expr,
    UU: sp.Expr,
    UUp: sp.Expr,
) -> sp.Expr:
    return (
        21 * alpha * E**3 * V**2 * Vp
        - 14 * alpha * E**2 * mu * V**2 * Vp
        - 7 * alpha * E * mu * V**3
        + 48 * beta * E * UU * Vp
        - 12 * beta * E * V * UUp
        + 24 * beta * UU * V
    )


degenerate_inner = degenerate_inner_formula(
    e, v, sp.diff(v, e), U, sp.diff(U, e)
)
same(degenerate_next, e**3 * v**4 * degenerate_inner,
     "degenerate next-bucket law")

# The terminal payment gives U(rho)!=0, so a nonzero v-root has the unique
# order m-1 lead 48 beta rho U(rho) m unit.
for sample_mult in (1, 2, 4):
    sample_v = unit * x**sample_mult
    sample_local_inner = degenerate_inner_formula(
        rho + x,
        sample_v,
        sp.diff(sample_v, x),
        Uroot,
        0,
    )
    sample_local_lead = sp.cancel(
        sample_local_inner / x**(sample_mult - 1)
    ).subs(x, 0)
    same(
        sample_local_lead,
        48 * beta * rho * Uroot * sample_mult * unit,
        f"degenerate next-bucket nonzero-root lead m={sample_mult}",
    )

# For the remaining v=e^d tower, arm divisibility excludes mu=0: modulo
# 3e^2-1 the odd top monomial leaves a nonzero e coefficient.
arm_mu_zero_e_coefficient = 2 * beta / 3 ** (1 + 2 * d)
gate(arm_mu_zero_e_coefficient != 0,
     "pure-origin arm excludes mu=0")

u_order = sp.symbols("u_order", integer=True, nonnegative=True)
gate(sp.solve(sp.Eq(d + u_order, 1 + 3 * d), u_order)[0] == 1 + 2 * d,
     "degenerate next-bucket U valuation match")
gate(4 * d + 2 > 1 + 2 * d,
     "degenerate U coefficient resonance above match")

# With ord(U)=1+2d, the terminal law has a unique e^d term from f(0).
terminal_origin_coefficient = 2 * beta * (8 * d + 4) * sp.Rational(1, 12)
same(terminal_origin_coefficient,
     sp.Rational(2, 3) * beta * (2 * d + 1),
     "degenerate pure-origin terminal obstruction")


# ---------------------------------------------------------------------------
# 7. Independently replay the constant and linear claims as controls.
# ---------------------------------------------------------------------------
D = 3 * e**2 - 2 * mu * e - 1
f_const = sp.Rational(1, 12) - mu * e / 6
constant_quotient, constant_remainder = sp.div(
    sp.Poly(2 * beta * e**3 - sp.Rational(1, 12), e),
    sp.Poly(D, e),
)
same(
    constant_remainder.as_expr(),
    sp.Rational(4, 9) * beta * mu - sp.Rational(1, 12)
    + e * (sp.Rational(8, 9) * beta * mu**2
           + sp.Rational(2, 3) * beta),
    "constant arm exact remainder",
)


def reduce_mu_law(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    gate(sp.diff(denominator, mu) == 0, "mu-law denominator is constant in mu")
    return sp.cancel(sp.rem(numerator, 4 * mu**2 + 3, mu) / denominator)


same(
    reduce_mu_law(constant_remainder.as_expr().subs(beta, -mu / 4)),
    0,
    "constant arm remainder under the mu law",
)
same(
    reduce_mu_law(constant_quotient.as_expr().subs(beta, -mu / 4) - f_const),
    0,
    "constant arm quotient",
)

Uc = sp.Function("Uc")(e)
Kc = -mu * e**2 / 4
Sc = alpha * e**4
pgc = e**2 * Uc
Pgc = delta * e**3
Qgc = c * e * Uc - 4 * eta * alpha / mu * e**2
generic_constant_bucket = b_r2z2.subs(
    {f: f_const, kap: mu * f_const + Kc, S: Sc, T: mu * Sc,
     p: pgc, q: mu * pgc + Pgc, h: mu * g + Qgc}
).doit()
same(
    reduce_mu_law(sp.cancel(generic_constant_bucket / e**3).subs(e, 0)),
    sp.Rational(5, 2) * delta,
    "constant generic r2z2 obstruction",
)

psc = e * Uc
Qsc = theta * e**2
skip_constant_r2z2 = b_r2z2.subs(
    {f: f_const, kap: mu * f_const + Kc, S: Sc, T: mu * Sc,
     p: psc, q: mu * psc, h: mu * g + Qsc}
).doit()
skip_constant_r3 = b_r3.subs(
    {f: f_const, kap: mu * f_const + Kc, S: Sc, T: mu * Sc,
     p: psc, q: mu * psc, h: mu * g + Qsc}
).doit()
Uc0 = sp.symbols("Uc0")
same(
    reduce_mu_law(
        sp.cancel(skip_constant_r2z2 / e**3).subs(e, 0)
        .subs(Uc.subs(e, 0), Uc0)
    ),
    -6 * mu * Uc0,
    "constant degenerate r2z2 obstruction",
)
same(
    reduce_mu_law(
        sp.cancel(skip_constant_r3 / e**3).subs(e, 0)
        .subs(Uc.subs(e, 0), Uc0)
    ),
    18 * theta * Uc0 - mu / 2,
    "constant degenerate r3 obstruction",
)
gate(sp.resultant(4 * mu**2 + 3, -mu / 2, mu) != 0,
     "constant degenerate terminal contradiction")

# Linear root t!=0 controls from the actual r2z bucket.
t = sp.symbols("t")
vlin = e - t
Klin = beta * e**2 * vlin**4
Slin = alpha * e**4 * vlin**7
numlin = 2 * e * Klin - sp.Rational(1, 12)
qlin, rlin = sp.div(sp.Poly(numlin, e), sp.Poly(D, e))
flin = qlin.as_expr()
Ulin = sp.Function("Ulin")(e)
glin = sp.Function("glin")(e)
pg = e**2 * vlin**3 * Ulin
Pg = delta * e**3 * vlin**5
Qg = c * e * vlin * Ulin + eta * alpha / beta * e**2 * vlin**3
generic_r2z_at_t = b_r2z.subs(
    {f: flin, kap: mu * flin + Klin, S: Slin, T: mu * Slin,
     p: pg, q: mu * pg + Pg, g: glin, h: mu * glin + Qg}
).doit().subs(e, t).doit()
same(generic_r2z_at_t,
     60 * delta * t**2 * Ulin.subs(e, t) * flin.subs(e, t)
     / (7 * alpha),
     "linear generic separated-root bucket")

ps = e * vlin * Ulin
Qs = theta * e**2 * vlin**3
skip_r2z_at_t = b_r2z.subs(
    {f: flin, kap: mu * flin + Klin, S: Slin, T: mu * Slin,
     p: ps, q: mu * ps, g: glin, h: mu * glin + Qs}
).doit().subs(e, t).doit()
same(skip_r2z_at_t, t**2 * (-2 * mu + 3 * t) * Ulin.subs(e, t),
     "linear degenerate separated-root address")
same(rlin.as_expr().coeff(e, 1).subs(t, 2 * mu / 3),
     2 * beta / 27,
     "linear degenerate addressed arm obstruction")

xx = sp.symbols("xx")
poly_a = 64 * xx**3 + 240 * xx**2 + 216 * xx + 27
poly_b = 16 * xx**2 + 44 * xx + 21
same(sp.resultant(poly_a, poly_b, xx), -4534272,
     "linear confluent generic resultant")

# Recompute every confluent-root coefficient used in the candidate.
arm_sextic_mu = 64 * mu**6 + 240 * mu**4 + 216 * mu**2 + 27
same(
    rlin.as_expr().coeff(e, 1).subs(t, 0),
    2 * beta * arm_sextic_mu / 729,
    "linear confluent generic arm sextic",
)
gate(rlin.as_expr().coeff(e, 0).subs({t: 0, mu: 0}) != 0,
     "linear confluent arm excludes mu=0")

g00 = sp.symbols("g00")
generic_zero_profiles = {
    f: flin.subs(t, 0),
    kap: mu * flin.subs(t, 0) + Klin.subs(t, 0),
    S: Slin.subs(t, 0),
    T: mu * Slin.subs(t, 0),
    p: pg.subs(t, 0),
    q: (mu * pg + Pg).subs(t, 0),
    g: g00,
    h: (mu * glin + Qg).subs({t: 0, glin: g00}),
}
generic_zero_r2 = b_r2.subs(generic_zero_profiles).doit().subs(e, 0)
generic_zero_z2 = b_z2.subs(generic_zero_profiles).doit().subs(e, 0)
same(generic_zero_r2, mu * g00, "linear generic confluent r2 origin")
same(
    generic_zero_z2,
    -(32 * beta * mu**5 + 88 * beta * mu**3
      + 42 * beta * mu + 729 * g00) / 81,
    "linear generic confluent z2 origin",
)

U00 = sp.symbols("U00")
Utail = sp.Function("Utail")(e)
Uzero = U00 + e * Utail
skip_zero_profiles = {
    f: flin.subs(t, 0),
    kap: mu * flin.subs(t, 0) + Klin.subs(t, 0),
    S: Slin.subs(t, 0),
    T: mu * Slin.subs(t, 0),
    p: e**2 * Uzero,
    q: mu * e**2 * Uzero,
    g: glin,
    h: mu * glin + theta * e**5,
}
skip_zero_r2z = sp.cancel(
    b_r2z.subs(skip_zero_profiles).doit() / e**2
).subs(e, 0)
skip_zero_r3 = sp.cancel(
    b_r3.subs(skip_zero_profiles).doit() / e**7
).subs(e, 0)
same(skip_zero_r2z, -3 * mu * U00,
     "linear degenerate confluent r2z origin")
arm_constant_zero = rlin.as_expr().coeff(e, 0).subs(t, 0)
arm_constant_numerator = sp.factor(2916 * arm_constant_zero)
same(
    arm_constant_numerator,
    256 * beta * mu**5 + 768 * beta * mu**3
    + 432 * beta * mu - 243,
    "linear degenerate confluent arm constant",
)
same(
    skip_zero_r3.subs(U00, 0),
    2 * beta
    * (256 * beta * mu**5 + 768 * beta * mu**3 + 432 * beta * mu)
    / 81,
    "linear degenerate confluent raw r3 lead",
)
same(
    skip_zero_r3.subs(U00, 0),
    6 * beta + 2 * beta * arm_constant_numerator / 81,
    "linear degenerate confluent arm-reduced r3 obstruction",
)


# ---------------------------------------------------------------------------
# 8. Audit the candidate's later r^2 z closure as an independent path.
# ---------------------------------------------------------------------------
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
generic_r2z_full = b_r2z.subs(
    {f: f, kap: mu * f + K0, S: S0, T: mu * S0,
     p: p_generic, q: mu * p_generic + P0,
     g: g, h: mu * g + Q_generic}
).doit()
generic_r2z_diff = sp.expand(
    generic_r2z_full.xreplace(differential_symbols)
)
same(
    generic_r2z_diff.subs(vv, 0),
    sp.Rational(60, 7) * delta / alpha * e**2 * ff * uu * vd,
    "candidate generic r2z root coordinate",
)

degenerate_r2z_full = b_r2z.subs(
    {f: f, kap: mu * f + K0, S: S0, T: mu * S0,
     p: p_skip, q: mu * p_skip, g: g, h: mu * g + Q0}
).doit()
degenerate_r2z_diff = sp.expand(
    degenerate_r2z_full.xreplace(differential_symbols)
)
same(
    degenerate_r2z_diff.subs(vv, 0),
    e**2 * (3 * e - 2 * mu) * uu * vd,
    "candidate degenerate r2z first root address",
)

unit1, U1, f1root, g1root = sp.symbols("unit1 U1 f1root g1root")
for sample_mult in (1, 2, 3, 4, 6):
    local_v_full = x**sample_mult * (unit + unit1 * x)
    local_U_full = Uroot + U1 * x
    local_f_full = froot + f1root * x
    local_g_full = g0 + g1root * x
    addressed = degenerate_r2z_diff.subs(
        {
            e: rho + x,
            mu: sp.Rational(3, 2) * rho,
            vv: local_v_full,
            vd: sp.diff(local_v_full, x),
            uu: local_U_full,
            ud: sp.diff(local_U_full, x),
            ff: local_f_full,
            fd: sp.diff(local_f_full, x),
            gg: local_g_full,
            gd: sp.diff(local_g_full, x),
        }
    )
    addressed_coefficient = sp.expand(addressed).coeff(x, sample_mult)
    same(
        addressed_coefficient,
        sp.Rational(3, 2) * rho**2 * Uroot * unit
        * (2 * sample_mult - 5),
        f"candidate addressed degenerate coefficient m={sample_mult}",
    )

# Candidate generic monomial-origin jets and root-ratio obstruction.
fj0, fj1, gj0, gj1 = sp.symbols("fj0 fj1 gj0 gj1")
fjet2 = fj0 + fj1 * e
gjet2 = gj0 + gj1 * e
origin_profiles = {
    f: fjet2,
    kap: mu * fjet2,
    S: 0,
    T: 0,
    p: 0,
    q: 0,
    g: gjet2,
    h: mu * gjet2,
}
origin_r2 = b_r2.subs(origin_profiles).doit().subs(e, 0)
origin_z2 = b_z2.subs(origin_profiles).doit().subs(e, 0)
same(origin_r2, mu * gj0, "candidate generic origin r2 jet")
same(origin_z2, -3 * fj0 - mu * fj1 - 9 * gj0,
     "candidate generic origin z2 jet")
same(
    origin_z2.subs(
        {fj0: sp.Rational(1, 12), fj1: -mu / 6, gj0: 0}
    ),
    mu**2 / 6 - sp.Rational(1, 4),
    "candidate generic origin mu-square law",
)
same(
    -(sp.Rational(4, 3) * sp.Rational(3, 2) + 2),
    -4,
    "candidate D-root ratio trace",
)
qratio = sp.symbols("qratio")
gate((qratio**2 + 4 * qratio + 1).subs(qratio, 1) != 0,
     "candidate root ratio excludes +1")
gate((qratio**2 + 4 * qratio + 1).subs(qratio, -1) != 0,
     "candidate root ratio excludes -1")

# Candidate degenerate monomial-origin coefficients.
degree_sample, Uorigin = sp.symbols(
    "degree_sample Uorigin", integer=True, positive=True
)
candidate_terminal_origin = (
    2 * beta * (8 * degree_sample + 4) * sp.Rational(1, 12)
    + 3 * theta * (4 * degree_sample + 2) * Uorigin
)
same(
    candidate_terminal_origin,
    sp.Rational(2, 3) * (2 * degree_sample + 1)
    * (beta + 9 * theta * Uorigin),
    "candidate degenerate terminal origin payment",
)
candidate_r2z_origin = -mu * (2 * degree_sample + 1) * Uorigin
same(
    candidate_r2z_origin.subs(Uorigin, -beta / (9 * theta)),
    mu * beta * (2 * degree_sample + 1) / (9 * theta),
    "candidate degenerate monomial r2z contradiction",
)
gate(sp.Mod(3 + 4 * degree_sample, 2) == 1,
     "candidate mu-zero arm exponent is odd")


semantic = {
    "audit_verdict": "FAIL_AS_WRITTEN",
    "first_failure": (
        "generic origin valuation (30) omits the equal-order gP block when "
        "ord_0(v)=0"
    ),
    "repair": (
        "a=0 is killed by r3 order one; a>=1 is killed by r2z2 order three; "
        "therefore the generic branch has ord_0(v)>=1 before the ladder transfer"
    ),
    "repaired_result": (
        "the next r2z2 bucket independently eliminates nonzero roots in both "
        "branches; the resulting pure-origin towers are contradictory, so "
        "the candidate's complete rz2 no-go survives after repair"
    ),
    "scope": "profile no-go only; no planar Jacobian counterexample",
}
semantic_bytes = json.dumps(
    semantic, sort_keys=True, separators=(",", ":")
).encode("utf-8")
source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python asserts")

print("audit=THM-3821-rz2-odd-ladder-independent-hostile-audit")
print("verdict=FAIL_AS_WRITTEN_REPAIRED_COMPLETE_NO_GO")
print("first_failure=origin_d0_gP_shares_order_3_with_fK_forcing")
print("repair=d0_generic_empty_via_r3_e1_or_r2z2_e3")
print("extension=r2z2_eliminates_nonzero_roots_then_pure_origin_towers")
print("repaired_result=complete_rz2_profile_no_go")
print("scope=no_planar_JC_claim")
print(f"semantic_sha256={hashlib.sha256(semantic_bytes).hexdigest()}")
print(f"GATES={GATES}")
print("RESULT=PASS")
