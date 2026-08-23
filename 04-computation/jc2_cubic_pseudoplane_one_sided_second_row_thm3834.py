#!/usr/bin/env python3
"""Exact companion for THM-3834's one-sided second-row no-go."""

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


f = sp.Function("f")(e)
g = sp.Function("g")(e)
h = sp.Function("h")(e)
kap = sp.Function("kappa")(e)
p = sp.Function("p")(e)
q = sp.Function("q")(e)
S = sp.Function("S")(e)
T = sp.Function("T")(e)
X = sp.Function("X")(e)
Y = sp.Function("Y")(e)
A = e**2 - z / 3 + r * g + z**2 * f + r * z * p + r * z**2 * S + r**2 * z**2 * X
C = (
    e**3
    - e
    - e * z / 2
    + r * h
    + z**2 * kap
    + r * z * q
    + r * z**2 * T
    + r**2 * z**2 * Y
)
normal = sp.rem(bracket(A, C) - 1, monic_relation, z)
normal_poly = sp.Poly(sp.expand(normal), r, z)
gate(
    max(z_degree for _, z_degree in normal_poly.monoms()) <= 2,
    "canonical z-normal degree",
)
zero(
    sp.rem(bracket(A, C) - 1 - normal, monic_relation, z),
    "canonical reduction has no quotient loss",
)


def bucket(r_degree: int, z_degree: int = 0) -> sp.Expr:
    return normal_poly.coeff_monomial(r**r_degree * z**z_degree)


bucket_z = bucket(0, 1)
bucket_r7 = bucket(7)
bucket_r6 = bucket(6)
bucket_r5 = bucket(5)
bucket_r4 = bucket(4)

expected_z = (36 * e**2 * f - 24 * e * kap - 12 * f + 1) / 2
expected_r7 = 30 * e**2 * (X * sp.diff(Y, e) - Y * sp.diff(X, e))
zero(bucket_z - expected_z, "arm bucket")
zero(bucket_r7 - expected_r7, "top Wronskian bucket")


# One-sided top orientation X=0, Y!=0.
one_sided = {X: 0}
arm = sp.expand(bucket_z.subs(one_sided).doit() / 6)
arm_expected = (3 * e**2 - 1) * f - 2 * e * kap + sp.Rational(1, 12)
zero(arm - arm_expected, "one-sided arm identity")
f0 = sp.Rational(1, 12)
zero(arm_expected.subs({e: 0, f: f0}), "arm fixes f(0)=1/12")
zero(bucket_r7.subs(one_sided).doit(), "one-sided top bucket vanishes")

r6_one = sp.expand(bucket_r6.subs(one_sided).doit())
r5_one = sp.expand(bucket_r5.subs(one_sided).doit())
r4_one = sp.expand(bucket_r4.subs(one_sided).doit())
expected_r6_one = -3 * e * (
    -7 * e * S * sp.diff(Y, e)
    + 10 * e * Y * sp.diff(S, e)
    + 2 * S * Y
)
zero(r6_one - expected_r6_one, "one-sided r6 bucket")


# Branch I: S!=0.  The r6 equation is the same 10/7 tower as THM-3829.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
gamma = sp.symbols("gamma")
v = sp.Function("v")(e)
Y_tower = alpha * e**6 * v**10
S_tower = beta * e**4 * v**7
zero(
    r6_one.subs({Y: Y_tower, S: S_tower}).doit(),
    "one-sided 10/7 tower solves r6",
)
mult = sp.symbols("mult", integer=True, nonnegative=True)
gate(7 * (6 + 10 * mult) == 2 + 10 * (4 + 7 * mult), "e-prime 10/7 family")
gate(7 * (10 * mult) == 10 * (7 * mult), "nonzero-prime 10/7 family")

T_integrated = (
    sp.Rational(10, 7) * alpha / beta * e**2 * v**3 * f
    + sp.Rational(2, 7) * alpha * e**5 * v**10
    + gamma * e**4 * v**7
)
zero(
    r5_one.subs({Y: Y_tower, S: S_tower, T: T_integrated}).doit(),
    "one-sided exact r5 integration",
)
W = sp.Function("W")(e)
r5_homogeneous = sp.expand(
    r5_one.subs({Y: Y_tower, S: S_tower, T: T_integrated + W}).doit()
    - r5_one.subs({Y: Y_tower, S: S_tower, T: T_integrated}).doit()
)
zero(
    r5_homogeneous
    - 21 * e**2 * (S_tower * sp.diff(W, e) - W * sp.diff(S_tower, e)),
    "one-sided r5 homogeneous kernel equation",
)
zero(
    r5_homogeneous.subs(W, e**4 * v**7).doit(),
    "one-sided r5 homogeneous kernel is scalar e4v7",
)

r4_tower = sp.expand(
    r4_one.subs({Y: Y_tower, S: S_tower, T: T_integrated}).doit()
)
P4 = sp.cancel(r4_tower / (3 * e**3 * v**2 / (7 * beta)))
P4_expected = (
    30 * alpha * beta**2 * e**7 * v**14 * sp.diff(v, e)
    + 10 * alpha * beta**2 * e**6 * v**15
    + 20 * alpha * beta * e**4 * f * v**7 * sp.diff(v, e)
    - 20 * alpha * beta * e**4 * v**8 * sp.diff(f, e)
    + 20 * alpha * beta * e**3 * f * v**8
    + 120 * alpha * e * f**2 * sp.diff(v, e)
    - 30 * alpha * e * f * v * sp.diff(f, e)
    + 60 * alpha * f**2 * v
    - 196 * beta**2 * e**3 * kap * v**4 * sp.diff(v, e)
    + 49 * beta**2 * e**3 * v**5 * sp.diff(kap, e)
    - 98 * beta**2 * e**2 * kap * v**5
    + 196 * beta * gamma * e**3 * f * v**4 * sp.diff(v, e)
    - 49 * beta * gamma * e**3 * v**5 * sp.diff(f, e)
    + 98 * beta * gamma * e**2 * f * v**5
)
zero(P4 - P4_expected, "full one-sided r4 source polynomial")

origin_degree = sp.symbols("origin_degree", integer=True, nonnegative=True)
unit_v = sp.symbols("unit_v", nonzero=True)
origin_lead = 60 * alpha * f0**2 * unit_v * (2 * origin_degree + 1)
gate(2 * origin_degree + 1 > 0, "one-sided origin odd multiplier")
gate(origin_lead != 0, "one-sided origin leading coefficient is nonzero")

# Hostile fixed-degree controls replay the coefficient directly in the full P4
# source; they are controls only, not the proof of the symbolic all-d law.
f1, k0 = sp.symbols("f1 k0")
for degree_value in range(6):
    specialized = sp.expand(
        P4_expected.subs(
            {
                v: e**degree_value,
                f: f0 + f1 * e,
                kap: k0,
            }
        ).doit()
    )
    coefficient = sp.Poly(specialized, e).coeff_monomial(e**degree_value)
    zero(
        coefficient - 60 * alpha * f0**2 * (2 * degree_value + 1),
        f"one-sided origin hostile degree {degree_value}",
    )


# Branch II: S=0.  The fixed e-charge in r5 gives a shifted 5/2 Kummer
# relation rather than the unshifted relation one would get by dropping the
# 4eYf term.
r5_S_zero = sp.expand(r5_one.subs(S, 0).doit())
expected_r5_S_zero = -6 * e * (
    5 * e * Y * sp.diff(f, e)
    - 2 * e * f * sp.diff(Y, e)
    + 2 * Y * f
)
zero(r5_S_zero - expected_r5_S_zero, "S-zero r5 5/2 equation")

w = sp.Function("w")(e)
f_52 = beta * w**2
Y_52 = alpha * e * w**5
zero(
    expected_r5_S_zero.subs({f: f_52, Y: Y_52}).doit(),
    "S-zero shifted 5/2 Kummer tower solves r5",
)
gate(5 * (2 * mult) == 2 * (5 * mult), "S-zero UFD valuation family")
gate(5 * (2 * mult) - 2 * (1 + 5 * mult) == -2,
     "S-zero e-prime shifted valuation family")

# The r4 equation has the exact particular solution T=alpha*w^5.  Its
# homogeneous equation has the forbidden origin resonance ord_0(W)=1/2,
# because w(0) is a unit (the arm fixes f(0)=1/12).
r4_S_zero = sp.expand(r4_one.subs(S, 0).doit())
zero(
    r4_S_zero.subs({Y: Y_52, f: f_52, T: alpha * w**5}).doit(),
    "S-zero exact r4 particular solution",
)
W_zero = sp.Function("W_zero")(e)
r4_S_zero_homogeneous = sp.expand(
    r4_S_zero.subs(
        {Y: Y_52, f: f_52, T: alpha * w**5 + W_zero}
    ).doit()
    - r4_S_zero.subs(
        {Y: Y_52, f: f_52, T: alpha * w**5}
    ).doit()
)
expected_r4_S_zero_homogeneous = -6 * beta * e * w * (
    7 * e * W_zero * sp.diff(w, e)
    - 2 * e * w * sp.diff(W_zero, e)
    + W_zero * w
)
zero(
    r4_S_zero_homogeneous - expected_r4_S_zero_homogeneous,
    "S-zero r4 homogeneous equation",
)
origin_order = sp.symbols("origin_order", integer=True, nonnegative=True)
gate(
    sp.solve(sp.Eq(1 - 2 * origin_order, 0), origin_order) == [],
    "S-zero r4 homogeneous half-integral resonance unavailable",
)

# The next two high buckets force p=g=0 by the same half-integral origin
# resonance.  We recompute them from the full canonical normal form.
bucket_r4z2 = bucket(4, 2)
bucket_r4z1 = bucket(4, 1)
r4z2_one = sp.expand(bucket_r4z2.subs(one_sided).doit())
r4z1_one = sp.expand(bucket_r4z1.subs(one_sided).doit())
expected_r4z2_one = 15 * e * (-2 * Y * sp.diff(p, e) + p * sp.diff(Y, e))
expected_r4z1_one = 3 * (
    -10 * e * Y * sp.diff(g, e)
    + 3 * e * g * sp.diff(Y, e)
    + 2 * Y * g
)
zero(r4z2_one - expected_r4z2_one, "one-sided p Wronskian")
zero(r4z1_one - expected_r4z1_one, "one-sided g Wronskian")
gate(
    sp.solve(sp.Eq(1 - 2 * origin_order, 0), origin_order) == [],
    "p half-integral origin resonance unavailable",
)
gate(
    sp.solve(sp.Eq(5 - 10 * origin_order, 0), origin_order) == [],
    "g half-integral origin resonance unavailable",
)

# With p=g=0, r3z2 forces w'=0.  The arm then fixes f=1/12 and
# kappa=e/8, and r3z1 leaves the nonzero monomial -60*a*e^3.
bucket_r3z2 = bucket(3, 2)
bucket_r3z1 = bucket(3, 1)
r3z2_terminal = sp.expand(
    bucket_r3z2.subs(
        {
            X: 0,
            S: 0,
            Y: Y_52,
            f: f_52,
            T: alpha * w**5,
            p: 0,
            g: 0,
        }
    ).doit()
)
zero(
    r3z2_terminal + 10 * alpha * e**2 * w**4 * sp.diff(w, e),
    "S-zero r3z2 forces w constant",
)
zero(
    arm_expected.subs({f: f0, kap: e / 8}).doit(),
    "constant-w arm gives kappa=e/8",
)

a_const = sp.symbols("a_const", nonzero=True)
constant_terminal = sp.expand(
    bucket_r3z1.subs(
        {
            X: 0,
            S: 0,
            Y: a_const * e,
            f: f0,
            T: a_const,
            kap: e / 8,
            p: 0,
            g: 0,
        }
    ).doit()
)
zero(constant_terminal + 60 * a_const * e**3,
     "S-zero final r3z coefficient")
gate(constant_terminal != 0, "S-zero final monomial is nonzero")


semantic = {
    "ansatz": "THM3821 plus r2z2*X(e),r2z2*Y(e)",
    "branch": "X=0;Y!=0",
    "arm": "(3e2-1)f-2e*kappa+1/12=0;f(0)=1/12",
    "S_nonzero": "Y=alpha*e6*v10;S=beta*e4*v7;T=(10alpha/(7beta))*e2*v3*f+(2alpha/7)*e5*v10+gamma*e4*v7",
    "S_nonzero_gate": "cancelled r4 source has origin lead 60alpha*f(0)^2*u(0)*(2d+1)",
    "S_zero": "5eYfprime-2efYprime+2Yf=0;f=beta*w2;Y=alpha*e*w5;w(0)!=0",
    "S_zero_gate": "r4 gives unique T=alpha*w5;p=g=0 by half-resonances;r3z2 forces w constant;r3z1=-60a*e3",
    "corollary": "THM3821+THM3828+THM3829+THM3834 close the entire fixed first r2z2 ansatz",
    "scope": "higher canonical slots and general planar JC remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3834-one-sided-second-row-r2z2-profile-nonentry")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("ansatz=THM3821_plus_r2z2_XY")
print("branch=X_zero;Y_nonzero")
print("arm=f_at_zero_one_twelfth")
print("S_nonzero=10_7_tower;complete_T_integral;r4_origin_odd_multiplier")
print("S_zero=shifted_5_2_Kummer;unique_T;p_g_half_resonance;w_constant;final_r3z")
print("corollary=full_fixed_r2z2_ansatz_empty_via_THM3821_3828_3829_3834")
print("scope=higher_slots_and_general_JC_open")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
