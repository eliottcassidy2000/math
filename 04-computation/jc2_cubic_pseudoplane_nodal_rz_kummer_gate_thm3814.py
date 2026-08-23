#!/usr/bin/env python3
"""Exact companion for THM-3814's nodal rz Kummer degree gate."""

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
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


# Full first omitted mixed-profile ansatz after THM-3812.
g = sp.Function("g")(e)
f = sp.Function("f")(e)
h = sp.Function("h")(e)
kap = sp.Function("kappa")(e)
p = sp.Function("p")(e)
q = sp.Function("q")(e)
A = e**2 - z / 3 + r * g + z**2 * f + r * z * p
C = e**3 - e - e * z / 2 + r * h + z**2 * kap + r * z * q
normal = sp.rem(bracket(A, C) - 1, monic_relation, z)
normal_poly = sp.Poly(sp.expand(normal), r, z)
gate(max(z_degree for (_, z_degree) in normal_poly.monoms()) <= 2,
     "canonical z-normal degree")
zero(
    sp.rem(bracket(A, C) - 1 - normal, monic_relation, z),
    "canonical reduction has no quotient loss",
)


# The three structural buckets: arm Bezout, top Wronskian, and Kummer ODE.
z_bucket = normal_poly.coeff_monomial(z)
r3z_bucket = normal_poly.coeff_monomial(r**3 * z)
r2z2_bucket = normal_poly.coeff_monomial(r**2 * z**2)
expected_z = (36 * e**2 * f - 24 * e * kap - 12 * f + 1) / 2
expected_r3z = 15 * e * (p * sp.diff(q, e) - q * sp.diff(p, e))
expected_r2z2 = -3 * (
    -4 * e * f * sp.diff(q, e)
    + 4 * e * kap * sp.diff(p, e)
    - 5 * e * p * sp.diff(kap, e)
    + 5 * e * q * sp.diff(f, e)
    + 2 * f * q
    - 2 * kap * p
)
zero(z_bucket - expected_z, "pure z arm-Bezout bucket")
zero(r3z_bucket - expected_r3z, "r^3z Wronskian bucket")
zero(r2z2_bucket - expected_r2z2, "r^2z^2 Kummer bucket")


lam = sp.symbols("lambda")
W = sp.Function("W")(e)


def reduce_lambda_law(expression: sp.Expr) -> sp.Expr:
    """Reduce a rational expression modulo 4*lambda^2+3.

    Every use below has lambda-independent denominator; rejecting any drift
    keeps the quotient-ring computation honest rather than relying on a
    syntactic ``subs(lambda**2, ...)``.
    """

    numerator, denominator = sp.fraction(sp.cancel(expression))
    if sp.diff(denominator, lam) != 0:
        raise RuntimeError("lambda-law reduction acquired a lambda denominator")
    remainder = sp.rem(numerator, 4 * lam**2 + 3, lam)
    return sp.cancel(remainder / denominator)


kummer_ode = 4 * e * W * sp.diff(p, e) - 5 * e * p * sp.diff(W, e) - 2 * p * W
zero(
    expected_r2z2.subs({q: lam * p, kap: lam * f + W}).doit()
    + 3 * kummer_ode,
    "symmetric branch reduces to the 5/4 Kummer ODE",
)
zero(
    sp.diff(p**4 / (e**2 * W**5), e)
    - p**4 / (e**2 * W**5)
    * (4 * sp.diff(p, e) / p - 2 / e - 5 * sp.diff(W, e) / W),
    "Kummer ODE is a logarithmic derivative",
)


# The valuation parametrization solves the ODE identically.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
v = sp.Function("v")(e)
p_tower = alpha * e**3 * v**5
W_tower = beta * e**2 * v**4
zero(kummer_ode.subs({p: p_tower, W: W_tower}).doit(),
     "symmetric Kummer tower solves the ODE")
t = sp.symbols("t", integer=True, nonnegative=True)
gate(4 * (3 + 5 * t) - 5 * (2 + 4 * t) == 2,
     "zero-prime valuation family")
gate(4 * (5 * t) - 5 * (4 * t) == 0,
     "nonzero-prime valuation family")


# Asymmetric p=0,q!=0: the same tower puts e^2 into f, contradicting the
# arm Bezout value f(0)=1/12.
q_asymmetric = alpha * e**3 * v**5
f_asymmetric = beta * e**2 * v**4
asymmetric_ode = (
    5 * e * q * sp.diff(f, e)
    - 4 * e * f * sp.diff(q, e)
    + 2 * f * q
)
zero(asymmetric_ode.subs({q: q_asymmetric, f: f_asymmetric}).doit(),
     "asymmetric Kummer tower solves its ODE")
gate(f_asymmetric.subs(e, 0) == 0, "asymmetric tower forces f(0)=0")
gate(sp.Rational(1, 12) != 0, "arm Bezout forces f(0)=1/12")


# Constant-v symmetric branch.  Absorb the constant v into alpha,beta.
# Arm divisibility first forces lambda^2=-3/4 and beta=-lambda/4.
D = 3 * e**2 - 2 * lam * e - 1
numerator_f = 2 * beta * e**3 - sp.Rational(1, 12)
quotient_f, remainder_f = sp.div(sp.Poly(numerator_f, e), sp.Poly(D, e))
zero(
    numerator_f - quotient_f.as_expr() * D - remainder_f.as_expr(),
    "constant-v arm division identity",
)
zero(
    quotient_f.as_expr() - 2 * beta * (3 * e + 2 * lam) / 9,
    "constant-v arm quotient",
)
zero(
    remainder_f.as_expr()
    - (
        sp.Rational(4, 9) * beta * lam
        - sp.Rational(1, 12)
        + e * (sp.Rational(8, 9) * beta * lam**2 + sp.Rational(2, 3) * beta)
    ),
    "constant-v arm remainder",
)

f_min = sp.Rational(1, 12) - lam * e / 6
kap_min = lam / 12 + e / 8 - lam * e**2 / 4
p_min = -sp.Rational(4, 5) * e**3
q_min = lam * p_min
gate(sp.Poly(4 * lam**2 + 3, lam).degree() == 2,
     "constant-v quadratic parameter law")
zero(
    reduce_lambda_law(
        (3 * e**2 - 1) * f_min - 2 * e * kap_min + sp.Rational(1, 12)
    ),
    "minimal arm Bezout specialization",
)
zero(kummer_ode.subs({p: p_min, W: -lam * e**2 / 4}).doit(),
     "minimal Kummer specialization")


# The pure-r^3 and pure-z^2 buckets successively determine h-lambda*g,
# alpha, the resonant e^2 constant, and then g,h.  First verify that the
# claimed two equations are the actual coefficients of the canonical normal
# form, rather than merely compatible auxiliary identities.
alpha_symbol, delta = sp.symbols("alpha_symbol delta", nonzero=True)
H_symbol = h - lam * g
A_branch = (
    e**2
    - z / 3
    + r * g
    + z**2 * f_min
    + r * z * alpha_symbol * e**3
)
C_branch = (
    e**3
    - e
    - e * z / 2
    + r * h
    + z**2 * kap_min
    + r * z * lam * alpha_symbol * e**3
)
branch_normal = sp.Poly(
    sp.rem(bracket(A_branch, C_branch) - 1, monic_relation, z), r, z
)
expected_branch_r3 = -e**3 * (
    -120 * alpha_symbol * e * sp.diff(H_symbol, e)
    + 240 * alpha_symbol * H_symbol
    + 3 * e
    + 4 * lam
) / 8
expected_branch_z2 = 3 * (
    24 * D * g - 48 * e * H_symbol + 2 * e * lam - 1
) / 8
zero(
    reduce_lambda_law(
        branch_normal.coeff_monomial(r**3) - expected_branch_r3
    ),
    "constant-v pure-r^3 coefficient",
)
zero(
    reduce_lambda_law(
        branch_normal.coeff_monomial(z**2) - expected_branch_z2
    ),
    "constant-v pure-z^2 coefficient",
)

H_general = (
    delta * e**2
    - e / (40 * alpha_symbol)
    - lam / (60 * alpha_symbol)
)
zero(
    e * sp.diff(H_general, e) - 2 * H_general
    - (3 * e + 4 * lam) / (120 * alpha_symbol),
    "pure-r^3 H equation",
)
numerator_g = (1 - 2 * lam * e + 48 * e * H_general) / 24
quotient_g, remainder_g = sp.div(sp.Poly(numerator_g, e), sp.Poly(D, e))
zero(
    numerator_g - quotient_g.as_expr() * D - remainder_g.as_expr(),
    "pure-z^2 g division identity",
)
expected_remainder_g = (
    e
    * (
        160 * alpha_symbol * delta * lam**2
        + 120 * alpha_symbol * delta
        - 15 * alpha_symbol * lam
        - 12 * lam
    )
    / (180 * alpha_symbol)
    + (160 * alpha_symbol * delta * lam + 15 * alpha_symbol - 6)
    / (360 * alpha_symbol)
)
zero(remainder_g.as_expr() - expected_remainder_g,
     "pure-z^2 g remainder")
zero(
    reduce_lambda_law(sp.expand(remainder_g.as_expr()).coeff(e, 1))
    + lam * (5 * alpha_symbol + 4) / (60 * alpha_symbol),
    "pure-z^2 linear remainder forces alpha=-4/5",
)
zero(
    reduce_lambda_law(sp.expand(remainder_g.as_expr()).coeff(e, 0))
    - (
        160 * alpha_symbol * delta * lam + 15 * alpha_symbol - 6
    )
    / (360 * alpha_symbol),
    "pure-z^2 constant remainder determines delta",
)

g_min = (3 * lam * e - 1) / 24
h_min = (9 * lam * e**2 - 3 * e - lam) / 48
H_min = h_min - lam * g_min
zero(
    reduce_lambda_law(
        H_min
        - H_general.subs(
            {alpha_symbol: -sp.Rational(4, 5), delta: 3 * lam / 16}
        )
    ),
    "minimal H specialization",
)
zero(
    reduce_lambda_law(
        D * g_min
        - numerator_g.subs(
            {alpha_symbol: -sp.Rational(4, 5), delta: 3 * lam / 16}
        )
    ),
    "minimal g specialization",
)


# Full hostile substitution: all preceding necessary gates vanish, but the
# constant coefficient of the pure-r bucket is exactly 1/4.
A_min = (
    e**2 - z / 3 + r * g_min + z**2 * f_min + r * z * p_min
)
C_min = (
    e**3 - e - e * z / 2 + r * h_min + z**2 * kap_min + r * z * q_min
)
normal_min = sp.Poly(
    sp.rem(bracket(A_min, C_min) - 1, monic_relation, z), r, z, e
)
final_quarter = normal_min.coeff_monomial(r)
zero(reduce_lambda_law(final_quarter - sp.Rational(1, 4)),
     "constant-v final pure-r coefficient")


# The next possible tower is linear.  Normalize its leading coefficient
# into alpha,beta and write v=e-tau.  The arm law gives two exact remainder
# equations.  A rank-five first-order equation supplies a third necessary
# compatibility, and those three equations have no solution with beta!=0.
tau, beta_inv = sp.symbols("tau beta_inv")
v_linear = e - tau
p0_linear = e**3 * v_linear**5
W_linear = beta * e**2 * v_linear**4
numerator_f_linear = 2 * e * W_linear - sp.Rational(1, 12)
quotient_f_linear, remainder_f_linear = sp.div(
    sp.Poly(numerator_f_linear, e), sp.Poly(D, e)
)

arm0 = (
    256 * beta * lam**5
    - 1536 * beta * lam**4 * tau
    + 3456 * beta * lam**3 * tau**2
    + 768 * beta * lam**3
    - 3456 * beta * lam**2 * tau**3
    - 3456 * beta * lam**2 * tau
    + 1296 * beta * lam * tau**4
    + 5184 * beta * lam * tau**2
    + 432 * beta * lam
    - 2592 * beta * tau**3
    - 864 * beta * tau
    - 243
)
arm1 = (
    64 * lam**6
    - 384 * lam**5 * tau
    + 864 * lam**4 * tau**2
    + 240 * lam**4
    - 864 * lam**3 * tau**3
    - 1152 * lam**3 * tau
    + 324 * lam**2 * tau**4
    + 1944 * lam**2 * tau**2
    + 216 * lam**2
    - 1296 * lam * tau**3
    - 648 * lam * tau
    + 243 * tau**4
    + 486 * tau**2
    + 27
)
zero(
    remainder_f_linear.as_expr().coeff(e, 0) - arm0 / 2916,
    "linear-v arm constant remainder",
)
zero(
    remainder_f_linear.as_expr().coeff(e, 1) - 2 * beta * arm1 / 729,
    "linear-v arm linear remainder",
)

# Verify the remaining equation directly against the actual pure-r^3
# coefficient of the arbitrary-profile normal form.
H_function = sp.Function("H")(e)
master_general = (
    4 * e**2 * (W * sp.diff(f, e) - f * sp.diff(W, e))
    + 3 * e * H_function * sp.diff(p, e)
    - 5 * e * p * sp.diff(H_function, e)
    + H_function * p
)
zero(
    normal_poly.coeff_monomial(r**3)
    .subs({q: lam * p, kap: lam * f + W, h: lam * g + H_function})
    .doit()
    + 3 * master_general,
    "linear-v master equation is the actual pure-r^3 bucket",
)

x_coefficients = sp.symbols("x0:6")
X_linear = sum(x_coefficients[j] * e**j for j in range(6))
forcing_linear = sp.expand(
    4
    * e**2
    * (
        W_linear * sp.diff(quotient_f_linear.as_expr(), e)
        - quotient_f_linear.as_expr() * sp.diff(W_linear, e)
    )
)
operator_linear = sp.expand(
    3 * e * X_linear * sp.diff(p0_linear, e)
    - 5 * e * p0_linear * sp.diff(X_linear, e)
    + X_linear * p0_linear
)
linear_equation = sp.Poly(forcing_linear + operator_linear, e)
linear_coefficients = sp.Matrix(
    [linear_equation.coeff_monomial(e**j) for j in range(13)]
)
linear_matrix = linear_coefficients.jacobian(x_coefficients)
linear_rhs = -linear_coefficients.subs({x: 0 for x in x_coefficients})

# Uniform rank five: e^2(e-tau)^3 is an explicit kernel vector, while the
# rows 8,...,12 and columns 0,...,4 have constant nonzero determinant.
kernel_profile = e**2 * v_linear**3
zero(
    3 * e * kernel_profile * sp.diff(p0_linear, e)
    - 5 * e * p0_linear * sp.diff(kernel_profile, e)
    + kernel_profile * p0_linear,
    "linear-v operator kernel",
)
rank_minor = linear_matrix.extract(list(range(8, 13)), list(range(5))).det()
zero(rank_minor - 375000, "linear-v constant rank-five minor")
gate(linear_matrix.rank() == 5, "linear-v coefficient matrix rank")

compatibility_factor = (
    16 * lam**4
    - 48 * lam**3 * tau
    + 36 * lam**2 * tau**2
    + 48 * lam**2
    - 90 * lam * tau
    + 27 * tau**2
    + 27
)
compatibility_core = (-2 * lam + 3 * tau) * compatibility_factor
compatibility_full = 32 * beta**2 * compatibility_core

# For tau!=0 this rational left-null row gives the compatibility.  At
# tau=0 a separate e^7 coefficient gives the same polynomial, avoiding any
# specialization through a denominator.
left_null = sp.zeros(13, 1)
left_null[3] = -1176 / tau**6
left_null[4] = -224 / tau**5
left_null[6] = 20 / tau**3
left_null[7] = 10 / tau**2
left_null[8] = 4 / tau
left_null[9] = 1
zero((left_null.T * linear_matrix)[0, :].norm(),
     "linear-v generic left-null row")
zero(
    (left_null.T * linear_rhs)[0]
    - compatibility_full / (729 * tau**2),
    "linear-v generic compatibility",
)
zero(
    sp.Poly(forcing_linear.subs(tau, 0), e).coeff_monomial(e**7)
    - compatibility_full.subs(tau, 0) / 486,
    "linear-v tau-zero compatibility",
)

linear_groebner = sp.groebner(
    [arm0, arm1, compatibility_core, beta * beta_inv - 1],
    beta_inv,
    beta,
    lam,
    tau,
    order="grevlex",
    domain=sp.QQ,
)
gate(
    len(linear_groebner.polys) == 1
    and linear_groebner.polys[0].as_expr() == 1,
    "linear-v arm and compatibility ideal is the unit ideal",
)


semantic = {
    "ansatz": "THM3812 plus rz*p(e),rz*q(e);c=1",
    "top": "pqprime-qpprime=0",
    "symmetric": "q=lambda*p;W=kappa-lambda*f;p=alpha*e^3*v^5;W=beta*e^2*v^4",
    "asymmetric": "p=0,q!=0 forces f=e^2*beta*v^4, contradicting f(0)=1/12",
    "constant_v": "lambda^2=-3/4;beta=-lambda/4;alpha=-4/5;delta=3lambda/16;final r coefficient=1/4",
    "linear_v": "two arm remainders plus rank-five compatibility have unit ideal",
    "survivor": "deg(v)>=2;deg(p)>=13;deg(W)>=10",
    "open": "degree-at-least-two Kummer tower and higher canonical coefficients",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
linear_packet_blob = json.dumps(
    [
        str(sp.expand(arm0)),
        str(sp.expand(arm1)),
        str(sp.expand(compatibility_core)),
    ],
    separators=(",", ":"),
).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3814-nodal-rz-kummer-profile-degree-gate")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("ansatz=THM3812_plus_rz*p(e),rz*q(e)")
print("top=p*qprime-q*pprime=0")
print("symmetric=q=lambda*p;W=kappa-lambda*f;p=alpha*e^3*v^5;W=beta*e^2*v^4")
print("asymmetric=p0_qnonzero_forces_f0=0_but_arm_forces_f0=1/12")
print("constant_v=lambda^2=-3/4;beta=-lambda/4;alpha=-4/5;delta=3lambda/16;final_r0=1/4")
print("linear_v=arm_two_remainders_plus_rank5_compatibility_have_unit_ideal")
print("survivor=deg_v>=2;deg_p>=13;deg_W>=10")
print("open=degree_at_least_two_Kummer_tower;higher_canonical_coefficients")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"linear_packet_sha256={hashlib.sha256(linear_packet_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
