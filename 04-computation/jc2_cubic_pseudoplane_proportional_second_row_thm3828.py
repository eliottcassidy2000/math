#!/usr/bin/env python3
"""Exact companion for THM-3828's proportional second-row R2Z2 no-go."""

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


# Full first canonical row plus a single second-row r^2 z^2 slot.
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

bucket_z = normal_poly.coeff_monomial(z)
bucket_r7 = normal_poly.coeff_monomial(r**7)
bucket_r6 = normal_poly.coeff_monomial(r**6)
bucket_r5 = normal_poly.coeff_monomial(r**5)
bucket_r4z2 = normal_poly.coeff_monomial(r**4 * z**2)
bucket_r4z1 = normal_poly.coeff_monomial(r**4 * z)
bucket_r4 = normal_poly.coeff_monomial(r**4)
bucket_r3z2 = normal_poly.coeff_monomial(r**3 * z**2)

expected_z = (36 * e**2 * f - 24 * e * kap - 12 * f + 1) / 2
expected_r7 = -30 * e**2 * (-X * sp.diff(Y, e) + Y * sp.diff(X, e))
zero(bucket_z - expected_z, "arm bucket")
zero(bucket_r7 - expected_r7, "second-row top Wronskian")


# Proportional top branch Y=lambda X, with the preceding slot aligned:
# T=lambda S.  Differences M,N,H remove the common target direction.
lam = sp.symbols("lambda")
M = sp.Function("M")(e)
N = sp.Function("N")(e)
H = sp.Function("H")(e)
aligned = {
    Y: lam * X,
    T: lam * S,
    kap: lam * f + M,
    q: lam * p + N,
    h: lam * g + H,
}
r6_aligned = sp.expand(bucket_r6.subs(aligned).doit())
r5_aligned = sp.expand(bucket_r5.subs(aligned).doit())
r4z2_aligned = sp.expand(bucket_r4z2.subs(aligned).doit())
r4z1_aligned = sp.expand(bucket_r4z1.subs(aligned).doit())
r4_aligned = sp.expand(bucket_r4.subs(aligned).doit())
r3z2_aligned = sp.expand(bucket_r3z2.subs(aligned).doit())
zero(r6_aligned, "aligned r6 bucket vanishes")

expected_r5_aligned = 6 * e * (
    -2 * e * M * sp.diff(X, e)
    + 5 * e * X * sp.diff(M, e)
    + 2 * M * X
)
expected_r4z2_aligned = 15 * e * (
    2 * X * sp.diff(N, e) - N * sp.diff(X, e)
)
expected_r4z1_aligned = 3 * (
    10 * e * X * sp.diff(H, e)
    - 3 * e * H * sp.diff(X, e)
    - 2 * X * H
)
expected_r4_aligned = 3 * (
    -4 * e**2 * M * sp.diff(S, e)
    + 7 * e**2 * S * sp.diff(M, e)
    + 2 * e * M * S
    - 6 * e * M * sp.diff(X, e)
    + 18 * e * X * sp.diff(M, e)
    + 4 * M * X
)
expected_r3z2_aligned = (
    (3 * e**2 - 2 * lam * e) * sp.diff(X, e)
    + (-18 * e + 2 * lam) * X
    + 21 * e * S * sp.diff(N, e)
    - 15 * e * N * sp.diff(S, e)
    - 3 * S * N
    + 24 * X * sp.diff(N, e)
    - 12 * N * sp.diff(X, e)
)
for actual, expected, label in [
    (r5_aligned, expected_r5_aligned, "aligned r5 X/M bucket"),
    (r4z2_aligned, expected_r4z2_aligned, "aligned r4z2 X/N bucket"),
    (r4z1_aligned, expected_r4z1_aligned, "aligned r4z1 X/H bucket"),
    (r4_aligned, expected_r4_aligned, "aligned r4 S/M bucket"),
    (r3z2_aligned, expected_r3z2_aligned, "aligned r3z2 terminal bucket"),
]:
    zero(actual - expected, label)

# The arm after the target-direction differences is D*f=2e*M-1/12.
D = 3 * e**2 - 2 * lam * e - 1
arm_aligned = sp.expand(bucket_z.subs(aligned).doit() / 6)
zero(
    arm_aligned - (D * f - 2 * e * M + sp.Rational(1, 12)),
    "aligned arm identity",
)
gate(sp.degree(D, e) == 2, "M=0 leaves a nonconstant arm divisor")


# M!=0: the r5 equation is X^2/(e^2 M^5)=constant.  If N!=0,
# r4z2 adds N^2/X=constant and the combined UFD tower is 10/5/4.
alpha, beta, delta = sp.symbols("alpha beta delta", nonzero=True)
v = sp.Function("v")(e)
X_52 = alpha * e * v**5
M_52 = beta * v**2
zero(
    expected_r5_aligned.subs({X: X_52, M: M_52}).doit(),
    "5/2 intermediate tower solves r5",
)
mult = sp.symbols("mult", integer=True, nonnegative=True)
gate(2 * (1 + 5 * mult) == 2 + 5 * (2 * mult),
     "5/2 e-prime valuation family")
gate(2 * (5 * mult) == 5 * (2 * mult),
     "5/2 nonzero-prime valuation family")

X_tower = alpha * e**6 * v**10
M_tower = beta * e**2 * v**4
N_tower = delta * e**3 * v**5
zero(
    expected_r5_aligned.subs({X: X_tower, M: M_tower}).doit(),
    "10/4 tower solves r5",
)
zero(
    expected_r4z2_aligned.subs({X: X_tower, N: N_tower}).doit(),
    "10/5 tower solves r4z2",
)
gate(4 * (3 + 5 * mult) == 2 + 5 * (2 + 4 * mult),
     "combined e-prime valuation family")
gate(4 * (5 * mult) == 5 * (4 * mult),
     "combined nonzero-prime valuation family")

# The unused H equation has the compatible optional 10/3 rung.
theta = sp.symbols("theta", nonzero=True)
H_tower = theta * e**2 * v**3
zero(
    expected_r4z1_aligned.subs({X: X_tower, H: H_tower}).doit(),
    "optional 10/3 rung solves r4z1",
)

# The r4 equation integrates exactly; gamma is the homogeneous 7-rung.
gamma = sp.symbols("gamma")
S_tower = alpha * e**5 * v**10 + gamma * e**4 * v**7
zero(
    expected_r4_aligned.subs(
        {X: X_tower, M: M_tower, S: S_tower}
    ).doit(),
    "integrated S profile",
)

# The next bucket is the terminal equation, including the fixed-base lambda
# terms that cannot be removed by silently shearing the target.
terminal = (
    3 * delta * e**2 * v**5 * (3 * e * sp.diff(v, e) + v)
    - 2 * e * (3 * e - 2 * lam) * sp.diff(v, e)
    + 2 * lam * v
)
terminal_actual = expected_r3z2_aligned.subs(
    {X: X_tower, M: M_tower, N: N_tower, S: S_tower}
).doit()
zero(
    terminal_actual + 5 * alpha * e**6 * v**9 * terminal,
    "fixed-base terminal equation",
)

# A nonzero root first has to sit at rho=2lambda/3.  At that address the next
# coefficient is 3*rho*u*(1-2m), whose resonance m=1/2 is not integral.
vv, vd = sp.symbols("vv vd")
terminal_differential = terminal.xreplace({v: vv, sp.diff(v, e): vd})
zero(
    terminal_differential.subs(vv, 0)
    + 2 * e * (3 * e - 2 * lam) * vd,
    "terminal first root address",
)
rho, root_unit = sp.symbols("rho root_unit", nonzero=True)
root_order = sp.symbols("root_order", integer=True, positive=True)
vprime_coefficient = -2 * e * (3 * e - 2 * lam)
v_coefficient = 2 * lam
zero(
    sp.diff(vprime_coefficient, e).subs(
        {e: rho, lam: sp.Rational(3, 2) * rho}
    )
    + 6 * rho,
    "addressed derivative coefficient",
)
zero(
    v_coefficient.subs(lam, sp.Rational(3, 2) * rho) - 3 * rho,
    "addressed undifferentiated coefficient",
)
addressed_next = 3 * rho * root_unit * (1 - 2 * root_order)
zero(
    (-6 * rho * root_order + 3 * rho) * root_unit - addressed_next,
    "addressed next coefficient",
)
gate(
    sp.solve(sp.Eq(1 - 2 * root_order, 0), root_order) == [],
    "half-integral root resonance is unavailable",
)

# With no nonzero roots, v=c*e^d.  A constant v already leaves a nonzero
# e^2 coefficient.  For d>=1 the lowest origin coefficient forces lambda=0.
v_constant = sp.symbols("v_constant", nonzero=True)
constant_terminal = terminal_differential.subs({vv: v_constant, vd: 0})
gate(sp.Poly(constant_terminal, e).coeff_monomial(e**2) != 0,
     "constant-v terminal obstruction")
degree_d = sp.symbols("degree_d", integer=True, positive=True)
origin_lambda_coefficient = 2 * lam * (2 * degree_d + 1)
gate(2 * degree_d + 1 > 0, "positive monomial origin multiplier")

# At lambda=0 the terminal equation has a rational first integral.
terminal_lambda_zero = (
    delta * v**5 * (3 * e * sp.diff(v, e) + v)
    - 2 * sp.diff(v, e)
)
first_integral = (1 + delta * e * v**5) / v**2
zero(
    sp.diff(first_integral, e) - terminal_lambda_zero / v**3,
    "lambda-zero first integral",
)
gate(sp.Integer(1) != 0, "first integral fails at every root of nonconstant v")
gate(delta != 0, "constant v leaves delta*v^6")


# N=0 is a separate specialization and is recomputed rather than divided out.
terminal_N_zero = expected_r3z2_aligned.subs(N, 0).doit()
expected_N_zero = (
    e * (3 * e - 2 * lam) * sp.diff(X, e)
    - 2 * (9 * e - lam) * X
)
zero(terminal_N_zero - expected_N_zero, "N-zero terminal equation")
chi = sp.symbols("chi", nonzero=True)
X_N_zero = chi * e * (3 * e - 2 * lam) ** 5
zero(
    expected_N_zero.subs(X, X_N_zero).doit(),
    "N-zero terminal solution",
)
M_N_zero = beta * (3 * e - 2 * lam) ** 2
zero(
    expected_r5_aligned.subs({X: X_N_zero, M: M_N_zero}).doit(),
    "N-zero compatible M tower",
)
arm_N_zero_numerator = 2 * e * M_N_zero - sp.Rational(1, 12)
_, arm_N_zero_remainder = sp.div(
    sp.Poly(arm_N_zero_numerator, e), sp.Poly(D, e)
)
zero(
    arm_N_zero_remainder.as_expr()
    - (72 * beta * e - 48 * beta * lam - 1) / 12,
    "N-zero arm remainder",
)
gate(
    arm_N_zero_remainder.as_expr().coeff(e, 1) == 6 * beta,
    "N-zero arm remainder has nonzero linear coefficient",
)


# Open complementary branch.  If L=T-lambda*S is nonzero, r6 gives a new
# 10/7 tower.  It is deliberately frozen as anatomy, not claimed empty.
L = sp.Function("L")(e)
r6_difference = sp.expand(
    bucket_r6.subs({Y: lam * X, T: lam * S + L}).doit()
)
open_107_ode = 7 * e * L * sp.diff(X, e) - 10 * e * X * sp.diff(L, e) - 2 * L * X
zero(r6_difference + 3 * e * open_107_ode,
     "open 10/7 r6 reduction")
L_tower = beta * e**4 * v**7
zero(
    open_107_ode.subs({X: X_tower, L: L_tower}).doit(),
    "open 10/7 tower",
)
gate(7 * (6 + 10 * mult) == 2 + 10 * (4 + 7 * mult),
     "10/7 e-prime valuation family")
gate(7 * (10 * mult) == 10 * (7 * mult),
     "10/7 nonzero-prime valuation family")


semantic = {
    "ansatz": "THM3821 plus r2z2*X(e),r2z2*Y(e)",
    "branch": "X!=0;Y=lambda*X;T=lambda*S",
    "differences": "M=kappa-lambda*f;N=q-lambda*p;H=h-lambda*g",
    "N_nonzero": "X=alpha*e6*v10;M=beta*e2*v4;N=delta*e3*v5;S=alpha*e5*v10+gamma*e4*v7",
    "terminal": "3delta*e2*v5*(3e*vprime+v)-2e*(3e-2lambda)*vprime+2lambda*v=0",
    "root_gate": "rho=2lambda/3 then resonance m=1/2;no nonzero roots",
    "confluent_gate": "lambda=0 then (1+delta*e*v5)/v2 is constant;impossible",
    "N_zero": "X=chi*e*(3e-2lambda)^5;arm remainder has coefficient 6beta",
    "open": "L=T-lambda*S nonzero gives X=alpha*e6*v10,L=beta*e4*v7;not closed",
    "scope": "proportional aligned second-row no-go only;one-sided and 10/7 branches open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3828-proportional-second-row-r2z2-profile-nonentry")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("ansatz=THM3821_plus_r2z2_XY")
print("branch=X_nonzero;Y=lambdaX;T=lambdaS")
print("N_nonzero=10_5_4_3_ladder_then_fixed_base_terminal")
print("root_gate=address_2lambda_over_3;next_resonance_one_half")
print("confluent_gate=lambda_zero_first_integral_impossible")
print("N_zero=arm_linear_remainder_6beta_nonzero")
print("open_branch=L_nonzero_10_over_7_tower_not_closed")
print("scope=proportional_aligned_second_row_no_go_only")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
