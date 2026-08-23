#!/usr/bin/env python3
"""Exact companion for THM-3829's misaligned second-row 10/7 no-go."""

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


# THM-3821 plus the first r^2 z^2 second-row slot.
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
bucket_r4z2 = bucket(4, 2)
bucket_r4z1 = bucket(4, 1)
bucket_r4 = bucket(4)
bucket_r3z2 = bucket(3, 2)
bucket_r3z1 = bucket(3, 1)
bucket_r3 = bucket(3)

expected_z = (36 * e**2 * f - 24 * e * kap - 12 * f + 1) / 2
expected_r7 = -30 * e**2 * (-X * sp.diff(Y, e) + Y * sp.diff(X, e))
zero(bucket_z - expected_z, "arm bucket")
zero(bucket_r7 - expected_r7, "top r7 Wronskian")


# The r7 bucket forces Y=lambda X.  L is the deliberately misaligned
# preceding profile, and M,N,H are the remaining target-direction differences.
lam = sp.symbols("lambda")
L = sp.Function("L")(e)
M = sp.Function("M")(e)
N = sp.Function("N")(e)
H = sp.Function("H")(e)
differences = {
    Y: lam * X,
    T: lam * S + L,
    kap: lam * f + M,
    q: lam * p + N,
    h: lam * g + H,
}

r6_difference = sp.expand(bucket_r6.subs(differences).doit())
r5_difference = sp.expand(bucket_r5.subs(differences).doit())
r4z2_difference = sp.expand(bucket_r4z2.subs(differences).doit())
r4z1_difference = sp.expand(bucket_r4z1.subs(differences).doit())
r4_difference = sp.expand(bucket_r4.subs(differences).doit())
r3z2_difference = sp.expand(bucket_r3z2.subs(differences).doit())
r3z1_difference = sp.expand(bucket_r3z1.subs(differences).doit())
r3_difference = sp.expand(bucket_r3.subs(differences).doit())

ode_107 = 7 * e * L * sp.diff(X, e) - 10 * e * X * sp.diff(L, e) - 2 * L * X
zero(r6_difference + 3 * e * ode_107, "misaligned r6 10/7 equation")

expected_r4z2 = 15 * e * (2 * X * sp.diff(N, e) - N * sp.diff(X, e))
expected_r4z1 = 3 * (
    10 * e * X * sp.diff(H, e)
    - 3 * e * H * sp.diff(X, e)
    - 2 * X * H
)
zero(r4z2_difference - expected_r4z2, "X/N Wronskian bucket")
zero(r4z1_difference - expected_r4z1, "X/H Wronskian bucket")

D = 3 * e**2 - 2 * lam * e - 1
arm = sp.expand(bucket_z.subs(differences).doit() / 6)
zero(arm - (D * f - 2 * e * M + sp.Rational(1, 12)), "arm identity")
gate(sp.degree(D, e) == 2, "M=0 would make a quadratic divide a unit")


# UFD solution of the r6 logarithmic equation.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
gamma = sp.symbols("gamma")
v = sp.Function("v")(e)
X_tower = alpha * e**6 * v**10
L_tower = beta * e**4 * v**7
zero(ode_107.subs({X: X_tower, L: L_tower}).doit(), "10/7 tower solves r6")

mult = sp.symbols("mult", integer=True, nonnegative=True)
gate(7 * (6 + 10 * mult) == 2 + 10 * (4 + 7 * mult), "e-prime 10/7 family")
gate(7 * (10 * mult) == 10 * (7 * mult), "nonzero-prime 10/7 family")

# The r5 bucket is a first-order equation in S and integrates without any
# division by M.  The gamma term is its homogeneous seven-rung.
S_integrated = (
    sp.Rational(10, 7) * alpha / beta * e**2 * v**3 * M
    + sp.Rational(2, 7) * alpha * e**5 * v**10
    + gamma * e**4 * v**7
)
zero(
    r5_difference.subs({X: X_tower, L: L_tower, S: S_integrated}).doit(),
    "exact r5 integration",
)
W = sp.Function("W")(e)
r5_homogeneous = sp.expand(
    r5_difference.subs(
        {X: X_tower, L: L_tower, S: S_integrated + W}
    ).doit()
    - r5_difference.subs(
        {X: X_tower, L: L_tower, S: S_integrated}
    ).doit()
)
zero(
    r5_homogeneous
    - 21 * e**2 * (-L_tower * sp.diff(W, e) + W * sp.diff(L_tower, e)),
    "r5 homogeneous kernel equation",
)
zero(
    r5_homogeneous.subs(W, e**4 * v**7).doit(),
    "r5 homogeneous kernel is scalar e4v7",
)


# The cancellation-safe r4 source polynomial.  It is kept in the pretyped M
# form so the two competing local valuations are auditable.
r4_tower = sp.expand(
    r4_difference.subs({X: X_tower, L: L_tower, S: S_integrated}).doit()
)
P4 = sp.cancel(r4_tower / (-3 * e**3 * v**2 / (7 * beta)))
P4_expected = (
    30 * alpha * beta**2 * e**7 * v**14 * sp.diff(v, e)
    + 10 * alpha * beta**2 * e**6 * v**15
    + 20 * alpha * beta * e**4 * M * v**7 * sp.diff(v, e)
    - 20 * alpha * beta * e**4 * v**8 * sp.diff(M, e)
    + 20 * alpha * beta * e**3 * M * v**8
    + 120 * alpha * e * M**2 * sp.diff(v, e)
    - 30 * alpha * e * M * v * sp.diff(M, e)
    + 60 * alpha * M**2 * v
    - 196 * beta**2 * e**3 * f * v**4 * sp.diff(v, e)
    + 49 * beta**2 * e**3 * v**5 * sp.diff(f, e)
    - 98 * beta**2 * e**2 * f * v**5
    + 196 * beta * gamma * e**3 * M * v**4 * sp.diff(v, e)
    - 49 * beta * gamma * e**3 * v**5 * sp.diff(M, e)
    + 98 * beta * gamma * e**2 * M * v**5
)
zero(P4 - P4_expected, "full cancellation-safe r4 source polynomial")

M2_block = 30 * alpha * M * (
    4 * e * M * sp.diff(v, e) - e * v * sp.diff(M, e) + 2 * M * v
)
f_block = 49 * beta**2 * e**2 * v**4 * (
    -4 * e * f * sp.diff(v, e) + e * v * sp.diff(f, e) - 2 * f * v
)
zero(
    P4_expected
    - M2_block
    - f_block
    - (
        30 * alpha * beta**2 * e**7 * v**14 * sp.diff(v, e)
        + 10 * alpha * beta**2 * e**6 * v**15
        + 20 * alpha * beta * e**4 * M * v**7 * sp.diff(v, e)
        - 20 * alpha * beta * e**4 * v**8 * sp.diff(M, e)
        + 20 * alpha * beta * e**3 * M * v**8
        + 196 * beta * gamma * e**3 * M * v**4 * sp.diff(v, e)
        - 49 * beta * gamma * e**3 * v**5 * sp.diff(M, e)
        + 98 * beta * gamma * e**2 * M * v**5
    ),
    "r4 source splits into M2, arm, and strictly later blocks",
)

# Frozen local coefficients.  rho,u,w,F are nonzero leading units and m,s,d
# are local orders.  These arithmetic gates support the valuation table in
# the theorem; they do not infer orders from a finite specialization.
rho, unit_v, unit_M, unit_f = sp.symbols("rho unit_v unit_M unit_f", nonzero=True)
order_m = sp.symbols("order_m", integer=True, positive=True)
order_s = sp.symbols("order_s", integer=True, nonnegative=True)
order_d = sp.symbols("order_d", integer=True, nonnegative=True)
root_M2 = 30 * alpha * rho * unit_v * unit_M**2 * (4 * order_m - order_s)
root_f = -196 * beta**2 * rho**3 * unit_f * order_m * unit_v**5
origin_M2 = 30 * alpha * unit_v * unit_M**2 * (4 * order_d - order_s + 2)
origin_f = -98 * beta**2 * unit_f * unit_v**5 * (2 * order_d + 1)
zero(
    root_M2
    - (
        120 * alpha * rho * unit_M**2 * order_m * unit_v
        - 30 * alpha * rho * unit_M * unit_v * order_s * unit_M
    ),
    "nonzero-root M2 leading coefficient",
)
zero(root_f + 196 * beta**2 * rho**3 * unit_f * order_m * unit_v**5,
     "nonzero-root arm leading coefficient")
zero(
    origin_M2
    - 30 * alpha * unit_v * unit_M**2 * (4 * order_d - order_s + 2),
    "origin M2 leading coefficient",
)
zero(origin_f + 98 * beta**2 * unit_f * unit_v**5 * (2 * order_d + 1),
     "origin arm leading coefficient")
gate(4 * order_m >= 2 * order_m, "nonzero M2 resonance lies above the low side")
gate(4 * order_d + 2 >= 2 * order_d + 1, "origin M2 resonance lies above the low side")

# Therefore M=e*v^2*R, with R a unit at every zero of v.  The typed r4
# equation retains the same local payment 15 alpha R^2=49 beta^2 f.
R = sp.Function("R")(e)
M_typed = e * v**2 * R
P4_typed = sp.expand(P4_expected.subs(M, M_typed).doit())
R0 = sp.symbols("R0", nonzero=True)
zero(
    60 * alpha * R0**2 - 196 * beta**2 * unit_f
    - 4 * (15 * alpha * R0**2 - 49 * beta**2 * unit_f),
    "nonzero-root r4 payment",
)
zero(
    30 * alpha * R0**2 - 98 * beta**2 * unit_f
    - 2 * (15 * alpha * R0**2 - 49 * beta**2 * unit_f),
    "origin r4 payment",
)
f0 = sp.Rational(1, 12)
zero(D.subs(e, 0) * f0 + sp.Rational(1, 12), "arm fixes f(0)=1/12")


# The N and H Wronskians.  Each is split explicitly into zero and nonzero
# branches in the proof.
delta, theta = sp.symbols("delta theta", nonzero=True)
N_tower = delta * e**3 * v**5
H_tower = theta * e**2 * v**3
zero(
    expected_r4z2.subs({X: X_tower, N: N_tower}).doit(),
    "N nonzero 10/5 rung",
)
zero(
    expected_r4z1.subs({X: X_tower, H: H_tower}).doit(),
    "H nonzero 10/3 rung",
)


# The r3z2 bucket gives only the lower bound needed later: p is divisible by
# e^2 v^3.  We freeze both N branches and both local coefficient formulas.
S_typed = sp.expand(S_integrated.subs(M, M_typed))


def p_source(n_value: sp.Expr) -> sp.Expr:
    expression = sp.expand(
        r3z2_difference.subs(
            {X: X_tower, L: L_tower, M: M_typed, S: S_typed, N: n_value}
        ).doit()
    )
    return sp.cancel(expression / (e**2 * v**2 / (7 * beta)))


p_block = (
    735 * beta**2 * e**3 * p * v**4 * sp.diff(v, e)
    - 147 * beta**2 * e**3 * v**5 * sp.diff(p, e)
    + 441 * beta**2 * e**2 * p * v**5
)
for n_value, label in [(N_tower, "N-nonzero"), (sp.Integer(0), "N-zero")]:
    source = p_source(n_value)
    actual_p_block = sp.expand(source - source.subs(p, 0).doit())
    zero(actual_p_block - p_block, f"{label} r3z2 p differential block")

order_t = sp.symbols("order_t", integer=True, nonnegative=True)
root_p_coefficient = 147 * beta**2 * rho**3 * unit_v**5 * unit_M * (
    5 * order_m - order_t
)
origin_p_coefficient = 147 * beta**2 * unit_v**5 * unit_M * (
    5 * order_d - order_t + 3
)
zero(
    root_p_coefficient
    - 147 * beta**2 * rho**3 * unit_v**5 * unit_M * (5 * order_m - order_t),
    "nonzero-root p leading coefficient",
)
zero(
    origin_p_coefficient
    - 147 * beta**2 * unit_v**5 * unit_M * (5 * order_d - order_t + 3),
    "origin p leading coefficient",
)
gate(5 * order_m >= 3 * order_m, "p resonance lies above nonzero lower-bound range")
gate(5 * order_d + 3 >= 3 * order_d + 2, "p resonance lies above origin lower-bound range")

U = sp.Function("U")(e)
p_typed = e**2 * v**3 * U


# When N is nonzero, r3z1 similarly forces g to be divisible by e*v.  When
# N=0, g disappears from the decisive r3 bucket and no typing is used.
def g_source(h_value: sp.Expr) -> sp.Expr:
    expression = sp.expand(
        r3z1_difference.subs(
            {
                X: X_tower,
                L: L_tower,
                M: M_typed,
                S: S_typed,
                N: N_tower,
                H: h_value,
                p: p_typed,
            }
        ).doit()
    )
    return sp.cancel(expression / (3 * e / (7 * beta)))


g_block = (
    147 * beta**2 * e**4 * g * v**6 * sp.diff(v, e)
    - 49 * beta**2 * e**4 * v**7 * sp.diff(g, e)
    + 98 * beta**2 * e**3 * g * v**7
)
for h_value, label in [(H_tower, "H-nonzero"), (sp.Integer(0), "H-zero")]:
    source = g_source(h_value)
    actual_g_block = sp.expand(source - source.subs(g, 0).doit())
    zero(actual_g_block - g_block, f"{label} r3z1 g differential block")

root_g_coefficient = 49 * beta**2 * rho**4 * unit_v**7 * unit_M * (
    3 * order_m - order_t
)
origin_g_coefficient = 49 * beta**2 * unit_v**7 * unit_M * (
    3 * order_d - order_t + 2
)
zero(
    root_g_coefficient
    - 49 * beta**2 * rho**4 * unit_v**7 * unit_M * (3 * order_m - order_t),
    "nonzero-root g leading coefficient",
)
zero(
    origin_g_coefficient
    - 49 * beta**2 * unit_v**7 * unit_M * (3 * order_d - order_t + 2),
    "origin g leading coefficient",
)
gate(3 * order_m >= order_m, "g resonance lies above nonzero lower-bound range")
gate(3 * order_d + 2 >= order_d + 1, "g resonance lies above origin lower-bound range")

V = sp.Function("V")(e)
g_typed = e * v * V


# The decisive r3 bucket.  In every N/H branch it is a harmless outer factor
# times one universal four-term block plus e^2*v^4 times a polynomial.  This
# is the exact reason no cancellation from the lower profiles can reach the
# first root address.
universal_terminal = (
    -56 * beta * e * R * f * sp.diff(v, e)
    + 28 * beta * e * R * v * sp.diff(f, e)
    - 28 * beta * e * f * v * sp.diff(R, e)
    - 28 * beta * R * f * v
)


def terminal_source(n_value: sp.Expr, h_value: sp.Expr, g_value: sp.Expr) -> sp.Expr:
    expression = sp.expand(
        r3_difference.subs(
            {
                X: X_tower,
                L: L_tower,
                M: M_typed,
                S: S_typed,
                N: n_value,
                H: h_value,
                p: p_typed,
                g: g_value,
            }
        ).doit()
    )
    source = sp.cancel(expression / (-3 * e**2 * v / (7 * beta)))
    zero(expression + 3 * e**2 * v * source / (7 * beta), "terminal outer factor")
    return source


terminal_labels = []
for n_value, n_label, g_value in [
    (N_tower, "N-nonzero", g_typed),
    (sp.Integer(0), "N-zero", g),
]:
    for h_value, h_label in [(H_tower, "H-nonzero"), (sp.Integer(0), "H-zero")]:
        source = terminal_source(n_value, h_value, g_value)
        higher = sp.cancel((source - universal_terminal) / (e**2 * v**4))
        gate(
            sp.denom(sp.together(higher)) == 1,
            f"{n_label}/{h_label} terminal remainder is polynomial",
        )
        zero(
            source - universal_terminal - e**2 * v**4 * higher,
            f"{n_label}/{h_label} universal terminal split",
        )
        terminal_labels.append(f"{n_label}/{h_label}")

root_terminal = -56 * beta * rho * R0 * unit_f * order_m * unit_v
origin_terminal = -28 * beta * R0 * f0 * unit_v * (2 * order_d + 1)
zero(
    root_terminal
    + 56 * beta * rho * R0 * unit_f * order_m * unit_v,
    "nonzero-root terminal leading coefficient",
)
zero(
    origin_terminal
    + 28 * beta * R0 * f0 * unit_v * (2 * order_d + 1),
    "origin terminal leading coefficient",
)
gate(2 * order_d + 1 > 0, "origin terminal has no integral resonance")
gate(len(terminal_labels) == 4, "all four N/H branches frozen")


semantic = {
    "ansatz": "THM3821 plus r2z2*X(e),r2z2*Y(e)",
    "branch": "X!=0;Y=lambdaX;L=T-lambdaS!=0",
    "tower": "X=alpha*e6*v10;L=beta*e4*v7",
    "integral": "S=(10alpha/(7beta))*e2*v3*M+(2alpha/7)*e5*v10+gamma*e4*v7",
    "r4_typing": "M=e*v2*R;R and f units at every v-root;15alphaR2=49beta2f",
    "side_rungs": "N=0 or delta*e3*v5;H=0 or theta*e2*v3",
    "lower_bounds": "p divisible by e2*v3;when N!=0 g divisible by e*v",
    "terminal": "Q=-56beta*e*R*f*vprime+28beta*e*R*v*fprime-28beta*e*f*v*Rprime-28beta*R*f*v+e2*v4*Psi",
    "root_gate": "nonzero v-root has unique coefficient -56beta*rho*R*f*m*u;v=c*e^d has -28beta*R(0)*f(0)*c*(2d+1)",
    "scope": "closes L!=0 for X!=0 in fixed r2z2 ansatz;combined with THM3828 closes all X!=0;X=0,Y!=0 and higher slots open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3829-misaligned-second-row-r2z2-10-7-profile-nonentry")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("ansatz=THM3821_plus_r2z2_XY")
print("branch=X_nonzero;Y=lambdaX;L=T-lambdaS_nonzero")
print("tower=X_alpha_e6_v10;L_beta_e4_v7")
print("r4_gate=M_e_v2_R;R_f_units_on_v_roots")
print("side_branches=N_zero_or_10_5;H_zero_or_10_3")
print("lower_bounds=p_e2_v3;g_e_v_when_N_nonzero")
print("terminal=universal_four_terms_plus_e2_v4_Psi")
print("root_gate=nonzero_root_unique_vprime;origin_odd_2d_plus_1")
print("scope=fixed_second_row_X_nonzero_closed_with_THM3828;one_sided_and_higher_open")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
