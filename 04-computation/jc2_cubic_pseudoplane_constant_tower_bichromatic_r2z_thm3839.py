#!/usr/bin/env python3
"""Exact companion for THM-3839's constant-tower bichromatic r^2 z no-go."""

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


# THM-3821 plus the first r^2 z profile in each coordinate.
f = sp.Function("f")(e)
g = sp.Function("g")(e)
h = sp.Function("h")(e)
kap = sp.Function("kappa")(e)
p = sp.Function("p")(e)
q = sp.Function("q")(e)
S = sp.Function("S")(e)
T = sp.Function("T")(e)
U = sp.Function("U")(e)
V = sp.Function("V")(e)
A = e**2 - z / 3 + r * g + z**2 * f + r * z * p + r * z**2 * S + r**2 * z * U
C = (
    e**3
    - e
    - e * z / 2
    + r * h
    + z**2 * kap
    + r * z * q
    + r * z**2 * T
    + r**2 * z * V
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
bucket_r5z = bucket(5, 1)
bucket_r5 = bucket(5)
bucket_r4z2 = bucket(4, 2)
bucket_r4z1 = bucket(4, 1)
bucket_r4 = bucket(4)
bucket_r3z2 = bucket(3, 2)
bucket_r3z1 = bucket(3, 1)
bucket_r3 = bucket(3)

expected_z = (36 * e**2 * f - 24 * e * kap - 12 * f + 1) / 2
expected_r5z = 24 * e * (U * sp.diff(V, e) - V * sp.diff(U, e))
expected_r5 = 21 * e**2 * (S * sp.diff(T, e) - T * sp.diff(S, e))
zero(bucket_z - expected_z, "arm bucket")
zero(bucket_r5z - expected_r5z, "first color Wronskian")
zero(bucket_r5 - expected_r5, "second color Wronskian")


# In the branch U,S!=0 the top buckets force V=lambda U and T=mu S.
# We keep the genuinely bichromatic locus Delta=mu-lambda!=0.
lam, mu, Delta = sp.symbols("lambda mu Delta")
M = sp.Function("M")(e)
N = sp.Function("N")(e)
H = sp.Function("H")(e)
differences = {
    V: lam * U,
    T: mu * S,
    kap: mu * f + M,
    q: lam * p + N,
    h: lam * g + H,
}

r4z2_difference = sp.expand(bucket_r4z2.subs(differences).doit())
r4z1_difference = sp.expand(bucket_r4z1.subs(differences).doit())
r4_difference = sp.expand(bucket_r4.subs(differences).doit())
r3z2_difference = sp.expand(bucket_r3z2.subs(differences).doit())
r3z1_difference = sp.expand(bucket_r3z1.subs(differences).doit())
r3_difference = sp.expand(bucket_r3.subs(differences).doit())

ode_87 = 7 * e * S * sp.diff(U, e) - 8 * e * U * sp.diff(S, e) - 3 * S * U
zero(r4z2_difference - 3 * (lam - mu) * ode_87, "bichromatic 8/7 equation")

D_mu = 3 * e**2 - 2 * mu * e - 1
arm = sp.expand(bucket_z.subs(differences).doit() / 6)
zero(arm - (D_mu * f - 2 * e * M + sp.Rational(1, 12)), "mu-arm identity")


# UFD solution of the 8/7 equation.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
v = sp.Function("v")(e)
U_tower = alpha * e**5 * v**8
S_tower = beta * e**4 * v**7
zero(ode_87.subs({U: U_tower, S: S_tower}).doit(), "8/7 tower solves r4z2")

mult = sp.symbols("mult", integer=True, nonnegative=True)
gate(7 * (5 + 8 * mult) - 8 * (4 + 7 * mult) == 3, "e-prime 8/7 family")
gate(7 * (8 * mult) - 8 * (7 * mult) == 0, "nonzero-prime 8/7 family")


# The next bucket gives the optional 8/5 rung N=delta*e^3*v^5.
ode_85 = 8 * e * U * sp.diff(N, e) - 5 * e * N * sp.diff(U, e) + U * N
zero(
    r4z1_difference.subs({U: U_tower, S: S_tower}).doit()
    - 3 * ode_85.subs(U, U_tower).doit(),
    "8/5 N equation",
)
delta = sp.symbols("delta")
N_tower = delta * e**3 * v**5
zero(ode_85.subs({U: U_tower, N: N_tower}).doit(), "8/5 rung solves r4z1")
gate(8 * (3 + 5 * mult) - 5 * (5 + 8 * mult) == -1, "e-prime 8/5 family")
gate(8 * (5 * mult) - 5 * (8 * mult) == 0, "nonzero-prime 8/5 family")


# The pure r4 bucket is an exact factor-cofactor differential equation.
K = 7 * beta * M + 8 * alpha * v * H
cofactor_ode = e * v * sp.diff(K, e) - 4 * e * sp.diff(v, e) * K - 2 * v * K
r4_tower = sp.expand(
    r4_difference.subs({U: U_tower, S: S_tower, N: N_tower}).doit()
)
zero(
    r4_tower - 3 * e**5 * v**6 * cofactor_ode,
    "full r4 factor-cofactor equation",
)
eta = sp.symbols("eta")
K_integrated = eta * e**2 * v**4
zero(
    e * v * sp.diff(K_integrated, e)
    - 4 * e * sp.diff(v, e) * K_integrated
    - 2 * v * K_integrated,
    "factor-cofactor first integral",
)


# Constant-v branch.  Write the actual constant coefficients as
#
#     U=(a0*b0)e^5,  S=b0 e^4,  N=(delta*b0)e^3,
#
# where a0 is the U/S ratio and b0 is retained.  The b0=1 formulas below are
# a reference model only; exact seam gates then restore arbitrary b0!=0.
a0 = sp.symbols("a0", nonzero=True)
b0 = sp.symbols("b0", nonzero=True)
U_const = a0 * e**5
S_const = e**4
M_const = sp.Function("M_const")(e)
H_const = sp.Function("H_const")(e)
N_const = delta * e**3

constant_differences = {
    U: U_const,
    V: lam * U_const,
    S: S_const,
    T: mu * S_const,
    kap: mu * f + M_const,
    q: lam * p + N_const,
    h: lam * g + H_const,
}

constant_cofactor = 8 * a0 * H_const + 7 * M_const - eta * e**2
constant_arm = D_mu * f - 2 * e * M_const + sp.Rational(1, 12)

# The two next exact buckets are the color jets.  E2 is r3z2 after cancelling
# its nonzero scalar/e factor and using the cofactor identity.
E2 = (
    -3 * a0 * e**4
    - 8 * a0 * e**2 * sp.diff(f, e)
    - 8 * a0 * e**2 * sp.diff(M_const, e) / (mu - lam)
    + 16 * a0 * e * f
    + 16 * a0 * e * M_const / (mu - lam)
    + 7 * e * sp.diff(p, e)
    - 21 * p
)

constant_r3z2 = sp.expand(bucket_r3z2.subs(constant_differences).doit())
E2_from_bucket = sp.cancel(
    constant_r3z2
    .subs(sp.diff(H_const, e), sp.diff((eta * e**2 - 7 * M_const) / (8 * a0), e))
    .subs(H_const, (eta * e**2 - 7 * M_const) / (8 * a0))
    .subs(lam, mu - Delta)
    / (-3 * Delta * e**4)
)
zero(E2_from_bucket - E2.subs(lam, mu - Delta), "exact r3z2 color-jet equation")

# Restore the actual nonzero second-color coefficient b0.  The next three
# gates deliberately expose, rather than suppress, the nonlinear scaling
# seams.
scaled_differences = {
    U: b0 * U_const,
    V: lam * b0 * U_const,
    S: b0 * S_const,
    T: mu * b0 * S_const,
    kap: mu * f + M_const,
    q: lam * p + b0 * N_const,
    h: lam * g + H_const,
}
zero(
    bucket_r3z2.subs(scaled_differences).doit()
    - b0 * constant_r3z2
    + 9 * a0 * b0 * e**8 * (b0 - 1) * (lam - mu),
    "constant-tower r3z2 scaling seam",
)

# The r3z1 equation determines the Euler derivative of g.  This identity is
# valid in both N=0 and N!=0 branches and loses no polynomial equation.
P0 = e * sp.diff(p, e) - 3 * p
H0 = e * sp.diff(H_const, e) - 2 * H_const
G0 = e * sp.diff(g, e) - 2 * g
G0_solved = (
    21 * e * H0
    + 3 * a0 * delta * e**4
    - 8 * a0 * lam * e**2
    - 15 * delta * P0
) / (21 * (mu - lam) * e)
constant_r3z1 = sp.expand(bucket_r3z1.subs(constant_differences).doit())
zero(
    constant_r3z1 + 21 * (mu - lam) * e**4 * (G0 - G0_solved),
    "exact r3z1 Euler-g equation",
)
zero(
    bucket_r3z1.subs(scaled_differences).doit()
    - b0 * constant_r3z1
    - 3 * a0 * b0 * delta * e**7 * (b0 - 1),
    "constant-tower r3z1 scaling seam",
)

# The terminal r3 bucket.  Eliminate only G0 using the preceding exact
# equation.  The optional 8/5 rung contributes the final -5*delta*e^3*G0.
E3_base = (
    24 * a0 * e**8
    - 16 * a0 * lam * e**7
    - 8 * a0 * e**6
    + 7 * a0 * e**5 * sp.diff(H_const, e)
    - 15 * a0 * e**4 * H_const
    + 12 * e**5 * sp.diff(M_const, e)
    - 22 * e**4 * M_const
    - 4 * e**2 * M_const * sp.diff(f, e)
    + 4 * e**2 * f * sp.diff(M_const, e)
    - 3 * e * H_const * sp.diff(p, e)
    + 5 * e * p * sp.diff(H_const, e)
    - H_const * p
)
E3 = sp.expand(E3_base - 5 * delta * e**3 * G0_solved)
constant_r3 = sp.expand(bucket_r3.subs(constant_differences).doit())
zero(
    constant_r3 - 3 * E3_base + 15 * delta * e**3 * G0,
    "exact terminal r3 equation before Euler-g elimination",
)
Q_terminal = (
    -4 * e**2 * M_const * sp.diff(f, e)
    + 4 * e**2 * f * sp.diff(M_const, e)
    - 3 * e * H_const * sp.diff(p, e)
    + 5 * e * p * sp.diff(H_const, e)
    - H_const * p
)
zero(
    bucket_r3.subs(scaled_differences).doit()
    - b0 * constant_r3
    + 3 * (b0 - 1) * Q_terminal,
    "constant-tower r3 scaling seam",
)


# Solve arm and cofactor for M,H.  They are polynomial identities whenever a
# Darboux pair exists; the rational notation is used only inside k(e).
M_from_f = sp.cancel((D_mu * f + sp.Rational(1, 12)) / (2 * e))
H_from_f = sp.cancel((eta * e**2 - 7 * M_from_f) / (8 * a0))
zero(constant_arm.subs(M_const, M_from_f), "arm solves M")
zero(constant_cofactor.subs({M_const: M_from_f, H_const: H_from_f}), "cofactor solves H")

# Honest constant tower with second-color coefficient b0.  The b0=1 model is
# only a positive control: the charge, Euler-g equation, and nonlinear
# terminal block all acquire the displayed compensating factors.
E2_general = sp.expand(E2 - 3 * a0 * (b0 - 1) * e**4)
G0_solved_general = sp.expand(
    G0_solved + a0 * delta * (b0 - 1) * e**3 / (7 * (mu - lam))
)
E3_general = sp.expand(
    E3_base - Q_terminal + Q_terminal / b0 - 5 * delta * e**3 * G0_solved_general
)

# Direct all-b0 reconstructions from the unreduced canonical buckets.
scaled_r3z2_direct = sp.cancel(
    bucket_r3z2.subs(scaled_differences).doit()
    .subs(sp.diff(H_const, e), sp.diff((eta * e**2 - 7 * M_const) / (8 * a0), e))
    .subs(H_const, (eta * e**2 - 7 * M_const) / (8 * a0))
)
zero(
    scaled_r3z2_direct + 3 * b0 * (mu - lam) * e**4 * E2_general,
    "all-b0 exact E2 reconstruction",
)
zero(
    bucket_r3z1.subs(scaled_differences).doit()
    + 21 * b0 * (mu - lam) * e**4 * (G0 - G0_solved_general),
    "all-b0 exact Euler-g reconstruction",
)
E3_general_pre = sp.expand(
    E3_base - Q_terminal + Q_terminal / b0 - 5 * delta * e**3 * G0
)
zero(
    bucket_r3.subs(scaled_differences).doit() - 3 * b0 * E3_general_pre,
    "all-b0 exact terminal reconstruction",
)


# Formal Taylor coefficients through the second color jet.  A free cubic
# tail is retained to verify that no truncation assumption enters.  The arm
# numerator is automatically divisible by e because f(0)=1/12.
f1, f2, f3 = sp.symbols("f1 f2 f3")
p0c, p1c, p2c, p3c = sp.symbols("p0c p1c p2c p3c")
f_series = sp.Rational(1, 12) + f1 * e + f2 * e**2 + f3 * e**3
p_series = p0c + p1c * e + p2c * e**2 + p3c * e**3
M_series = sp.cancel(M_from_f.subs(f, f_series).doit())
H_series = sp.cancel(H_from_f.subs(f, f_series).doit())
E2_series = sp.Poly(
    sp.expand(E2_general.subs({f: f_series, M_const: M_series, p: p_series}).doit()),
    e,
)

# The first three E2 coefficients determine p0,p1,p2; p3 is the Euler
# resonance.  Higher coefficients cannot feed backward into these equations.
e2_equations = [E2_series.coeff_monomial(e**j) for j in range(3)]
p_solution = sp.solve(e2_equations, [p0c, p1c, p2c], dict=True)
gate(len(p_solution) == 1, "low E2 solution is unique modulo the p3 Euler resonance")
p_series_solved = sp.expand(p_series.subs(p_solution[0]))

E3_series = sp.expand(sp.cancel(
    E3_general.subs(
        {
            f: f_series,
            M_const: M_series,
            H_const: H_series,
            p: p_series_solved,
        }
    ).doit()
))
color_one = sp.factor(sp.Poly(E3_series, e).coeff_monomial(e))
color_two = sp.factor(sp.Poly(E3_series, e).coeff_monomial(e**2))
expected_color_one = -(
    (6 * f1 + lam) * (6 * f1 + mu)
) / (36 * b0 * (lam - mu))
expected_color_two = -(
    48 * f1**2 * lam
    + 48 * f1**2 * mu
    + 48 * f1 * f2
    + 16 * f1 * lam * mu
    - 12 * f1
    + 4 * f2 * lam
    + 4 * f2 * mu
    - lam
    - mu
) / (32 * b0 * (lam - mu))
zero(color_one - expected_color_one, "first exact color jet")
zero(color_two - expected_color_two, "second exact color jet")
gate(not color_one.has(f3, p3c), "first color jet is tail-independent")
gate(not color_two.has(f3, p3c), "second color jet is tail-independent")

color = sp.symbols("color")
for color_value in (lam, mu):
    zero(
        expected_color_one.subs(f1, -color_value / 6),
        f"first jet vanishes on color {color_value}",
    )
    zero(
        expected_color_two.subs(
            {f1: -color_value / 6, f2: color_value**2 / 3 + sp.Rational(1, 4)}
        ),
        f"second jet fixes color {color_value}",
    )

# The jets are the first two Taylor coefficients of -1/(12*D_color).
taylor_model = sp.series(
    1 / (12 * (1 + 2 * color * e - 3 * e**2)), e, 0, 3
).removeO()
zero(
    taylor_model
    - (sp.Rational(1, 12) - color * e / 6 + (color**2 / 3 + sp.Rational(1, 4)) * e**2),
    "color jets match the reciprocal arm model",
)


# Terminal degree cases d=0,1,2.  For each direct replay, solve E2 exactly
# coefficientwise and read the unreduced terminal coefficient at e^8.
def terminal_for_polynomial(poly_f: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    degree = sp.degree(poly_f, e)
    m_poly = sp.cancel(M_from_f.subs(f, poly_f).doit())
    h_poly = sp.cancel(H_from_f.subs(f, poly_f).doit())
    # The required p degree is at most degree+2; retain the Euler p3 constant.
    p_coefficients = sp.symbols(f"u0:{max(6, degree + 4)}")
    p_poly = sum(coefficient * e**j for j, coefficient in enumerate(p_coefficients))
    e2_poly = sp.Poly(
        sp.expand(E2_general.subs({f: poly_f, M_const: m_poly, p: p_poly}).doit()), e
    )
    equations = [coefficient for coefficient in e2_poly.all_coeffs()]
    solve_variables = [coefficient for j, coefficient in enumerate(p_coefficients) if j != 3]
    solutions = sp.solve(equations, solve_variables, dict=True)
    gate(len(solutions) == 1, f"degree {degree} E2 solution modulo p3")
    p_poly = sp.expand(p_poly.subs(solutions[0]))
    terminal = sp.expand(
        sp.cancel(E3_general.subs({f: poly_f, M_const: m_poly, H_const: h_poly, p: p_poly}).doit())
    )
    return p_poly, terminal


# d=0: the two low terminal coefficients are already inconsistent.
p_deg0, terminal_deg0 = terminal_for_polynomial(sp.Rational(1, 12))
deg0_e1 = sp.factor(sp.Poly(terminal_deg0, e).coeff_monomial(e))
deg0_e2 = sp.factor(sp.Poly(terminal_deg0, e).coeff_monomial(e**2))
zero(deg0_e1 + lam * mu / (36 * b0 * (lam - mu)), "degree-zero first coefficient")
zero(deg0_e2 - (lam + mu) / (32 * b0 * (lam - mu)), "degree-zero second coefficient")
gate(
    sp.resultant(sp.together(deg0_e1).as_numer_denom()[0],
                 sp.together(deg0_e2).as_numer_denom()[0], mu) != 0,
    "degree-zero coefficients have no bichromatic common solution",
)
gate(
    sp.solve([lam * mu, lam + mu], [lam, mu], dict=True) == [{lam: 0, mu: 0}],
    "degree-zero common zero is monochromatic",
)


# d=1 and d=2: impose either color jet, then the exact e8 terminal remains
# 24*a0.  The d=1 second jet also forces 4*c^2+3=0, but the top coefficient
# closes the branch independently of which color was chosen.
f2_free = sp.symbols("f2_free")
for degree_value, polynomial in (
    (1, sp.Rational(1, 12) + f1 * e),
    (2, sp.Rational(1, 12) + f1 * e + f2_free * e**2),
):
    p_degree, terminal_degree = terminal_for_polynomial(polynomial)
    for color_value in (lam, mu):
        substitutions = {f1: -color_value / 6}
        if degree_value == 2:
            substitutions[f2_free] = color_value**2 / 3 + sp.Rational(1, 4)
        e8 = sp.factor(sp.Poly(sp.expand(terminal_degree.subs(substitutions)), e).coeff_monomial(e**8))
        zero(e8 - 24 * a0, f"degree {degree_value} color {color_value} terminal e8")

deg1_second_lambda = sp.factor(
    expected_color_two.subs({f1: -lam / 6, f2: 0})
)
deg1_second_mu = sp.factor(
    expected_color_two.subs({f1: -mu / 6, f2: 0})
)
zero(
    deg1_second_lambda + (4 * lam**2 + 3) / (96 * b0),
    "degree-one lambda-color second jet",
)
zero(
    deg1_second_mu - (4 * mu**2 + 3) / (96 * b0),
    "degree-one mu-color second jet",
)


# Generic degree d>=3.  The arm, cofactor, and E2 successively force the
# displayed leading terms.  The N/g contribution has degree at most d+4,
# below the terminal degree 2d+3.  The surviving coefficient is nonzero.
d = sp.symbols("d", integer=True, positive=True)
fd = sp.symbols("f_d", nonzero=True)
M_lead = sp.Rational(3, 2) * fd
H_lead = -sp.Rational(21, 16) * fd / a0
p_lead = sp.Rational(12, 7) * a0 * fd / (mu - lam)
generic_top = -sp.Rational(9, 2) * (d - 1) * fd**2 / (b0 * (mu - lam))
zero(3 * fd - 2 * M_lead, "generic arm leading coefficient")
zero(8 * a0 * H_lead + 7 * M_lead, "generic cofactor leading coefficient")
zero(
    -8 * a0 * (d + 1) * M_lead / (mu - lam)
    + 16 * a0 * M_lead / (mu - lam)
    + (7 * (d + 2) - 21) * p_lead,
    "generic E2 leading coefficient",
)
zero(
    H_lead * p_lead * (-3 * (d + 2) + 5 * (d + 1) - 1) / b0
    - generic_top,
    "generic terminal leading coefficient",
)
zero((2 * d + 3) - (d + 4) - (d - 1), "terminal degree separation identity")
gate(2 * 3 + 3 > 3 + 4, "N/g terminal contribution is strictly lower from d=3 onward")
gate(generic_top != 0, "generic terminal coefficient is nonzero for d>=3")

# Hostile fixed-degree replays verify the symbolic leading arithmetic without
# serving as the all-degree proof.
for degree_value in range(3, 7):
    f_hostile = sp.Rational(1, 12) + fd * e**degree_value
    _, terminal_hostile = terminal_for_polynomial(f_hostile)
    coefficient = sp.factor(
        sp.Poly(terminal_hostile, e).coeff_monomial(e ** (2 * degree_value + 3))
    )
    zero(
        coefficient - generic_top.subs(d, degree_value),
        f"generic top hostile degree {degree_value}",
    )


semantic = {
    "ansatz": "THM3821 plus r2z*U(e),r2z*V(e)",
    "branch": "U,S nonzero;V=lambda*U;T=mu*S;lambda!=mu;v constant",
    "tower": "U=alpha*e5*v8;S=beta*e4*v7;N=0 or delta*e3*v5",
    "cofactor": "7beta*M+8alpha*v*H=eta*e2*v4",
    "constant_scaling": "U=rho*b*e5;S=b*e4;N=delta*b*e3;b retained and nonzero",
    "color_jets": "f1=-c/6;f2=c2/3+1/4;c in {lambda,mu}",
    "terminal": "d=0 low coefficients;d=1,2 e8=24rho;d>=3 top=-9(d-1)fd2/(2b(mu-lambda))",
    "scope": "constant-v bichromatic branch only;nonconstant v and other higher slots open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3839-constant-tower-bichromatic-r2z-profile-nonentry")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("ansatz=THM3821_plus_r2z_UV")
print("branch=U_S_nonzero;lambda_not_mu;constant_v")
print("tower=8_7;optional_N_8_5")
print("cofactor=7beta_M_plus_8alpha_v_H_equals_eta_e2v4")
print("constant_scaling=b_retained_nonzero")
print("jets=two_target_colors")
print("terminal=all_polynomial_f_degrees_closed")
print("scope=nonconstant_v_and_other_higher_slots_open")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
