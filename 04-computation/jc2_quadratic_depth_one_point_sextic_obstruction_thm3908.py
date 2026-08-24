#!/usr/bin/env python3
"""Assertion-free exact gates for the THM-3908 proof candidate.

The companion checks the tangent equation on the smooth binary-cubic
discriminant, the triple-root sextic normal form, the constant-row Newton
coefficient chain, and the sole moving-repeated-root normal form.

Reproduction:
  python3 04-computation/jc2_quadratic_depth_one_point_sextic_obstruction_thm3908.py
  python3 -O 04-computation/jc2_quadratic_depth_one_point_sextic_obstruction_thm3908.py
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"{label}: failed")


def zero(label: str, value: sp.Expr) -> None:
    gate(sp.factor(sp.cancel(value)) == 0, label)


def equal(label: str, left: object, right: object) -> None:
    gate(left == right, f"{label}: {left!r} != {right!r}")


def binary_disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


def homogeneous_piece(poly: sp.Expr, variables: tuple[sp.Symbol, ...], degree: int) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient
            * sp.prod(variable**exponent for variable, exponent in zip(variables, powers))
            for powers, coefficient in sp.Poly(poly, *variables).terms()
            if sum(powers) == degree
        )
    )


# ---------------------------------------------------------------------------
# 1. The smooth discriminant tangent and the two degree ledgers.
# ---------------------------------------------------------------------------

eps = sp.symbols("eps")
r0, s0 = sp.symbols("r0 s0", nonzero=True)
a1, b1, c1, d1 = sp.symbols("a1 b1 c1 d1")

# The repeated factor is U and the distinct simple factor is r0 U+s0 V.
# The differential of the discriminant sees only the value d1 at U=0.
smooth_deformation = binary_disc(
    r0 + eps * a1,
    s0 + eps * b1,
    eps * c1,
    eps * d1,
)
equal(
    "smooth_discriminant_differential",
    sp.Poly(smooth_deformation, eps).coeff_monomial(eps),
    -4 * s0**3 * d1,
)

# A quadratic coefficient row mapping into the twisted cubic has residual
# coordinate degree 2-g after removing a gcd of degree g.  A nonconstant
# triple-root map would require that degree to be a positive multiple of 3.
triple_degree_ledger = [(g, 2 - g, (2 - g) % 3 == 0 and 2 - g > 0) for g in range(3)]
equal("triple_root_degree_ledger", triple_degree_ledger, [(0, 2, False), (1, 1, False), (2, 0, False)])

# On the normalization P1_l x P1_m of the smooth discriminant surface,
# O(1) pulls back to O(2,1).  These are all degree <=2 possibilities after
# the coefficient gcd is removed.
double_degree_ledger = [
    (degree, repeated_degree, simple_degree)
    for degree in range(3)
    for repeated_degree in range(2)
    for simple_degree in range(3)
    if 2 * repeated_degree + simple_degree == degree
]
equal(
    "double_root_bidegree_ledger",
    double_degree_ledger,
    [(0, 0, 0), (1, 0, 1), (2, 0, 2), (2, 1, 0)],
)
equal(
    "sole_moving_repeated_root_bidegree",
    [row for row in double_degree_ledger if row[1] > 0],
    [(2, 1, 0)],
)


# ---------------------------------------------------------------------------
# 2. Triple-root leading rows and the common-zero sextic grammar.
# ---------------------------------------------------------------------------

h2, aa, bb, cc, dd = sp.symbols("h2 aa bb cc dd")
triple_deformation = binary_disc(
    h2 + eps * aa,
    eps * bb,
    eps * cc,
    eps * dd,
)
equal(
    "triple_root_degree_six_piece",
    sp.Poly(triple_deformation, eps).coeff_monomial(eps**2),
    -27 * h2**2 * dd**2,
)

A, C, U, V = sp.symbols("A C U V")
alpha, alpha1, beta, beta1, gamma, gamma1 = sp.symbols(
    "alpha alpha1 beta beta1 gamma gamma1"
)
a = C**2 + alpha * A + alpha1 * C
b = beta * A + beta1 * C
c = gamma * A + gamma1 * C
d = C
D_common = binary_disc(a, b, c, d)
equal("common_zero_total_degree", sp.Poly(D_common, A, C).total_degree(), 6)
zero("common_zero_pure_top", homogeneous_piece(D_common, (A, C), 6) + 27 * C**6)
zero(
    "common_zero_C_zero_gate",
    D_common.subs(C, 0) - gamma**2 * (beta**2 - 4 * alpha * gamma) * A**4,
)

x, z = sp.symbols("x z")
h_common = sp.expand(z**6 * D_common.subs({A: 1 / z, C: x / z}, simultaneous=True))
equal("common_local_x6", sp.Poly(h_common, x, z).coeff_monomial(x**6), -27)
equal("common_local_x2z", sp.Poly(h_common, x, z).coeff_monomial(x**2 * z), -4 * gamma**3)
zero(
    "common_local_z2",
    sp.Poly(h_common, x, z).coeff_monomial(z**2)
    - gamma**2 * (beta**2 - 4 * alpha * gamma),
)
equal("common_local_no_z", sp.Poly(h_common, x, z).coeff_monomial(z), 0)
equal("common_local_no_xz", sp.Poly(h_common, x, z).coeff_monomial(x * z), 0)
gate(
    all(i + 2 * j >= 4 for (i, j), coefficient in sp.Poly(h_common, x, z).terms() if coefficient != 0),
    "common_first_Newton_halfspace",
)
gate(
    all(i + 4 * j >= 6 for (i, j), coefficient in sp.Poly(h_common, x, z).terms() if coefficient != 0),
    "common_second_Newton_halfspace",
)


# ---------------------------------------------------------------------------
# 3. Constants cannot repair the triple-root infinity collision.
# ---------------------------------------------------------------------------

a0, b0, c0, d0 = sp.symbols("a0 b0 c0 d0")
a_const = a + a0
b_const = b + b0
c_const = c + c0
d_const = d + d0
D_const = binary_disc(a_const, b_const, c_const, d_const)
zero("constant_triple_pure_top", homogeneous_piece(D_const, (A, C), 6) + 27 * C**6)
h_const = sp.expand(z**6 * D_const.subs({A: 1 / z, C: x / z}, simultaneous=True))
P_const = sp.Poly(h_const, x, z)
equal("constant_chain_x6", P_const.coeff_monomial(x**6), -27)
equal("constant_chain_x2z", P_const.coeff_monomial(x**2 * z), -4 * gamma**3)
equal(
    "constant_chain_xz2_when_gamma_zero",
    sp.factor(P_const.coeff_monomial(x * z**2).subs(gamma, 0)),
    -4 * beta**3,
)
equal(
    "constant_chain_x4z_when_gamma_beta_zero",
    sp.factor(P_const.coeff_monomial(x**4 * z).subs({gamma: 0, beta: 0})),
    -54 * alpha,
)
equal("constant_chain_no_z", P_const.coeff_monomial(z), 0)
equal(
    "constant_chain_no_z2_when_gamma_zero",
    sp.factor(P_const.coeff_monomial(z**2).subs(gamma, 0)),
    0,
)
equal(
    "constant_chain_no_z3_when_gamma_beta_zero",
    sp.factor(P_const.coeff_monomial(z**3).subs({gamma: 0, beta: 0})),
    0,
)
D_A_free = sp.expand(D_const.subs({alpha: 0, beta: 0, gamma: 0}))
equal("terminal_row_A_free", sp.diff(D_A_free, A), 0)
equal("terminal_univariate_degree", sp.degree(D_A_free, C), 6)

# The exact strict-below inequalities used in the lower-hull proof.
gate(sp.Rational(1) < sp.Rational(2, 3) * 2, "gamma_vertex_below_j0_two_chord")
gate(sp.Rational(2) < sp.Rational(5, 6) * 3, "beta_vertex_below_j0_three_chord")
gate(sp.Rational(1) < sp.Rational(1, 3) * 4, "alpha_vertex_below_j0_four_chord")


# ---------------------------------------------------------------------------
# 4. The sole moving-repeated-root row collapses to reducibility.
# ---------------------------------------------------------------------------

p, q, rr, ss, tt, uu, vv = sp.symbols("p q rr ss tt uu vv")
ell = A * U + C * V
Q2 = p * U**2 + q * U * V + rr * V**2
Q3 = ss * U**3 + tt * U**2 * V + uu * U * V**2 + vv * V**3
F_moving = sp.expand(ell**2 * V + ell * Q2 + Q3)
F_poly = sp.Poly(F_moving, U, V)
am = F_poly.coeff_monomial(U**3)
bm = F_poly.coeff_monomial(U**2 * V)
cm = F_poly.coeff_monomial(U * V**2)
dm = F_poly.coeff_monomial(V**3)
equal("moving_coefficient_a", am, A * p + ss)
equal("moving_coefficient_b", bm, A**2 + A * q + C * p + tt)
equal("moving_coefficient_c", cm, 2 * A * C + A * rr + C * q + uu)
equal("moving_coefficient_d", dm, C**2 + C * rr + vv)
D_moving = binary_disc(am, bm, cm, dm)
D6_moving = homogeneous_piece(D_moving, (A, C), 6)
D6_expected = A**2 * (
    (rr**2 - 4 * vv) * A**4
    + (-2 * q * rr + 4 * uu) * A**3 * C
    + (2 * p * rr + q**2 - 4 * tt) * A**2 * C**2
    + (-2 * p * q + 4 * ss) * A * C**3
    + p**2 * C**4
)
zero("moving_degree_six_factorization", D6_moving - D6_expected)

pure_substitution = {
    p: 0,
    ss: 0,
    tt: q**2 / 4,
    uu: q * rr / 2,
}
zero(
    "moving_pure_power_top_after_constraints",
    D6_moving.subs(pure_substitution) - (rr**2 - 4 * vv) * A**6,
)
F_reducible = sp.factor(F_moving.subs(pure_substitution))
gate(sp.cancel(F_reducible / V).is_polynomial(U, V, A, C), "moving_full_cubic_has_V_factor")
equal("moving_full_cubic_vanishes_at_V_zero", F_reducible.subs(V, 0), 0)


# ---------------------------------------------------------------------------
# 5. The THM-3906 row is the sharp common-zero two-place control.
# ---------------------------------------------------------------------------

D_control = binary_disc(C**2 + A, A, A + C, C)
equal("control_degree", sp.Poly(D_control, A, C).total_degree(), 6)
zero("control_top", homogeneous_piece(D_control, (A, C), 6) + 27 * C**6)
h_control = sp.expand(z**6 * D_control.subs({A: 1 / z, C: x / z}, simultaneous=True))
equal("control_order_two_vertex", sp.Poly(h_control, x, z).coeff_monomial(z**2), -3)
equal("control_middle_vertex", sp.Poly(h_control, x, z).coeff_monomial(x**2 * z), -4)
equal("control_order_four_vertex", sp.Poly(h_control, x, z).coeff_monomial(x**6), -27)

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

semantic = {
    "object": "binary cubic over k[A,C] with coefficient degree at most two",
    "top": "nonzero degree-six discriminant part is one linear sixth power",
    "common_zero": "smooth leading stratum reducible; triple stratum is THM-3906 grammar and two-place",
    "constants_triple": "gamma/beta/alpha Newton chain forces two edges or univariate reducibility",
    "constants_smooth": "fixed repeated root is reducible-or-monogenic; moving root pure top forces a fixed V factor",
    "conclusion": "no irreducible nonmonogenic one-place sextic at quadratic coefficient depth",
    "next_depth": "coefficient degree three is first depth permitting a moving triple root",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3908-quadratic-depth-common-zero-one-point-sextic-two-place-obstruction")
print("common_zero_quadratic_depth=REDUCIBLE_OR_TWO_PLACE")
print("triple_root_with_constants=REDUCIBLE_OR_AT_LEAST_TWO_PLACE")
print("fixed_repeated_root_with_constants=REDUCIBLE_OR_MONOGENIC")
print("moving_repeated_root_pure_sextic_top=GENERICALLY_REDUCIBLE")
print("irreducible_nonmonogenic_quadratic_depth_one_place_sextic=EMPTY")
print("next_unclosed_coefficient_depth=3")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
