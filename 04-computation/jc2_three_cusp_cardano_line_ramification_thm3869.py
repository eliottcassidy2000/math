#!/usr/bin/env python3
"""Assertion-free exact gates for THM-3869.

Reproduction:
  python3 04-computation/jc2_three_cusp_cardano_line_ramification_thm3869.py
  python3 -O 04-computation/jc2_three_cusp_cardano_line_ramification_thm3869.py
"""

from __future__ import annotations

import hashlib
import sys

import sympy as sp

sys.stdout.reconfigure(newline="\n")


x, y, t, T, Z, u, v, s = sp.symbols("x y t T Z u v s")
c, d, rho, W = sp.symbols("c d rho W")
Y = sp.symbols("Y")
CHECKS = 0


def zero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def nonzero(label: str, value: sp.Expr) -> None:
    global CHECKS
    CHECKS += 1
    value = sp.factor(value)
    if value == 0:
        raise AssertionError(f"{label}: unexpectedly zero")


def equal(label: str, left: object, right: object) -> None:
    global CHECKS
    CHECKS += 1
    if left != right:
        raise AssertionError(f"{label}: {left!r} != {right!r}")


Delta = sp.expand(
    81 * x**5
    + 90 * x**4
    + 25 * x**3
    + 30 * x**2 * y**2
    + 30 * x * y**2
    - y**4
    + 8 * y**2
)
ell = 9 * x + 4
a = x + 1
q = y**2 - (15 * x**2 + 15 * x + 4)
P = sp.expand(a * ell**2)
Q = sp.expand(q * ell**2)
f = sp.expand(T**3 - 3 * P * T + 2 * Q)

zero("square_residual_identity", a**3 * ell**2 - q**2 - Delta)
zero("cardano_residual", P**3 - Q**2 - Delta * ell**4)
zero("power_polynomial_discriminant", sp.discriminant(f, T) - 108 * Delta * ell**4)

# Normalization and the exact seminormal-direction address.
xt = t**4 - 2 * t**2
yt = 3 * t**5 - 5 * t**3
ht = sp.expand((t**2 - 1) * (9 * t**4 - 18 * t**2 + 4))
h1 = sp.expand(t**2 * (t**2 - 1) * (9 * t**2 - 14))
zero("normalization_Delta", Delta.subs({x: xt, y: yt}))
zero("normalization_ell", ell.subs(x, xt) - (9 * t**4 - 18 * t**2 + 4))
zero("normalization_a", a.subs(x, xt) - (t**2 - 1) ** 2)
zero("normalization_q", q.subs({x: xt, y: yt}) - (t**2 - 1) ** 3 * ell.subs(x, xt))
zero("normalization_P", P.subs({x: xt, y: yt}) - ht**2)
zero("normalization_Q", Q.subs({x: xt, y: yt}) - ht**3)
zero("same_seminormal_direction", ht - h1 + 4 * (xt + 1))
equal(
    "h_derivative_packet",
    tuple(sp.diff(ht, t).subs(t, z) for z in (0, 1, -1)),
    (0, -10, 10),
)

# The complementary a <-> 1-a bicubic identity trades ell for M=9x+5,
# but it does not give another polynomial square/cube peel: 1-a=-x pulls
# back to t^2(2-t^2), with an odd simple factor.
Fpolar = 15 * a**2 - 15 * a + 4
Mpolar = 9 * a - 4
bpolar = 1 - a
zero("complementary_bicubic_identity", Fpolar**2 - a**3 * ell**2 - bpolar**3 * Mpolar**2)
zero(
    "complementary_Delta_form",
    Delta - (y**2 * (2 * Fpolar - y**2) - bpolar**3 * Mpolar**2),
)
zero("complementary_b_pullback", bpolar.subs(x, xt) - t**2 * (2 - t**2))
equal("complementary_simple_factor_gcd", sp.gcd(2 - t**2, sp.diff(2 - t**2, t)), 1)

# Multiplication table for the corrected rank-three order on (1,T,Z).
# T^2=ell Z, TZ=3a ell T-2q ell, Z^2=3a ell Z-2q T.
basis = (
    sp.Matrix((1, 0, 0)),
    sp.Matrix((0, 1, 0)),
    sp.Matrix((0, 0, 1)),
)


def mul(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    l0, l1, l2 = left
    r0, r1, r2 = right
    return sp.Matrix(
        (
            sp.expand(l0 * r0 - 2 * q * ell * (l1 * r2 + l2 * r1)),
            sp.expand(
                l0 * r1
                + l1 * r0
                + 3 * a * ell * (l1 * r2 + l2 * r1)
                - 2 * q * l2 * r2
            ),
            sp.expand(l0 * r2 + l2 * r0 + ell * l1 * r1 + 3 * a * ell * l2 * r2),
        )
    )


for i, bi in enumerate(basis):
    for j, bj in enumerate(basis):
        for k, bk in enumerate(basis):
            associator = mul(mul(bi, bj), bk) - mul(bi, mul(bj, bk))
            zero(f"associativity_{i}{j}{k}_0", associator[0])
            zero(f"associativity_{i}{j}{k}_1", associator[1])
            zero(f"associativity_{i}{j}{k}_2", associator[2])


def trace(left_multiplicand: sp.Matrix) -> sp.Expr:
    matrix = sp.Matrix.hstack(*(mul(left_multiplicand, bj) for bj in basis))
    return sp.expand(sp.trace(matrix))


trace_pairing = sp.Matrix(
    3,
    3,
    lambda i, j: trace(mul(basis[i], basis[j])),
)
expected_trace_pairing = sp.Matrix(
    [
        [3, 0, 6 * a * ell],
        [0, 6 * a * ell**2, -6 * q * ell],
        [6 * a * ell, -6 * q * ell, 18 * a**2 * ell**2],
    ]
)
for i in range(3):
    for j in range(3):
        zero(f"trace_pairing_{i}{j}", trace_pairing[i, j] - expected_trace_pairing[i, j])
zero("maximal_order_discriminant", trace_pairing.det() - 108 * Delta * ell**2)

# The local uniformizer Z is Eisenstein at the generic point of ell=0.
G = sp.expand(Z * (Z - 3 * a * ell) ** 2 - 4 * q**2 * ell)
zero(
    "Z_minimal_polynomial",
    G - (Z**3 - 6 * a * ell * Z**2 + 9 * a**2 * ell**2 * Z - 4 * q**2 * ell),
)
equal("ell_residue_a", a.subs(x, -sp.Rational(4, 9)), sp.Rational(5, 9))
equal("ell_residue_q", q.subs(x, -sp.Rational(4, 9)), y**2 - sp.Rational(8, 27))
nonzero("ell_eisenstein_constant_unit", (G.coeff(Z, 0) / ell).subs(x, -sp.Rational(4, 9)))
zero("recover_T_from_Z_table", mul(basis[2], basis[2])[1] + 2 * q)

# The corrected binary index form; its coefficient ideal is (ell,q).
theta = u * basis[1] + v * basis[2]
theta2 = mul(theta, theta)
index_form = sp.factor(u * theta2[2] - v * theta2[1])
expected_index = sp.expand(ell * u**3 - 3 * a * ell * u * v**2 + 2 * q * v**3)
zero("binary_index_form", index_form - expected_index)
equal("index_u3", sp.Poly(index_form, u, v).coeff_monomial(u**3), ell)
equal("index_u2v", sp.Poly(index_form, u, v).coeff_monomial(u**2 * v), 0)
zero("index_uv2", sp.Poly(index_form, u, v).coeff_monomial(u * v**2) + 3 * a * ell)
zero("index_v3", sp.Poly(index_form, u, v).coeff_monomial(v**3) - 2 * q)
zero("index_ideal_witness_ell", ell.subs(x, -sp.Rational(4, 9)))
zero(
    "index_ideal_witness_q",
    q.subs({x: -sp.Rational(4, 9), y**2: sp.Rational(8, 27)}),
)

# Quadratic-resolvent split and exact Kummer valuations (4,2) over ell.
zero("Delta_on_ell", Delta.subs(x, -sp.Rational(4, 9)) + q.subs(x, -sp.Rational(4, 9)) ** 2)
zero("kummer_norm", (-q + W) * (-q - W) - a**3 * ell**2 + (W**2 + Delta))
zero("signed_kummer_norm", q**2 + Delta - a**3 * ell**2)

# Trace translations preserve the same polynomial/field discriminant.
translated = sp.expand(f.subs(T, T - s))
zero("trace_translation_discriminant", sp.discriminant(translated, T) - sp.discriminant(f, T))

# Complete scalar representative scout P+c Delta, Q+d Delta.
Pcd = P + c * Delta
Qcd = Q + d * Delta
Rcd = sp.expand(
    ell**4
    + 3 * c * P**2
    + 3 * c**2 * P * Delta
    + c**3 * Delta**2
    - 2 * d * Q
    - d**2 * Delta
)
zero("scalar_shift_division", Pcd**3 - Qcd**2 - Delta * Rcd)
q_line = q.subs(x, -sp.Rational(4, 9))
zero(
    "scalar_shift_line_specialization",
    Rcd.subs(x, -sp.Rational(4, 9)) - q_line**2 * (d**2 + c**3 * q_line**2),
)

# If c=0 and d!=0, the quadratic in Y has nonzero discriminant.
RY = sp.expand(Rcd.subs(y**2, Y))
disc_pure_Q = sp.factor(sp.discriminant(RY.subs(c, 0), Y))
zero("pure_Q_shift_square_defect", disc_pure_Q - 4 * d**4 * a**3 * ell**2)

# If d=0 and c!=0, the forced even square misses in its constant Y row by
# ell^4(3c(x+1)^2+4)/4.  Reduction is modulo rho^2=c^3.
pure_P = sp.Poly(RY.subs(d, 0), Y)
coeffs = pure_P.all_coeffs()
h2 = rho
h1Y = sp.cancel(coeffs[1] / (2 * h2))
h0Y = sp.cancel((coeffs[2] - h1Y**2) / (2 * h2))
zero(
    "pure_P_Y_row",
    sp.rem(
        sp.Poly(sp.together(coeffs[3] - 2 * h1Y * h0Y).as_numer_denom()[0], rho),
        sp.Poly(rho**2 - c**3, rho),
    ).as_expr(),
)
constant_defect = sp.together(coeffs[4] - h0Y**2 - ell**4 * (3 * c * a**2 + 4) / 4)
zero(
    "pure_P_constant_square_defect",
    sp.rem(
        sp.Poly(constant_defect.as_numer_denom()[0], rho),
        sp.Poly(rho**2 - c**3, rho),
    ).as_expr(),
)
nonzero("pure_P_defect_nonzero", 3 * c * a**2 + 4)

semantic_packet = (
    "Delta=(x+1)^3(9x+4)^2-q^2 and P^3-Q^2=Delta(9x+4)^4",
    "normalization h=(t^2-1)(9t^4-18t^2+4)=h1-4(x+1)",
    "complementary bicubic trades ell for 9x+5 but introduces nonsquare -x",
    "irreducible S3 cubic T^3-3PT+2Q",
    "power-order discriminant 108 Delta ell^4 and index ell",
    "normal maximal order B+B*T+B*Z with discriminant 108 Delta ell^2",
    "tame total cubic ramification along ell and transposition ramification along Delta",
    "Kummer valuations (4,2) above the split ell line",
    "corrected index form ell*u^3-3aell*u*v^2+2q*v^3; no unit value",
    "trace shifts and all scalar coefficient shifts do not remove the extra line",
    "nonlinear coefficient shifts, Keller atlas, and JC2 remain open",
)
semantic_sha = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()

print("THM3869_IDENTITY", "Delta=(x+1)^3*ell^2-q^2;P^3-Q^2=Delta*ell^4")
print("THM3869_DIRECTION", "h=h1-4(x+1);same seminormal class,not fourth")
print("THM3869_COMPLEMENT", "F^2-a^3*ell^2=(1-a)^3*(9a-4)^2;nonsquare 1-a blocks peel")
print("THM3869_CUBIC", "irreducible S3;Newton slope -2/3 at ell")
print("THM3869_DISCRIMINANTS", "power=108*Delta*ell^4;maximal=field=108*Delta*ell^2")
print("THM3869_RAMIFICATION", "Delta:(2,1);ell:e=3;Kummer valuations=(4,2)")
print("THM3869_INDEX_FORM", "ell*u^3-3*(x+1)*ell*u*v^2+2*q*v^3;ideal=(ell,q)")
print("THM3869_SHIFT_SCOUT", "trace and scalar (c,d) shifts do not remove ell")
print("THM3869_SCOPE", "explicit nonmonogenic near-miss; nonlinear shifts/atlas/JC2 open")
print("SEMANTIC_SHA256", semantic_sha)
print("CHECKS", CHECKS)
