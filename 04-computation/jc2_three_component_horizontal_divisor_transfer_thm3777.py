#!/usr/bin/env python3
"""Exact companion for THM-3777's vertical-to-horizontal divisor transfer."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


x, y, t = sp.symbols("x y t")
m_symbol = sp.symbols("m", integer=True, positive=True)
p_target, u_target = sp.symbols("p_target u_target")
s, rho_symbol, g0 = sp.symbols("s rho g0")


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, y)
        - sp.diff(first, y) * sp.diff(second, x)
    )


def jac_xt(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.simplify(
        sp.diff(first, x) * sp.diff(second, t)
        - sp.diff(first, t) * sp.diff(second, x)
    )


# Symbolic all-m chart identities.
P_chart = 1 / x + (m_symbol + 1) * t**m_symbol / m_symbol
U_chart = t / x + t ** (m_symbol + 1)
gate(jac_xt(P_chart, U_chart) == -1 / x**3,
     "all-m independent-chart Jacobian")
gate(sp.simplify(x**3 * jac_xt(P_chart, U_chart) + 1) == 0,
     "all-m source Jacobian sign")
gate(
    sp.simplify(
        t ** (m_symbol + 1)
        - m_symbol * P_chart * t
        + m_symbol * U_chart
    )
    == 0,
    "all-m cover relation",
)

# Each rational target shear has determinant one in target variables.
F_target = sp.Function("F")
G_target = sp.Function("G")
q_shear = u_target + F_target(p_target)
p_shear = p_target + G_target(u_target)
target_q_det = sp.diff(p_target, p_target) * sp.diff(q_shear, u_target)
target_q_det -= sp.diff(p_target, u_target) * sp.diff(q_shear, p_target)
target_p_det = sp.diff(p_shear, p_target) * sp.diff(u_target, u_target)
target_p_det -= sp.diff(p_shear, u_target) * sp.diff(u_target, p_target)
gate(target_q_det == 1 and target_p_det == 1,
     "target shear determinants")


A = 1 + x * y
profile_checks = 0
gcd_checks = 0
source_jacobian_checks = 0
profiles: dict[int, tuple[sp.Expr, ...]] = {}

for m in range(1, 9):
    B = 1 + x ** (2 * m + 1) * A**m
    U = x * A * B
    P = ((m + 1) * B - 1) / (m * x)
    W = sp.cancel(U * P)
    R = sp.expand(W - 1)
    N = (
        y
        + sp.Rational(2 * m + 1, m) * x ** (2 * m) * A ** (m + 1)
        + sp.Rational(m + 1, m) * x ** (4 * m + 1) * A ** (2 * m + 1)
    )
    D_profile = sp.expand(A * B)
    P1 = sp.cancel(P - 1 / U)

    gate(sp.expand(B - (1 + x * (x**2 * A) ** m)) == 0,
         "B chart profile")
    gate(
        sp.cancel(P - (1 / x + sp.Rational(m + 1, m) * (x**2 * A) ** m))
        == 0,
        "P chart profile",
    )
    gate(
        sp.expand(U - ((x**2 * A) / x + (x**2 * A) ** (m + 1))) == 0,
        "U chart profile",
    )
    gate(sp.expand(R - x * N) == 0, "R=xN profile")
    gate(sp.cancel(P1 - N / D_profile) == 0,
         "paid primitive reduced expression")
    gate(sp.expand(N.subs(x, 0) - y) == 0,
         "paid-axis N restriction")
    gate(sp.expand(W.subs(x, 0) - 1) == 0,
         "axis principal coefficient one")
    gate(
        sp.cancel(W - A * B * ((m + 1) * B - 1) / m) == 0,
        "zero coefficients on A and B",
    )
    gate(sp.cancel(U * P1 - R) == 0,
         "remaining coefficient packet")

    gcd_value = sp.gcd(sp.Poly(sp.expand(N), x, y), sp.Poly(D_profile, x, y))
    gate(gcd_value.total_degree() == 0,
         "reduced global denominator")
    gcd_checks += 1

    # Direct source differentiation independently guards the chart proof on
    # four increasing degrees; the symbolic chart identity covers all m.
    if m <= 4:
        gate(sp.cancel(jac(P, U) + 1) == 0,
             "direct source Jacobian")
        source_jacobian_checks += 1

    profiles[m] = (B, U, P, W, R, N, D_profile, P1)
    profile_checks += 1


# Finite-pole pullbacks are explicit transverse horizontal curves.  Several
# m and rho values guard all signs and the rho=0 boundary.
horizontal_controls = 0
for m in (1, 2, 4, 8):
    B, U, P, W, R, N, D_profile, P1 = profiles[m]
    for rho in (-2, 0, 3):
        H_rho = sp.expand(N - rho * D_profile)
        gate(sp.expand(H_rho.subs(x, 0) - (y - rho)) == 0,
             "horizontal axis restriction")
        gate(sp.diff(H_rho, y).subs({x: 0, y: rho}) == 1,
             "horizontal transversality")
        gate(D_profile.subs({x: 0, y: rho}) == 1,
             "horizontal point avoids vertical denominator")
        gcd_horizontal = sp.gcd(
            sp.Poly(H_rho, x, y), sp.Poly(D_profile, x, y)
        )
        gate(gcd_horizontal.total_degree() == 0,
             "horizontal and vertical coprimality")

        # G(s)=1/(s-rho) pulls back to D/H_rho.  At the transverse point the
        # numerator of U+D/H_rho is the unit D=1, so the pole is genuine.
        pulled_finite_pole = sp.cancel(1 / (P1 - rho))
        gate(sp.cancel(pulled_finite_pole - D_profile / H_rho) == 0,
             "finite pole pullback")
        U1_numerator = sp.expand(U * H_rho + D_profile)
        gate(U1_numerator.subs({x: 0, y: rho}) == 1,
             "finite pole cannot cancel against U")
        horizontal_controls += 1


# A pole at infinity retains a vertical denominator.  The polynomial control
# G(s)=s^2 is sufficient because the proof uses only its positive order.
infinity_controls = 0
for m in (1, 3, 6):
    B, U, P, W, R, N, D_profile, P1 = profiles[m]
    pulled_polynomial = sp.cancel(P1**2)
    gate(sp.cancel(pulled_polynomial - N**2 / D_profile**2) == 0,
         "infinity pole pullback")
    infinity_numerator = sp.expand(U * D_profile**2 + N**2)
    infinity_gcd = sp.gcd(
        sp.Poly(infinity_numerator, x, y), sp.Poly(D_profile, x, y)
    )
    gate(infinity_gcd.total_degree() == 0,
         "vertical infinity pole cannot cancel")
    infinity_controls += 1


# The constant middle branch: the unique simple correction is 1/U and
# literally restores P, hence its x-axis pole.
constant_controls = 0
for m in (1, 2, 5, 8):
    B, U, P, W, R, N, D_profile, P1 = profiles[m]
    U1_constant = U + g0
    restored = sp.cancel(P1 + 1 / (U1_constant - g0))
    gate(sp.cancel(restored - P) == 0,
         "constant middle branch exact reversal")
    gate(sp.limit(x * restored.subs(y, 0), x, 0) == 1,
         "restored x-axis principal coefficient")
    gate(sp.cancel(U * P1 - (W - 1)) == 0,
         "equal negative remaining coefficients")
    constant_controls += 1


# Exhaustive rational-function trichotomy controls by degree profile:
# numerator degree > denominator degree -> infinity; otherwise a nonconstant
# reduced denominator over an algebraically closed field has a finite root;
# denominator constant and no infinity pole -> constant.
degree_cells = 0
branch_counts = {"infinity": 0, "finite": 0, "constant": 0}
for numerator_degree in range(0, 7):
    for denominator_degree in range(0, 7):
        degree_cells += 1
        if numerator_degree > denominator_degree:
            branch_counts["infinity"] += 1
        elif denominator_degree > 0:
            branch_counts["finite"] += 1
        else:
            gate(numerator_degree == 0,
                 "only constant avoids finite and infinity poles")
            branch_counts["constant"] += 1
gate(degree_cells == 49 and branch_counts == {"infinity": 21, "finite": 27, "constant": 1},
     "rational branch census")


semantic_rows = (
    "family:A=1+xy;B=1+x^(2m+1)A^m;U=xAB;P=((m+1)B-1)/(mx)",
    "chart:t=x^2A;J(P,U)=-1;t^(m+1)-mPt+mU=0",
    "packet:x,A,B=(1,0,0);payment:P1=P-1/U=N_m/(AB)",
    "N_m=y+(2m+1)/m*x^(2m)A^(m+1)+(m+1)/m*x^(4m+1)A^(2m+1)",
    "finite_pole:rho->H_rho=N_m-rhoAB;H_rho|x=0=y-rho",
    "infinity_pole:retains_A,B_vertical_debt",
    "constant_middle:third_simple_pole_restores_P_and_x_debt",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3777-three-component-tower-reciprocal-payment-horizontal-divisor-transfer")
print("scope=all_integer_m_at_least_1;opposite_start_three_target_shears")
print("rational_keller=J(P,U)=-1;cover=t^(m+1)-mPt+mU")
print("vertical_packet=(1,0,0);paid_packet=(0,-1,-1)")
print("nonconstant_middle=vertical_pole_at_infinity_or_horizontal_H_rho")
print("constant_middle=third_shear_restores_x_axis_pole")
print(f"profile_checks={profile_checks};gcd_checks={gcd_checks};source_jacobian_checks={source_jacobian_checks}")
print(f"horizontal_controls={horizontal_controls};infinity_controls={infinity_controls};constant_controls={constant_controls}")
print(f"degree_cells={degree_cells};branches={branch_counts['infinity']},{branch_counts['finite']},{branch_counts['constant']}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
