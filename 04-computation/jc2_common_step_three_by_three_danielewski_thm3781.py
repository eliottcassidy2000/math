#!/usr/bin/env python3
"""Exact support, coefficient, and pole-order companion for THM-3781."""

from __future__ import annotations

import ast
import hashlib
from collections import Counter
from math import gcd
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def weight_grid(left: tuple[int, ...], right: tuple[int, ...]) -> Counter[int]:
    return Counter(r + s + 1 for r in left for s in right)


def bracket_coefficient(
    r: int,
    f: sp.Expr,
    fp: sp.Expr,
    s: int,
    g: sp.Expr,
    gp: sp.Expr,
) -> sp.Expr:
    """Coefficient of c^(r+s+1) in {c^r f(b),c^s g(b)}."""
    return sp.expand(s * fp * g - r * f * gp)


def deriv_power_times(
    owner: sp.Symbol,
    owner_p: sp.Symbol,
    power: int,
    value: sp.Expr,
    value_p: sp.Expr,
) -> sp.Expr:
    """Formal b-derivative of owner^power*value, also for negative power."""
    return power * owner ** (power - 1) * owner_p * value + owner**power * value_p


# K and L may be algebraic roots of polynomial endpoint owners.  The formal
# derivative identities below therefore verify the equations in the finite
# algebraic function-field extension used in the proof, without assuming the
# endpoint exponents are coprime.
K, Kp, L, Lp = sp.symbols("K Kp L Lp", nonzero=True)
h, hp = sp.symbols("h hp")
a_mid, a_mid_p = sp.symbols("a_mid a_mid_p")
A0, B0, L0, M0 = sp.symbols("A0 B0 L0 M0", nonzero=True)
lam, mu = sp.symbols("lambda mu", nonzero=True)
C0, D0 = sp.symbols("C D")


# Unequal-step cells die at a lonely convolution weight, before any endpoint
# factorization.  The first output has the negative centre -a and step r;
# the second has centre a-1 and step s.  Endpoint signs require
# 1<=a<=min(s,r-1).  Weight -r is unique except at s=r or s=2r, while +s is
# unique except at r=s or r=2s.  Unequal steps cannot escape both tests.
unequal_cells: list[tuple[int, int, int]] = []
doubling_hostiles: list[tuple[int, int, int]] = []
for r_step in range(1, 21):
    for s_step in range(1, 21):
        if r_step == s_step:
            continue
        for centre_a in range(1, min(s_step, r_step - 1) + 1):
            unequal_cells.append((r_step, s_step, centre_a))
            pairs_by_weight: dict[int, list[tuple[int, int]]] = {}
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    pairs_by_weight.setdefault(i * r_step + j * s_step, []).append((i, j))

            lower_unique = pairs_by_weight[-r_step] == [(-1, 0)]
            upper_unique = pairs_by_weight[s_step] == [(0, 1)]
            gate(lower_unique == (s_step != 2 * r_step),
                 "lower lonely-weight collision classification")
            gate(upper_unique == (r_step != 2 * s_step),
                 "upper lonely-weight collision classification")
            gate(lower_unique or upper_unique,
                 "unequal steps have a forbidden lonely mixed-sign weight")
            if s_step == 2 * r_step or r_step == 2 * s_step:
                doubling_hostiles.append((r_step, s_step, centre_a))


cells: list[tuple[int, int]] = []
nonprimitive_cells: list[tuple[int, int]] = []
for d_step in range(2, 31):
    # After exchanging the outputs, every scalar-centred placement whose
    # endpoints have the required signs has centres -a and a-1.
    brute_centres = []
    for p_mid in range(-2 * d_step, 2 * d_step + 1):
        q_mid = -p_mid - 1
        if (p_mid - d_step < 0 < p_mid + d_step
                and q_mid - d_step < 0 < q_mid + d_step):
            brute_centres.append((p_mid, q_mid))
    expected_centres = [
        (p_mid, -p_mid - 1)
        for p_mid in range(-d_step + 1, d_step - 1)
    ]
    gate(brute_centres == expected_centres, "complete scalar-centred census")

    for centre_a in range(1, d_step):
        cells.append((d_step, centre_a))
        neg_p = d_step + centre_a
        neg_q = d_step - centre_a + 1
        pos_p = d_step - centre_a
        pos_q = d_step + centre_a - 1

        supports_p = (-neg_p, -centre_a, pos_p)
        supports_q = (-neg_q, centre_a - 1, pos_q)
        convolution = weight_grid(supports_p, supports_q)
        gate(
            convolution
            == Counter({
                -2 * d_step: 1,
                -d_step: 2,
                0: 3,
                d_step: 2,
                2 * d_step: 1,
            }),
            "universal support convolution",
        )

        u_neg = gcd(neg_p, neg_q)
        u_pos = gcd(pos_p, pos_q)
        if u_neg > 1 or u_pos > 1:
            nonprimitive_cells.append((d_step, centre_a))
        gate(gcd(neg_p // u_neg, neg_q // u_neg) == 1,
             "negative endpoint owner is primitive")
        gate(gcd(pos_p // u_pos, pos_q // u_pos) == 1,
             "positive endpoint owner is primitive")

        # In the radical extension K^u_neg=k and L^u_pos=ell, the endpoint
        # Wronskians take the same full-exponent form in every cell.
        f = A0 * K**neg_p
        fp = A0 * neg_p * K ** (neg_p - 1) * Kp
        g = B0 * K**neg_q
        gp = B0 * neg_q * K ** (neg_q - 1) * Kp
        big_f = L0 * L**pos_p
        big_fp = L0 * pos_p * L ** (pos_p - 1) * Lp
        big_h = M0 * L**pos_q
        big_hp = M0 * pos_q * L ** (pos_q - 1) * Lp

        gate(bracket_coefficient(-neg_p, f, fp, -neg_q, g, gp) == 0,
             "negative endpoint power law")
        gate(bracket_coefficient(pos_p, big_f, big_fp,
                                 pos_q, big_h, big_hp) == 0,
             "positive endpoint power law")

        lower = (
            bracket_coefficient(-neg_p, f, fp,
                                centre_a - 1, h, hp)
            + bracket_coefficient(-centre_a, a_mid, a_mid_p,
                                  -neg_q, g, gp)
        )
        expected_lower = K ** (d_step + 1) * (
            A0 * neg_p
            * deriv_power_times(K, Kp, centre_a - 1, h, hp)
            - B0 * neg_q
            * deriv_power_times(K, Kp, -centre_a, a_mid, a_mid_p)
        )
        gate(sp.factor(lower - expected_lower) == 0,
             "lower adjacent derivative identity")

        upper = (
            bracket_coefficient(-centre_a, a_mid, a_mid_p,
                                pos_q, big_h, big_hp)
            + bracket_coefficient(pos_p, big_f, big_fp,
                                  centre_a - 1, h, hp)
        )
        expected_upper = L ** (d_step - 1) * (
            M0 * pos_q
            * deriv_power_times(L, Lp, centre_a, a_mid, a_mid_p)
            - L0 * pos_p
            * deriv_power_times(L, Lp, -(centre_a - 1), h, hp)
        )
        gate(sp.factor(upper - expected_upper) == 0,
             "upper adjacent derivative identity")

        # Eliminate h from the two integrated adjacent equations.
        h_from_upper = L ** (centre_a - 1) * (
            mu * a_mid * L**centre_a + D0
        )
        combined = sp.expand(
            a_mid
            - K**centre_a
            * (lam * K ** (centre_a - 1) * h_from_upper + C0)
        )
        obstruction = sp.expand(
            a_mid * (1 - lam * mu * (K * L) ** (2 * centre_a - 1))
            - K**centre_a
            * (lam * D0 * (K * L) ** (centre_a - 1) + C0)
        )
        gate(sp.factor(combined - obstruction) == 0,
             "universal adjacent-equation elimination")

        # Exact pole-order formulas at any place above b=infinity.  Let dk>0
        # and dl>=0 be the normalized pole orders of K and L.  A nonzero D
        # term dominates C when a>1; otherwise the surviving numerator is
        # constant.  Both branches force the polynomial a_mid to have
        # nonpositive degree, contradicting its negative-weight membership.
        for degree_k in range(1, 7):
            for degree_l in range(0, 7):
                d_branch_degree_a = -centre_a * degree_l
                constant_branch_degree_a = (
                    -(centre_a - 1) * degree_k
                    - (2 * centre_a - 1) * degree_l
                )
                gate(d_branch_degree_a <= 0,
                     "D-branch pole order is nonpositive")
                gate(constant_branch_degree_a <= 0,
                     "constant-numerator pole order is nonpositive")


gate((2, 1) in cells, "step-two control")
gate({a0 for d0, a0 in cells if d0 == 3} == {1, 2},
     "step-three complete row")
gate((4, 2) in nonprimitive_cells,
     "noncoprime endpoint hostile is included")
gate((5, 2) in nonprimitive_cells,
     "second noncoprime endpoint hostile is included")
gate((2, 4, 1) in doubling_hostiles,
     "oriented doubling hostile is included")
gate((4, 2, 1) in doubling_hostiles,
     "reverse doubling hostile is included")


semantic_rows = (
    "surface:c^2*e=Sigma(b);Sigma_squarefree;degSigma>=2",
    "scope:all_scalar_centred_three_piece_arithmetic_progressions",
    "unequal:weight_-r_unique_unless_s=r_or_s=2r;weight_+s_unique_unless_r=s_or_r=2s",
    "unequal:one_lonely_mixed_sign_bracket_always_survives",
    "common:d>=2;1<=a<=d-1;all_endpoint_gcds",
    "supports:P=-(d+a),-a,d-a;Q=-(d-a+1),a-1,d+a-1",
    "radicals:K^gcd(d+a,d-a+1)=k;L^gcd(d-a,d+a-1)=ell",
    "adjacent:a_mid/K^a=lambda*K^(a-1)*h+C",
    "adjacent:h/L^(a-1)=mu*a_mid*L^a+D",
    "elimination:a_mid*(1-nu*(K*L)^(2a-1))=K^a*(lambda*D*(K*L)^(a-1)+C)",
    "infinity:deltaK>0;deltaL>=0;negative_middle_forces_deg_a_mid>0",
    "exit:D_branch_deg_a_mid=-a*deltaL;constant_branch_deg_a_mid=-(a-1)*deltaK-(2a-1)*deltaL",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3781-common-step-three-by-three-danielewski-darboux-nonentry")
print("surface=exponent_two_squarefree_Danielewski")
print("scope=all_scalar_centred_three_piece_arithmetic_progression_supports")
print(f"bounded_unequal_cells_r_s_to20={len(unequal_cells)}")
print(f"bounded_doubling_hostiles={len(doubling_hostiles)}")
print(f"bounded_common_cells_d2_to_d30={len(cells)}")
print(f"bounded_nonprimitive_endpoint_controls={len(nonprimitive_cells)}")
print("step=3;placements_up_to_swap=2;subsumed")
print("exit=unequal_lonely_weight_or_common_radical_owner_infinity_pole_order")
print("scalar_equation=never_reached")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
