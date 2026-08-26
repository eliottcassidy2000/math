#!/usr/bin/env python3
"""Exact symbolic audit for the uniform A5 arbitrary-context bound.

This script uses only the exact THM-4208 response jet and ordinal-transfer
identities.  It derives the eleven coordinates of

    D(B) = lambda((C3 ▷ P5) ▷ B) - lambda(B),

checks the right-context lower-bound decomposition, and checks the unary
terminal estimate.  All arithmetic is exact SymPy integer/rational algebra.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    need(sp.expand(left - right) == 0, message)


def left_jet(
    h: sp.Expr,
    w: sp.Expr,
    s0: sp.Expr,
    s1: sp.Expr,
    q00: sp.Expr,
    q01: sp.Expr,
    q11: sp.Expr,
    ell0: sp.Expr,
    ell1: sp.Expr,
) -> tuple[sp.Expr, ...]:
    return (
        h * w,
        2 * ell0,
        2 * ell1,
        2 * h * s0,
        2 * h * s1,
        2 * s0**2 + 6 * q00,
        4 * s0 * s1 + 12 * q01,
        2 * s1**2 + 6 * q11,
        2 * q00 - 10 * s0**2,
        4 * q01 - 20 * s0 * s1,
        2 * q11 - 10 * s1**2,
    )


def main() -> None:
    h, w = sp.symbols("h w")
    ell0, ell1 = sp.symbols("ell0 ell1")
    q00, q01, q11 = sp.symbols("q00 q01 q11")
    n00, n01, n10, n11 = sp.symbols("n00 n01 n10 n11")

    # B totals.  Here q_ef=sum_i U_i^e U_i^f and
    # n_ef=sum_i U_i^e V_i^f.
    b = w + 2 * h
    s0 = w / 2
    s1 = w / 2 + h

    # Exact A5=C3 ▷ P5 data from THM-4208, with 2^5=32.
    hp, wp = sp.Integer(3), sp.Integer(378)
    ps0, ps1 = sp.Integer(189), sp.Integer(192)
    pq00, pq01, pq11 = sp.Integer(7677), sp.Integer(7677), sp.Integer(7686)
    pell0, pell1 = sp.Integer(116064), sp.Integer(116622)

    # Exact ordinal transfer T=A5 ▷ B.
    ht = 3 * h
    wt = 384 * w + 762 * h
    ts0 = wt / 2
    ts1 = wt / 2 + ht

    # Old A5 start states become balanced after convolution with E(B),
    # while the B start states scale by H(A5)=3.
    prefix_uplus_square = pq00 + 2 * pq01 + pq11
    common_q = b**2 * prefix_uplus_square / 4
    tq00 = 9 * q00 + common_q
    tq01 = 9 * q01 + common_q
    tq11 = 9 * q11 + common_q

    # For an A5 vertex i, p_i=2<U_i(A5),s(B)>.  For a B vertex j,
    # q_j=378 V_j^0+384 V_j^1.  Expanding the two blockwise linear
    # moments gives the following exact common and B-dependent parts.
    prefix_uplus_linear = pell0 + pell1
    prefix_weighted_cross = 2 * ((pq00 + pq01) * s0 + (pq01 + pq11) * s1)
    common_ell = (
        b * h * prefix_uplus_linear / 2
        + b * 384 * (w + h) * (ps0 + ps1) / 2
        + 3 * b * prefix_weighted_cross / 2
    )
    tell0 = (
        common_ell
        + 9 * ell0
        + 3 * (381 * w + 762 * h) * s0
        - 3 * (378 * n00 + 384 * n01)
    )
    tell1 = (
        common_ell
        + 9 * ell1
        + 3 * (381 * w + 762 * h) * s1
        - 3 * (378 * n10 + 384 * n11)
    )

    jet_b = left_jet(h, w, s0, s1, q00, q01, q11, ell0, ell1)
    jet_t = left_jet(ht, wt, ts0, ts1, tq00, tq01, tq11, tell0, tell1)
    derived = tuple(sp.factor(x - y) for x, y in zip(jet_t, jet_b))

    expected = (
        h * (2286 * h + 1151 * w),
        2
        * (
            471168 * h**2
            + 475182 * h * w
            + 8 * ell0
            - 1134 * n00
            - 1152 * n01
            + 119799 * w**2
        ),
        2
        * (
            473454 * h**2
            + 476325 * h * w
            + 8 * ell1
            - 1134 * n10
            - 1152 * n11
            + 119799 * w**2
        ),
        h * (2286 * h + 1151 * w),
        1151 * h * (2 * h + w),
        474624 * h**2 + 476910 * h * w + 48 * q00 + 119803 * w**2,
        2 * (476910 * h**2 + 478061 * h * w + 48 * q01 + 119803 * w**2),
        479212 * h**2 + 479212 * h * w + 48 * q11 + 119803 * w**2,
        -1390176 * h**2 - 1401606 * h * w + 16 * q00 - 353279 * w**2,
        2 * (-1401606 * h**2 - 1407361 * h * w + 16 * q01 - 353279 * w**2),
        -1413116 * h**2 - 1413116 * h * w + 16 * q11 - 353279 * w**2,
    )
    for index, (actual, target) in enumerate(zip(derived, expected)):
        same(actual, target, f"D coordinate {index}")

    d = expected

    # For C write k=H(C), z=W(C), v_i=V_i^0+V_i^1, and use THM-4208's
    # coordinatewise identity t_i=Start_i=V_i^1-V_i^0.  The aggregate symbols below are
    # m=sum v_i^2, p=sum v_i t_i, and tau=sum t_i^2.
    k, z = sp.symbols("k z")
    m, p, tau = sp.symbols("m p tau")
    cell0, cell1 = sp.symbols("cell0 cell1")
    cm00 = (m - 2 * p + tau) / 4
    cm01 = (m - tau) / 4
    cm11 = (m + 2 * p + tau) / 4

    # mu(C)-k^2 mu(P1).  The ell coordinates are left literal for now.
    mu_difference = (
        k * z,
        k * z / 2,
        k * z / 2,
        cell0,
        cell1,
        z**2 / 4,
        z**2 / 4 + z * k / 2,
        z**2 / 4 + z * k,
        cm00,
        cm01,
        cm11 - k**2,
    )
    context_slack = sp.expand(sum(x * y for x, y in zip(d, mu_difference)))

    # Canonical coefficients after changing from (V0,V1) to (v,t).
    coefficient_kz = sp.factor(d[0] + (d[1] + d[2]) / 2)
    coefficient_z2 = sp.factor((d[5] + d[6] + d[7]) / 4)
    coefficient_zk = sp.factor(d[6] / 2 + d[7])
    coefficient_m = sp.factor((d[8] + d[9] + d[10]) / 4)
    coefficient_p = sp.factor((d[10] - d[8]) / 2)
    coefficient_tau = sp.factor((d[8] - d[9] + d[10]) / 4)
    coefficient_k2 = sp.factor(-d[10])
    weighted_ell = d[3] * cell0 + d[4] * cell1

    canonical_slack = (
        coefficient_kz * z * k
        + weighted_ell
        + coefficient_z2 * z**2
        + coefficient_zk * z * k
        + coefficient_m * m
        + coefficient_p * p
        + coefficient_tau * tau
        + coefficient_k2 * k**2
    )
    same(context_slack, canonical_slack, "canonical right-context slack")

    # Lower coefficients obtained only from nonnegative-coordinate bounds on B.
    min_kz = 945756 * h**2 + 950363 * h * w + 238455 * w**2
    min_z2 = 476914 * h**2 + 478061 * h * w + 119803 * w**2
    min_zk = 956122 * h**2 + 957273 * h * w + 239606 * w**2
    min_m = -(1401626 * h**2 + 1407361 * h * w + 353279 * w**2)
    min_p = -(11470 * h**2 + 5755 * h * w)
    min_tau = -20 * h**2
    min_k2 = 1413100 * h**2 + 1413100 * h * w + 353275 * w**2
    ell_weight_v = h * (2294 * h + 1151 * w)

    # Each difference below has the displayed nonnegative-debt form.  On B,
    # U_i^1-U_i^0=End_i coordinatewise.  In particular
    # q11-q00=sum_i End_i(U_i^0+U_i^1)>=0, and
    # q00-2q01+q11=sum_i End_i^2>=0.
    a = w + h
    deficit_n0 = a * w / 2 - (n00 + n10)
    deficit_n1 = a * b / 2 - (n01 + n11)
    same(
        coefficient_kz - min_kz,
        8 * (ell0 + ell1) + 1134 * deficit_n0 + 1152 * deficit_n1,
        "kz coefficient debt",
    )
    start_mass_square = q00 + 2 * q01 + q11
    same(coefficient_z2 - min_z2, 12 * start_mass_square, "z2 coefficient debt")
    same(coefficient_zk - min_zk, 48 * (q01 + q11), "zk coefficient debt")
    same(coefficient_m - min_m, 4 * start_mass_square, "m coefficient debt")
    same(
        coefficient_p - min_p,
        8 * (q11 - q00),
        "p coefficient debt",
    )
    endpoint_end_square = q00 - 2 * q01 + q11
    same(coefficient_tau - min_tau, 4 * endpoint_end_square, "tau coefficient debt")
    same(coefficient_k2 - min_k2, 4 * (b**2 - 4 * q11), "k2 coefficient debt")
    same(
        weighted_ell,
        ell_weight_v * (cell0 + cell1) + 8 * h**2 * (cell1 - cell0),
        "weighted ell endpoint change",
    )

    # Apply the C-side inequalities:
    #   weighted_ell >= -4z[ell_weight_v(z+k)+8h^2 k],
    #   m <= (z^2+4zk+3k^2)/3, p <= (z+k)k, tau <= k^2.
    lower_bound = sp.expand(
        min_kz * z * k
        - 4 * z * (ell_weight_v * (z + k) + 8 * h**2 * k)
        + min_z2 * z**2
        + min_zk * z * k
        + min_m * (z**2 + 4 * z * k + 3 * k**2) / 3
        + min_p * (z + k) * k
        + min_tau * k**2
        + min_k2 * k**2
    )
    advertised_lower_bound = (
        2 * (794 * h**2 + 6505 * h * w + 3065 * w**2) * z**2 / 3
        + (37096 * h**2 + 62387 * h * w + 21067 * w**2) * z * k / 3
        - 4 * (2 * h + w) ** 2 * k**2
    )
    same(lower_bound, advertised_lower_bound, "advertised right-context lower bound")

    # For a non-singleton C, W(C)>=H(C), so z>=k.  Both positive-z
    # coefficients are strictly positive and the minimum occurs at z=k.
    terminal_lower_bound = (
        k**2 * (38636 * h**2 + 75349 * h * w + 27185 * w**2) / 3
    )
    same(advertised_lower_bound.subs(z, k), terminal_lower_bound, "z=k boundary")

    # The singleton right jet has nonzero coordinates 2,7,10 only.
    unary = sp.expand(d[2] + d[7] + d[10])
    unary_identity = (
        13004 * h**2
        + 18746 * h * w
        + 6122 * w**2
        + 16 * ell1
        + 64 * q11
        - 2268 * n10
        - 2304 * n11
    )
    same(unary, unary_identity, "unary exact identity")

    unary_slack = sp.expand(unary - 10764 * h**2)
    unary_floor = -64 * h**2 + 15308 * h * w + 4979 * w**2 + 64 * q11
    weighted_mixed_deficit = b * (1143 * w + 1152 * h) - 2268 * n10 - 2304 * n11
    same(
        unary_slack - unary_floor,
        16 * ell1 + weighted_mixed_deficit,
        "unary nonnegative debt",
    )
    nontrivial_unary_floor = 20223 * h**2
    same(unary_floor.subs({w: h, q11: 0}), nontrivial_unary_floor, "w=h unary boundary")

    lines = [
        "theorem=THM-4212-candidate",
        "symbolic_D_coordinates=11",
        "D_coordinate_audit=PASS",
        "right_context_canonicalization=PASS",
        "right_context_debt_decomposition=PASS",
        "right_context_LB_numerator=2*(794h^2+6505hw+3065w^2)z^2+(37096h^2+62387hw+21067w^2)zk-12*(2h+w)^2k^2",
        "right_context_nontrivial_terminal_numerator=38636h^2+75349hw+27185w^2",
        "unary_identity_audit=PASS",
        "unary_floor=-64h^2+15308hw+4979w^2+64q11",
        "nontrivial_unary_floor_at_w=h_q11=0=20223h^2",
        "equality_boundary=B=C=P1",
    ]
    payload = "\n".join(lines) + "\n"
    print(payload, end="")
    print(f"semantic_sha256={hashlib.sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
