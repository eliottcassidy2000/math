#!/usr/bin/env python3
"""Exact primary audit for THM-4215.

The theorem concerns

    Z_n = C3 ▷ C3 ▷ P_n,
    F_n(B,C) = R_+(Z_n ▷ B,C) - R_+(B,C).

This audit has two deliberately visible layers.  The inherited exact ordinal
engine reconstructs the fixed prefix Z_7 and the singleton obstruction rows.
SymPy then starts only from the finite response jet and exact ordinal-transfer
identities, derives all eleven coordinates of

    lambda(Z_7 ▷ B) - lambda(B),

and checks the coefficient-debt proof of the sharp uniform floor 967788.
"""

from __future__ import annotations

import hashlib
import itertools

import sympy as sp

import tournament_ordinal_cocycle_parity_thm4184 as base


FLOOR = 967_788


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    need(sp.expand(left - right) == 0, message)


def ordinal(*factors: base.TournamentData) -> base.TournamentData:
    need(bool(factors), "ordinal product needs a factor")
    value = factors[0]
    for factor in factors[1:]:
        value = base.ordinal_data(value, factor)
    return value


def response(
    prefix: base.TournamentData,
    middle: base.TournamentData,
    right: base.TournamentData,
) -> int:
    return base.remainder(base.ordinal_data(prefix, middle), right) - base.remainder(
        middle, right
    )


def directed_masses(
    data: base.TournamentData,
) -> tuple[list[int], list[int], list[int]]:
    degree = [sum(row) for row in data.capacities]
    outgoing = [
        sum(
            data.capacities[i][j]
            for j in range(len(data.out))
            if data.out[i] & (1 << j)
        )
        for i in range(len(data.out))
    ]
    incoming = [degree[i] - outgoing[i] for i in range(len(data.out))]
    return degree, outgoing, incoming


def prefix_coordinates(
    data: base.TournamentData,
) -> tuple[int, int, int, int, int, int, int, int, int]:
    degree, outgoing, _ = directed_masses(data)
    linear = [
        data.mass - degree[i] + 4 * outgoing[i] for i in range(len(data.out))
    ]
    q00 = sum(row[0] * row[0] for row in data.starts)
    q01 = sum(row[0] * row[1] for row in data.starts)
    q11 = sum(row[1] * row[1] for row in data.starts)
    ell0 = sum(data.starts[i][0] * linear[i] for i in range(len(data.out)))
    ell1 = sum(data.starts[i][1] * linear[i] for i in range(len(data.out)))
    return (
        data.hamilton,
        data.mass,
        data.mass // 2,
        data.mass // 2 + data.hamilton,
        q00,
        q01,
        q11,
        ell0,
        ell1,
    )


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
    one = base.tournament_data(base.transitive(1))
    cycle = base.tournament_data(base.parse("101", 3))
    p7 = one
    for _ in range(6):
        p7 = base.ordinal_data(p7, one)
    prefix = ordinal(cycle, cycle, p7)

    fixed = prefix_coordinates(prefix)
    expected_fixed = (
        9,
        18414,
        9207,
        9216,
        17031141,
        17031141,
        17031222,
        271372032,
        271454814,
    )
    need(fixed == expected_fixed, "exact Z7 prefix coordinates")
    ph, pw, ps0, ps1, pq00, pq01, pq11, pell0, pell1 = map(
        sp.Integer, expected_fixed
    )

    h, w = sp.symbols("h w")
    ell0, ell1 = sp.symbols("ell0 ell1")
    q00, q01, q11 = sp.symbols("q00 q01 q11")
    n00, n01, n10, n11 = sp.symbols("n00 n01 n10 n11")
    a = w + h
    b = w + 2 * h
    s0 = w / 2
    s1 = w / 2 + h

    # Exact transfer from B to T=Z7 ▷ B.  The mixed n_ef are the same-vertex
    # moments sum_i U_i^e(B)V_i^f(B).
    ht = ph * h
    wt = pw * b + 2 * ph * a
    ts0 = wt / 2
    ts1 = wt / 2 + ht

    prefix_u_square = pq00 + 2 * pq01 + pq11
    common_q = b**2 * prefix_u_square / 4
    tq00 = ph**2 * q00 + common_q
    tq01 = ph**2 * q01 + common_q
    tq11 = ph**2 * q11 + common_q

    prefix_weighted_cross = 2 * (
        (pq00 + pq01) * s0 + (pq01 + pq11) * s1
    )
    common_ell = b * (
        h * (pell0 + pell1)
        + (pw + 2 * ph) * a * (ps0 + ps1)
        + 3 * prefix_weighted_cross
    ) / 2
    tell0 = (
        common_ell
        + ph**2 * ell0
        + ph * (wt - ph * w) * s0
        - 2 * ph * (ps0 * n00 + ps1 * n01)
    )
    tell1 = (
        common_ell
        + ph**2 * ell1
        + ph * (wt - ph * w) * s1
        - 2 * ph * (ps0 * n10 + ps1 * n11)
    )

    jet_b = left_jet(h, w, s0, s1, q00, q01, q11, ell0, ell1)
    jet_t = left_jet(ht, wt, ts0, ts1, tq00, tq01, tq11, tell0, tell1)
    derived = tuple(sp.factor(x - y) for x, y in zip(jet_t, jet_b))
    expected = (
        h * (331614 * h + 165887 * w),
        2
        * (
            1086773760 * h**2
            + 1087499358 * h * w
            + 80 * ell0
            - 165726 * n00
            - 165888 * n01
            + 272056239 * w**2
        ),
        2
        * (
            1087105374 * h**2
            + 1087665165 * h * w
            + 80 * ell1
            - 165726 * n10
            - 165888 * n11
            + 272056239 * w**2
        ),
        h * (331614 * h + 165887 * w),
        165887 * h * (2 * h + w),
        1087561728 * h**2
        + 1087893342 * h * w
        + 480 * q00
        + 272056279 * w**2,
        2
        * (
            1087893342 * h**2
            + 1088059229 * h * w
            + 480 * q01
            + 272056279 * w**2
        ),
        1088225116 * h**2
        + 1088225116 * h * w
        + 480 * q11
        + 272056279 * w**2,
        5
        * (
            -651564000 * h**2
            - 651895614 * h * w
            + 32 * q00
            - 163056847 * w**2
        ),
        10
        * (
            -651895614 * h**2
            - 652061501 * h * w
            + 32 * q01
            - 163056847 * w**2
        ),
        5
        * (
            -652227388 * h**2
            - 652227388 * h * w
            + 32 * q11
            - 163056847 * w**2
        ),
    )
    for index, (actual, target) in enumerate(zip(derived, expected)):
        same(actual, target, f"D coordinate {index}")
    d = expected

    # Canonicalize C through k=H(C), z=W(C), v_i=V_i^0+V_i^1,
    # t_i=Start_i(C), and m=sum v_i^2, p=sum v_i t_i, tau=sum t_i^2.
    k, z = sp.symbols("k z")
    m, p, tau = sp.symbols("m p tau")
    cell0, cell1 = sp.symbols("cell0 cell1")
    cm00 = (m - 2 * p + tau) / 4
    cm01 = (m - tau) / 4
    cm11 = (m + 2 * p + tau) / 4
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
    coefficient_kz = sp.factor(d[0] + (d[1] + d[2]) / 2)
    coefficient_z2 = sp.factor((d[5] + d[6] + d[7]) / 4)
    coefficient_zk = sp.factor(d[6] / 2 + d[7])
    coefficient_m = sp.factor((d[8] + d[9] + d[10]) / 4)
    coefficient_p = sp.factor((d[10] - d[8]) / 2)
    coefficient_tau = sp.factor((d[8] - d[9] + d[10]) / 4)
    coefficient_k2 = sp.factor(-d[10])
    canonical_slack = (
        coefficient_kz * z * k
        + d[3] * cell0
        + d[4] * cell1
        + coefficient_z2 * z**2
        + coefficient_zk * z * k
        + coefficient_m * m
        + coefficient_p * p
        + coefficient_tau * tau
        + coefficient_k2 * k**2
    )
    same(context_slack, canonical_slack, "canonical right-context slack")

    # B-side coefficient floors and their exact nonnegative debts.
    min_kz = 2174044860 * h**2 + 2174998715 * h * w + 543946671 * w**2
    min_z2 = 1087893382 * h**2 + 1088059229 * h * w + 272056279 * w**2
    min_zk = (2 * h + w) * (1088059229 * h + 544112558 * w)
    min_m = -5 * (651895654 * h**2 + 652061501 * h * w + 163056847 * w**2)
    min_p = -5 * h * (331694 * h + 165887 * w)
    min_tau = -200 * h**2
    min_k2 = 815284195 * (2 * h + w) ** 2

    n0_defect = a * w / 2 - (n00 + n10)
    n1_defect = a * b / 2 - (n01 + n11)
    same(
        coefficient_kz - min_kz,
        80 * (ell0 + ell1) + 165726 * n0_defect + 165888 * n1_defect,
        "kz coefficient debt",
    )
    qsum = q00 + 2 * q01 + q11
    same(coefficient_z2 - min_z2, 120 * qsum, "z2 coefficient debt")
    same(coefficient_zk - min_zk, 480 * (q01 + q11), "zk coefficient debt")
    same(coefficient_m - min_m, 40 * qsum, "m coefficient debt")
    same(coefficient_p - min_p, 80 * (q11 - q00), "p coefficient debt")
    same(
        coefficient_tau - min_tau,
        40 * (q00 - 2 * q01 + q11),
        "tau coefficient debt",
    )
    same(
        coefficient_k2 - min_k2,
        40 * (b**2 - 4 * q11),
        "k2 coefficient debt",
    )

    ell_average = h * (331694 * h + 165887 * w)
    ell_difference = 80 * h**2
    same(
        d[3] * cell0 + d[4] * cell1,
        ell_average * (cell0 + cell1)
        + ell_difference * (cell1 - cell0),
        "right linear endpoint split",
    )

    # Use g_i=z-d_i-4r_i>=-4z and the endpoint-energy bounds
    # m<=(z^2+4zk+3k^2)/3, p<=(z+k)k, tau<=k^2.
    lower_bound = sp.expand(
        min_kz * z * k
        - 4 * z * (ell_average * (z + k) + ell_difference * k)
        + min_z2 * z**2
        + min_zk * z * k
        + min_m * (z**2 + 4 * z * k + 3 * k**2) / 3
        + min_p * (z + k) * k
        + min_tau * k**2
        + min_k2 * k**2
    )
    advertised = (
        (221548 * h**2 + 1879538 * h * w + 884602 * w**2) * z**2 / 3
        + (3620176 * h**2 + 8140211 * h * w + 3040747 * w**2) * z * k / 3
        - 40 * (2 * h + w) ** 2 * k**2
    )
    same(lower_bound, advertised, "advertised right-context lower bound")
    terminal = (
        k**2
        * (3841244 * h**2 + 10019269 * h * w + 3925229 * w**2)
        / 3
    )
    same(advertised.subs(z, k), terminal, "non-singleton z=k boundary")

    # The singleton right jet has nonzero entries 2,7,10.
    unary = sp.expand(d[2] + d[7] + d[10])
    unary_identity = 2 * (
        649462 * h**2
        + 1209253 * h * w
        + 80 * ell1
        - 165726 * n10
        - 165888 * n11
        + 320 * q11
        + 442261 * w**2
    )
    same(unary, unary_identity, "singleton-right exact identity")
    unary_floor = 967148 * h**2 + 1921004 * h * w + 718715 * w**2 + 640 * q11
    mixed_defect = (
        ph * b * (ps0 * w + ps1 * b)
        - 4 * ph * ps0 * n10
        - 4 * ph * ps1 * n11
    )
    same(
        unary - unary_floor,
        160 * ell1 + mixed_defect,
        "singleton-right nonnegative debt",
    )
    unary_slack = sp.expand(unary_floor - FLOOR * h**2)
    expected_unary_slack = (
        -640 * h**2 + 1921004 * h * w + 718715 * w**2 + 640 * q11
    )
    same(unary_slack, expected_unary_slack, "sharp unary floor")
    same(
        unary_slack.subs({w: h, q11: 0}),
        2639079 * h**2,
        "non-singleton middle boundary",
    )

    # Exact singleton obstruction and crossing formula.
    singleton_formula = lambda n: 108 * (2 * 4**n - (12 * n + 102) * 2**n + 1)
    singleton_values: list[int] = []
    p_tail: base.TournamentData | None = None
    for n in range(8):
        if n == 0:
            z_n = ordinal(cycle, cycle)
        else:
            p_tail = one if p_tail is None else base.ordinal_data(p_tail, one)
            z_n = ordinal(cycle, cycle, p_tail)
        value = response(z_n, one, one)
        need(value == singleton_formula(n), f"singleton formula n={n}")
        singleton_values.append(value)
    need(
        tuple(singleton_values)
        == (-10692, -23652, -50868, -105300, -203796, -338580, -317844, FLOOR),
        "sharp singleton transition",
    )

    # Small exact controls of the consequence object, not assumptions.
    contexts: list[base.TournamentData] = []
    for order in range(1, 4):
        width = order * (order - 1) // 2
        for integer in range(1 << width):
            bits = "" if width == 0 else format(integer, f"0{width}b")
            contexts.append(base.tournament_data(base.parse(bits, order)))
    control_rows = 0
    equality_rows = 0
    for middle, right in itertools.product(contexts, repeat=2):
        value = response(prefix, middle, right)
        target = FLOOR * middle.hamilton**2 * right.hamilton**2
        need(value >= target, "small-context uniform floor")
        equality = value == target
        if equality:
            need(len(middle.out) == len(right.out) == 1, "small equality boundary")
        equality_rows += int(equality)
        control_rows += 1
    need(control_rows == 121 and equality_rows == 1, "small-context summary")

    lines = [
        "theorem=THM-4215",
        "prefix=Z7=C3_tri_C3_tri_P7",
        "prefix_order=13",
        "prefix_coordinates=" + ",".join(str(value) for value in expected_fixed),
        "symbolic_D_coordinates=11",
        "D_coordinate_audit=PASS",
        "right_context_canonicalization=PASS",
        "right_context_debt_decomposition=PASS",
        "right_context_LB_numerator=(221548h2+1879538hw+884602w2)z2+(3620176h2+8140211hw+3040747w2)zk-120(2h+w)^2k2",
        "right_context_z_eq_k_numerator=3841244h2+10019269hw+3925229w2",
        "unary_identity_audit=PASS",
        "unary_floor=967148h2+1921004hw+718715w2+640q11",
        "sharp_uniform_floor=967788",
        "singleton_formula=108*(2*4^n-(12n+102)*2^n+1)",
        "singleton_values_n0_to_n7=" + ",".join(str(value) for value in singleton_values),
        f"small_labelled_context_rows={control_rows}",
        f"small_labelled_equalities={equality_rows}",
        "equality_boundary=n7,B=P1,C=P1",
    ]
    payload = "\n".join(lines) + "\n"
    print(payload, end="")
    print(f"semantic_sha256={hashlib.sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
