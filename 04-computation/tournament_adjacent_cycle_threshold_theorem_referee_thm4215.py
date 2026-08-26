#!/usr/bin/env python3
"""Clean-room theorem referee for THM-4215.

This audit imports no tournament engine from the repository.  It rebuilds the
rooted Z_0 and Z_7 data from labelled adjacency, Hamilton-subset DP, and
directed-simple-path subset DP.  Starting from the proved THM-4208 response
jet, it then independently checks the eleven displayed response coordinates,
every coefficient debt, the signs used in the moment substitutions, the
strict equality boundary, and the all-n singleton formula via universal-sink
recurrences.

The already frozen literal C++ audit independently evaluates G_+ itself.  This
script is complementary: it audits the rooted coordinates and analytic cone
argument on which the all-context quantifiers depend.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


FLOOR = 967_788


def need(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    difference = sp.expand(left - right)
    need(difference == 0 or sp.simplify(difference) == 0, label)


@dataclass(frozen=True)
class Tournament:
    out: tuple[int, ...]

    @property
    def order(self) -> int:
        return len(self.out)

    def arc(self, source: int, target: int) -> bool:
        return bool(self.out[source] & (1 << target))


def transitive(order: int) -> Tournament:
    rows = [0] * order
    for source in range(order):
        for target in range(source + 1, order):
            rows[source] |= 1 << target
    return Tournament(tuple(rows))


def cycle3() -> Tournament:
    return Tournament((1 << 1, 1 << 2, 1 << 0))


def ordinal(left: Tournament, right: Tournament) -> Tournament:
    shift = left.order
    right_mask = ((1 << right.order) - 1) << shift
    rows = [row | right_mask for row in left.out]
    rows.extend(row << shift for row in right.out)
    return Tournament(tuple(rows))


def converse(tournament: Tournament) -> Tournament:
    rows = [0] * tournament.order
    for source in range(tournament.order):
        for target in range(tournament.order):
            if source != target and tournament.arc(target, source):
                rows[source] |= 1 << target
    return Tournament(tuple(rows))


def hamilton_table(tournament: Tournament) -> list[int]:
    n = tournament.order
    states = 1 << n
    end = [[0] * n for _ in range(states)]
    hamilton = [0] * states
    hamilton[0] = 1
    for mask in range(1, states):
        total = 0
        for last in range(n):
            bit = 1 << last
            if not mask & bit:
                continue
            if mask == bit:
                count = 1
            else:
                previous = mask ^ bit
                count = sum(
                    end[previous][before]
                    for before in range(n)
                    if previous & (1 << before) and tournament.arc(before, last)
                )
            end[mask][last] = count
            total += count
        hamilton[mask] = total
    return hamilton


def rooted_states_and_capacity(
    tournament: Tournament, hamilton: list[int]
) -> tuple[list[list[int]], list[list[int]]]:
    n = tournament.order
    states = 1 << n
    full = states - 1
    rooted = [[0, 0] for _ in range(n)]
    capacity = [[0] * n for _ in range(n)]
    for start in range(n):
        path = [[0] * n for _ in range(states)]
        for mask in range(1, states):
            if not mask & (1 << start):
                continue
            parity = mask.bit_count() & 1
            complement_paths = hamilton[full ^ mask]
            for last in range(n):
                bit = 1 << last
                if not mask & bit:
                    continue
                if mask == (1 << start) and last == start:
                    count = 1
                else:
                    previous = mask ^ bit
                    count = sum(
                        path[previous][before]
                        for before in range(n)
                        if previous & (1 << before)
                        and tournament.arc(before, last)
                    )
                path[mask][last] = count
                weighted = count * complement_paths
                rooted[start][parity] += weighted
                if last != start and parity == 0:
                    low, high = sorted((start, last))
                    capacity[low][high] += 2 * weighted
                    capacity[high][low] += 2 * weighted
    return rooted, capacity


@dataclass(frozen=True)
class RootedData:
    h: int
    w: int
    s0: int
    s1: int
    q00: int
    q01: int
    q11: int
    ell0: int
    ell1: int
    mass_square: int
    q01_plus_q11: int
    ell_sum: int


def analyze(tournament: Tournament) -> RootedData:
    hamilton = hamilton_table(tournament)
    rooted_u, capacity = rooted_states_and_capacity(tournament, hamilton)
    rooted_v, _ = rooted_states_and_capacity(converse(tournament), hamilton)
    n = tournament.order
    h = hamilton[-1]
    w = sum(capacity[i][j] for i in range(n) for j in range(i + 1, n))
    degree = [sum(capacity[i]) for i in range(n)]
    outgoing = [
        sum(capacity[i][j] for j in range(n) if i != j and tournament.arc(i, j))
        for i in range(n)
    ]
    ell = [
        sum(
            rooted_u[i][parity] * (w - degree[i] + 4 * outgoing[i])
            for i in range(n)
        )
        for parity in range(2)
    ]
    s0 = sum(row[0] for row in rooted_u)
    s1 = sum(row[1] for row in rooted_u)
    need((s0, s1) == (w // 2, w // 2 + h), "rooted total identity")
    need(
        [sum(row[p] for row in rooted_v) for p in range(2)] == [s0, s1],
        "start/end total agreement",
    )
    q00 = sum(row[0] ** 2 for row in rooted_u)
    q01 = sum(row[0] * row[1] for row in rooted_u)
    q11 = sum(row[1] ** 2 for row in rooted_u)
    return RootedData(
        h,
        w,
        s0,
        s1,
        q00,
        q01,
        q11,
        ell[0],
        ell[1],
        sum((row[0] + row[1]) ** 2 for row in rooted_u),
        q01 + q11,
        ell[0] + ell[1],
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


def nonnegative_coefficients(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool:
    return all(coefficient >= 0 for coefficient in sp.Poly(sp.expand(expr), *variables).coeffs())


def main() -> None:
    c3 = cycle3()
    z0 = ordinal(c3, c3)
    z7 = ordinal(z0, transitive(7))
    raw0 = analyze(z0)
    raw7 = analyze(z7)
    need(
        (
            raw7.h,
            raw7.w,
            raw7.s0,
            raw7.s1,
            raw7.q00,
            raw7.q01,
            raw7.q11,
            raw7.ell0,
            raw7.ell1,
        )
        == (
            9,
            18414,
            9207,
            9216,
            17031141,
            17031141,
            17031222,
            271372032,
            271454814,
        ),
        "clean-room Z7 rooted coordinates",
    )
    need(
        (raw0.w, raw0.mass_square, raw0.q01_plus_q11, raw0.ell_sum)
        == (126, 4131, 2106, 27540),
        "clean-room Z0 sink-recurrence seed",
    )

    h, w = sp.symbols("h w", positive=True)
    ell0, ell1 = sp.symbols("ell0 ell1", nonnegative=True)
    q00, q01, q11 = sp.symbols("q00 q01 q11", nonnegative=True)
    n00, n01, n10, n11 = sp.symbols("n00 n01 n10 n11", nonnegative=True)
    a, b = w + h, w + 2 * h
    s0, s1 = w / 2, w / 2 + h
    fixed = tuple(
        map(
            sp.Integer,
            (
                raw7.h,
                raw7.w,
                raw7.s0,
                raw7.s1,
                raw7.q00,
                raw7.q01,
                raw7.q11,
                raw7.ell0,
                raw7.ell1,
            ),
        )
    )
    ah, aw, as0, as1, aq00, aq01, aq11, aell0, aell1 = fixed

    # Direct block-transfer derivation for T=Z7 ▷ B.
    th = ah * h
    tw = aw * b + 2 * ah * a
    ts0, ts1 = tw / 2, tw / 2 + th
    gram_common = b**2 * (aq00 + 2 * aq01 + aq11) / 4
    tq00 = ah**2 * q00 + gram_common
    tq01 = ah**2 * q01 + gram_common
    tq11 = ah**2 * q11 + gram_common
    cross_rooted = 2 * (
        (aq00 + aq01) * s0 + (aq01 + aq11) * s1
    )
    common_linear = b * (
        h * (aell0 + aell1)
        + (aw + 2 * ah) * a * (as0 + as1)
        + 3 * cross_rooted
    ) / 2
    tell0 = (
        common_linear
        + ah**2 * ell0
        + ah * (tw - ah * w) * s0
        - 2 * ah * (as0 * n00 + as1 * n01)
    )
    tell1 = (
        common_linear
        + ah**2 * ell1
        + ah * (tw - ah * w) * s1
        - 2 * ah * (as0 * n10 + as1 * n11)
    )
    d = tuple(
        sp.expand(left - right)
        for left, right in zip(
            left_jet(th, tw, ts0, ts1, tq00, tq01, tq11, tell0, tell1),
            left_jet(h, w, s0, s1, q00, q01, q11, ell0, ell1),
        )
    )
    advertised_d = (
        h * (331614 * h + 165887 * w),
        2 * (1086773760 * h**2 + 1087499358 * h * w + 80 * ell0 - 165726 * n00 - 165888 * n01 + 272056239 * w**2),
        2 * (1087105374 * h**2 + 1087665165 * h * w + 80 * ell1 - 165726 * n10 - 165888 * n11 + 272056239 * w**2),
        h * (331614 * h + 165887 * w),
        165887 * h * (2 * h + w),
        1087561728 * h**2 + 1087893342 * h * w + 480 * q00 + 272056279 * w**2,
        2 * (1087893342 * h**2 + 1088059229 * h * w + 480 * q01 + 272056279 * w**2),
        1088225116 * h**2 + 1088225116 * h * w + 480 * q11 + 272056279 * w**2,
        5 * (-651564000 * h**2 - 651895614 * h * w + 32 * q00 - 163056847 * w**2),
        10 * (-651895614 * h**2 - 652061501 * h * w + 32 * q01 - 163056847 * w**2),
        5 * (-652227388 * h**2 - 652227388 * h * w + 32 * q11 - 163056847 * w**2),
    )
    for index, (actual, expected) in enumerate(zip(d, advertised_d)):
        same(actual, expected, f"response coordinate D_{index}")

    # Extract the right-context coefficients directly from D dot mu.
    k, z, m, p, tau = sp.symbols("k z m p tau", nonnegative=True)
    lm0, lm1 = sp.symbols("lm0 lm1")
    mu_delta = (
        k * z,
        k * z / 2,
        k * z / 2,
        lm0,
        lm1,
        z**2 / 4,
        z**2 / 4 + z * k / 2,
        z**2 / 4 + z * k,
        (m - 2 * p + tau) / 4,
        (m - tau) / 4,
        (m + 2 * p + tau) / 4 - k**2,
    )
    slack = sp.expand(sum(left * right for left, right in zip(d, mu_delta)))
    kz0 = sp.expand(d[0] + (d[1] + d[2]) / 2)
    z2 = sp.expand((d[5] + d[6] + d[7]) / 4)
    kz1 = sp.expand(d[6] / 2 + d[7])
    am = sp.expand((d[8] + d[9] + d[10]) / 4)
    ap = sp.expand((d[10] - d[8]) / 2)
    atau = sp.expand((d[8] - d[9] + d[10]) / 4)
    ak = sp.expand(-d[10])
    same(
        slack,
        kz0 * z * k + d[3] * lm0 + d[4] * lm1 + z2 * z**2
        + kz1 * z * k + am * m + ap * p + atau * tau + ak * k**2,
        "right-context canonicalization",
    )

    floors = (
        2174044860 * h**2 + 2174998715 * h * w + 543946671 * w**2,
        1087893382 * h**2 + 1088059229 * h * w + 272056279 * w**2,
        (2 * h + w) * (1088059229 * h + 544112558 * w),
        -5 * (651895654 * h**2 + 652061501 * h * w + 163056847 * w**2),
        -5 * h * (331694 * h + 165887 * w),
        -200 * h**2,
        815284195 * (2 * h + w) ** 2,
    )
    f_kz0, f_z2, f_kz1, f_m, f_p, f_tau, f_k = floors
    debts = (
        80 * (ell0 + ell1) + 165726 * (a * w / 2 - n00 - n10) + 165888 * (a * b / 2 - n01 - n11),
        120 * (q00 + 2 * q01 + q11),
        480 * (q01 + q11),
        40 * (q00 + 2 * q01 + q11),
        80 * (q11 - q00),
        40 * (q00 - 2 * q01 + q11),
        40 * (b**2 - 4 * q11),
    )
    for actual, floor, debt, label in zip(
        (kz0, z2, kz1, am, ap, atau, ak),
        floors,
        debts,
        ("K0", "K2", "K1", "Am", "Ap", "Atau", "Ak"),
    ):
        same(actual - floor, debt, f"{label} debt")
    need(all(nonnegative_coefficients(expr, (h, w)) for expr in floors[:3]), "positive K floors")
    need(all(nonnegative_coefficients(-expr, (h, w)) for expr in floors[3:6]), "negative moment floors")
    need(nonnegative_coefficients(f_k, (h, w)), "positive k2 floor")

    average_weight = h * (331694 * h + 165887 * w)
    endpoint_weight = 80 * h**2
    same(
        d[3] * lm0 + d[4] * lm1,
        average_weight * (lm0 + lm1) + endpoint_weight * (lm1 - lm0),
        "linear endpoint split",
    )
    moment_cap = (z**2 + 4 * z * k + 3 * k**2) / 3
    cone_lower = sp.expand(
        f_kz0 * z * k
        - 4 * z * (average_weight * (z + k) + endpoint_weight * k)
        + f_z2 * z**2
        + f_kz1 * z * k
        + f_m * moment_cap
        + f_p * (z + k) * k
        + f_tau * k**2
        + f_k * k**2
    )
    advertised_lower = sp.expand(
        (
            (221548 * h**2 + 1879538 * h * w + 884602 * w**2) * z**2
            + (3620176 * h**2 + 8140211 * h * w + 3040747 * w**2) * z * k
            - 120 * (2 * h + w) ** 2 * k**2
        )
        / 3
    )
    same(cone_lower, advertised_lower, "collected cone lower bound")
    boundary = sp.expand(3 * advertised_lower.subs(z, k) / k**2)
    same(
        boundary,
        3841244 * h**2 + 10019269 * h * w + 3925229 * w**2,
        "strict non-singleton right boundary",
    )
    need(nonnegative_coefficients(boundary, (h, w)), "positive right boundary coefficients")

    # Singleton-right face and unique middle equality.
    unary = sp.expand(d[2] + d[7] + d[10])
    unary_floor = 967148 * h**2 + 1921004 * h * w + 718715 * w**2 + 640 * q11
    mixed_debt = 9 * b * (9207 * w + 9216 * b) - 36 * 9207 * n10 - 36 * 9216 * n11
    same(unary - unary_floor, 160 * ell1 + mixed_debt, "unary nonnegative debt")
    same(
        unary_floor - FLOOR * h**2,
        -640 * h**2 + 1921004 * h * w + 718715 * w**2 + 640 * q11,
        "unary normalized slack",
    )
    same(
        (unary_floor - FLOOR * h**2).subs({w: h, q11: 0}),
        2639079 * h**2,
        "strict non-singleton middle boundary",
    )

    # Independent all-n singleton derivation.  T_m=Z_m and T_(m+1)=T_m ▷ P1.
    index = sp.symbols("index", integer=True, nonnegative=True)
    wm = 144 * 2**index - 18
    massm = 4158 * 4**index - 27
    rm = 2079 * 4**index + 27
    km = 33210 * 4**index - (648 * index + 5508) * 2**index - 162
    same(wm.subs(index, 0), raw0.w, "sink recurrence w seed")
    same(massm.subs(index, 0), raw0.mass_square, "sink recurrence M seed")
    same(rm.subs(index, 0), raw0.q01_plus_q11, "sink recurrence r seed")
    same(km.subs(index, 0), raw0.ell_sum, "sink recurrence K seed")
    same(wm.subs(index, index + 1), 2 * wm + 18, "sink recurrence w")
    same(massm.subs(index, index + 1), 4 * massm + 81, "sink recurrence M")
    same(rm.subs(index, index + 1), 2 * massm + 81, "sink recurrence r")
    same(
        km.subs(index, index + 1),
        2 * km + 2 * (wm + 18) * (wm + 9) + 12 * rm + 9 * wm,
        "sink recurrence K",
    )
    n = sp.symbols("n", integer=True, nonnegative=True)
    ell1_next = 66420 * 4**n - (648 * n + 5508) * 2**n - 162
    q11_next = 4158 * 4**n + 54
    s1_next = 144 * 2**n
    singleton = sp.expand(2 * ell1_next + 8 * q11_next - 8 * s1_next**2)
    advertised_singleton = 108 * (2 * 4**n - (12 * n + 102) * 2**n + 1)
    same(singleton, advertised_singleton, "all-n singleton formula")
    boundary_values = tuple(int(advertised_singleton.subs(n, value)) for value in range(8))
    need(
        boundary_values
        == (-10692, -23652, -50868, -105300, -203796, -338580, -317844, FLOOR),
        "n=0,...,7 singleton boundary",
    )

    print("theorem=THM-4215")
    print("audit_role=external_theorem_referee")
    print("repo_tournament_imports=NONE")
    print("raw_engine=labelled_adjacency+Hamilton_subset_DP+directed_path_subset_DP")
    print("z0_sink_seed=126,4131,2106,27540")
    print("z7_rooted_coordinates=9,18414,9207,9216,17031141,17031141,17031222,271372032,271454814")
    print("response_coordinates=PASS:11")
    print("coefficient_debts=PASS:7")
    print("moment_substitution_signs=PASS:K_positive;Am_Ap_Atau_negative;Ak_positive")
    print("right_boundary_numerator=3841244h2+10019269hw+3925229w2")
    print("middle_boundary_gap=2639079h2")
    print("all_n_singleton_sink_recurrence=PASS")
    print("singleton_values_n0_to_n7=" + ",".join(map(str, boundary_values)))
    print("equality_ledger=n7_and_B=P1_and_C=P1_only;later_tails_strict")
    print("propagation_ledger=exact_telescope+THM4187_nonnegative_transitive_prefix")
    print("osplus_ledger=n>=8_peel_last_P1+THM4187_no_sink_strictness")
    print("OVERALL=ACCEPT")


if __name__ == "__main__":
    main()
