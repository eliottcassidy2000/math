#!/usr/bin/env python3
"""Primary exact audit for THM-4449's dyadic seventh-rounding theorem.

The proof is infinite and appears in the theorem/result note.  This script
checks its finite residue tables, independently reconstructs every pair comb
through denominator 151 from strict quotient walls, checks the compatibility
calculation behind the all-height energy/physical bounds, and replays the
transfer and hostile controls.  Standard library only; checks remain active
under ``python -O``.
"""

from fractions import Fraction as Q
from itertools import combinations
from math import gcd


DELTA = Q(1, 14)
CHECKS = 0


def need(condition, detail):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(detail)


def gap(z):
    z %= 1
    return min(z, 1 - z)


def odd_three_unit(n):
    return n > 0 and n % 2 == 1 and n % 3 != 0


def seventh_residue(n):
    """Unique representative of n modulo 7 in {-3,...,3}."""
    return n - 7 * ((n + 3) // 7)


def pair_formula(a, b):
    """Haar mass of the two-lift cross-comb for distinct positive odd tails."""
    common = gcd(a, b)
    p, q = sorted((a // common, b // common))
    need(p > 0 and p < q and p % 2 and q % 2 and gcd(p, q) == 1, (a, b, "pair type"))
    alpha = (p + q) // 2
    beta = (q - p) // 2
    d = seventh_residue(alpha)
    e = seventh_residue(beta)
    need(-3 <= d <= 3 and -3 <= e <= 3, (p, q, d, e, "seventh residues"))
    return Q(2, 49) * (1 + Q(e * e - d * d, p * q))


def quotient_walls(tails):
    """Every strict-danger wall for lifts (y+j)/2, j=0,1, on 0<=y<=1."""
    walls = {Q(0), Q(1)}
    for tail in tails:
        need(tail > 0 and tail % 2 == 1, (tails, "odd positive tails"))
        for label in (0, 1):
            # tail*(y+label)/2 = k +/- 1/14.
            for k in range(tail + 1):
                for sign in (-1, 1):
                    y = (Q(2 * k) + Q(sign, 7)) / tail - label
                    if 0 <= y <= 1:
                        walls.add(y)
    return sorted(walls)


def quotient_failure(tails, y):
    return all(
        any(gap(Q(tail) * (y + label) / 2) < DELTA for tail in tails)
        for label in (0, 1)
    )


def quotient_stats(tails):
    """Exact open-set mass, longest component, and circular component count."""
    walls = quotient_walls(tuple(tails))
    live = []
    total = Q(0)
    for left, right in zip(walls, walls[1:]):
        value = quotient_failure(tails, (left + right) / 2)
        live.append(value)
        if value:
            total += right - left

    live_indices = [i for i, value in enumerate(live) if value]
    if not live_indices:
        return Q(0), Q(0), 0

    parent = list(range(len(live)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        i, j = find(i), find(j)
        if i != j:
            parent[j] = i

    for i in range(1, len(walls) - 1):
        if live[i - 1] and live[i] and quotient_failure(tails, walls[i]):
            union(i - 1, i)
    if live[-1] and live[0] and quotient_failure(tails, Q(0)):
        union(len(live) - 1, 0)

    lengths = {}
    for i in live_indices:
        root = find(i)
        lengths[root] = lengths.get(root, Q(0)) + walls[i + 1] - walls[i]
    return total, max(lengths.values()), len(lengths)


def ratio_measure(ratio):
    if ratio < 1:
        ratio = 1 / ratio
    return pair_formula(ratio.denominator, ratio.numerator)


def energy(labels):
    total = Q(0)
    for a, b in combinations(labels, 2):
        ratio = Q(a) / Q(b)
        total += ratio_measure(ratio)
    return total


def merge(intervals):
    merged = []
    for left, right in sorted(intervals):
        if left == right:
            continue
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        else:
            merged[-1][1] = max(merged[-1][1], right)
    return tuple((left, right) for left, right in merged)


def danger_union(speeds):
    intervals = []
    for speed in set(speeds):
        radius = Q(1, 14 * speed)
        for k in range(speed):
            left = Q(k, speed) - radius
            right = Q(k, speed) + radius
            if left < 0:
                intervals.extend(((Q(0), right), (1 + left, Q(1))))
            elif right > 1:
                intervals.extend(((left, Q(1)), (Q(0), right - 1)))
            else:
                intervals.append((left, right))
    return merge(intervals)


def safe_measure(speeds):
    return 1 - sum((right - left for left, right in danger_union(speeds)), Q(0))


def clearance(speeds, x):
    return min(gap(Q(speed) * x) for speed in set(speeds))


def lift_masks(sheet_count, tails, y):
    return tuple(
        tuple(
            label
            for label in range(sheet_count)
            if gap(Q(tail) * (y + label) / sheet_count) < DELTA
        )
        for tail in tails
    )


def check_pair_formula_and_residue_tables():
    literal_checks = 0
    for q in range(3, 152, 2):
        for p in range(1, q, 2):
            if gcd(p, q) != 1:
                continue
            literal = quotient_stats((p, q))[0]
            formula = pair_formula(p, q)
            need(literal == formula, (p, q, literal, formula, "literal/formula"))
            literal_checks += 1
    need(literal_checks == 2350, (literal_checks, "literal check count"))

    low33 = []
    for p in range(1, 33, 2):
        for q in range(p + 2, 33, 2):
            if p * q < 33 and gcd(p, q) == 1 and odd_three_unit(p) and odd_three_unit(q):
                low33.append(((p, q), pair_formula(p, q)))
    expected_low33 = [
        ((1, 5), Q(0)),
        ((1, 7), Q(2, 49)),
        ((1, 11), Q(4, 77)),
        ((1, 13), Q(4, 91)),
        ((1, 17), Q(4, 119)),
        ((1, 19), Q(4, 133)),
        ((1, 23), Q(8, 161)),
        ((1, 25), Q(8, 175)),
        ((1, 29), Q(8, 203)),
        ((1, 31), Q(8, 217)),
    ]
    need(low33 == expected_low33, (low33, "product below 33"))
    need(max(value for _, value in low33) == Q(4, 77), (low33, "ternary pair cap"))
    need([ratio for ratio, value in low33 if value == Q(4, 77)] == [(1, 11)], (low33, "pair equality"))

    high = []
    boundary = []
    for p in range(1, 117, 2):
        for q in range(p + 2, 117, 2):
            if p * q >= 117 or gcd(p, q) != 1 or not (odd_three_unit(p) and odd_three_unit(q)):
                continue
            value = pair_formula(p, q)
            if value > Q(4, 91):
                high.append(((p, q), value))
            elif value == Q(4, 91):
                boundary.append((p, q))
    expected_high = [
        ((1, 11), Q(4, 77)),
        ((1, 23), Q(8, 161)),
        ((1, 25), Q(8, 175)),
        ((1, 37), Q(12, 259)),
        ((5, 11), Q(18, 385)),
    ]
    need(high == expected_high, (high, "high alphabet"))
    need(boundary == [(1, 13), (1, 65), (5, 13)], (boundary, "4/91 boundary"))
    return literal_checks, high, boundary


def compatibility_tables():
    exceptional = (("23", Q(23)), ("11/5", Q(11, 5)), ("37", Q(37)), ("25", Q(25)))
    one_eleven = []
    for name, ratio in exceptional:
        values = []
        for endpoint in (Q(1), Q(11)):
            for oriented in (ratio, 1 / ratio):
                third = endpoint * oriented
                if third not in (Q(1), Q(11)):
                    values.append(energy((Q(1), Q(11), third)))
        one_eleven.append((name, max(values)))

    no_eleven = []
    for _, first in exceptional:
        row = []
        for _, second in exceptional:
            values = []
            for x in (first, 1 / first):
                for y in (second, 1 / second):
                    if x == y:
                        continue
                    ratios = tuple(z if z >= 1 else 1 / z for z in (x, y, x / y))
                    if Q(11) not in ratios:
                        values.append(sum((ratio_measure(z) for z in ratios), Q(0)))
            row.append(max(values))
        no_eleven.append(row)
    return one_eleven, no_eleven


def check_energy_and_physical_bounds():
    one_eleven, no_eleven = compatibility_tables()
    expected_one = [
        ("23", Q(36, 253)),
        ("11/5", Q(592, 4235)),
        ("37", Q(36, 259)),
        ("25", Q(268, 1925)),
    ]
    expected_matrix = [
        [Q(520, 3703), Q(174, 1265), Q(816, 5957), Q(548, 4025)],
        [Q(174, 1265), Q(406, 3025), Q(382, 2849), Q(258, 1925)],
        [Q(816, 5957), Q(382, 2849), Q(1280, 9583), Q(172, 1295)],
        [Q(548, 4025), Q(258, 1925), Q(172, 1295), Q(116, 875)],
    ]
    need(one_eleven == expected_one, (one_eleven, "one 1:11 table"))
    need(no_eleven == expected_matrix, (no_eleven, "two-exception table"))

    gp = (1, 11, 121)
    gp_energy = energy(tuple(Q(x) for x in gp))
    need(gp_energy == Q(124, 847), (gp_energy, "sharp energy"))
    case_bounds = (
        Q(8, 161) + 2 * Q(4, 91),
        max(max(row) for row in no_eleven),
        Q(4, 77) + 2 * Q(4, 91),
        max(value for _, value in one_eleven),
    )
    need(max(case_bounds) == Q(36, 253), (case_bounds, "non-GP cap"))

    gp_stats = quotient_stats(gp)
    need(gp_stats == (Q(108, 847), Q(2, 77), 34), (gp_stats, "GP physical union"))
    need(gp_stats[0] < Q(36, 253) < gp_energy, (gp_stats, gp_energy, "energy/union distinction"))

    leader = (1, 7, 11)
    leader_stats = quotient_stats(leader)
    need(leader_stats == (Q(72, 539), Q(2, 77), 6), (leader_stats, "leader quotient geometry"))
    # Doubling pullback preserves mass, halves widths, and doubles components.
    physical_stats = (leader_stats[0], leader_stats[1] / 2, 2 * leader_stats[2])
    need(physical_stats == (Q(72, 539), Q(1, 77), 12), (physical_stats, "leader physical geometry"))

    for dilation in (5, 7, 11):
        need(odd_three_unit(dilation), (dilation, "equality-class dilation"))
        need(energy(tuple(Q(dilation * x) for x in gp)) == gp_energy, (dilation, "energy dilation"))
        need(quotient_stats(tuple(dilation * x for x in gp))[0] == gp_stats[0], (dilation, "union dilation"))
    return one_eleven, no_eleven, gp_energy, gp_stats, leader_stats, physical_stats


def check_transfer_and_hostiles():
    transfer_samples = []
    samples = (
        ((1, 2, 3), 5),
        ((1, 4, 7, 10), 11),
        ((2, 3, 5, 8, 13), 17),
        (tuple(range(1, 11)), 13),
    )
    for body, r in samples:
        mu_body = safe_measure(body)
        mu_q4 = safe_measure(tuple(2 * c for c in body) + (r,))
        mu_q2 = safe_measure(tuple(body) + (r,))
        need(mu_q4 >= mu_body / 2, (body, r, mu_body, mu_q4, "q4 half-mass"))
        need(mu_q2 >= mu_body - Q(1, 7), (body, r, mu_body, mu_q2, "q2 deletion loss"))
        transfer_samples.append((len(body), r, mu_body, mu_q4, mu_q2))
    need(Q(8, 77) / 2 == Q(4, 77), "q4 entry arithmetic")
    need(Q(15, 77) - Q(1, 7) == Q(4, 77), "q2 entry arithmetic")
    need(Q(8, 91) / 2 == Q(4, 91), "q4 localization arithmetic")
    need(Q(17, 91) - Q(1, 7) == Q(4, 91), "q2 localization arithmetic")

    beta_expected = {
        (1, 11): Q(2, 77),
        (1, 23): Q(2, 161),
        (5, 11): Q(9, 385),
        (1, 37): Q(2, 259),
        (1, 25): Q(2, 175),
    }
    beta = {}
    for ratio, expected in beta_expected.items():
        stats = quotient_stats(ratio)
        need(stats[0] == pair_formula(*ratio) and stats[1] == expected, (ratio, stats, expected, "beta"))
        beta[ratio] = stats[1]

    body0 = tuple(range(1, 11))
    y0 = Q(1, 11)
    masks0 = lift_masks(2, (1, 7, 11), y0)
    need(clearance(body0, y0) == Q(1, 11), "zero-even hostile body")
    need(masks0 == ((0,), (), (1,)), (masks0, "zero-even masks"))
    row0 = tuple(2 * c for c in body0) + (1, 7, 11)
    need(clearance(row0, Q(181, 2352)) == Q(15, 196), "zero-even positive row")

    # The original even tail is 26=2*13; absorbed r=13 is an odd 3-unit.
    r1 = 13
    need(odd_three_unit(r1), (r1, "absorbed speed"))
    h1 = body0 + (r1,)
    masks1 = lift_masks(2, (1, 11), y0)
    need(clearance(h1, y0) == Q(1, 11), "one-even absorbed body")
    need(masks1 == ((0,), (1,)), (masks1, "one-even masks"))
    row1 = tuple(2 * c for c in body0) + (26, 1, 11)
    need(clearance(row1, Q(229, 560)) == Q(5, 56), "one-even positive row")

    body2 = (1, 2, 3, 4, 5, 7, 8, 9, 10, 11)
    tails2 = (1, 11, 14)
    y2 = Q(86, 539)
    masks2 = lift_masks(4, tails2, y2)
    need(clearance(body2, y2) == Q(9, 77), "q4 hostile body")
    need(masks2 == ((0,), (2,), (1, 3)), (masks2, "q4 masks"))
    row2 = tuple(4 * c for c in body2) + tails2
    need(clearance(row2, Q(37, 480)) == Q(3, 40), "q4 positive row")
    hostiles = ((masks0, Q(15, 196)), (masks1, Q(5, 56)), (masks2, Q(3, 40)))
    return transfer_samples, beta, hostiles


def main():
    literal_checks, high, boundary = check_pair_formula_and_residue_tables()
    one_eleven, matrix, gp_energy, gp_stats, leader_q, leader_x = check_energy_and_physical_bounds()
    transfers, beta, hostiles = check_transfer_and_hostiles()

    print("LRC14_DYADIC_SEVENTH_ROUNDING_THM4449_PRIMARY")
    print("STATUS=EXACT_SUPPORT_FOR_PROVED_RESULTS;LRC14_OPEN")
    print(f"formula_literal_checks={literal_checks}")
    print("ternary_pair_cap=4/77 equality_primitive=(1,11) equality_dilations=odd_3_units")
    print("high_edge_alphabet=" + ",".join(f"{p}:{q}:{value}" for (p, q), value in high))
    print("four_over_91_equalities=" + ",".join(f"{p}:{q}" for p, q in boundary))
    print(f"pair_energy_max={gp_energy} equality_shape=(1,11,121) equality_dilations=odd_3_units")
    print("one_1_11_edge_maxima=" + ",".join(str(value) for _, value in one_eleven))
    print("compatibility_max=" + str(max(max(row) for row in matrix)))
    print(f"energy_shape_quotient_union=mass:{gp_stats[0]},longest:{gp_stats[1]},components:{gp_stats[2]}")
    print("proved_physical_owner_cut_cap=36/253 not_claimed_sharp")
    print(
        f"height199_leader_control=(1,7,11) mass:{leader_x[0]} "
        f"physical_x_longest:{leader_x[1]} physical_x_components:{leader_x[2]} "
        f"quotient_y_longest:{leader_q[1]} quotient_y_components:{leader_q[2]}"
    )
    print("beta_table=" + ",".join(f"{p}:{q}:{value}" for (p, q), value in beta.items()))
    print("transfer_samples=" + ",".join(f"{n}:{r}:{mc}:{m4}:{m2}" for n, r, mc, m4, m2 in transfers))
    print("entry_gates=q2_zero:36/253,q2_one_even:15/77,q4_one_v2_tail:8/77")
    print("localization_gates=q2_one_even:17/91,q4_one_v2_tail:8/91")
    print("q2_one_even_control=tail26_absorbed_r13_odd_3_unit")
    print("hostile_masks=" + ",".join(f"{masks}:{safe_gap}" for masks, safe_gap in hostiles))
    print(f"checks={CHECKS}")
    print("PASS")


if __name__ == "__main__":
    main()
