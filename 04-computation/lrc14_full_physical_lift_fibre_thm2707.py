#!/usr/bin/env python3
"""Exact full THM-2657 lift-fibre scan above the THM-2694 terminal vertex.

The terminal private address is n0=110 in the slope-seven coordinate
q_n={13x}+7n/13^6.  A physical THM-2657 lift has numerator k prime to 13;
it sends n0 to n0+7^{-1}k modulo 13^6.  This script enumerates every such
lift, verifies the entire inherited open cylinder rather than its midpoint,
and emits one exact transcript row for every surviving lift.

The retained packet is exactly the fixed THM-2694 typed packet: literal rail,
present factor, delayed Boolean prefix, predecessor carry, future/private
half-tooth, and primitive coefficient-unit flag.  It does not include the
whole following THM-2680 atom and does not manufacture a THM-2365 target
action or semantic endpoint current.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_mixed_slope_long_word_probe as old


m = old.m


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def strict_interval_index(value, starts, intervals):
    index = bisect_right(starts, value) - 1
    return index >= 0 and intervals[index][0] < value < intervals[index][1]


def packet_address_data(n, z, rows):
    q = frac(z + Fraction(7 * n, m.R))
    carry = (2 + 7 * n) % m.P
    root = (6 + n) % m.P
    unit = bool(root) and m.is_unit(rows[0][0][carry][1][6], root, 26)
    return q, carry, root, unit


def main():
    p = 13
    R = p**6
    S = p**5
    inv7 = pow(7, -1, R)
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    I = (
        Fraction(960117507257, 1930018885886),
        Fraction(324519717452867, 652346383429468),
    )
    radius = Fraction(1, 1304692766858936)
    n0 = 110
    require(R == m.R and inv7 == 4137265, "base lift constants changed")
    require(I == (x - radius, x + radius), "THM-2694 cylinder changed")
    require(I[1] - I[0] == 2 * radius, "cylinder length changed")
    require(frac(13 * x) == z, "frozen D endpoint changed")

    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    pair_prefixes = m.build_pair_prefixes(module)
    rows = m.shard((0, 1))[6][0]

    # Fixed THM-2694 skeleton (rail,sector,edge,kappa,h,shallow).
    skeleton = (0, 0, 0, 1, 6, 1)
    rail_support = old.merge_support(rails[0][3])
    present_support = tuple(present[1, (-6) % p])
    rail_starts = tuple(left for left, _ in rail_support)
    present_starts = tuple(left for left, _ in present_support)

    good_residues = []
    for residue in range(p):
        carry = (2 + 7 * residue) % p
        root = (6 + residue) % p
        if root and m.is_unit(rows[0][0][carry][1][6], root, 26):
            good_residues.append(residue)
    good_residues = tuple(good_residues)
    require(good_residues == (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
            "private root/unit residue bank changed")

    # Exact integer midpoint scan, inherited independently from THM-2694.
    denominator = (z * m.T).denominator
    base_numerator = (z * m.T).numerator
    modulus = m.T * denominator
    step = (7 * m.T // R) * denominator
    require(m.T % R == 0 and step * R == 7 * m.T * denominator,
            "orbit translation grid changed")
    scaled_rail = tuple((left * denominator, right * denominator)
                        for left, right in rail_support)
    scaled_present = tuple((left * denominator, right * denominator)
                           for left, right in present_support)
    scaled_rail_starts = tuple(left for left, _ in scaled_rail)
    scaled_present_starts = tuple(left for left, _ in scaled_present)

    good = []
    candidate_count = 0
    point = base_numerator
    for n in range(R):
        # The affine map k=7(n-n0) is a permutation of Z/RZ.  Hence this
        # single complete address scan simultaneously enumerates every
        # THM-2657 numerator k: the nonzero quotient fibres are exactly the
        # addresses with n mod13 different from n0 mod13.
        if n % p != n0 % p:
            k = (7 * (n - n0)) % R
            require(k % p != 0
                    and (n0 + inv7 * k) % R == n,
                    "address/lift affine bijection changed")
            candidate_count += 1
        if (n % p in good_residues
                and strict_interval_index(
                    point, scaled_rail_starts, scaled_rail)
                and strict_interval_index(
                    point, scaled_present_starts, scaled_present)):
            good.append(n)
        point = (point + step) % modulus
    good = tuple(good)
    good_set = set(good)
    require(len(good) == 3346 and n0 in good_set,
            "complete midpoint packet census changed")

    residue_counts = Counter(n % p for n in good)
    expected_counts = {
        0: 304, 1: 305, 2: 304, 3: 305, 4: 304, 5: 305,
        6: 304, 9: 301, 10: 304, 11: 305, 12: 305,
    }
    require(dict(sorted(residue_counts.items())) == expected_counts,
            "good-address residue census changed")

    # Whole-cylinder verification.  The delayed coordinate is independent of
    # n; the carry fractional part is also independent of n.  Rail, present,
    # private root, and unit are checked at every one of the 3346 addresses.
    y = frac(R * z) * m.T
    delayed_fraction = frac(R * z)
    delayed_left = delayed_fraction - 13 * R * radius
    delayed_right = delayed_fraction + 13 * R * radius
    require(0 < delayed_left < delayed_right < 1,
            "common cylinder crossed a delayed integer boundary")
    require(int(13 * delayed_left) == int(13 * delayed_fraction)
            == int(13 * delayed_right) == 6,
            "future digit changed across the common cylinder")
    require(int(26 * delayed_left) - 12
            == int(26 * delayed_fraction) - 12
            == int(26 * delayed_right) - 12 == 1,
            "future half-digit changed across the common cylinder")
    prefix = pair_prefixes[0][1][6][1]
    delayed = tuple((left, left + length)
                    for left, length in zip(prefix[0], prefix[1]))
    delayed_radius = 13 * R * m.T * radius
    require(old.open_arc_is_contained(
        y, delayed_radius, delayed, m.T
    ), "common cylinder escaped the delayed Boolean prefix")

    carry_fraction = frac(R * z)
    carry_radius = 13 * R * radius
    require(carry_radius <= min(carry_fraction, 1 - carry_fraction),
            "common cylinder crossed a predecessor-carry boundary")

    base_radius = 13 * m.T * radius
    deep_radius = 13 * module.C3 * 182 * radius
    rail_lower_slacks = []
    rail_upper_slacks = []
    present_lower_slacks = []
    present_upper_slacks = []
    full_I_targets = 0
    for n in good:
        q, carry, root, unit = packet_address_data(n, z, rows)
        require(unit, "midpoint census retained a non-unit private root")
        base = q * m.T
        require(old.open_arc_is_contained(
            base, base_radius, rail_support, m.T
        ), "common cylinder escaped a literal rail")
        require(old.open_arc_is_contained(
            base, base_radius, present_support, m.T
        ), "common cylinder escaped a present factor")
        deep = frac(module.C3 * q) * 182
        deep_support = ((Fraction(14 * root - 13), Fraction(14 * root)),)
        require(old.open_arc_is_contained(
            deep, deep_radius, deep_support, 182
        ), "common cylinder escaped a private half-tooth")
        require(carry == (2 + 7 * n) % p and root == (6 + n) % p,
                "private digit covariance changed")

        # Record the closest exact rail/present walls around the midpoint.
        ri = bisect_right(rail_starts, base) - 1
        pi = bisect_right(present_starts, base) - 1
        require(ri >= 0 and pi >= 0, "containing support index vanished")
        rail_lower_slacks.append(base - rail_support[ri][0] - base_radius)
        rail_upper_slacks.append(rail_support[ri][1] - base - base_radius)
        present_lower_slacks.append(
            base - present_support[pi][0] - base_radius)
        present_upper_slacks.append(
            present_support[pi][1] - base - base_radius)
        full_I_targets += 1
    require(full_I_targets == len(good), "not every good midpoint retained I")
    require(min(rail_lower_slacks) >= 0 and min(rail_upper_slacks) >= 0
            and min(present_lower_slacks) >= 0
            and min(present_upper_slacks) >= 0,
            "a rail/present whole-cylinder slack became negative")

    # Enumerate every nonzero THM-2657 quotient lift from n0.  Because
    # n'=n0+7^{-1}k, k prime to 13 is equivalent to n' mod13 != n0 mod13.
    survivor_rows = []
    delta_counts = Counter()
    for target in good:
        if target % p == n0 % p:
            continue
        k = (7 * (target - n0)) % R
        delta = (target - n0) % p
        require((n0 + inv7 * k) % R == target
                and delta == 2 * k % p
                and k % p == 7 * delta % p,
                "THM-2657 root/carry congruence changed")
        q, carry, root, unit = packet_address_data(target, z, rows)
        require(unit, "a surviving lift lost the private unit")
        next_count = len(good) - residue_counts[target % p]
        require(next_count > 0, "a surviving lift became terminal")
        survivor_rows.append(
            (k, delta, target, target % p, carry, root, next_count)
        )
        delta_counts[delta] += 1
    survivor_rows.sort()
    require(candidate_count == 12 * S,
            "full lift fibre did not biject onto the twelve residue cosets")
    require(len(survivor_rows) == 3042,
            "full physical lift survivor count changed")
    require(dict(sorted(delta_counts.items())) == {
        3: 301, 4: 304, 5: 305, 6: 305, 7: 304,
        8: 305, 9: 304, 10: 305, 11: 304, 12: 305,
    }, "survivor delta census changed")

    # The twelve least positive representatives are precisely the old failed
    # n=111,...,122 tests.  This is the sharp boundary repaired here.
    canonical_failures = []
    for delta in range(1, p):
        k = 7 * delta
        target = (n0 + inv7 * k) % R
        require(target == n0 + delta, "least representative changed")
        q = frac(z + Fraction(7 * target, R))
        base = q * m.T
        rail_meets = old.open_arc_meets_support(
            base, base_radius, rail_support, m.T)
        present_meets = old.open_arc_meets_support(
            base, base_radius, present_support, m.T)
        canonical_failures.append(
            (delta, k, target, rail_meets, present_meets)
        )
    canonical_failures = tuple(canonical_failures)
    require(all(row[3] and not row[4] for row in canonical_failures),
            "a canonical THM-2694 next lift unexpectedly survived")

    # The full lift graph has one vertex for every packet address and an edge
    # exactly when the source and target residues differ.  It is therefore a
    # complete directed 11-partite graph, one SCC, with diameter at most two.
    directed_edges = len(good) ** 2 - sum(
        count**2 for count in residue_counts.values()
    )
    require(directed_edges == 10177910,
            "complete multipartite edge count changed")
    require(len(residue_counts) == 11 and min(residue_counts.values()) > 0,
            "an active residue part vanished")
    require(set(row[-1] for row in survivor_rows) == {3041, 3042, 3045},
            "continuation outdegree census changed")

    def address_from_cumulative_lift(K):
        return (n0 + inv7 * (K % R)) % R

    def verify_cycle(cumulative, steps):
        require(len(cumulative) == len(steps) + 1
                and cumulative[0] == cumulative[-1] == 0,
                "cycle endpoint data changed")
        require(tuple(b - a for a, b in zip(
            cumulative, cumulative[1:])) == steps,
            "cycle step spelling changed")
        nodes = tuple(address_from_cumulative_lift(K) for K in cumulative)
        require(nodes[0] == nodes[-1] == n0,
                "cycle failed to close physically")
        for K, node in zip(cumulative, nodes):
            require(node in good_set, "cycle left the common packet simplex")
            require((node - n0) % p == 2 * (K % p) % p,
                    "cycle root increment changed")
        for step_value in steps:
            require(step_value % p != 0,
                    "cycle used a quotient-zero lift")
        return nodes

    triangle_K = (0, -1, -2, 0)
    triangle_steps = (-1, -1, 2)
    triangle_nodes = verify_cycle(triangle_K, triangle_steps)
    require(tuple(n % p for n in triangle_nodes) == (6, 4, 2, 6),
            "nonbacktracking triangle residues changed")

    cycle11_K = (0, 2, 4, 5, 8, 10, 11, 12, 9, 6, 3, 0)
    cycle11_steps = (2, 2, 1, 3, 2, 1, 1, -3, -3, -3, -3)
    cycle11_nodes = verify_cycle(cycle11_K, cycle11_steps)
    cycle11_data = tuple(
        (K, node, node % p, (2 + 7 * node) % p, (6 + node) % p)
        for K, node in zip(cycle11_K[:-1], cycle11_nodes[:-1])
    )
    require(cycle11_data == (
        (0, 110, 6, 5, 12),
        (2, 3447831, 10, 7, 3),
        (4, 2068743, 1, 9, 7),
        (5, 1379199, 3, 10, 9),
        (8, 4137376, 9, 0, 2),
        (10, 2758288, 0, 2, 6),
        (11, 2068744, 2, 3, 8),
        (12, 1379200, 4, 4, 10),
        (9, 3447832, 11, 1, 4),
        (6, 689655, 5, 11, 11),
        (3, 2758287, 12, 8, 5),
    ), "eleven-residue cycle data changed")
    require(set(row[2] for row in cycle11_data) == set(good_residues),
            "eleven-cycle stopped visiting every active residue")

    # Hostile scope control: the old weighted following atom is not one of the
    # retained packet factors.  It is already absent at the terminal source
    # and at ten of the eleven distinct nodes in the displayed cycle.
    following = old.d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    following_support = old.merge_support(following["pieces"])
    following_cycle_hits = tuple(
        any(left < packet_address_data(node, z, rows)[0] * m.T < right
            for left, right in following_support)
        for node in cycle11_nodes[:-1]
    )
    require(following_cycle_hits == (
        False, False, False, False, False, True,
        False, False, False, False, False,
    ), "frozen following-atom hostile changed")

    print("LRC14 FULL PHYSICAL LIFT-FIBRE / COMMON-SIMPLEX AUDIT")
    print(f"p={p} R={R} S={S} inv7={inv7}")
    print(f"source_address={n0} source_residue={n0 % p}")
    print(f"fixed_skeleton={skeleton}")
    print(f"common_I=({I[0]},{I[1]}) length={I[1]-I[0]}")
    print(f"active_residues={good_residues}")
    print(f"good_addresses={len(good)} residue_counts={tuple(sorted(residue_counts.items()))}")
    print(
        "whole_I_targets="
        f"{full_I_targets} min_wall_slacks="
        f"(rail_left={min(rail_lower_slacks)},"
        f"rail_right={min(rail_upper_slacks)},"
        f"present_left={min(present_lower_slacks)},"
        f"present_right={min(present_upper_slacks)})"
    )
    print(f"all_nonzero_quotient_lifts={candidate_count}")
    print(f"terminal_survivors={len(survivor_rows)} delta_counts={tuple(sorted(delta_counts.items()))}")
    print(f"canonical_least_lift_failures={canonical_failures}")
    print(
        "restricted_packet_nerve="
        "Delta^3345;terminal_nonzero_lift_cone=Delta^3042"
    )
    print(
        "lift_graph=complete_directed_11_partite "
        f"vertices={len(good)} edges={directed_edges} scc=1 diameter_le_2=True"
    )
    print(f"triangle_K={triangle_K} steps={triangle_steps} nodes={triangle_nodes}")
    print(
        f"cycle11_K={cycle11_K} steps={cycle11_steps} "
        f"data={cycle11_data} repeatable=True"
    )
    print(
        "following_atom_excluded=True "
        f"cycle_midpoint_hits={following_cycle_hits}"
    )
    print("SURVIVOR_ROWS_BEGIN")
    for k, delta, target, residue, carry, root, next_count in survivor_rows:
        print(
            f"k={k} delta={delta} target={target} residue={residue} "
            f"carry={carry} root={root} common_I=({I[0]},{I[1]}) "
            f"next_full_lifts={next_count}"
        )
    print("SURVIVOR_ROWS_END")
    print("scope=fixed_skeleton_support_only;no_target_action;no_semantic_endpoint")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
