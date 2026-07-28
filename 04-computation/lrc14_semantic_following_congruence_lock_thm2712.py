#!/usr/bin/env python3
"""Exact semantic-following refinement of the THM-2707 lift graph.

The fixed THM-2680 following atom forces one outer base-13 address residue.
This kills every nonzero-root THM-2657 edge inside the semantic locus, while
leaving a complete root-neutral physical graph one odometer layer lower.
All interval tests use exact rational arithmetic on the full inherited open
cylinder, not only midpoint samples.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_full_physical_lift_fibre_thm2707 as lift


m = lift.m
old = lift.old


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def strict_contains(value, starts, intervals):
    index = bisect_right(starts, value) - 1
    return index >= 0 and intervals[index][0] < value < intervals[index][1]


def containing_index(value, intervals):
    starts = tuple(left for left, _ in intervals)
    index = bisect_right(starts, value) - 1
    require(index >= 0 and intervals[index][0] < value < intervals[index][1],
            "strict point has no containing component")
    return index


def wall_slacks(value, radius, intervals):
    index = containing_index(value, intervals)
    left, right = intervals[index]
    left_slack = value - left - radius
    right_slack = right - value - radius
    require(left_slack >= 0 and right_slack >= 0,
            "open cylinder escaped a claimed support component")
    return left_slack, right_slack, index


def prefix_intervals(prefix):
    starts, lengths, _ = prefix
    return tuple((left, left + length)
                 for left, length in zip(starts, lengths))


def atom_whole_I(atom, base, base_radius, delayed, delayed_radius,
                 prefixes):
    support = old.merge_support(atom["pieces"])
    if not support:
        return False
    delayed_support = prefix_intervals(prefixes[atom["future"]][atom["h"]])
    return (old.open_arc_is_contained(base, base_radius, support, m.T)
            and old.open_arc_is_contained(
                delayed, delayed_radius, delayed_support, m.T))


def valuation_13(value):
    require(value != 0, "valuation requested at zero")
    value = abs(value)
    exponent = 0
    while value % 13 == 0:
        value //= 13
        exponent += 1
    return exponent


def main():
    p = 13
    R = p**6
    S = p**5
    n0 = 110
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    I = (
        Fraction(960117507257, 1930018885886),
        Fraction(324519717452867, 652346383429468),
    )
    radius = Fraction(1, 1304692766858936)
    require(I == (x - radius, x + radius), "inherited I changed")

    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    rows = m.shard((0, 1))[6][0]
    following = old.d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    current = old.d.build_atom(
        module, prefixes, present, starts, rails[2], 3, 2, 0, 0
    )
    require(tuple(following[key] for key in (
        "future", "j", "h", "epsilon", "kappa", "value"
    )) == (1, 2, 6, 1, 1, 867661831383617727737280),
            "frozen following atom changed")
    require(tuple(current[key] for key in (
        "future", "j", "h", "epsilon", "kappa"
    )) == (3, 5, 2, 0, 0), "frozen current atom changed")

    following_root = (-following["h"] - 1) % p
    require(following["j"] == 2 and following_root == 6,
            "frozen semantic labels changed")

    rail_support = old.merge_support(rails[0][3])
    present_support = tuple(present[1, (-6) % p])
    rail_starts = tuple(left for left, _ in rail_support)
    present_starts = tuple(left for left, _ in present_support)
    following_support = old.merge_support(following["pieces"])
    following_starts = tuple(left for left, _ in following_support)
    require(len(following["pieces"]) == len(following_support) == 304,
            "frozen following component count changed")

    denominator = (z * m.T).denominator
    modulus = m.T * denominator
    base_numerator = (z * m.T).numerator
    orbit_step = (7 * m.T // R) * denominator
    step13 = 13 * orbit_step
    require(orbit_step * R == 7 * m.T * denominator,
            "address scan grid changed")
    scaled_rail = tuple((left * denominator, right * denominator)
                        for left, right in rail_support)
    scaled_present = tuple((left * denominator, right * denominator)
                           for left, right in present_support)
    scaled_following = tuple((left * denominator, right * denominator)
                             for left, right in following_support)
    scaled_rail_starts = tuple(left for left, _ in scaled_rail)
    scaled_present_starts = tuple(left for left, _ in scaled_present)
    scaled_following_starts = tuple(left for left, _ in scaled_following)

    # The frozen predecessor and root labels force n=0 modulo 13:
    # 2+7n=2 and 6+n=6.  Therefore scanning this one coset is exhaustive for
    # semantic candidates, while THM-2707 supplies the full packet universe.
    require(all(
        ((2 + 7 * residue) % p == following["j"]
         and (6 + residue) % p == following_root)
        == (residue == 0)
        for residue in range(p)
    ), "semantic label gate is not exactly one residue")

    packet_nodes = []
    following_mid_nodes = []
    point = base_numerator
    for j in range(S):
        n = 13 * j
        packet_hit = (
            strict_contains(point, scaled_rail_starts, scaled_rail)
            and strict_contains(point, scaled_present_starts, scaled_present)
            and m.is_unit(rows[0][0][2][1][6], 6, 26)
        )
        following_hit = strict_contains(
            point, scaled_following_starts, scaled_following)
        if packet_hit:
            packet_nodes.append(n)
        if following_hit:
            following_mid_nodes.append(n)
        point = (point + step13) % modulus
    packet_nodes = tuple(packet_nodes)
    following_mid_nodes = tuple(following_mid_nodes)
    require(len(packet_nodes) == 304
            and following_mid_nodes == packet_nodes,
            "semantic midpoint locus changed")

    base_radius = 13 * m.T * radius
    delayed = frac(R * z) * m.T
    delayed_radius = 13 * R * m.T * radius
    following_prefix = prefix_intervals(
        prefixes[following["future"]][following["h"]])
    dleft, dright, dindex = wall_slacks(
        delayed, delayed_radius, following_prefix)
    require((dleft, dright) == (0, 0)
            and following_prefix[dindex]
            == (160364030987160, 160364059635780),
            "sharp delayed pullback changed")

    whole_nodes = []
    component_indices = []
    weights = Counter()
    left_slacks = []
    right_slacks = []
    for n in packet_nodes:
        q = frac(z + Fraction(7 * n, R))
        base = q * m.T
        require(atom_whole_I(
            following, base, base_radius,
            delayed, delayed_radius, prefixes
        ), "semantic midpoint did not retain the whole cylinder")
        left_slack, right_slack, component_index = wall_slacks(
            base, base_radius, following_support)
        left_slacks.append(left_slack)
        right_slacks.append(right_slack)
        component_indices.append(component_index)
        weights[old.strict_weighted_hit(base, following["pieces"])[2]] += 1
        whole_nodes.append(n)
    whole_nodes = tuple(whole_nodes)
    require(whole_nodes == packet_nodes and len(whole_nodes) == 304,
            "whole-I semantic locus changed")
    require(tuple(sorted(component_indices)) == tuple(range(304)),
            "semantic nodes stopped occupying following components bijectively")
    require(weights == Counter({27582102210: 266, 27581135604: 38}),
            "semantic weight census changed")
    require(min(left_slacks) == Fraction(368024580, 2197)
            and min(right_slacks)
            == Fraction(10574833707900, 371293),
            "semantic base-support slack changed")

    # The fixed current atom contains the same inherited cylinder, so these
    # are genuine current/following support overlaps rather than endpoint-only
    # following witnesses.
    current_support = old.merge_support(current["pieces"])
    current_base = x * m.T
    current_base_radius = m.T * radius
    current_delayed = frac(R * x) * m.T
    current_delayed_radius = R * m.T * radius
    require(atom_whole_I(
        current, current_base, current_base_radius,
        current_delayed, current_delayed_radius, prefixes
    ), "frozen current stopped containing I")
    cbl, cbr, _ = wall_slacks(
        current_base, current_base_radius, current_support)
    cpl, cpr, _ = wall_slacks(
        current_delayed, current_delayed_radius,
        prefix_intervals(prefixes[current["future"]][current["h"]]))
    require((cbl, cbr, cpl, cpr) == (
        Fraction(344122815960, 28561),
        Fraction(80124658752600, 4826809),
        28648620,
        312931080,
    ), "frozen current whole-I slack changed")

    # Outer nonzero-root dynamics.  Every semantic node has one residue, so
    # the induced THM-2657 delta!=0 graph is edgeless.  The nonsemantic source
    # n0 has 304 edges into and 304 edges out of the locus, but those cycles
    # lose the frozen following atom at n0.
    outer_edge_count = sum(
        1 for a in whole_nodes for b in whole_nodes
        if a % p != b % p
    )
    require(outer_edge_count == 0,
            "semantic outer lift graph acquired an edge")
    terminal_rows = []
    for n in whole_nodes:
        k = 7 * (n - n0) % R
        delta = (n - n0) % p
        require(delta == 7 and k % p == 10,
                "terminal semantic lift congruence changed")
        terminal_rows.append((k, n))
    terminal_rows.sort()
    require(terminal_rows[0] == (10, 2758288),
            "shortest positive terminal semantic lift changed")
    signed_rows = sorted(
        ((k if k <= R // 2 else k - R), n)
        for k, n in terminal_rows
    )
    require(min(signed_rows, key=lambda row: abs(row[0])) == (-3, 2068742),
            "shortest signed terminal semantic lift changed")

    # Root-neutral inner dynamics.  Write n=13j and divide the forced lift
    # numerator k=91(j'-j) by 13.  The resulting gain 7(j'-j) is the same
    # exact address coboundary one level lower.
    inner_counts = Counter((n // 13) % p for n in whole_nodes)
    expected_inner_counts = (
        (0, 23), (1, 27), (2, 22), (3, 27), (4, 22),
        (5, 23), (6, 22), (7, 23), (8, 23), (9, 23),
        (10, 23), (11, 23), (12, 23),
    )
    require(tuple(sorted(inner_counts.items())) == expected_inner_counts,
            "inner residue part sizes changed")
    inner_edges = len(whole_nodes)**2 - sum(
        count**2 for count in inner_counts.values())
    require(inner_edges == 85274,
            "exact-valuation-one inner edge count changed")

    valuation_census = Counter()
    for index, a in enumerate(whole_nodes):
        for b in whole_nodes[index + 1:]:
            exponent = valuation_13(7 * (b - a))
            valuation_census[exponent] += 1
            ja, jb = a // 13, b // 13
            k_over_13 = 7 * (jb - ja)
            require(k_over_13 == 7 * jb - 7 * ja,
                    "inner address coboundary changed")
    require(valuation_census == Counter({
        1: 42637, 2: 3222, 3: 159, 4: 38,
    }), "semantic pair valuation census changed")
    require(sum(valuation_census.values()) == 304 * 303 // 2
            and 2 * valuation_census[1] == inner_edges,
            "semantic edge/pair invoice changed")

    inner_cycle_nodes = (0, 13, 26, 0)
    inner_cycle_steps = (91, 91, -182)
    require(all(n in set(whole_nodes) for n in inner_cycle_nodes)
            and sum(inner_cycle_steps) == 0
            and tuple(7 * (b - a) for a, b in zip(
                inner_cycle_nodes, inner_cycle_nodes[1:]
            )) == inner_cycle_steps,
            "inner semantic triangle changed")

    # Upstairs the full address gain is already a coboundary.  It cannot
    # descend to the outer residue quotient because moving within one residue
    # fibre by 13 changes the potential by 91, whose divided kernel class is
    # 7 modulo 13: exactly the THM-2657 obstruction.
    require((7 * (n0 + 13) - 7 * n0) % R == 91
            and 2 * 91 % p == 0
            and (91 // p) % p == 7,
            "address-coboundary descent obstruction changed")

    # Same-current-compatible sibling atoms on this rail/future clock.
    siblings = []
    for h in range(1, 12):
        for epsilon in (0, 1):
            for kappa in (0, 1):
                atom = old.d.build_atom(
                    module, prefixes, present, starts,
                    rails[0], 1, h, epsilon, kappa
                )
                if atom["j"] != current["h"]:
                    continue
                root = (-h - 1) % p
                count = 0
                for n in packet_nodes:
                    q = frac(z + Fraction(7 * n, R))
                    if atom_whole_I(
                        atom, q * m.T, base_radius,
                        delayed, delayed_radius, prefixes
                    ):
                        count += 1
                siblings.append((h, epsilon, kappa, root,
                                 atom["value"], count))
    siblings = tuple(siblings)
    require(tuple(row[:4] for row in siblings) == (
        (6, 1, 1, 6), (7, 0, 1, 5),
        (7, 1, 0, 5), (8, 0, 0, 4),
    ) and tuple(row[-1] for row in siblings) == (304, 0, 0, 0),
            "same-current semantic sibling census changed")

    print("LRC14 SEMANTIC FOLLOWING CONGRUENCE-LOCK AUDIT")
    print(f"p={p} R={R} S={S} source_address={n0}")
    print(f"common_I=({I[0]},{I[1]}) length={I[1]-I[0]}")
    print(
        "following_labels="
        f"(future={following['future']},j={following['j']},h={following['h']},"
        f"epsilon={following['epsilon']},kappa={following['kappa']},"
        f"root={following_root})"
    )
    print(
        f"outer_label_compatible={S} outer_label_incompatible={11*S} "
        f"same_congruence_support_failures={S-len(whole_nodes)}"
    )
    print(
        f"packet_survivors=3042 semantic_whole_I={len(whole_nodes)} "
        f"skeleton_only={3042-len(whole_nodes)} residue=0"
    )
    print(
        f"following_components={len(following_support)} occupancy=one_each "
        f"weights={tuple(sorted(weights.items()))}"
    )
    print(
        "following_base_min_slacks="
        f"({min(left_slacks)},{min(right_slacks)}) "
        f"delayed_slacks=({dleft},{dright})"
    )
    print(
        "current_slacks="
        f"(base_left={cbl},base_right={cbr},"
        f"prefix_left={cpl},prefix_right={cpr})"
    )
    print(
        "outer_nonzero_root_graph="
        f"vertices={len(whole_nodes)} edges={outer_edge_count} "
        "sccs=304 cycle=False"
    )
    print(
        "terminal_packet_star="
        f"bidirected_edges={2*len(whole_nodes)} "
        f"shortest_positive={terminal_rows[0]} shortest_signed=(-3,2068742) "
        "semantic_at_terminal=False"
    )
    print(
        "inner_root_neutral_graph="
        f"vertices={len(whole_nodes)} all_nonzero_edges={304*303} "
        f"valuation1_edges={inner_edges} parts={expected_inner_counts} "
        "scc=1 diameter_le_2=True"
    )
    print(f"unordered_pair_valuation_census={tuple(sorted(valuation_census.items()))}")
    print(f"inner_triangle_nodes={inner_cycle_nodes} steps={inner_cycle_steps}")
    print("address_gain=k_ab=7(b-a)=phi(b)-phi(a);outer_descent_class=7")
    print(f"same_current_siblings={siblings}")
    print("SEMANTIC_ROWS_BEGIN")
    component_by_node = dict(zip(whole_nodes, component_indices))
    weight_by_node = {}
    for n in whole_nodes:
        q = frac(z + Fraction(7 * n, R))
        weight_by_node[n] = old.strict_weighted_hit(
            q * m.T, following["pieces"])[2]
    for n in whole_nodes:
        k = 7 * (n - n0) % R
        signed_k = k if k <= R // 2 else k - R
        print(
            f"n={n} j={n//13} inner_residue={(n//13)%13} "
            f"terminal_k={k} signed_k={signed_k} delta=7 "
            f"component={component_by_node[n]} weight={weight_by_node[n]} "
            f"common_I=({I[0]},{I[1]})"
        )
    print("SEMANTIC_ROWS_END")
    print("scope=frozen_current/rail/future-clock;outer_semantic_current_open")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
