#!/usr/bin/env python3
"""Exact THM-2707 packet macros based at the frozen THM-2680 atom.

THM-2707 proves that 3,346 lift addresses share one open packet cylinder and
form a complete directed eleven-partite graph.  Its displayed eleven-cycle
meets the old following-atom support once.  This companion determines the
whole atom-bearing locus, verifies the full open cylinder (including the
delayed prefix), and counts the resulting two-step based macros.

The output is a support statement.  A physical lift residue is not thereby a
THM-2334 relation index or a THM-2350 target-dipole residue, and no target
action, endpoint current, row exclusion, or LRC(14) conclusion is asserted.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_full_physical_lift_fibre_thm2707 as fibre
import lrc14_mixed_slope_long_word_probe as old


m = old.m


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def strict_interval_index(value, starts, intervals):
    index = bisect_right(starts, value) - 1
    return index >= 0 and intervals[index][0] < value < intervals[index][1]


def main():
    p = 13
    R = p**6
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    radius = Fraction(1, 1304692766858936)
    interval = (x - radius, x + radius)
    require(R == m.R, "lift modulus changed")

    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    pair_prefixes = m.build_pair_prefixes(module)
    rows = m.shard((0, 1))[6][0]

    require(rails[0][:3] == (1, 0, 12), "following rail metadata changed")
    following = old.d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    require(
        {key: following[key]
         for key in ("future", "j", "h", "epsilon", "kappa")}
        == {"future": 1, "j": 2, "h": 6, "epsilon": 1, "kappa": 1},
        "following atom metadata changed",
    )

    rail_support = old.merge_support(rails[0][3])
    present_support = tuple(present[1, (-6) % p])
    following_support = old.merge_support(following["pieces"])
    rail_starts = tuple(left for left, _ in rail_support)
    present_starts = tuple(left for left, _ in present_support)
    following_starts = tuple(left for left, _ in following_support)

    good_residues = tuple(
        residue
        for residue in range(p)
        if (6 + residue) % p
        and m.is_unit(
            rows[0][0][(2 + 7 * residue) % p][1][6],
            (6 + residue) % p,
            26,
        )
    )
    require(
        good_residues == (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
        "active residue bank changed",
    )

    # Repeat THM-2707's independent integer-grid packet scan.
    denominator = (z * m.T).denominator
    point = (z * m.T).numerator
    modulus = m.T * denominator
    step = (7 * m.T // R) * denominator
    scaled_rail = tuple(
        (left * denominator, right * denominator)
        for left, right in rail_support
    )
    scaled_present = tuple(
        (left * denominator, right * denominator)
        for left, right in present_support
    )
    scaled_following = tuple(
        (left * denominator, right * denominator)
        for left, right in following_support
    )
    scaled_rail_starts = tuple(left for left, _ in scaled_rail)
    scaled_present_starts = tuple(left for left, _ in scaled_present)
    scaled_following_starts = tuple(left for left, _ in scaled_following)

    good = []
    midpoint_atom = []
    for n in range(R):
        if (
            n % p in good_residues
            and strict_interval_index(point, scaled_rail_starts, scaled_rail)
            and strict_interval_index(
                point, scaled_present_starts, scaled_present
            )
        ):
            good.append(n)
            if strict_interval_index(
                point, scaled_following_starts, scaled_following
            ):
                midpoint_atom.append(n)
        point = (point + step) % modulus

    good = tuple(good)
    midpoint_atom = tuple(midpoint_atom)
    require(len(good) == 3346, "packet census changed")
    require(
        dict(sorted(Counter(n % p for n in good).items()))
        == {
            0: 304,
            1: 305,
            2: 304,
            3: 305,
            4: 304,
            5: 305,
            6: 304,
            9: 301,
            10: 304,
            11: 305,
            12: 305,
        },
        "packet residue census changed",
    )

    # The delayed following-word coordinate is independent of the lift
    # address because R*(z+7n/R) differs from R*z by an integer.
    prefix = prefixes[following["future"]][following["h"]]
    delayed = tuple(
        (left, left + length) for left, length in zip(prefix[0], prefix[1])
    )
    delayed_center = old.frac(R * z) * m.T
    delayed_radius = 13 * R * m.T * radius
    require(
        old.open_arc_is_contained(
            delayed_center, delayed_radius, delayed, m.T
        ),
        "common cylinder escaped the following delayed prefix",
    )

    # A perturbation of initial x by e moves every q_n by 13e.  Test the
    # whole open initial cylinder against the base support of the atom.
    base_radius = 13 * m.T * radius
    whole_atom = []
    for n in good:
        q, _, _, unit = fibre.packet_address_data(n, z, rows)
        require(unit, "packet census retained a nonunit root")
        if old.open_arc_is_contained(
            q * m.T, base_radius, following_support, m.T
        ):
            whole_atom.append(n)
    whole_atom = tuple(whole_atom)

    residue_zero_packets = tuple(n for n in good if n % p == 0)
    require(
        midpoint_atom == whole_atom == residue_zero_packets,
        "following atom is not exactly the residue-zero packet part",
    )
    require(len(whole_atom) == 304, "atom endpoint bank changed")
    require(
        all(
            old.atom_margin(following, q, prefixes, 13) >= radius
            for q, _, _, _ in (
                fibre.packet_address_data(n, z, rows) for n in whole_atom
            )
        ),
        "an atom endpoint lost the whole open cylinder",
    )

    atom_set = set(whole_atom)
    transit = tuple(n for n in good if n not in atom_set)
    require(len(transit) == 3042, "transit packet bank changed")
    require(all(n % p for n in transit), "transit bank met residue zero")

    # Edges in THM-2707 exist exactly across distinct residue parts.  Hence
    # the atom/transit cut is complete in both directions, and every ordered
    # atom endpoint pair has every transit packet as a physical midpoint.
    one_way_edges = len(whole_atom) * len(transit)
    based_macro_paths = len(whole_atom) ** 2 * len(transit)
    require(one_way_edges == 924768, "atom/transit edge count changed")
    require(based_macro_paths == 281129472, "based macro count changed")

    a = 0
    b = 1
    require(a in atom_set and b in transit, "canonical based loop disappeared")
    signed_steps = (7 * (b - a), 7 * (a - b))
    positive_steps = tuple(step % R for step in signed_steps)
    require(
        signed_steps == (7, -7)
        and sum(signed_steps) == 0
        and all(step % p for step in positive_steps),
        "canonical physical two-step loop changed",
    )

    print("LRC14 FOLLOWING-ATOM-BASED PHYSICAL MACRO AUDIT")
    print(f"p={p} R={R} I=({interval[0]},{interval[1]}) length={2*radius}")
    print(
        "following_atom="
        f"rail={rails[0][:3]} future={following['future']} j={following['j']} "
        f"h={following['h']} epsilon={following['epsilon']} "
        f"kappa={following['kappa']} mass={following['value']}"
    )
    print(
        f"packet_addresses={len(good)} atom_endpoints={len(whole_atom)} "
        f"transit_packets={len(transit)}"
    )
    print(
        "atom_endpoint_residue_counts="
        f"{tuple(sorted(Counter(n % p for n in whole_atom).items()))}"
    )
    print(
        "atom_endpoint_characterization="
        "midpoint_atom=whole_open_I_atom=packet_residue_0"
    )
    print(
        f"atom_to_transit_edges={one_way_edges} "
        f"transit_to_atom_edges={one_way_edges} "
        f"two_step_atom_macros={based_macro_paths}"
    )
    print(
        f"canonical_loop={a}->{b}->{a} signed_lift_numerators={signed_steps} "
        f"positive_representatives={positive_steps} phase_sum=0"
    )
    print(
        "scope=whole-I fixed-atom support macro;"
        "no_target_action;no_relation-index-identification;"
        "no_semantic_endpoint_current"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
