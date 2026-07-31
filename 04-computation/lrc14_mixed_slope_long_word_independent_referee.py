#!/usr/bin/env python3
"""Independent exact referee for the THM-2694 mixed long word.

Unlike the primary complete-orbit scanner, this companion reconstructs the
initial D component from local affine inequalities and then verifies the
claimed finite lift directly on that entire open interval.  It does not
import the primary probe or use its gap/component routines.
"""

from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_predecessor_carry_private_root_atlas_thm2640 as m
import lrc14_dilation_reversed_clock_fibre_product_probe as d


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def merge_weighted(pieces):
    raw = sorted((left, right) for left, right, weight in pieces if weight)
    out = []
    for left, right in raw:
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return tuple(out)


def containing_interval(value, intervals):
    starts = tuple(left for left, _ in intervals)
    index = bisect_right(starts, value) - 1
    require(index >= 0 and intervals[index][0] < value < intervals[index][1],
            "strict point left its claimed interval")
    return intervals[index]


def perturbation_bounds(value, intervals, multiplier):
    """Bounds on initial perturbation e when value changes by multiplier*e."""
    scaled = frac(value) * m.T
    left, right = containing_interval(scaled, intervals)
    return ((left - scaled) / (multiplier * m.T),
            (right - scaled) / (multiplier * m.T))


def interval_contained(left, right, support):
    return any(a <= left and right <= b for a, b in support)


def interval_disjoint(left, right, support):
    return all(min(right, b) <= max(left, a) for a, b in support)


def weighted_hit(value, pieces):
    hits = tuple(piece for piece in pieces if piece[0] < value < piece[1])
    require(len(hits) == 1, "expected one strict weighted hit")
    return hits[0]


def main():
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    claimed = (
        Fraction(960117507257, 1930018885886),
        Fraction(324519717452867, 652346383429468),
    )
    require(frac(13 * x) == z, "D endpoint changed")

    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    pair_prefixes = m.build_pair_prefixes(module)
    current = d.build_atom(
        module, prefixes, present, starts, rails[2], 3, 2, 0, 0
    )
    following = d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )

    # Rebuild the local D-fibre component from four independent affine
    # restrictions: current/following base atoms and their delayed words.
    current_support = merge_weighted(current["pieces"])
    following_support = merge_weighted(following["pieces"])
    current_delayed = tuple(
        (left, left + length)
        for left, length in zip(prefixes[3][2][0], prefixes[3][2][1])
    )
    following_delayed = tuple(
        (left, left + length)
        for left, length in zip(prefixes[1][6][0], prefixes[1][6][1])
    )
    bounds = (
        perturbation_bounds(x, current_support, 1),
        perturbation_bounds(z, following_support, 13),
        perturbation_bounds(m.R * x, current_delayed, m.R),
        perturbation_bounds(m.R * z, following_delayed, 13 * m.R),
    )
    lower = max(row[0] for row in bounds)
    upper = min(row[1] for row in bounds)
    rebuilt = (x + lower, x + upper)
    require(rebuilt == claimed and lower == -upper,
            "independent D-component reconstruction changed")
    radius = upper

    # Rebuild the fixed THM-2640 packet and its allowed residue line.
    rail_index, sector, edge, kappa, h, shallow = 0, 0, 0, 1, 6, 1
    carry0, root0 = 2, 6
    rows = m.shard((0, 1))[6][0]
    allowed = tuple(
        residue for residue in range(13)
        if (root0 + residue) % 13
        and m.is_unit(
            rows[sector][edge][(carry0 + 7 * residue) % 13][kappa][h],
            (root0 + residue) % 13,
            26,
        )
    )
    require(allowed == (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
            "private unit residue bank changed")
    path = tuple(
        13 * height + residue
        for height in range(8)
        for residue in allowed
    ) + tuple(range(104, 111))
    require(len(path) == 95 and path[0] == 0 and path[-1] == 110,
            "direct lifted path census changed")
    steps = tuple(right - left for left, right in zip(path, path[1:]))
    require(set(steps) == {1, 3}
            and steps.count(1) == 86 and steps.count(3) == 8,
            "direct lifted slope word changed")
    require(len(allowed) == 11 and len(allowed) % 2 == 1,
            "residue quotient stopped being an odd cycle")
    path_labels = tuple(
        index // 2 if index % 2 == 0 else 94 - index // 2
        for index in range(95)
    )
    cycle_labels = (0, 4, 6, 5, 8, 3, 9, 2, 10, 1, 11)
    require(tuple(abs(a - b) for a, b in zip(
        path_labels, path_labels[1:]
    )) == tuple(range(94, 0, -1)), "P95 graceful control failed")
    require(tuple(abs(a - b) for a, b in zip(
        cycle_labels, cycle_labels[1:] + cycle_labels[:1]
    )) == (4, 2, 1, 3, 5, 6, 7, 8, 9, 10, 11),
            "C11 graceful control failed")

    rail_support = merge_weighted(rails[rail_index][3])
    present_support = tuple(present[shallow, (-h) % 13])
    prefix = pair_prefixes[sector][shallow][h][kappa]
    delayed_support = tuple(
        (left, left + length)
        for left, length in zip(prefix[0], prefix[1])
    )
    e_left, e_right = claimed[0] - x, claimed[1] - x

    metadata = {}
    for n in path:
        q = frac(z + Fraction(7 * n, m.R))
        q_left = q + 13 * e_left
        q_right = q + 13 * e_right
        require(0 < q_left < q_right < 1, "local q interval wrapped")
        require(interval_contained(
            q_left * m.T, q_right * m.T, rail_support
        ), "lifted cylinder escaped the literal rail")
        require(interval_contained(
            q_left * m.T, q_right * m.T, present_support
        ), "lifted cylinder escaped the present factor")

        delayed_center = frac(m.R * q)
        delayed_left = delayed_center + 13 * m.R * e_left
        delayed_right = delayed_center + 13 * m.R * e_right
        require(0 < delayed_left < delayed_right < 1,
                "delayed cylinder crossed an integer boundary")
        require(interval_contained(
            delayed_left * m.T, delayed_right * m.T, delayed_support
        ), "lifted cylinder escaped the delayed Boolean word")
        require(int(13 * delayed_left) == int(13 * delayed_center)
                == int(13 * delayed_right) == h,
                "future digit changed across the cylinder")
        require(int(26 * delayed_left) - 2 * h
                == int(26 * delayed_center) - 2 * h
                == int(26 * delayed_right) - 2 * h == kappa,
                "future half-digit changed across the cylinder")

        carry = (carry0 + 7 * n) % 13
        root = (root0 + n) % 13
        raw_center = m.R * q
        carry_fraction = frac(raw_center)
        require(0 < carry_fraction + 13 * m.R * e_left
                < carry_fraction + 13 * m.R * e_right < 1,
                "predecessor carry changed across the cylinder")
        require((raw_center.numerator // raw_center.denominator) % 13 == carry,
                "predecessor carry covariance changed")
        require(root and m.is_unit(
            rows[sector][edge][carry][kappa][h], root, 26
        ), "private unit row disappeared")

        deep_center = frac(module.C3 * q) * 182
        deep_left = deep_center + 13 * module.C3 * 182 * e_left
        deep_right = deep_center + 13 * module.C3 * 182 * e_right
        tooth = (Fraction(14 * root - 13), Fraction(14 * root))
        require(tooth[0] <= deep_left < deep_right <= tooth[1],
                "private deep half-tooth changed across the cylinder")

        metadata[n] = (
            rail_index, sector, edge, kappa, h, shallow, carry, root,
            weighted_hit(q * m.T, rails[rail_index][3])[2],
        )

    # All twelve legal next increments from n=110 fall in a present-free
    # strip on this entire D component.  The first private-unit candidate is
    # n=113, the failed lift of the projected residue edge 6 -> 9.
    continuation = []
    for n in range(111, 123):
        q = frac(z + Fraction(7 * n, m.R))
        left = (q + 13 * e_left) * m.T
        right = (q + 13 * e_right) * m.T
        rail_meets = not interval_disjoint(left, right, rail_support)
        present_meets = not interval_disjoint(left, right, present_support)
        unit = n % 13 in allowed
        continuation.append((n, unit, rail_meets, present_meets))
    require(all(rail and not present
                for _, _, rail, present in continuation),
            "terminal all-increment present-gap check failed")
    require(next(n for n, unit, _, _ in continuation if unit) == 113,
            "first admissible post-path address changed")

    # A full projected C7 x C13 return occurs inside the path.  It preserves
    # every retained metadata field, while the physical point drifts.  A
    # Fourier endpoint character at integer frequency X sees the same phase
    # exactly when 13^5 divides X.
    require(metadata[0] == metadata[91],
            "C91 return changed retained metadata")
    return_atom_weights = tuple(
        weighted_hit(
            frac(z + Fraction(7 * n, m.R)) * m.T,
            following["pieces"],
        )[2]
        for n in (0, 91)
    )
    require(return_atom_weights[0] == return_atom_weights[1],
            "C91 return changed the inherited atom weight")
    drift = Fraction(7 * 91, m.R)
    require(drift == Fraction(49, 13**5), "return drift changed")

    print("LRC14 MIXED D/SLOPE-SEVEN LONG-WORD INDEPENDENT REFEREE")
    print(f"rebuilt_D_component={rebuilt} radius={radius}")
    print(
        f"direct_path=(vertices={len(path)},slope_edges={len(steps)},"
        f"step1={steps.count(1)},step3={steps.count(3)},end={path[-1]})"
    )
    print("all_95_vertices_preserve=rail,present,delayed,carry,half-digit,root,unit")
    print(f"terminal_all_delta_1_to_12={tuple(continuation)}")
    print(
        "projected_cycle=(0,1,2,3,4,5,6,9,10,11,12,0);"
        "completed_lifts=8;first_failure=110->113"
    )
    print(
        "selected_chronology=P_95:partial_cube=True;"
        "residue_quotient=C_11:bipartite=False:partial_cube=False"
    )
    print(
        f"graceful_control=P_95:True;C_11:True;"
        f"C_11_labels={cycle_labels}"
    )
    print(
        f"C91_return=(0,91):drift={drift}:metadata_equal=True:"
        "endpoint_phase_equal_iff_13^5_divides_X"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
