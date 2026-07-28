#!/usr/bin/env python3
"""Exact long-word probe on the mixed D/slope-seven witness configuration.

The frozen THM-2680 endpoint at ``z`` refines to one THM-2640 packet.  Moving
by ``7*n/13^6`` preserves its delayed coordinate and advances carry/root by
``(7*n,n)`` modulo thirteen.  This script scans the complete translation
orbit ``n in Z/(13^6)`` using exact integer arithmetic.  Consecutive retained
orbit points at gaps 1..12 form lawful nonzero slope-seven steps.  The scan
therefore computes the longest composable same-configuration slope word after
the initial physical D edge and tests whether a closed all-depth word exists.

This is a support computation only.  It preserves the literal rail, present
factor, safe-sector delayed word, carry, future half-digit, private root, and
primitive-unit flag; it does not create a THM-2365 target action or semantic
endpoint current.
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


def strict_interval_index(value, starts, intervals):
    index = bisect_right(starts, value) - 1
    return index >= 0 and intervals[index][0] < value < intervals[index][1]


def strict_margin(value, intervals):
    """Distance to the boundary of the containing open interval union."""
    starts = tuple(left for left, _ in intervals)
    index = bisect_right(starts, value) - 1
    require(index >= 0 and intervals[index][0] < value < intervals[index][1],
            "point left a claimed strict interval")
    return min(value - intervals[index][0], intervals[index][1] - value)


def merge_support(pieces):
    raw = sorted((left, right) for left, right, weight in pieces if weight)
    out = []
    for left, right in raw:
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return tuple(out)


def strict_weighted_hit(value, pieces):
    """Return the unique weighted interval containing a strict point."""
    hits = tuple(piece for piece in pieces if piece[0] < value < piece[1])
    require(len(hits) == 1, "weighted point did not have a unique strict hit")
    return hits[0]


def open_arc_segments(center, radius, period):
    """Split an open circular interval into ordinary open intervals."""
    require(0 <= center < period and 0 < 2 * radius < period,
            "invalid circular test interval")
    left = center - radius
    right = center + radius
    if left < 0:
        return ((0, right), (left + period, period))
    if right > period:
        return ((left, period), (0, right - period))
    return ((left, right),)


def open_arc_meets_support(center, radius, support, period):
    """Whether an open circular interval has positive overlap with support."""
    return any(
        max(left, a) < min(right, b)
        for left, right in open_arc_segments(center, radius, period)
        for a, b in support
    )


def open_arc_is_contained(center, radius, support, period):
    """Whether every point of an open circular interval lies in support."""
    return all(
        any(a <= left and right <= b for a, b in support)
        for left, right in open_arc_segments(center, radius, period)
    )


def atom_margin(atom, value, prefixes, derivative):
    """Symmetric initial-x radius preserving one strict THM-2680 atom."""
    support = merge_support(atom["pieces"])
    base_margin = strict_margin(frac(value) * m.T, support)
    prefix = prefixes[atom["future"]][atom["h"]]
    delayed = tuple((left, left + length)
                    for left, length in zip(prefix[0], prefix[1]))
    delayed_margin = strict_margin(frac(m.R * value) * m.T, delayed)
    return min(base_margin / (derivative * m.T),
               delayed_margin / (derivative * m.R * m.T))


def cyclic_gap_components(good, modulus):
    """Return circular good-point components cut at gaps greater than twelve."""
    require(good, "empty orbit support")
    gaps = tuple(
        ((good[(i + 1) % len(good)] - good[i]) % modulus)
        for i in range(len(good))
    )
    breaks = tuple(i for i, gap in enumerate(gaps) if gap > 12)
    if not breaks:
        return gaps, (tuple(good),)
    components = []
    for j, cut in enumerate(breaks):
        next_cut = breaks[(j + 1) % len(breaks)]
        start = (cut + 1) % len(good)
        stop = next_cut
        row = []
        i = start
        while True:
            row.append(good[i])
            if i == stop:
                break
            i = (i + 1) % len(good)
        components.append(tuple(row))
    return gaps, tuple(components)


def forward_component_from(good, gaps, start):
    """Increasing directed component beginning at a specified good address."""
    index = good.index(start)
    row = [start]
    while gaps[index] <= 12:
        index = (index + 1) % len(good)
        if good[index] <= row[-1]:
            # A wrap by R is a lawful final edge only if the orbit then closes.
            # Starting at n=0, wrapping would revisit the initial state and is
            # recorded separately rather than appended to a finite prefix.
            break
        row.append(good[index])
    return tuple(row)


def main():
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    component = (
        Fraction(960117507257, 1930018885886),
        Fraction(324519717452867, 652346383429468),
    )
    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    pair_prefixes = m.build_pair_prefixes(module)
    require(component[0] < x < component[1] and frac(13 * x) == z,
            "frozen D-edge point changed")
    current = d.build_atom(
        module, prefixes, present, starts, rails[2], 3, 2, 0, 0
    )
    following = d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    require((current["j"], current["h"], following["j"], following["h"])
            == (5, 2, 2, 6), "frozen D labels changed")

    rail_index, sector, edge, kappa, h, shallow = 0, 0, 0, 1, 6, 1
    carry0, root0 = 2, 6
    shard = m.shard((rail_index, rail_index + 1))
    require(shard[1] == 26 and shard[5] == ((1, 0, 12),),
            "frozen private rail changed")
    rows = shard[6][0]
    good_residues = []
    for residue in range(m.P):
        carry = (carry0 + 7 * residue) % m.P
        root = (root0 + residue) % m.P
        if root and m.is_unit(rows[sector][edge][carry][kappa][h], root, 26):
            good_residues.append(residue)
    good_residues = tuple(good_residues)
    require(good_residues == (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
            "frozen unit/root residue bank changed")

    # The delayed y={Rz}, future half-digit, and deep-root membership are
    # invariant/covariant along z_n=z+7n/R.  Check the invariant once.
    y = frac(m.R * z) * m.T
    prefix = pair_prefixes[sector][shallow][h][kappa]
    delayed = tuple((left, left + length)
                    for left, length in zip(prefix[0], prefix[1]))
    delayed_starts = tuple(left for left, _ in delayed)
    require(strict_interval_index(y, delayed_starts, delayed),
            "frozen delayed word is absent")
    require(((m.R * z).numerator // (m.R * z).denominator) % m.P == carry0,
            "frozen predecessor carry changed")
    deep0 = frac(module.C3 * z) * 182
    require(Fraction(14 * root0 - 13) < deep0 < Fraction(14 * root0),
            "frozen private half-tooth is absent")

    rail_support = merge_support(rails[rail_index][3])
    present_support = tuple(present[shallow, (-h) % m.P])
    rail_starts = tuple(left for left, _ in rail_support)
    present_starts = tuple(left for left, _ in present_support)

    denominator = (z * m.T).denominator
    base_numerator = (z * m.T).numerator
    modulus = m.T * denominator
    step = (7 * m.T // m.R) * denominator
    require(m.T % m.R == 0 and step * m.R == 7 * m.T * denominator,
            "translation grid changed")
    scaled_rail = tuple((left * denominator, right * denominator)
                        for left, right in rail_support)
    scaled_present = tuple((left * denominator, right * denominator)
                           for left, right in present_support)
    scaled_rail_starts = tuple(left for left, _ in scaled_rail)
    scaled_present_starts = tuple(left for left, _ in scaled_present)

    good = []
    point = base_numerator
    for n in range(m.R):
        if (n % m.P in good_residues
                and strict_interval_index(
                    point, scaled_rail_starts, scaled_rail)
                and strict_interval_index(
                    point, scaled_present_starts, scaled_present)):
            good.append(n)
        point = (point + step) % modulus

    gaps, components = cyclic_gap_components(good, m.R)
    gap_hist = tuple(sorted(
        (gap, gaps.count(gap)) for gap in set(gaps)
    ))
    largest = max(components, key=len)
    start_component = next(row for row in components if 0 in row)
    start_forward = forward_component_from(good, gaps, 0)
    expected_forward = tuple(
        n for n in range(111) if n % m.P in good_residues
    )
    require(start_forward == expected_forward and len(start_forward) == 95,
            "mixed D endpoint forward component changed")
    forward_steps = tuple(b - a for a, b in zip(
        start_forward, start_forward[1:]
    ))
    require(len(forward_steps) == 94
            and forward_steps.count(1) == 86
            and forward_steps.count(3) == 8,
            "mixed long-word step census changed")

    # The already-audited local word visits all eleven active residues before
    # any cocycle wrap.  Record its exact slope-step spelling explicitly.
    local_vertices = good_residues
    local_steps = tuple(b - a for a, b in zip(local_vertices, local_vertices[1:]))
    require(local_steps == (1, 1, 1, 1, 1, 1, 3, 1, 1, 1),
            "local mixed-word spelling changed")
    projected_cycle = local_vertices + (0,)
    projected_cycle_steps = local_steps + (1,)
    require(projected_cycle_steps == (1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1),
            "projected residue cycle changed")
    require(len(local_vertices) == 11 and len(local_vertices) % 2 == 1,
            "projected quotient stopped being an odd cycle")
    path_labels = tuple(
        index // 2 if index % 2 == 0 else 94 - index // 2
        for index in range(95)
    )
    cycle_labels = (0, 4, 6, 5, 8, 3, 9, 2, 10, 1, 11)
    require(set(path_labels) == set(range(95))
            and tuple(abs(a - b) for a, b in zip(
                path_labels, path_labels[1:]
            )) == tuple(range(94, 0, -1)),
            "standard graceful path labelling changed")
    require(len(set(cycle_labels)) == 11
            and set(cycle_labels).issubset(set(range(12)))
            and tuple(abs(a - b) for a, b in zip(
                cycle_labels, cycle_labels[1:] + cycle_labels[:1]
            )) == (4, 2, 1, 3, 5, 6, 7, 8, 9, 10, 11),
            "explicit graceful C11 certificate changed")
    completed_cycles = 8
    for height in range(completed_cycles):
        lifted_cycle = tuple(13 * height + residue
                             for residue in local_vertices) + (13 * (height + 1),)
        require(all(address in start_forward for address in lifted_cycle),
                "a claimed projected cycle did not lift through the word")
    require(start_forward[-7:] == tuple(range(104, 111)),
            "terminal lifted-cycle prefix changed")

    # The lifted path is longer than one C7 x C13 residue period.  Addresses
    # n=0 and n=91 have identical fixed configuration, carry/root, unit,
    # delayed, present, and rail-weight metadata, but remain distinct physical
    # points.  This is a concrete pumping obstruction: quotient recurrence
    # forgets the odometer height floor(n/13) (and ultimately n itself).
    return_address = 91
    require(return_address in start_forward
            and return_address % 7 == return_address % 13 == 0,
            "coarse C91 return left the lifted component")
    return_drift = Fraction(7 * return_address, m.R)
    require(return_drift == Fraction(49, 13**5),
            "coarse C91 return drift changed")

    # Every retained factor is strict at the same initial x.  Compute an
    # explicit common symmetric cylinder rather than infer positivity from 95
    # unrelated point tests.  A perturbation e of x perturbs z_n by 13e.
    radius = min(
        atom_margin(current, x, prefixes, 1),
        atom_margin(following, z, prefixes, 13),
    )
    delayed_margin = strict_margin(y, delayed) / (13 * m.R * m.T)
    carry_fraction = frac(m.R * z)
    carry_margin = min(carry_fraction, 1 - carry_fraction) / (13 * m.R)
    for n in start_forward:
        q = frac(z + Fraction(7 * n, m.R))
        carry = (carry0 + 7 * n) % m.P
        root = (root0 + n) % m.P
        require(root and m.is_unit(
            rows[sector][edge][carry][kappa][h], root, 26
        ), "forward word lost its primitive unit")
        base = frac(q) * m.T
        rail_margin = strict_margin(base, rail_support) / (13 * m.T)
        present_margin = strict_margin(base, present_support) / (13 * m.T)
        deep = frac(module.C3 * q) * 182
        deep_interval = (Fraction(14 * root - 13), Fraction(14 * root))
        deep_margin = min(
            deep - deep_interval[0], deep_interval[1] - deep
        ) / (13 * module.C3 * 182)
        require(deep_margin > 0, "forward word lost its private half-tooth")
        radius = min(radius, rail_margin, present_margin, deep_margin,
                     delayed_margin, carry_margin)
    require(radius > 0, "mixed long word has no common positive cylinder")

    return_metadata = []
    for n in (0, return_address):
        q = frac(z + Fraction(7 * n, m.R))
        base = q * m.T
        carry = (carry0 + 7 * n) % m.P
        root = (root0 + n) % m.P
        rail_piece = strict_weighted_hit(base, rails[rail_index][3])
        atom_piece = strict_weighted_hit(base, following["pieces"])
        return_metadata.append((
            rail_index, sector, edge, kappa, h, shallow, carry, root,
            m.is_unit(rows[sector][edge][carry][kappa][h], root, 26),
            rail_piece[2], atom_piece[2],
            strict_interval_index(base, present_starts, present_support),
            strict_interval_index(y, delayed_starts, delayed),
        ))
    require(return_metadata[0] == return_metadata[1],
            "coarse C91 return did not preserve retained packet metadata")

    # Recheck the final cylinder as a set, rather than only through its
    # centre and separately computed margins.  The two endpoints may land on
    # support boundaries; the open interval between them is the carrier.
    q_half_width = 13 * radius * m.T
    for n in start_forward:
        base = frac(z + Fraction(7 * n, m.R)) * m.T
        require(open_arc_is_contained(
            base, q_half_width, rail_support, m.T
        ), "common cylinder escaped the rail support")
        require(open_arc_is_contained(
            base, q_half_width, present_support, m.T
        ), "common cylinder escaped the present support")

    first_post = 113
    post_point = (base_numerator + first_post * step) % modulus
    post_rail = strict_interval_index(
        post_point, scaled_rail_starts, scaled_rail
    )
    post_present = strict_interval_index(
        post_point, scaled_present_starts, scaled_present
    )
    require(first_post % m.P in good_residues
            and not (post_rail and post_present),
            "first post-component admissible residue unexpectedly survived")

    # From the terminal address n=110, every allowed next slope increment
    # delta=1,...,12 lands at one of n=111,...,122.  The rail still meets the
    # entire pulled-forward D cylinder in every case, but the present factor
    # misses it completely.  This proves maximality only for this fixed
    # configuration and inherited D component; it says nothing about a
    # configuration switch.
    continuation = []
    for n in range(111, 123):
        base = frac(z + Fraction(7 * n, m.R)) * m.T
        root_unit = n % m.P in good_residues
        rail_meets = open_arc_meets_support(
            base, q_half_width, rail_support, m.T
        )
        present_meets = open_arc_meets_support(
            base, q_half_width, present_support, m.T
        )
        continuation.append((n, n % m.P, root_unit,
                             rail_meets, present_meets))
    continuation = tuple(continuation)
    require(all(row[3] and not row[4] for row in continuation),
            "a next slope increment escaped the present-gap obstruction")
    require(next(row[0] for row in continuation if row[2]) == first_post,
            "first post-component root/unit-admissible address changed")
    require((start_forward[-1] % m.P, first_post % m.P,
             first_post - start_forward[-1]) == (6, 9, 3),
            "first failed lift of the projected cycle changed")

    print("LRC14 MIXED D/SLOPE-SEVEN LONG-WORD ORBIT PROBE")
    print(f"R={m.R} z={z} orbit_step=7/R")
    print(f"D_component={component} x={x} z=D(x)")
    print(f"unit_root_residues_mod13={good_residues}")
    print(f"good_orbit_points={len(good)}")
    print(f"gap_hist={gap_hist} max_gap={max(gaps)}")
    print(f"cyclic_components={len(components)} sizes={tuple(sorted(map(len, components)))}")
    print(
        f"largest_component=(size={len(largest)},start={largest[0]},"
        f"end={largest[-1]})"
    )
    print(
        f"component_through_mixed_D_endpoint=(size={len(start_component)},"
        f"forward_size={len(start_forward)},"
        f"forward_prefix={start_forward[:20]},end={start_forward[-1]})"
    )
    print(
        f"heterogeneous_word=(D_edges=1,slope_edges={len(forward_steps)},"
        f"total_edges={1 + len(forward_steps)},vertices={1 + len(start_forward)})"
    )
    print(
        f"forward_step_hist=((1,{forward_steps.count(1)}),"
        f"(3,{forward_steps.count(3)}))"
    )
    print(
        f"common_initial_cylinder=({x-radius},{x+radius}) "
        f"radius={radius} length={2*radius}"
    )
    print(
        f"first_post_component_admissible_address={first_post}:"
        f"rail={post_rail}:present={post_present}:fails="
        f"{'rail' if not post_rail else 'present'}"
    )
    print(
        "terminal_continuation_table=(n,residue,root_unit,rail_meets,"
        f"present_meets)={continuation}"
    )
    print(f"local_vertices={local_vertices} local_steps={local_steps}")
    print(
        f"projected_residue_cycle={projected_cycle} "
        f"steps={projected_cycle_steps} completed_lifts={completed_cycles}"
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
        f"coarse_C91_return=(addresses=0,{return_address},"
        f"heights=0,{return_address // 13},drift={return_drift},"
        "retained_metadata_equal=True,base_points_equal=False)"
    )
    print(
        "endpoint_character_boundary=same_phase_iff_13^5_divides_"
        "integer_frequency"
    )
    print(
        "first_projected_cycle_lift_failure="
        f"(height={start_forward[-1] // 13},address={start_forward[-1]}"
        f"->{first_post},residue=6->9,step=3,reason=present)"
    )
    print(
        "closed_all_depth_word="
        f"{max(gaps) <= 12}; criterion=max circular gap <=12"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
