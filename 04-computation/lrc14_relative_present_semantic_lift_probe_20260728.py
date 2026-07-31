#!/usr/bin/env python3
"""Exact relative-present probe on the THM-2707 source-one lift fibre.

The full THM-2707 packet SCC and its THM-2712 following-atom part have no
THM-2305 semantic handoff vertex.  This companion first verifies that sharp
negative statement.  It then removes the dynamically typed present factor
from the fourteen source-one rails and recomputes, rather than inherits, the
seven-coordinate primitive coefficient vector.  The difference

    (present-free vector) - (present vector)

is the exact coefficient on the present complement.  This is a bounded
mapping-cone scout on one canonical fibre; it is not a target current or a
row exclusion.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_central_half_odometer_full_local_cycle_thm2698 as half
import lrc14_full_physical_lift_fibre_thm2707 as lift
import lrc14_half_c221_semantic_source_fibre_census_20260728 as semantic
import lrc14_predecessor_carry_private_root_atlas_thm2640 as private


P = 13
R = P**6
Z = Fraction(46873542509301, 100360982066072)
I_RADIUS = Fraction(1, 1304692766858936)
Q_RADIUS = P * I_RADIUS
ACTIVE_RESIDUES = (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12)
SOURCE_ONE_RAILS = tuple(range(14))
CANONICAL_CONTENT = 26


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def strict_contains(value, starts, intervals):
    index = bisect_right(starts, value) - 1
    return (
        index >= 0
        and intervals[index][0] < value < intervals[index][1]
    )


def merge_open(intervals):
    """Union overlapping open intervals, but retain a touching-point gap."""
    rows = []
    for left, right in sorted(intervals):
        if rows and left < rows[-1][1]:
            rows[-1] = (rows[-1][0], max(rows[-1][1], right))
        else:
            rows.append((left, right))
    return tuple(rows)


def distance(value):
    return semantic.distance(value)


def predicate_boundary_radius(value, speed, danger_denominator):
    """Physical-x distance to the nearest boundary of one danger comb."""
    return abs(
        distance(speed * value) - Fraction(1, danger_denominator)
    ) / speed


def shifted_danger(value, speed, centre):
    phase = semantic.centered(speed * value - centre)
    radius = Fraction(1, 14)
    return -radius <= phase < radius


def full_target_labels(value):
    """Lawful THM-2365/2407 (s,t) labels at one E3 point."""
    s_labels = tuple(
        s for s in range(P)
        if not shifted_danger(value, semantic.UNITS[0], -Fraction(s, P))
        and not shifted_danger(value, semantic.BLOCKERS[1], Fraction(s, P))
    )
    t_labels = tuple(
        t for t in range(P)
        if not shifted_danger(value, semantic.UNITS[1], -Fraction(t, P))
        and not shifted_danger(value, semantic.BLOCKERS[2], Fraction(t, P))
    )
    return s_labels, t_labels


def semantic_stability_radius(value):
    """Radius on which the complete semantic record is unchanged."""
    source, _word = semantic.semantic_record(value)
    rows = [predicate_boundary_radius(value, 1, 7)]
    rows.extend(
        predicate_boundary_radius(value, speed, 14)
        for speed in semantic.UNITS + semantic.BLOCKERS
    )
    if source:
        steps = semantic.CLOCKS[source - 1]
        endpoint = semantic.ordinary_dilate(value, steps)
        expansion = P**steps
        rows.append(
            predicate_boundary_radius(endpoint, 1, 7) / expansion
        )
        rows.extend(
            predicate_boundary_radius(endpoint, speed, 14) / expansion
            for speed in semantic.UNITS + semantic.BLOCKERS
        )
    require(min(rows) > 0, "semantic record lies on a wall")
    return min(rows)


def complement_slack(value, intervals, starts, period):
    """Distance from an exterior point to a periodic interval union."""
    index = bisect_right(starts, value) - 1
    candidates = []
    if index >= 0:
        left, right = intervals[index]
        require(not (left < value < right),
                "complement slack requested inside present")
        candidates.append(value - right)
    if index + 1 < len(intervals):
        candidates.append(intervals[index + 1][0] - value)
    if intervals:
        candidates.extend((
            value + period - intervals[-1][1],
            intervals[0][0] + period - value,
        ))
    candidates = tuple(row for row in candidates if row >= 0)
    require(candidates, "empty complement-boundary ledger")
    return min(candidates)


def local_relative_radius(module, pair_prefixes, rails, present,
                          present_starts, row):
    """Strict q-radius preserving the relative-present endpoint state."""
    q = row["q"]
    rail_index = row["rail"]
    ell = row["shallow"]
    edge = row["edge"]
    root = row["root"]
    coordinate = q * lift.m.T
    rail_starts = tuple(left for left, *_ in rails[rail_index][3])
    rail_radius = half.interval_slack(
        coordinate, rails[rail_index][3], rail_starts
    ) / lift.m.T
    present_radius = complement_slack(
        coordinate,
        present[ell, 7],
        present_starts[ell, 7],
        lift.m.T,
    ) / lift.m.T

    y = semantic.frac(R * q)
    delayed_radius = half.interval_slack(
        y, half.prefix_intervals(pair_prefixes[0][ell][6][1])
    ) / R
    half_digit_radius = min(
        y - Fraction(13, 26), Fraction(14, 26) - y
    ) / R
    carry_radius = min(y, 1 - y) / R
    deep = semantic.frac(module.C3 * q) * 182
    low = 14 * root - 13 if edge == 0 else 14 * root
    deep_radius = min(deep - low, low + 13 - deep) / (182 * module.C3)
    rows = (
        semantic.semantic_cospan_radius(q, 3, (1, 2)),
        rail_radius,
        present_radius,
        delayed_radius,
        half_digit_radius,
        carry_radius,
        deep_radius,
        half.clock_slack(q, P),
        half.clock_slack(q, P**2),
    )
    require(all(radius > 0 for radius in rows),
            "relative-present state is not strict")
    return min(rows)


def coefficient_rows(module, pair_prefixes, rails, present, starts):
    """Old, present-free, and complement vectors on the two E3 residues."""
    # Address residue 7 has carry 12 and only right root 12.  Address
    # residue 8 has carry 6 and only left root 1.
    states = ((7, 1, 12, 12), (8, 0, 6, 1))
    caches = [[{} for _h in range(P)] for _ell in range(7)]
    result = []
    for rail_index in SOURCE_ONE_RAILS:
        pieces = rails[rail_index][3]
        for residue, edge, carry, root in states:
            low = 14 * root - 13 if edge == 0 else 14 * root
            free_half = private.old.intersect_weighted_comb(
                pieces, module.C3, 182, low, low + 13
            )
            old_vector = []
            free_vector = []
            for ell in range(7):
                overlap = private.old.intersect_weighted_union(
                    pieces, present[ell, 7], starts[ell, 7]
                )
                old_half = private.old.intersect_weighted_comb(
                    overlap, module.C3, 182, low, low + 13
                )
                old_values = private.delayed_carry_pair(
                    old_half,
                    pair_prefixes[0][ell][6],
                    caches[ell][6],
                )
                free_values = private.delayed_carry_pair(
                    free_half,
                    pair_prefixes[0][ell][6],
                    caches[ell][6],
                )
                old_vector.append(old_values[carry][1])
                free_vector.append(free_values[carry][1])
            old_vector = tuple(old_vector)
            free_vector = tuple(free_vector)
            relative_vector = tuple(
                free - old
                for old, free in zip(old_vector, free_vector)
            )
            require(min(relative_vector) >= 0
                    and tuple(old + relative for old, relative in zip(
                        old_vector, relative_vector
                    )) == free_vector,
                    "present/free coefficient additivity failed")
            result.append({
                "rail": rail_index,
                "residue": residue,
                "edge": edge,
                "carry": carry,
                "root": root,
                "old": old_vector,
                "free": free_vector,
                "relative": relative_vector,
            })
    return tuple(result)


def content(rows, key):
    value = 0
    for row in rows:
        for coefficient in row[key]:
            value = gcd(value, coefficient)
    return value


def unit_rails(rows, key):
    return {
        residue: tuple(
            row["rail"] for row in rows
            if row["residue"] == residue
            and private.is_unit(
                row[key], row["root"], CANONICAL_CONTENT
            )
        )
        for residue in (7, 8)
    }


def main():
    require((lift.m.P, lift.m.R, private.P, private.R) == (P, R, P, R),
            "inherited carrier scales changed")
    module, _prefixes, _, _, rails, present, starts = (
        lift.m.core.build_carrier_data()
    )
    pair_prefixes = private.build_pair_prefixes(module)
    require(tuple(rail[:3] for rail in rails[:14]) == tuple(
        (1, owner, deep)
        for owner in range(7) for deep in (12, 0)
    ), "source-one rail ordering changed")

    delayed = semantic.frac(R * Z)
    h = half.floor_fraction(P * delayed)
    kappa = half.floor_fraction(2 * P * delayed) - 2 * h
    sector_words = half.prior.sector_words(module)
    sector_starts = tuple(tuple(left for left, _ in word)
                          for word in sector_words)
    active_sectors = tuple(
        sector for sector in range(2)
        if half.strict_interval_member(
            delayed * lift.m.T,
            sector_words[sector], sector_starts[sector]
        )
    )
    active_clocks = tuple(
        ell for ell in range(7)
        if half.strict_interval_member(
            delayed,
            half.prefix_intervals(pair_prefixes[0][ell][h][kappa]),
        )
    )
    require((h, kappa, active_sectors, active_clocks)
            == (6, 1, (0,), (1, 2, 3, 4, 5, 6)),
            "delayed source-one state changed")

    denominator = (Z * lift.m.T).denominator
    modulus = lift.m.T * denominator
    point = (Z * lift.m.T).numerator
    orbit_step = (7 * lift.m.T // R) * denominator
    require(orbit_step * R == 7 * lift.m.T * denominator,
            "address orbit grid changed")

    old = lift.old
    fixed_rail = old.merge_support(rails[0][3])
    fixed_present = tuple(present[1, 7])
    scaled_fixed_rail = tuple(
        (left * denominator, right * denominator)
        for left, right in fixed_rail
    )
    scaled_fixed_present = tuple(
        (left * denominator, right * denominator)
        for left, right in fixed_present
    )
    fixed_rail_starts = tuple(left for left, _ in scaled_fixed_rail)
    fixed_present_starts = tuple(left for left, _ in scaled_fixed_present)

    rail_intervals = []
    per_rail = []
    for rail_index in SOURCE_ONE_RAILS:
        intervals = tuple(
            (left * denominator, right * denominator)
            for left, right, weight in rails[rail_index][3]
            if weight
        )
        per_rail.append((tuple(left for left, _ in intervals), intervals))
        rail_intervals.extend(intervals)
    rail_union = merge_open(rail_intervals)
    rail_union_starts = tuple(left for left, _ in rail_union)

    q_denominator = Z.denominator
    require(q_denominator % R == 0, "q orbit lost its fixed denominator")
    q_numerator = Z.numerator
    q_step = 7 * (q_denominator // R)

    packet_nodes = []
    packet_semantics = Counter()
    expanded_semantics = Counter()
    e3_rows = []
    for n in range(R):
        record_numerator = q_numerator % q_denominator
        if (n % P in ACTIVE_RESIDUES
                and strict_contains(
                    point, fixed_rail_starts, scaled_fixed_rail
                )
                and strict_contains(
                    point, fixed_present_starts, scaled_fixed_present
                )):
            record = semantic.semantic_record_numerator(
                record_numerator, q_denominator
            )
            packet_nodes.append(n)
            packet_semantics[record] += 1

        if strict_contains(point, rail_union_starts, rail_union):
            record = semantic.semantic_record_numerator(
                record_numerator, q_denominator
            )
            expanded_semantics[record] += 1
            if record == (3, (1, 2)):
                q = Fraction(record_numerator, q_denominator)
                require(semantic.semantic_record(q) == record,
                        "fraction/integer semantic paths disagree")
                rail_indices = tuple(
                    rail_index
                    for rail_index, (rail_starts, intervals)
                    in enumerate(per_rail)
                    if strict_contains(point, rail_starts, intervals)
                )
                require(len(rail_indices) == 1,
                        "E3 point does not have a unique source-one rail")
                rail_index = rail_indices[0]
                ell = half.shallow(q)
                owner = half.owner(q)
                require(rails[rail_index][1] == owner and ell in active_clocks,
                        "E3 rail/clock typing changed")
                present_labels = tuple(
                    label for label in range(P)
                    if half.strict_interval_member(
                        q * lift.m.T,
                        present[ell, label], starts[ell, label]
                    )
                )
                require(not present_labels,
                        "E3 point entered a present label")
                residue = n % P
                require(residue in (7, 8),
                        "E3 rail bank acquired another address residue")
                edge, carry, root = (
                    (1, 12, 12) if residue == 7 else (0, 6, 1)
                )
                require(carry == (2 + 7 * n) % P,
                        "E3 predecessor carry changed")
                forced_root = (
                    2 * carry + (2 * h + kappa) // P
                    + int(edge == 0)
                ) % P
                require(forced_root == root
                        and half.half_membership(
                            module, q, edge, root
                        ), "E3 private-root typing changed")
                e3_rows.append({
                    "n": n,
                    "q": q,
                    "rail": rail_index,
                    "shallow": ell,
                    "owner": owner,
                    "residue": residue,
                    "edge": edge,
                    "carry": carry,
                    "root": root,
                    "target_labels": full_target_labels(q),
                })
        point = (point + orbit_step) % modulus
        q_numerator = (q_numerator + q_step) % q_denominator

    packet_nodes = tuple(packet_nodes)
    e3_rows = tuple(e3_rows)
    require(len(packet_nodes) == 3346
            and packet_semantics == Counter({
                (0, None): 3081, (2, None): 265,
            }), "full packet semantic census changed")
    atom_nodes = tuple(n for n in packet_nodes if n % P == 0)
    atom_semantics = Counter(
        semantic.semantic_record(
            semantic.frac(Z + Fraction(7 * n, R))
        ) for n in atom_nodes
    )
    require(len(atom_nodes) == 304
            and atom_semantics == Counter({
                (0, None): 280, (2, None): 24,
            }), "following-atom semantic census changed")
    packet_radii = tuple(
        semantic_stability_radius(
            semantic.frac(Z + Fraction(7 * n, R))
        ) for n in packet_nodes
    )
    require(min(packet_radii) > Q_RADIUS,
            "packet semantic record changes inside common I")

    expected_expanded = Counter({
        (0, None): 301390,
        (1, None): 13639,
        (1, (3,)): 4996,
        (1, (2, 3)): 493,
        (2, None): 11406,
        (3, (1, 2)): 12848,
    })
    require(expanded_semantics == expected_expanded
            and len(e3_rows) == 12848,
            "expanded source-one semantic census changed")
    residue_counts = Counter(row["residue"] for row in e3_rows)
    require(residue_counts == Counter({7: 6404, 8: 6444}),
            "E3 residue polarization changed")

    coefficient_bank = coefficient_rows(
        module, pair_prefixes, rails, present, starts
    )
    contents = {
        key: content(coefficient_bank, key)
        for key in ("old", "free", "relative")
    }
    require(contents == {
        "old": 2122120,
        "free": 2122120,
        "relative": 2122120,
    } and all(
        coefficient % CANONICAL_CONTENT == 0
        for row in coefficient_bank
        for key in ("old", "free", "relative")
        for coefficient in row[key]
    ), "relative coefficient content changed")
    unit_banks = {
        key: unit_rails(coefficient_bank, key)
        for key in ("old", "free", "relative")
    }
    all_rails = tuple(range(14))
    except_3 = tuple(index for index in all_rails if index != 3)
    except_13 = tuple(index for index in all_rails if index != 13)
    require(unit_banks == {
        "old": {7: all_rails, 8: ()},
        "free": {7: except_13, 8: except_3},
        "relative": {7: all_rails, 8: except_3},
    }, "relative primitive-unit atlas changed")

    relative_rows = tuple(
        row for row in e3_rows
        if row["rail"] in unit_banks["relative"][row["residue"]]
    )
    relative_residues = Counter(row["residue"] for row in relative_rows)
    rail_residue_counts = Counter(
        (row["rail"], row["residue"]) for row in e3_rows
    )
    relative_rail_residue_counts = Counter(
        (row["rail"], row["residue"]) for row in relative_rows
    )
    require(relative_residues[7] == 6404
            and relative_residues[8] > 0,
            "relative carrier lost one residue side")
    directed_edges = 2 * relative_residues[7] * relative_residues[8]
    target_profile_counts = Counter(
        (len(row["target_labels"][0]), len(row["target_labels"][1]))
        for row in relative_rows
    )
    require(all(row["target_labels"][1]
                and 0 not in row["target_labels"][1]
                for row in relative_rows),
            "a relative E3 vertex lost every nonzero t label")

    by_residue_clock = defaultdict(Counter)
    by_residue_clock_rows = defaultdict(list)
    relative_radii = []
    for row in relative_rows:
        clock = (row["shallow"], row["owner"])
        by_residue_clock[row["residue"]][clock] += 1
        by_residue_clock_rows[row["residue"], clock].append(row)
        radius = local_relative_radius(
            module, pair_prefixes, rails, present, starts, row
        )
        row["relative_radius"] = radius
        relative_radii.append(radius)
    require(min(relative_radii) > 0,
            "relative semantic carrier lost local positive width")
    # Equality is sufficient because the inherited cylinder is open: a
    # binding wall may occur at either excluded endpoint.
    common_I_count = sum(radius >= Q_RADIUS for radius in relative_radii)
    require(common_I_count == len(relative_radii),
            "a relative state lost the inherited open cylinder")

    reverse_clock_pairs = sum(
        count * by_residue_clock[8][(owner, shallow)]
        for (shallow, owner), count in by_residue_clock[7].items()
    )
    require(reverse_clock_pairs > 0,
            "relative carrier has no reverse-clock pair")
    witness = None
    for first in sorted(
            (row for row in relative_rows if row["residue"] == 7),
            key=lambda row: row["n"]):
        targets = by_residue_clock_rows[
            8, (first["owner"], first["shallow"])
        ]
        if targets:
            second = min(targets, key=lambda row: row["n"])
            witness = (first, second)
            break
    require(witness is not None, "reverse-clock witness vanished")
    first, second = witness
    k_forward = 7 * (second["n"] - first["n"]) % R
    k_reverse = 7 * (first["n"] - second["n"]) % R
    require(
        semantic.frac(first["q"] + Fraction(k_forward, R))
        == second["q"]
        and semantic.frac(second["q"] + Fraction(k_reverse, R))
        == first["q"]
        and (k_forward + k_reverse) % R == 0
        and 2 * k_forward % P == 1
        and 2 * k_reverse % P == P - 1
        and (first["root"] + 2 * k_forward) % P == 0
        and (second["root"] + 2 * k_reverse) % P == 0
        and first["edge"] != second["edge"],
        "relative zero-root clutching obstruction changed",
    )
    common_s = tuple(sorted(set(first["target_labels"][0]).intersection(
        second["target_labels"][0]
    )))
    common_t = tuple(sorted(set(first["target_labels"][1]).intersection(
        second["target_labels"][1]
    )))

    print("LRC14 RELATIVE-PRESENT SEMANTIC LIFT PROBE")
    print("status=FINITE-EXACT relative-complement carrier and zero-root clutch no-go")
    print(f"R={R} delayed={delayed} state=(sector0,h6,kappa1) "
          f"active_clocks={active_clocks}")
    print(f"fixed_packet_vertices={len(packet_nodes)} "
          f"semantics={tuple(sorted(packet_semantics.items(), key=repr))}")
    print(f"atom_vertices={len(atom_nodes)} "
          f"semantics={tuple(sorted(atom_semantics.items(), key=repr))}")
    print(f"fixed_packet_min_semantic_stability={min(packet_radii)} "
          f"q_common_radius={Q_RADIUS}")
    print("existing_carrier_semantic_handoff_vertices=0")
    print("expanded_source_one_semantics="
          f"{tuple(sorted(expanded_semantics.items(), key=repr))}")
    print(f"E3_fork_residue_counts={tuple(sorted(residue_counts.items()))}")
    print("E3_present_label_union=EMPTY; unique_active_rail=True; "
          "private_root=True")
    print(f"coefficient_contents={tuple(sorted(contents.items()))} "
          f"canonical_content={CANONICAL_CONTENT}")
    print("unit_rail_banks="
          f"{tuple((key, tuple(sorted(value.items()))) for key, value in unit_banks.items())}")
    print("E3_rail_residue_counts="
          f"{tuple(sorted(rail_residue_counts.items()))}")
    print("relative_survivor_rail_residue_counts="
          f"{tuple(sorted(relative_rail_residue_counts.items()))}")
    print(f"relative_residue_counts={tuple(sorted(relative_residues.items()))} "
          f"directed_bipartite_edges={directed_edges}")
    print(f"full_target_label_size_profiles={tuple(sorted(target_profile_counts.items()))} "
          "every_vertex_has_nonzero_t=True")
    print("relative_clock_counts="
          f"{tuple((residue, tuple(sorted(rows.items()))) for residue, rows in sorted(by_residue_clock.items()))}")
    print(f"reverse_clock_unordered_pairs={reverse_clock_pairs} "
          f"directed_pairs={2 * reverse_clock_pairs}")
    print(f"relative_min_radius={min(relative_radii)} "
          f"common_I_survivors={common_I_count}/{len(relative_radii)}")
    print("relative_scalar_pair="
          f"((n={first['n']},q={first['q']},rail={first['rail']},"
          f"clock={(first['shallow'], first['owner'])},edge={first['edge']},root={first['root']}),"
          f"(n={second['n']},q={second['q']},rail={second['rail']},"
          f"clock={(second['shallow'], second['owner'])},edge={second['edge']},root={second['root']}),"
          f"k=({k_forward},{k_reverse}),common_s={common_s},common_t={common_t})")
    print("private_root_transport=EMPTY: every residue7<->8 scalar lift "
          "preserves edge and lands at root0; the lawful endpoint lies on "
          "the opposite half-tooth edge")
    print("mechanism=present-free=present+relative-complement coefficientwise; "
          "the complement repairs the residue8 determinant on 13/14 rails; "
          "full U_(s,t) labels survive but do not repair the zero-root clutch")
    print("SCOPE: relative-present endpoint/unit support and scalar lifts; "
          "private-root lift attachment is empty; no canonical transition "
          "amplitude, spectral endpoint pair, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
