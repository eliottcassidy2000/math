#!/usr/bin/env python3
"""Exact full-bank orbit/action audit for the THM-2859 q0 endpoint mask.

The distinguished q0 mask H is compared with all thirteen masks
Z^k H, where Z is translation by (0,1) on F_13^2.  Every labelled vertex
of the complete 685-component THM-2825 collar forest is inspected in both
the source and target endpoint frames.  The 587 cofiber-rooted components
and 98 common-only components remain separately typed.  Every hit is then
typed by its physical and
labelled coordinates, semantic value, six native factors in both frames,
13-twist carrier masks, and exact ancestry chamber.

There are two deliberately separate outputs:

* orbit membership is a statement about equality of 169-bit masks;
* a physical C13 action would additionally need a total, single-valued
  generator with identities, inverses, composable powers, and order 13.

The second condition is tested on the complete directed forest relation and on the
distinguished +66h physical displacement.  Merely seeing generator-valued
labels is never promoted to an action.

No executable Python ``assert`` statement is used.
"""

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_nearest_half_step_common_right_collar_thm2825.py":
        "bd9ffe7f6815b5c563bd483c300118fbdd683f3d9303babbab7912e031747c9a",
    "lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.py":
        "af924103c30a45eeab88520d6b569c00a5a8c9cafbc8052b3585dd86e07b5dac",
    "lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825.py":
        "5a74cfa0b9a5f30e15460e7e1bee17214193d1f307e79f1eeb68b8781b7d23d6",
    "lrc14_horn_collar_endpoint_coboundary_thm2859.py":
        "3d7702641a2df258b829538d8fcf1d066cdf5f426cceef5781bbcfb37747bc15",
    "lrc14_nearest_half_step_common_right_collar_path_operator_thm2825.py":
        "7f0780e70161cbfafc4499d02a7d1c5aa40366b6dfa9935b00dc652e3d54c8e0",
}
for name, expected in PINNED.items():
    payload = (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"pinned dependency changed: {name}",
    )


import lrc14_nearest_half_step_common_right_collar_thm2825 as collar
import lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825 as physical
import lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825 as independent
import lrc14_horn_collar_endpoint_coboundary_thm2859 as endpoint


P = 13
HSTEP = collar.H
SHIFT = collar.copies.SHIFT
I = (142004992589460, 142005019034340)


def safe_comb_contains(interval, module, speed, denominator, low, high):
    require(
        module.T % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = module.T // (denominator * speed)
    step = denominator * unit
    width = (high - low) * unit
    base = (low % denominator) * unit
    left, right = interval
    first = (left - base) // step
    for number in range(first - 1, first + 3):
        danger_left = base + number * step
        danger_right = danger_left + width
        if max(left, danger_left) < min(right, danger_right):
            return False
    return True


def factor_signature(interval, full_module, e3, clocks, cell):
    clock, sigma, target = cell
    return (
        collar.copies.contained(interval, e3),
        collar.copies.contained(interval, clocks[clock]),
        safe_comb_contains(
            interval, full_module, full_module.W[1], 182,
            -14 * sigma - 13, -14 * sigma + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.W[2], 182,
            -14 * target - 13, -14 * target + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.C2, 182,
            14 * sigma - 13, 14 * sigma + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.C3, 182,
            14 * target - 13, 14 * target + 13,
        ),
    )


def true_indices(row):
    return tuple(index for index, value in enumerate(row) if value)


def compact_carrier(row):
    return tuple(true_indices(mask) for mask in row)


def run_length_encode(labels):
    if not labels:
        return ()
    result = []
    start, value = labels[0]
    previous = start
    for index, label in labels[1:]:
        if index == previous + 1 and label == value:
            previous = index
            continue
        result.append((start, previous, value))
        start, previous, value = index, index, label
    result.append((start, previous, value))
    return tuple(result)


def build_paths():
    (
        _module,
        _rails,
        _present,
        details,
        full_module,
        e3,
        clocks,
        q_pairs,
        _delayed,
        _source_weight,
        _target_weight,
        _rail_common,
    ) = collar.copies.physical_setup()
    paths = []
    cells = {}
    common_atoms = 0
    right_roots = 0
    component_kinds = Counter()
    for clock in range(7):
        for sigma in collar.COMMON_S:
            for target in collar.COMMON_T:
                cell = (clock, sigma, target)
                source, target_physical, common, right = collar.cell_objects(
                    details,
                    full_module,
                    e3,
                    clocks,
                    clock,
                    sigma,
                    target,
                )
                if not right:
                    continue
                common_atoms += len(common)
                right_roots += len(right)
                cells[cell] = {
                    "source": source,
                    "target": target_physical,
                    "q_pair": q_pairs[clock],
                }
                common_by_left = {piece[0]: piece for piece in common}
                right_by_left = {piece[0]: piece for piece in right}
                root_index_by_left = {
                    piece[0]: index for index, piece in enumerate(right)
                }
                require(
                    len(common_by_left) == len(common)
                    and len(right_by_left) == len(right),
                    "one labelled cell acquired colliding left endpoints",
                )
                seen = set()
                starts = tuple(sorted(
                    left for left in common_by_left
                    if left - HSTEP not in common_by_left
                ))
                for start in starts:
                    common_vertices = []
                    cursor = start
                    while cursor in common_by_left:
                        require(
                            cursor not in seen,
                            "complete-forest components overlapped",
                        )
                        seen.add(cursor)
                        common_vertices.append(common_by_left[cursor])
                        cursor += HSTEP
                    predecessor = start - HSTEP
                    if predecessor in right_by_left:
                        kind = "cofiber_rooted"
                        root_index = root_index_by_left[predecessor]
                        vertices = (
                            right_by_left[predecessor],
                            *common_vertices,
                        )
                    else:
                        kind = "common_only"
                        root_index = None
                        vertices = tuple(common_vertices)
                    require(
                        (
                            kind == "cofiber_rooted"
                            and len(vertices) >= 3
                        )
                        or (
                            kind == "common_only"
                            and len(vertices) % P == 0
                            and (len(vertices) // P) % 2 == 1
                        ),
                        "a complete-forest component changed type",
                    )
                    component_kinds[kind] += 1
                    paths.append(
                        {
                            "number": len(paths),
                            "kind": kind,
                            "cell": cell,
                            "root_index": root_index,
                            "vertices": tuple(vertices),
                        }
                    )
                require(
                    seen == set(common_by_left),
                    "complete-forest decomposition missed common atoms",
                )
    require(
        len(paths) == 685
        and component_kinds
        == Counter({"cofiber_rooted": 587, "common_only": 98})
        and common_atoms == 63_308
        and right_roots == 587
        and sum(len(path["vertices"]) for path in paths) == 63_895,
        "complete THM-2825 forest census changed",
    )
    return (
        tuple(paths),
        cells,
        full_module,
        e3,
        clocks,
        component_kinds,
    )


def frame_action_audit(path_labels, frame):
    """Audit the directed forest relation on orbit-labelled vertices."""
    label_sequences = Counter()
    path_label_sets = Counter()
    directed_labels = Counter()
    directed_gaps = Counter()
    identity_nonreflexive = 0
    generator_outdegree = Counter()
    generator_indegree = Counter()
    inverse_outdegree = Counter()
    inverse_indegree = Counter()
    generator_square_triples = 0
    generator_inverse_triangles = 0
    first_generator = None
    first_inverse = None
    first_square = None

    total_hits = 0
    paths_with_hits = 0
    component_kind_census = Counter()
    for path_number, kind, cell, root_index, hits in path_labels:
        labelled = tuple(
            (index, source_label if frame == "source" else target_label)
            for index, source_label, target_label in hits
            if (source_label if frame == "source" else target_label) is not None
        )
        if not labelled:
            continue
        paths_with_hits += 1
        component_kind_census[kind] += 1
        total_hits += len(labelled)
        path_label_sets[tuple(sorted(set(label for _index, label in labelled)))] += 1
        label_sequences[run_length_encode(labelled)] += 1

        outgoing8 = Counter()
        incoming8 = Counter()
        outgoing5 = Counter()
        incoming5 = Counter()
        by_label = defaultdict(list)
        for index, label in labelled:
            by_label[label].append(index)

        identity_nonreflexive += sum(
            len(indices) * (len(indices) - 1) // 2
            for indices in by_label.values()
        )
        for source_position, (source_index, source_label) in enumerate(labelled):
            for target_index, target_label in labelled[source_position + 1:]:
                delta = (target_label - source_label) % P
                gap = target_index - source_index
                directed_labels[delta] += 1
                directed_gaps[(delta, gap)] += 1
                if delta == 8:
                    outgoing8[source_index] += 1
                    incoming8[target_index] += 1
                    if first_generator is None or gap < first_generator[0]:
                        first_generator = (
                            gap,
                            path_number,
                            kind,
                            cell,
                            root_index,
                            source_index,
                            target_index,
                            source_label,
                            target_label,
                        )
                if delta == 5:
                    outgoing5[source_index] += 1
                    incoming5[target_index] += 1
                    if first_inverse is None or gap < first_inverse[0]:
                        first_inverse = (
                            gap,
                            path_number,
                            kind,
                            cell,
                            root_index,
                            source_index,
                            target_index,
                            source_label,
                            target_label,
                        )

        for index, _label in labelled:
            generator_outdegree[outgoing8[index]] += 1
            generator_indegree[incoming8[index]] += 1
            inverse_outdegree[outgoing5[index]] += 1
            inverse_indegree[incoming5[index]] += 1

        # Z^8 followed by Z^8 must have labels k,k+8,k+3.
        index_to_label = dict(labelled)
        for middle_index, middle_label in labelled:
            lefts = tuple(
                index
                for index in by_label[(middle_label - 8) % P]
                if index < middle_index
            )
            rights = tuple(
                index
                for index in by_label[(middle_label + 8) % P]
                if index > middle_index
            )
            count = len(lefts) * len(rights)
            generator_square_triples += count
            if count and first_square is None:
                first_square = (
                    path_number,
                    kind,
                    cell,
                    root_index,
                    lefts[0],
                    middle_index,
                    rights[0],
                    index_to_label[lefts[0]],
                    middle_label,
                    index_to_label[rights[0]],
                )

            # A forward Z^8 followed by its inverse Z^5 closes to identity.
            inverse_rights = tuple(
                index
                for index in by_label[(middle_label + 5) % P]
                if index > middle_index
            )
            generator_inverse_triangles += len(lefts) * len(inverse_rights)

    require(
        sum(generator_outdegree.values()) == total_hits
        and sum(generator_indegree.values()) == total_hits
        and sum(inverse_outdegree.values()) == total_hits
        and sum(inverse_indegree.values()) == total_hits,
        f"{frame} action degree census lost vertices",
    )
    return {
        "frame": frame,
        "total_hits": total_hits,
        "paths_with_hits": paths_with_hits,
        "component_kind_census": component_kind_census,
        "path_label_sets": path_label_sets,
        "label_sequences": label_sequences,
        "identity_reflexive": total_hits,
        "identity_nonreflexive": identity_nonreflexive,
        "directed_labels": directed_labels,
        "directed_gaps": directed_gaps,
        "generator_outdegree": generator_outdegree,
        "generator_indegree": generator_indegree,
        "inverse_outdegree": inverse_outdegree,
        "inverse_indegree": inverse_indegree,
        "generator_square_triples": generator_square_triples,
        "generator_inverse_triangles": generator_inverse_triangles,
        "first_generator": first_generator,
        "first_inverse": first_inverse,
        "first_square": first_square,
    }


def fixed_displacement_audit(
    all_intervals,
    source_orbit_intervals,
    target_orbit_intervals,
):
    """Test powers of the +66h reference on unique physical intervals."""
    interval_universe = set(all_intervals)
    frames = (
        ("source", source_orbit_intervals),
        ("target", target_orbit_intervals),
    )
    results = {}
    for frame, orbit_intervals in frames:
        rows = []
        for power in (-13, -2, -1, 0, 1, 2, 13):
            physical_steps = power * 66
            expected_delta = (power * 8) % P
            census = Counter()
            first_failure = None
            for interval, label in sorted(orbit_intervals.items()):
                shifted = (
                    interval[0] + physical_steps * HSTEP,
                    interval[1] + physical_steps * HSTEP,
                )
                if shifted not in interval_universe:
                    outcome = ("outside_complete_forest", None)
                elif shifted not in orbit_intervals:
                    outcome = ("inside_complete_forest_outside_orbit", None)
                else:
                    observed_delta = (
                        orbit_intervals[shifted] - label
                    ) % P
                    outcome = (
                        "orbit",
                        observed_delta,
                        observed_delta == expected_delta,
                    )
                census[outcome] += 1
                if (
                    first_failure is None
                    and outcome
                    != ("orbit", expected_delta, True)
                ):
                    first_failure = (
                        interval,
                        label,
                        shifted,
                        expected_delta,
                        outcome,
                    )
            rows.append(
                (
                    power,
                    physical_steps,
                    expected_delta,
                    tuple(sorted(census.items(), key=repr)),
                    first_failure,
                )
            )
        results[frame] = tuple(rows)
    return results


def decoration_relation_audit(path_labels, typed_lookup, frame):
    """Type every rooted Z8 arrow by preservation of all decorations."""
    frame_index = 0 if frame == "source" else 1
    preservation = Counter()
    preservation_by_labels = Counter()
    preserving_gap_census = Counter()
    first_preserving = None
    first_failure = None
    total = 0
    component_kind_census = Counter()
    for path_number, kind, cell, root_index, hits in path_labels:
        labelled = tuple(
            (index, row[frame_index])
            for index, *row in hits
            if row[frame_index] is not None
        )
        for source_position, (source_index, source_label) in enumerate(labelled):
            for target_index, target_label in labelled[source_position + 1:]:
                if (target_label - source_label) % P != 8:
                    continue
                source = typed_lookup[path_number, source_index]
                target = typed_lookup[path_number, target_index]
                flags = (
                    source[7] == target[7],
                    tuple(value != 0 for value in source[7])
                    == tuple(value != 0 for value in target[7]),
                    source[8] == target[8],
                    source[9] == target[9],
                    source[10] == target[10],
                )
                total += 1
                component_kind_census[kind] += 1
                preservation[flags] += 1
                preservation_by_labels[
                    (source_label, target_label, flags)
                ] += 1
                witness = (
                    path_number,
                    kind,
                    cell,
                    root_index,
                    source_index,
                    target_index,
                    target_index - source_index,
                    source_label,
                    target_label,
                    flags,
                    source[7:11],
                    target[7:11],
                )
                if all(flags):
                    preserving_gap_census[target_index - source_index] += 1
                    if (
                        first_preserving is None
                        or witness[6] < first_preserving[6]
                    ):
                        first_preserving = witness
                elif (
                    first_failure is None
                    or witness[6] < first_failure[6]
                ):
                    first_failure = witness
    return {
        "frame": frame,
        "total": total,
        "component_kind_census": component_kind_census,
        "preservation": preservation,
        "by_labels": preservation_by_labels,
        "preserving_gaps": preserving_gap_census,
        "first_preserving": first_preserving,
        "first_failure": first_failure,
    }


def main():
    (
        paths,
        cells,
        full_module,
        e3,
        clocks,
        component_kinds,
    ) = build_paths()
    base = endpoint.endpoint_mask(I)
    first_axis, second_axis, product = endpoint.cartesian_factor(base)
    require(
        product
        and tuple(sorted(first_axis))
        == (0, 1, 2, 3, 4, 5, 6, 7, 12)
        and tuple(sorted(second_axis))
        == (0, 1, 2, 3, 4, 5, 8, 9, 10)
        and sum(base) == 81,
        "distinguished q0 endpoint rectangle changed",
    )
    vertical_masks = tuple(
        endpoint.translate_mask(base, (0, exponent))
        for exponent in range(P)
    )
    require(
        len(set(vertical_masks)) == P,
        "q0 endpoint rectangle acquired a vertical stabilizer",
    )
    vertical_lookup = {
        mask: exponent for exponent, mask in enumerate(vertical_masks)
    }

    source_census = Counter()
    target_census = Counter()
    pair_census = Counter()
    rooted_source_census = Counter()
    rooted_target_census = Counter()
    rooted_pair_census = Counter()
    path_labels = []
    raw_occurrences = []
    all_intervals = set()
    source_orbit_intervals = {}
    target_orbit_intervals = {}
    unique_source_masks = set()
    unique_target_masks = set()

    for path in paths:
        hits = []
        for index, piece in enumerate(path["vertices"]):
            interval = piece[:2]
            target_interval = (
                interval[0] + SHIFT,
                interval[1] + SHIFT,
            )
            source_mask = endpoint.endpoint_mask(interval)
            target_mask = endpoint.endpoint_mask(target_interval)
            unique_source_masks.add(source_mask)
            unique_target_masks.add(target_mask)
            all_intervals.add(interval)
            source_label = vertical_lookup.get(source_mask)
            target_label = vertical_lookup.get(target_mask)
            if source_label is not None:
                source_census[source_label] += 1
                if path["kind"] == "cofiber_rooted":
                    rooted_source_census[source_label] += 1
                prior = source_orbit_intervals.setdefault(
                    interval, source_label
                )
                require(
                    prior == source_label,
                    "one physical interval acquired two source labels",
                )
            if target_label is not None:
                target_census[target_label] += 1
                if path["kind"] == "cofiber_rooted":
                    rooted_target_census[target_label] += 1
                prior = target_orbit_intervals.setdefault(
                    interval, target_label
                )
                require(
                    prior == target_label,
                    "one physical interval acquired two target labels",
                )
            if source_label is None and target_label is None:
                continue
            pair_census[(source_label, target_label)] += 1
            if path["kind"] == "cofiber_rooted":
                rooted_pair_census[(source_label, target_label)] += 1
            hits.append((index, source_label, target_label))
            raw_occurrences.append(
                (
                    path["number"],
                    path["kind"],
                    path["cell"],
                    path["root_index"],
                    index,
                    piece,
                    source_label,
                    target_label,
                )
            )
        path_labels.append(
            (
                path["number"],
                path["kind"],
                path["cell"],
                path["root_index"],
                tuple(hits),
            )
        )

    # Complete typing is intentionally deferred until orbit hits are known.
    u_events, v_events = independent.ancestry_setup()
    supports_by_cell = {
        cell: physical.carrier_supports(
            cells[cell]["source"], cells[cell]["target"]
        )
        for cell in {row[2] for row in raw_occurrences}
    }
    semantic_caches = tuple({} for _clock in range(7))
    typed_rows = []
    typing_census = Counter()
    for (
        path_number,
        kind,
        cell,
        root_index,
        index,
        piece,
        source_label,
        target_label,
    ) in raw_occurrences:
        interval = piece[:2]
        target_interval = (
            interval[0] + SHIFT,
            interval[1] + SHIFT,
        )
        semantic = collar.semantic_value(
            piece,
            cells[cell]["q_pair"],
            semantic_caches[cell[0]],
        )
        factors = (
            factor_signature(interval, full_module, e3, clocks, cell),
            factor_signature(
                target_interval, full_module, e3, clocks, cell
            ),
        )
        carrier = compact_carrier(
            physical.carrier_mask(piece, supports_by_cell[cell])
        )
        ancestry = independent.ancestry_chamber(
            piece, u_events, v_events
        )
        row = (
            path_number,
            cell,
            root_index,
            index,
            interval,
            source_label,
            target_label,
            semantic,
            factors,
            carrier,
            ancestry,
            kind,
        )
        typed_rows.append(row)
        typing_census[
            (
                source_label,
                target_label,
                tuple(value != 0 for value in semantic),
                factors,
                carrier,
            )
        ] += 1

    source_action = frame_action_audit(tuple(path_labels), "source")
    target_action = frame_action_audit(tuple(path_labels), "target")
    rooted_path_labels = tuple(
        row for row in path_labels if row[1] == "cofiber_rooted"
    )
    rooted_source_action = frame_action_audit(
        rooted_path_labels, "source"
    )
    rooted_target_action = frame_action_audit(
        rooted_path_labels, "target"
    )
    displacement = fixed_displacement_audit(
        all_intervals,
        source_orbit_intervals,
        target_orbit_intervals,
    )

    present_source_labels = tuple(sorted(source_census))
    present_target_labels = tuple(sorted(target_census))
    missing_source_labels = tuple(
        exponent for exponent in range(P) if exponent not in source_census
    )
    missing_target_labels = tuple(
        exponent for exponent in range(P) if exponent not in target_census
    )
    typed_rows = tuple(typed_rows)
    typed_lookup = {
        (row[0], row[3]): row for row in typed_rows
    }
    require(
        len(typed_lookup) == len(typed_rows),
        "one labelled path vertex was typed twice",
    )
    source_decorations = decoration_relation_audit(
        tuple(path_labels), typed_lookup, "source"
    )
    target_decorations = decoration_relation_audit(
        tuple(path_labels), typed_lookup, "target"
    )
    typed_digest = sha256(repr(typed_rows).encode("ascii")).hexdigest()
    physical_digest = sha256(
        repr(
            (
                tuple(sorted(source_orbit_intervals.items())),
                tuple(sorted(target_orbit_intervals.items())),
            )
        ).encode("ascii")
    ).hexdigest()
    unique_source_census = Counter(source_orbit_intervals.values())
    unique_target_census = Counter(target_orbit_intervals.values())
    ancestry_census = Counter(row[10] for row in typed_rows)
    first_typed_by_pair = {}
    for row in typed_rows:
        first_typed_by_pair.setdefault((row[5], row[6]), row)

    require(
        present_source_labels == present_target_labels == (0, 4, 8, 9)
        and missing_source_labels
        == missing_target_labels
        == (1, 2, 3, 5, 6, 7, 10, 11, 12)
        and source_census == rooted_source_census
        and target_census == rooted_target_census
        and pair_census == rooted_pair_census
        and unique_source_census
        == unique_target_census
        == Counter({0: 12, 4: 23, 8: 23, 9: 23})
        and len(all_intervals) == 1_275
        and len(unique_source_masks) == 77
        and len(unique_target_masks) == 72
        and not any(row[11] == "common_only" for row in typed_rows),
        "complete-forest four-label orbit boundary changed",
    )
    require(
        source_action["directed_labels"]
        == target_action["directed_labels"]
        == Counter({0: 46_079, 8: 40_825})
        and source_action["first_inverse"] is None
        and target_action["first_inverse"] is None
        and source_action["generator_square_triples"] == 0
        and target_action["generator_square_triples"] == 0
        and source_action["generator_inverse_triangles"] == 0
        and target_action["generator_inverse_triangles"] == 0
        and source_action["component_kind_census"]
        == Counter({"cofiber_rooted": 148})
        and target_action["component_kind_census"]
        == Counter({"cofiber_rooted": 139}),
        "complete-forest partial-Z8 action obstruction changed",
    )
    require(
        rooted_source_action["directed_labels"]
        == rooted_target_action["directed_labels"]
        == Counter({0: 46_079, 8: 40_825})
        and rooted_source_action["first_inverse"] is None
        and rooted_target_action["first_inverse"] is None
        and rooted_source_action["generator_square_triples"] == 0
        and rooted_target_action["generator_square_triples"] == 0
        and rooted_source_action["generator_inverse_triangles"] == 0
        and rooted_target_action["generator_inverse_triangles"] == 0,
        "rooted partial-Z8 action obstruction changed",
    )
    expected_square_failure = (
        (("inside_complete_forest_outside_orbit", None), 81),
    )
    plus66_source = dict(next(
        row[3] for row in displacement["source"] if row[0] == 1
    ))
    plus66_target = dict(next(
        row[3] for row in displacement["target"] if row[0] == 1
    ))
    require(
        plus66_source.get(("orbit", 8, True)) == 22
        and plus66_target.get(("orbit", 8, True)) == 22
        and sum(plus66_source.values()) == 81
        and sum(plus66_target.values()) == 81
        and next(
            row[3] for row in displacement["source"] if row[0] == 2
        )
        == expected_square_failure
        and next(
            row[3] for row in displacement["target"] if row[0] == 2
        )
        == expected_square_failure,
        "fixed +132h Z3 obstruction changed",
    )

    print("THM-2859 COMPLETE ENDPOINT VERTICAL-ORBIT / ACTION AUDIT")
    print(
        f"components={len(paths)};"
        f"component_kinds={tuple(sorted(component_kinds.items()))};"
        f"path_vertices={sum(len(path['vertices']) for path in paths)};"
        f"unique_physical_intervals={len(all_intervals)};"
        f"source_masks={len(unique_source_masks)};"
        f"target_masks={len(unique_target_masks)};"
        f"endpoint_masks_evaluated={len(endpoint.mask_cache)}"
    )
    print(
        "base_rectangle="
        f"first_axis={tuple(sorted(first_axis))};"
        f"second_axis={tuple(sorted(second_axis))};mass={sum(base)};"
        "vertical_orbit_size=13"
    )
    print(
        "source_vertical_orbit="
        f"present={present_source_labels};missing={missing_source_labels};"
        f"labelled_census={tuple(sorted(source_census.items()))};"
        f"unique_interval_census="
        f"{tuple(sorted(unique_source_census.items()))};"
        f"unique_intervals={len(source_orbit_intervals)}"
    )
    print(
        "target_vertical_orbit="
        f"present={present_target_labels};missing={missing_target_labels};"
        f"labelled_census={tuple(sorted(target_census.items()))};"
        f"unique_interval_census="
        f"{tuple(sorted(unique_target_census.items()))};"
        f"unique_intervals={len(target_orbit_intervals)}"
    )
    print(
        "cofiber_rooted_subset="
        "components=587;vertices=55341;"
        f"source_labels={tuple(sorted(rooted_source_census.items()))};"
        f"target_labels={tuple(sorted(rooted_target_census.items()))};"
        f"pair_labels={tuple(sorted(rooted_pair_census.items(), key=repr))};"
        "common_only_orbit_hits=0;"
        "complete_relation_equals_rooted=1"
    )
    print(
        "two_frame_occurrences="
        f"count={len(typed_rows)};"
        f"label_pairs={tuple(sorted(pair_census.items(), key=repr))};"
        f"typed_census={tuple(sorted(typing_census.items(), key=repr))};"
        f"ancestry_census={tuple(sorted(ancestry_census.items()))};"
        f"representatives="
        f"{tuple(sorted(first_typed_by_pair.items(), key=repr))};"
        f"typed_rows_sha256={typed_digest};"
        f"unique_physical_orbit_sha256={physical_digest}"
    )
    for scope, action in (
        ("complete", source_action),
        ("complete", target_action),
    ):
        gap_rows = tuple(sorted(action["directed_gaps"].items()))
        gap_digest = sha256(
            repr(gap_rows).encode("ascii")
        ).hexdigest()
        nonzero_gaps = tuple(
            gap for (delta, gap), count in gap_rows if delta and count
        )
        print(
            f"{scope}_{action['frame']}_directed_relation="
            f"hits={action['total_hits']};"
            f"paths={action['paths_with_hits']};"
            f"component_kinds="
            f"{tuple(sorted(action['component_kind_census'].items()))};"
            f"path_label_sets="
            f"{tuple(sorted(action['path_label_sets'].items(), key=repr))};"
            f"sequence_census="
            f"{tuple(sorted(action['label_sequences'].items(), key=repr))};"
            f"identity_reflexive={action['identity_reflexive']};"
            f"identity_nonreflexive={action['identity_nonreflexive']};"
            f"directed_labels="
            f"{tuple(sorted(action['directed_labels'].items()))};"
            f"directed_gap_sha256={gap_digest};"
            f"Z8_gap_range={(min(nonzero_gaps), max(nonzero_gaps))};"
            f"Z8_gap66={action['directed_gaps'][(8, 66)]};"
            f"Z8_outdegree="
            f"{tuple(sorted(action['generator_outdegree'].items()))};"
            f"Z8_indegree="
            f"{tuple(sorted(action['generator_indegree'].items()))};"
            f"Z5_inverse_outdegree="
            f"{tuple(sorted(action['inverse_outdegree'].items()))};"
            f"Z5_inverse_indegree="
            f"{tuple(sorted(action['inverse_indegree'].items()))};"
            f"Z8_Z8_square_triples={action['generator_square_triples']};"
            f"Z8_Z5_identity_triangles="
            f"{action['generator_inverse_triangles']};"
            f"first_Z8={action['first_generator']};"
            f"first_Z5={action['first_inverse']};"
            f"first_Z8_square={action['first_square']}"
        )
    for decorations in (source_decorations, target_decorations):
        print(
            f"{decorations['frame']}_Z8_decoration_preservation="
            f"arrows={decorations['total']};"
            f"component_kinds="
            f"{tuple(sorted(decorations['component_kind_census'].items()))};"
            f"flags=(exact_semantic,live_pattern,native_factor_signatures,"
            "carrier,ancestry);"
            f"census="
            f"{tuple(sorted(decorations['preservation'].items()))};"
            f"by_labels="
            f"{tuple(sorted(decorations['by_labels'].items()))};"
            f"all_preserving_gaps="
            f"{tuple(sorted(decorations['preserving_gaps'].items()))};"
            f"first_all_preserving={decorations['first_preserving']};"
            f"first_failure={decorations['first_failure']}"
        )
    print(
        "fixed_plus66h_powers="
        f"{tuple(sorted(displacement.items()))}"
    )
    print(
        "TYPE_VERDICT="
        "the complete forest realizes exactly the labels {0,4,8,9}; "
        "its only nonidentity directed label is Z8, through the two "
        "disconnected arrow types 0->8 and 9->4.  The q0 mask has trivial "
        "vertical stabilizer, so nine of thirteen orbit objects are absent.  "
        "Orbit labels record exact equality of endpoint masks only.  "
        "A physical C13 action additionally requires a total single-valued "
        "Z8 generator, a directed Z5 inverse, a composable Z8-square "
        "producing Z3, and thirteenfold closure; the directed forest "
        "relation and fixed +66h "
        "power audit test those obligations separately."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
