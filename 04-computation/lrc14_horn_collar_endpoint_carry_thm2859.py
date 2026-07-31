#!/usr/bin/env python3
"""Exact endpoint/carry audit for the THM-2825/2847/2851 hinge.

Canonical THM-2859 companion.  This script verifies the distinguished horn root whose
even collar is the fixed source atom I, then compares the nonsplit C_169 carry
with the split F_13^2 endpoint-address translation law.

No executable Python ``assert`` statement is used.
"""

from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import ast
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
    "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    "lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.py":
        "af924103c30a45eeab88520d6b569c00a5a8c9cafbc8052b3585dd86e07b5dac",
    "lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825.py":
        "5a74cfa0b9a5f30e15460e7e1bee17214193d1f307e79f1eeb68b8781b7d23d6",
}
THM2857_CANDIDATE = (
    "lrc14_endpoint_galois_carry_torsor_thm2857.py",
    "0bae59c9b1460f37e1879a81746154593cb0699ee13b3e5e800ba0af95ea5e4c",
)
for name, expected in PINNED.items():
    payload = (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"pinned dependency changed: {name}",
    )
candidate_payload = (COMP / THM2857_CANDIDATE[0]).read_bytes().replace(
    b"\r\n", b"\n"
)
require(
    sha256(candidate_payload).hexdigest() == THM2857_CANDIDATE[1],
    "THM-2857 comparison fixture changed",
)


import lrc14_nearest_half_step_common_right_collar_thm2825 as collar
import lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825 as physical
import lrc14_nearest_half_step_common_right_collar_independent_audit_thm2825 as independent


P = 13
H = collar.H
T = collar.copies.T
LENGTH = collar.copies.LENGTH
SHIFT = collar.copies.SHIFT
I = (142004992589460, 142005019034340)
R = (I[0] - 2 * H, I[1] - 2 * H)
M1 = (I[0] - H, I[1] - H)
M2 = I

HORN_20 = {
    (1, s, t)
    for s in (0, 3, 8, 9, 12)
    for t in (5, 6, 9, 10)
}
COMMON_42 = {
    (1, s, t)
    for s in (0, 3, 8, 9, 10, 11, 12)
    for t in (5, 6, 7, 8, 9, 10)
}


# Direct endpoint evaluator.  This avoids materializing the 169 large present
# sets and evaluates the defining periodic inequalities on one small interval.
ENDPOINT = collar.copies.endpoint_base
KEYS = tuple(ENDPOINT.KEYS)
KEY_INDEX = {key: index for index, key in enumerate(KEYS)}
W = tuple(ENDPOINT.endpoint.W)
PATTERN = dict(ENDPOINT.PAT_E3)

require(
    KEYS == tuple((a, b) for a in range(P) for b in range(P)),
    "endpoint enumeration changed",
)


def periodic_parameters(index, ell_value):
    if PATTERN[index] == "gout":
        denominator = 91
        low = -13 - 7 * ell_value
        high = 13 - 7 * ell_value
    else:
        denominator = 182
        low = -13 - 14 * ell_value
        high = 13 - 14 * ell_value
    unit = T // (denominator * W[index])
    return (
        denominator * unit,
        (low % denominator) * unit,
        (high - low) * unit,
    )


def nearby_boundaries(left, right, period, base, width):
    first = (left - base) // period - 2
    last = (right - base) // period + 2
    rows = []
    for number in range(first, last + 1):
        start = base + number * period
        stop = start + width
        if left < start < right:
            rows.append(start)
        if left < stop < right:
            rows.append(stop)
    return tuple(rows)


def midpoint_in_comb(left, right, period, base, width):
    twice_x = left + right
    twice_period = 2 * period
    twice_base = 2 * base
    quotient = (twice_x - twice_base) // twice_period
    residue = twice_x - twice_base - quotient * twice_period
    return 0 <= residue < 2 * width


def address_indicator(interval, address):
    left, right = interval
    ell = ENDPOINT.REPS[address]
    factor_data = tuple(
        periodic_parameters(index, ell[index]) for index in range(len(W))
    )
    cuts = {left, right}
    for period, base, width in factor_data:
        cuts.update(nearby_boundaries(left, right, period, base, width))
    cuts = tuple(sorted(cuts))
    statuses = []
    for start, stop in zip(cuts, cuts[1:]):
        period, base, width = factor_data[8]
        live = midpoint_in_comb(start, stop, period, base, width)
        if live:
            for index in range(8):
                period, base, width = factor_data[index]
                if midpoint_in_comb(start, stop, period, base, width):
                    live = False
                    break
        statuses.append(live)
    if all(statuses):
        return 1
    if not any(statuses):
        return 0
    return 2


MASK_CACHE = {}


def endpoint_mask(interval):
    if interval not in MASK_CACHE:
        row = tuple(address_indicator(interval, address) for address in KEYS)
        require(2 not in row, f"endpoint mask cuts {interval}")
        MASK_CACHE[interval] = row
    return MASK_CACHE[interval]


def translate_mask(mask, delta):
    da, db = delta
    return tuple(
        mask[KEY_INDEX[((a + da) % P, (b + db) % P)]]
        for a, b in KEYS
    )


def translation_deltas(source, target):
    return tuple(
        delta
        for delta in KEYS
        if translate_mask(source, delta) == target
    )


def support(mask):
    return frozenset(
        address for address, value in zip(KEYS, mask) if value
    )


def cartesian_description(mask):
    points = support(mask)
    if not points:
        return (), (), 0, True
    first = frozenset(a for a, _b in points)
    second = frozenset(b for _a, b in points)
    product = frozenset(
        (a, b) for a in first for b in second
    )
    return tuple(sorted(first)), tuple(sorted(second)), len(points), points == product


def rank_mod(matrix, prime):
    rows = [list(value % prime for value in row) for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(rows))
             if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for index in range(len(rows)):
            if index == rank or not rows[index][column]:
                continue
            scalar = rows[index][column]
            rows[index] = [
                (value - scalar * pivot_value) % prime
                for value, pivot_value in zip(rows[index], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def convolution_matrix(mask):
    values = {key: mask[index] for index, key in enumerate(KEYS)}
    return tuple(
        tuple(
            values[((a - c) % P, (b - d) % P)]
            for c, d in KEYS
        )
        for a, b in KEYS
    )


def cyclic_matrix(values):
    return tuple(
        tuple(values[(row - column) % P] for column in range(P))
        for row in range(P)
    )


def safe_comb_contains(interval, module, speed, denominator, lo, hi):
    require(
        module.T % (denominator * speed) == 0,
        "comb grid ceased to resolve",
    )
    unit = module.T // (denominator * speed)
    step = denominator * unit
    width = (hi - lo) * unit
    base = (lo % denominator) * unit
    left, right = interval
    first = (left - base) // step
    for number in range(first - 1, first + 3):
        danger_left = base + number * step
        danger_right = danger_left + width
        if max(left, danger_left) < min(right, danger_right):
            return False
    return True


def factor_signature(interval, full_module, e3, clocks, clock, s, t):
    return (
        collar.copies.contained(interval, e3),
        collar.copies.contained(interval, clocks[clock]),
        safe_comb_contains(
            interval, full_module, full_module.W[1], 182,
            -14 * s - 13, -14 * s + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.C2, 182,
            14 * s - 13, 14 * s + 13,
        ),
        safe_comb_contains(
            interval, full_module, full_module.C3, 182,
            14 * t - 13, 14 * t + 13,
        ),
    )


def endpoint_address_carry_audit():
    source_masks = {
        "R": endpoint_mask(R),
        "M1": endpoint_mask(M1),
        "M2": endpoint_mask(M2),
    }
    target_masks = {
        name: endpoint_mask(
            (interval[0] + SHIFT, interval[1] + SHIFT)
        )
        for name, interval in (("R", R), ("M1", M1), ("M2", M2))
    }
    require(
        not any(source_masks["R"])
        and source_masks["M1"] == source_masks["M2"]
        and target_masks["R"] == target_masks["M1"] == target_masks["M2"]
        == source_masks["M2"],
        "distinguished endpoint ladder changed",
    )

    base = source_masks["M2"]
    base_description = cartesian_description(base)
    unit = T // P

    def shifted_atom(q, extra=0):
        left = (I[0] + q * unit + extra) % T
        require(
            left + LENGTH <= T,
            "tested pulled allocation atom crossed the seam",
        )
        return left, left + LENGTH

    allocation_masks = {
        q: (
            endpoint_mask(shifted_atom(q)),
            endpoint_mask(shifted_atom(q, SHIFT)),
        )
        for q in (0, 3, 7, 11)
    }
    allocation_descriptions = {
        q: tuple(cartesian_description(mask) for mask in masks)
        for q, masks in allocation_masks.items()
    }
    expected_first = (0, 1, 2, 3, 4, 5, 6, 7, 12)
    expected_second = (0, 1, 2, 3, 4, 5, 8, 9, 10)
    require(
        base_description
        == (expected_first, expected_second, 81, True),
        "selected endpoint rectangle changed",
    )
    require(
        allocation_descriptions[0] == (base_description, base_description)
        and allocation_descriptions[3][0][2]
        == allocation_descriptions[3][1][2] == 90
        and allocation_descriptions[7][0][2]
        == allocation_descriptions[7][1][2] == 0
        and allocation_descriptions[11][0][2]
        == allocation_descriptions[11][1][2] == 81
        and not translation_deltas(
            allocation_masks[0][0], allocation_masks[11][0]
        )
        and not translation_deltas(
            allocation_masks[0][1], allocation_masks[11][1]
        ),
        "q0/q3/q7/q11 endpoint allocation boundary changed",
    )

    z_mask = translate_mask(base, (0, 1))
    z8_mask = translate_mask(base, (0, 8))
    derivative = tuple(b - a for a, b in zip(base, z_mask))
    derivative8 = tuple(b - a for a, b in zip(base, z8_mask))
    require(
        translation_deltas(base, z_mask) == ((0, 1),)
        and translation_deltas(base, z8_mask) == ((0, 8),)
        and Counter(derivative) == Counter({0: 133, 1: 18, -1: 18}),
        "normalized central endpoint derivative changed",
    )
    require(
        Counter(derivative8) == Counter({0: 133, 1: 18, -1: 18}),
        "physical Z8 endpoint derivative changed",
    )

    prime = 1_000_003
    base_rank = rank_mod(convolution_matrix(base), prime)
    derivative_rank = rank_mod(convolution_matrix(derivative), prime)
    derivative8_rank = rank_mod(convolution_matrix(derivative8), prime)
    allocation_ranks = {
        role: {
            q: rank_mod(
                convolution_matrix(allocation_masks[q][role]), prime
            )
            for q in (0, 3, 7, 11)
        }
        for role in (0, 1)
    }
    allocation_difference_ranks = {}
    for role in (0, 1):
        masks = {
            q: allocation_masks[q][role] for q in (0, 3, 7, 11)
        }
        for source_q, target_q in ((0, 3), (3, 11), (11, 7), (0, 11)):
            difference = tuple(
                masks[target_q][index] - masks[source_q][index]
                for index in range(P * P)
            )
            allocation_difference_ranks[
                (role, source_q, target_q)
            ] = (
                rank_mod(convolution_matrix(difference), prime),
                sum(difference),
            )
    expected_allocation_ranks = {0: 169, 3: 169, 7: 0, 11: 169}
    expected_difference = {
        (0, 3): (169, 9),
        (3, 11): (169, -9),
        (11, 7): (169, -81),
        (0, 11): (168, 0),
    }
    require(
        allocation_ranks
        == {0: expected_allocation_ranks, 1: expected_allocation_ranks}
        and all(
            allocation_difference_ranks[(role, *edge)] == expected
            for role in (0, 1)
            for edge, expected in expected_difference.items()
        ),
        "allocation endpoint convolution ranks changed",
    )
    vertical = tuple(
        sum(base[KEY_INDEX[(a, b)]] for a in range(P))
        for b in range(P)
    )
    vertical_derivative = tuple(
        vertical[(b + 1) % P] - vertical[b] for b in range(P)
    )
    vertical_rank = rank_mod(cyclic_matrix(vertical_derivative), prime)
    require(
        (
            base_rank,
            derivative_rank,
            derivative8_rank,
            vertical_rank,
        ) == (169, 156, 156, 12),
        "endpoint carry ranks changed",
    )

    # The C_169 digit section and the split address plane are the same set,
    # but not the same additive group.  The high carry T is the central Z.
    def iota(point):
        return point[0] + P * point[1]

    def y(point):
        return ((point[0] + 1) % P, point[1])

    def z(point):
        return (point[0], (point[1] + 1) % P)

    digit_rows = []
    for point in KEYS:
        successor = (iota(point) + 1) % (P * P)
        split_successor = y(point) if point[0] != P - 1 else z(y(point))
        require(
            successor == iota(split_successor),
            "cyclic successor/digit carry comparison changed",
        )
        require(
            (iota(point) + P) % (P * P) == iota(z(point)),
            "high carry ceased to be vertical endpoint translation",
        )
        digit_rows.append((point, split_successor))
    require(
        len({iota(point) for point in KEYS}) == P * P,
        "digit section ceased to be bijective",
    )

    # Search the complete physical interval ladder available to every
    # distinguished horn root.  No normalized transverse endpoint direction
    # (a,1) occurs.  A later unnormalized (0,8) translate occurs only on long
    # paths and is retained as a hostile control rather than suppressed.
    path_masks = {
        step: endpoint_mask(
            (R[0] + step * H, R[1] + step * H)
        )
        for step in range(145)
    }
    translated_steps = {
        step: translation_deltas(base, mask)
        for step, mask in path_masks.items()
        if translation_deltas(base, mask)
    }
    normalized_hits = tuple(
        (step, deltas)
        for step, deltas in translated_steps.items()
        if any(delta[1] == 1 for delta in deltas)
    )
    require(
        not normalized_hits
        and tuple(
            step for step, deltas in translated_steps.items()
            if deltas == ((0, 0),)
        ) == tuple(range(1, 13))
        and tuple(
            step for step, deltas in translated_steps.items()
            if deltas == ((0, 8),)
        ) == tuple(range(68, 91)),
        "endpoint-translate hostile path profile changed",
    )

    return {
        "source_descriptions": {
            name: cartesian_description(mask)
            for name, mask in source_masks.items()
        },
        "target_description": cartesian_description(target_masks["R"]),
        "base_description": base_description,
        "z_hamming": sum(a != b for a, b in zip(base, z_mask)),
        "derivative_counts": Counter(derivative),
        "ranks": (
            base_rank,
            derivative_rank,
            derivative8_rank,
            vertical_rank,
        ),
        "translated_steps": translated_steps,
        "normalized_hits": normalized_hits,
        "mask_cache": len(MASK_CACHE),
        "allocation_descriptions": allocation_descriptions,
        "q0_q11_translation_deltas": (
            translation_deltas(
                allocation_masks[0][0], allocation_masks[11][0]
            ),
            translation_deltas(
                allocation_masks[0][1], allocation_masks[11][1]
            ),
        ),
        "allocation_ranks": allocation_ranks,
        "allocation_difference_ranks": allocation_difference_ranks,
    }


def common42_path_geometry():
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

    cache = {}
    paths = []
    horn_rows = []
    horn_common_atoms = 0
    for cell in sorted(COMMON_42):
        clock, s, t = cell
        source, target, common, right = collar.cell_objects(
            details, full_module, e3, clocks, clock, s, t
        )
        require(right, f"common42 cell {cell} left the collar bank")
        if cell in HORN_20:
            horn_common_atoms += len(common)
        common_by_left = {piece[0]: piece for piece in common}
        distinguished = []
        for root_index, root in enumerate(right):
            vertices = [root]
            cursor = root[0] + H
            while cursor in common_by_left:
                vertices.append(common_by_left[cursor])
                cursor += H
            require(len(vertices) >= 3, "root lost collar ladder")
            is_distinguished = vertices[2][:2] == I
            if is_distinguished:
                distinguished.append(root_index)
            paths.append({
                "cell": cell,
                "root_index": root_index,
                "vertices": tuple(vertices),
                "source": source,
                "target": target,
                "full_module": full_module,
                "e3": e3,
                "clocks": clocks,
                "q_pair": q_pairs[clock],
                "distinguished": is_distinguished,
            })
        require(
            distinguished == [0],
            f"cell {cell} lost its unique index-zero I root",
        )
        if cell in HORN_20:
            record = next(
                row for row in paths
                if row["cell"] == cell and row["distinguished"]
            )
            vertices = record["vertices"]
            values = tuple(
                collar.semantic_value(piece, q_pairs[clock], cache) != (0, 0)
                for piece in vertices[:3]
            )
            factor_rows = tuple(
                factor_signature(
                    piece[:2], full_module, e3, clocks, clock, s, t
                )
                for piece in vertices[:3]
            )
            require(
                vertices[0][:2] == R
                and vertices[1][:2] == M1
                and vertices[2][:2] == M2
                and values == (True, False, True)
                and factor_rows[0]
                == (False, True, True, True, True, True)
                and factor_rows[1] == factor_rows[2] == (True,) * 6
                and len(vertices) - 1 >= 14,
                f"distinguished ladder changed in {cell}",
            )
            horn_rows.append(
                (cell, record["root_index"], len(vertices) - 1, factor_rows)
            )

    require(
        len(paths) == 114
        and sum(row["cell"] in HORN_20 for row in paths) == 52
        and sum(row["distinguished"] for row in paths) == 42
        and horn_common_atoms == 4076
        and len(horn_rows) == 20
        and {row[1] for row in horn_rows} == {0}
        and Counter(row[2] for row in horn_rows)
        == Counter({14: 4, 40: 4, 118: 4, 144: 8}),
        "common42/distinguished horn census changed",
    )
    return tuple(paths), tuple(horn_rows), horn_common_atoms


def endpoint_leg_counts(paths):
    """Count exact h=8,9,4 endpoint translations on directed path pairs."""
    scopes = {
        "common42_all": paths,
        "common42_distinguished": tuple(
            row for row in paths if row["distinguished"]
        ),
        "horn20_all": tuple(
            row for row in paths if row["cell"] in HORN_20
        ),
        "horn20_distinguished": tuple(
            row for row in paths
            if row["cell"] in HORN_20 and row["distinguished"]
        ),
    }
    results = {}
    for scope_name, records in scopes.items():
        leg_counts = Counter()
        leg_witnesses = {}
        leg_paths = defaultdict(set)
        triangle_counts = Counter()
        triangle_witnesses = {}
        for record in records:
            masks = tuple(
                endpoint_mask(piece[:2]) for piece in record["vertices"]
            )
            index_by_mask = {}
            for index, mask in enumerate(masks):
                index_by_mask.setdefault(mask, []).append(index)
            for axis in (0, 1):
                deltas = {
                    h: (h, 0) if axis == 0 else (0, h)
                    for h in (8, 9, 4)
                }
                for h, delta in deltas.items():
                    for source_index, source_mask in enumerate(masks):
                        if not any(source_mask):
                            continue
                        target_mask = translate_mask(source_mask, delta)
                        for target_index in index_by_mask.get(target_mask, ()):
                            if target_index <= source_index:
                                continue
                            key = (axis, h)
                            leg_counts[key] += 1
                            leg_paths[key].add(
                                (record["cell"], record["root_index"])
                            )
                            leg_witnesses.setdefault(
                                key,
                                (
                                    record["cell"],
                                    record["root_index"],
                                    source_index,
                                    target_index,
                                    delta,
                                ),
                            )
                for source_index, source_mask in enumerate(masks):
                    if not any(source_mask):
                        continue
                    after8 = translate_mask(source_mask, deltas[8])
                    direct4 = translate_mask(source_mask, deltas[4])
                    for middle_index in index_by_mask.get(after8, ()):
                        if middle_index <= source_index:
                            continue
                        after9 = translate_mask(masks[middle_index], deltas[9])
                        if after9 != direct4:
                            continue
                        for target_index in index_by_mask.get(after9, ()):
                            if target_index <= middle_index:
                                continue
                            triangle_counts[axis] += 1
                            triangle_witnesses.setdefault(
                                axis,
                                (
                                    record["cell"],
                                    record["root_index"],
                                    source_index,
                                    middle_index,
                                    target_index,
                                ),
                            )
        results[scope_name] = {
            "paths": len(records),
            "leg_counts": leg_counts,
            "leg_witnesses": leg_witnesses,
            "leg_path_counts": {
                key: len(value) for key, value in leg_paths.items()
            },
            "leg_cells": {
                key: tuple(sorted({path[0] for path in value}))
                for key, value in leg_paths.items()
            },
            "triangle_counts": triangle_counts,
            "triangle_witnesses": triangle_witnesses,
        }
    return results


def step68_physical_audit(paths):
    candidates = tuple(
        row for row in paths
        if (
            row["cell"] in HORN_20
            and row["distinguished"]
            and len(row["vertices"]) > 68
        )
    )
    expected_cells = {
        (1, s, t)
        for s in (0, 3, 12)
        for t in (5, 6, 9, 10)
    }
    require(
        len(candidates) == 12
        and {row["cell"] for row in candidates} == expected_cells,
        "step-68 horn scope changed",
    )

    u_events, v_events = independent.ancestry_setup()
    ancestry = independent.ancestry
    e_intervals = tuple(
        ancestry.base.build_set(ancestry.base.PAT_E, ancestry.base.ZELL)
    )
    q_intervals = tuple(
        ancestry.base.build_set(ancestry.host.PAT_QB, ancestry.base.ZELL)
    )
    q_depth, q_starts = ancestry.scaled_intervals(
        q_intervals, ancestry.DEPTH
    )
    e_depth_pack, e_depth_pack_starts = ancestry.scaled_intervals(
        e_intervals, ancestry.DEPTH * ancestry.PACK
    )
    e_depth, e_depth_starts = ancestry.scaled_intervals(
        e_intervals, ancestry.DEPTH
    )
    ancestry_arguments = (
        q_depth,
        q_starts,
        e_depth_pack,
        e_depth_pack_starts,
        e_depth,
        e_depth_starts,
    )
    representative_source = candidates[0]["vertices"][2]
    representative_target = candidates[0]["vertices"][68]
    ancestry_coordinate_source = (
        representative_source[0] + representative_source[1]
    ) // 2
    ancestry_coordinate_target = (
        representative_target[0] + representative_target[1]
    ) // 2
    literal_source = ancestry.contributor_sets(
        ancestry_coordinate_source, *ancestry_arguments
    )
    literal_target = ancestry.contributor_sets(
        ancestry_coordinate_target, *ancestry_arguments
    )
    require(
        literal_source == literal_target
        and tuple(map(len, literal_source)) == (966606, 28534)
        and (
            ancestry.SUPPLIED_PATH[0] * ancestry.PACK
            + ancestry.SUPPLIED_PATH[1]
        ) in literal_source[0]
        and ancestry.SUPPLIED_PATH[2] in literal_source[1],
        "literal ancestry sets changed across step2-to-step68",
    )
    literal_ancestry_digest = ancestry.path_digest(*literal_source)
    source_endpoint_deltas = Counter()
    target_endpoint_deltas = Counter()
    factor_pairs = Counter()
    carrier_pairs = Counter()
    semantic_pairs = Counter()
    ancestry_pairs = Counter()
    rows = []

    for record in candidates:
        cell = record["cell"]
        clock, s, t = cell
        source_piece = record["vertices"][2]
        target_piece = record["vertices"][68]
        source_endpoint = endpoint_mask(source_piece[:2])
        target_endpoint = endpoint_mask(target_piece[:2])
        source_target_endpoint = endpoint_mask(
            (
                source_piece[0] + SHIFT,
                source_piece[1] + SHIFT,
            )
        )
        target_target_endpoint = endpoint_mask(
            (
                target_piece[0] + SHIFT,
                target_piece[1] + SHIFT,
            )
        )
        source_delta = translation_deltas(source_endpoint, target_endpoint)
        target_delta = translation_deltas(
            source_target_endpoint, target_target_endpoint
        )
        source_endpoint_deltas[source_delta] += 1
        target_endpoint_deltas[target_delta] += 1

        factor_pair = tuple(
            (
                factor_signature(
                    piece[:2],
                    record["full_module"],
                    record["e3"],
                    record["clocks"],
                    clock,
                    s,
                    t,
                ),
                factor_signature(
                    (
                        piece[0] + SHIFT,
                        piece[1] + SHIFT,
                    ),
                    record["full_module"],
                    record["e3"],
                    record["clocks"],
                    clock,
                    s,
                    t,
                ),
            )
            for piece in (source_piece, target_piece)
        )
        factor_pairs[factor_pair] += 1

        supports = physical.carrier_supports(
            record["source"], record["target"]
        )
        carrier_pair = tuple(
            physical.carrier_mask(piece, supports)
            for piece in (source_piece, target_piece)
        )
        carrier_pairs[carrier_pair] += 1

        semantic_pair = tuple(
            collar.semantic_value(piece, record["q_pair"], {}) != (0, 0)
            for piece in (source_piece, target_piece)
        )
        semantic_pairs[semantic_pair] += 1

        ancestry_pair = tuple(
            independent.ancestry_chamber(piece, u_events, v_events)
            for piece in (source_piece, target_piece)
        )
        ancestry_pairs[
            (
                ancestry_pair[0] == ancestry_pair[1],
                ancestry_pair,
            )
        ] += 1
        rows.append(
            (
                cell,
                source_piece[:2],
                target_piece[:2],
                source_delta,
                target_delta,
                semantic_pair,
                ancestry_pair[0] == ancestry_pair[1],
            )
        )

    all_six_pair = (
        (((True,) * 6, (True,) * 6)),
        (((True,) * 6, (True,) * 6)),
    )
    delta0 = (True,) + (False,) * 12
    expected_carrier = ((delta0, delta0), (delta0, delta0))
    require(
        source_endpoint_deltas == Counter({((0, 8),): 12})
        and target_endpoint_deltas == Counter({((0, 8),): 12})
        and factor_pairs == Counter({all_six_pair: 12})
        and carrier_pairs == Counter({expected_carrier: 12})
        and semantic_pairs == Counter({(True, True): 12})
        and all(key[0] for key in ancestry_pairs),
        "step2-to-step68 physical sidecar changed",
    )

    # The same endpoint translate persists through steps 68,...,90.  Audit
    # every labelled occurrence, not just the first semantic-preserving step.
    plateau_steps = tuple(range(68, 91))
    plateau_source_endpoint_deltas = Counter()
    plateau_target_endpoint_deltas = Counter()
    plateau_target_endpoint_by_step = {}
    plateau_semantic_pairs = Counter()
    plateau_semantic_by_step = {}
    plateau_factor_count = 0
    plateau_carrier_count = 0
    plateau_ancestry_count = 0
    for record in candidates:
        cell = record["cell"]
        clock, s, t = cell
        source_piece = record["vertices"][2]
        source_endpoint = endpoint_mask(source_piece[:2])
        source_target_endpoint = endpoint_mask(
            (source_piece[0] + SHIFT, source_piece[1] + SHIFT)
        )
        supports = physical.carrier_supports(
            record["source"], record["target"]
        )
        source_semantic = (
            collar.semantic_value(source_piece, record["q_pair"], {})
            != (0, 0)
        )
        source_ancestry = independent.ancestry_chamber(
            source_piece, u_events, v_events
        )
        require(source_semantic, f"plateau source ceased to be live in {cell}")
        for step in plateau_steps:
            target_piece = record["vertices"][step]
            target_endpoint = endpoint_mask(target_piece[:2])
            target_target_endpoint = endpoint_mask(
                (target_piece[0] + SHIFT, target_piece[1] + SHIFT)
            )
            plateau_source_endpoint_deltas[
                translation_deltas(source_endpoint, target_endpoint)
            ] += 1
            target_frame_delta = translation_deltas(
                source_target_endpoint, target_target_endpoint
            )
            plateau_target_endpoint_deltas[target_frame_delta] += 1
            plateau_target_endpoint_by_step.setdefault(
                step, Counter()
            )[target_frame_delta] += 1

            target_semantic = (
                collar.semantic_value(target_piece, record["q_pair"], {})
                != (0, 0)
            )
            plateau_semantic_pairs[(source_semantic, target_semantic)] += 1
            plateau_semantic_by_step.setdefault(step, Counter())[
                (source_semantic, target_semantic)
            ] += 1

            target_factor_pair = (
                factor_signature(
                    target_piece[:2],
                    record["full_module"],
                    record["e3"],
                    record["clocks"],
                    clock,
                    s,
                    t,
                ),
                factor_signature(
                    (
                        target_piece[0] + SHIFT,
                        target_piece[1] + SHIFT,
                    ),
                    record["full_module"],
                    record["e3"],
                    record["clocks"],
                    clock,
                    s,
                    t,
                ),
            )
            if target_factor_pair == ((True,) * 6, (True,) * 6):
                plateau_factor_count += 1
            if physical.carrier_mask(target_piece, supports) == (
                delta0,
                delta0,
            ):
                plateau_carrier_count += 1
            if independent.ancestry_chamber(
                target_piece, u_events, v_events
            ) == source_ancestry:
                plateau_ancestry_count += 1

    expected_semantic_by_step = {
        step: Counter({(True, step % 2 == 0): 12})
        for step in plateau_steps
    }
    require(
        plateau_source_endpoint_deltas
        == Counter({((0, 8),): 12 * len(plateau_steps)})
        and plateau_target_endpoint_deltas
        == Counter({((0, 8),): 264, (): 12})
        and plateau_target_endpoint_by_step
        == {
            step: Counter(
                {((0, 8),) if step < 90 else (): 12}
            )
            for step in plateau_steps
        }
        and plateau_semantic_pairs
        == Counter({(True, True): 144, (True, False): 132})
        and plateau_semantic_by_step == expected_semantic_by_step
        and plateau_factor_count == 12 * len(plateau_steps)
        and plateau_carrier_count == 12 * len(plateau_steps)
        and plateau_ancestry_count == 12 * len(plateau_steps),
        "full step68-to-step90 physical plateau changed: "
        f"source={plateau_source_endpoint_deltas};"
        f"target={plateau_target_endpoint_deltas};"
        f"semantic={plateau_semantic_pairs};"
        f"by_step={plateau_semantic_by_step};"
        f"factors={plateau_factor_count};carrier={plateau_carrier_count};"
        f"ancestry={plateau_ancestry_count}",
    )
    return {
        "rows": tuple(rows),
        "source_endpoint_deltas": source_endpoint_deltas,
        "target_endpoint_deltas": target_endpoint_deltas,
        "factor_pairs": factor_pairs,
        "carrier_pairs": carrier_pairs,
        "semantic_pairs": semantic_pairs,
        "ancestry_pairs": ancestry_pairs,
        "ancestry_event_counts": (len(u_events), len(v_events)),
        "literal_ancestry_counts": tuple(map(len, literal_source)),
        "literal_ancestry_digest": literal_ancestry_digest,
        "literal_ancestry_set_differences": (
            len(set(literal_source[0]) ^ set(literal_target[0])),
            len(set(literal_source[1]) ^ set(literal_target[1])),
        ),
        "plateau_steps": plateau_steps,
        "plateau_occurrences": 12 * len(plateau_steps),
        "plateau_source_endpoint_deltas": plateau_source_endpoint_deltas,
        "plateau_target_endpoint_deltas": plateau_target_endpoint_deltas,
        "plateau_target_endpoint_by_step": plateau_target_endpoint_by_step,
        "plateau_semantic_pairs": plateau_semantic_pairs,
        "plateau_semantic_by_step": plateau_semantic_by_step,
        "plateau_factor_count": plateau_factor_count,
        "plateau_carrier_count": plateau_carrier_count,
        "plateau_ancestry_count": plateau_ancestry_count,
    }


def semilinear_character_audit():
    """Compare the endpoint Z^8 germ with the THM-2857 character-3 fixture."""
    endpoint = ENDPOINT.endpoint
    delta = 66 * H
    cyclotomic_13_step = endpoint.NN // P
    x_power_exact = ENDPOINT.X0 * endpoint.RDIL * delta
    y_power_exact = ENDPOINT.Y0 * endpoint.RDIL * delta
    require(
        x_power_exact % cyclotomic_13_step == 0
        and y_power_exact % cyclotomic_13_step == 0,
        "physical plateau phase is not a 13th-root phase",
    )
    x_power = x_power_exact // cyclotomic_13_step
    y_power = y_power_exact // cyclotomic_13_step
    source_character = (-x_power) % P
    target_character = y_power % P

    terminal = tuple(endpoint.build_set(ENDPOINT.PAT_Q12, ENDPOINT.ZERO))
    q_starts = tuple(left for left, _right in terminal)
    tables = endpoint.make_tabs(
        terminal, ENDPOINT.X0, endpoint.MODS
    )
    translated_I = (I[0] + delta, I[1] + delta)
    source_before, overlap_before = endpoint.x_sweep(
        (I,),
        terminal,
        q_starts,
        ENDPOINT.X0,
        endpoint.MODS,
        tables,
    )
    source_after, overlap_after = endpoint.x_sweep(
        (translated_I,),
        terminal,
        q_starts,
        ENDPOINT.X0,
        endpoint.MODS,
        tables,
    )
    target_before = endpoint.endpoint_sum(
        (I,), -ENDPOINT.Y0, endpoint.MODS
    )
    target_after = endpoint.endpoint_sum(
        (translated_I,), -ENDPOINT.Y0, endpoint.MODS
    )
    specialization_rows = []
    for index, (prime, root) in enumerate(endpoint.MODS):
        omega = pow(root, cyclotomic_13_step, prime)
        require(
            source_before[index] != 0
            and target_before[index] != 0
            and source_after[index]
            == source_before[index]
            * pow(omega, source_character, prime)
            % prime
            and target_after[index]
            == target_before[index]
            * pow(omega, target_character, prime)
            % prime,
            f"physical phase ratio changed in embedding {index}",
        )
        specialization_rows.append(
            (
                prime,
                source_before[index],
                source_after[index],
                source_character,
                target_before[index],
                target_after[index],
                target_character,
            )
        )

    physical_b_shift = 8
    torsor_parameter_scale = 10
    torsor_parameter_shift = (
        torsor_parameter_scale * physical_b_shift
    ) % P
    candidate_character = (3 * torsor_parameter_shift) % P
    require(
        (x_power, y_power) == (396, 1254)
        and (source_character, target_character) == (7, 6)
        and overlap_before == overlap_after == 1568941674300
        and torsor_parameter_shift == 2
        and candidate_character == target_character
        and source_character == (-candidate_character) % P,
        "THM-2857 character comparison changed",
    )
    return {
        "delta": delta,
        "exact_powers": (x_power, y_power),
        "characters": (source_character, target_character),
        "overlap": overlap_before,
        "torsor_parameter_shift": torsor_parameter_shift,
        "candidate_character": candidate_character,
        "specializations": tuple(specialization_rows),
    }


def main():
    paths, horn_rows, common_atoms = common42_path_geometry()
    endpoint = endpoint_address_carry_audit()
    endpoint_legs = endpoint_leg_counts(paths)
    step68 = step68_physical_audit(paths)
    semilinear = semilinear_character_audit()

    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "scratch audit contains an executable assert",
    )

    print("THM-2859 Q3/Q11 COLLAR ENDPOINT-CARRY INTERTWINER AUDIT")
    print(f"pinned={tuple(PINNED.items())}")
    print(
        f"horn_cells={len(horn_rows)};common_atoms={common_atoms};"
        "unique_labelled_I_roots=20;physical_R_M1_M2="
        f"{(R, M1, M2)}"
    )
    print(
        f"horn_root_index_path_length="
        f"{tuple(row[:3] for row in horn_rows)};"
        "all_survive_14h=1;semantic=(live,dead,live)"
    )
    print(
        "endpoint_ladder_source="
        f"{endpoint['source_descriptions']};"
        f"target_all_equal={endpoint['target_description']}"
    )
    print(
        "C169_vs_F13xF13="
        "digit_section=iota(v,w)=v+13w;"
        "cyclic_successor=Y off v=12 and ZY on v=12;"
        "formal_coordinate_comparison=ancestry_T_label_and_endpoint_Z_"
        "both_vertical_(0,1);intertwiner_not_proved"
    )
    print(
        f"selected_I_endpoint_rectangle={endpoint['base_description']};"
        f"Z_derivative_hamming_l1={endpoint['z_hamming']};"
        f"coefficient_counts={tuple(sorted(endpoint['derivative_counts'].items()))};"
        "augmentation=0"
    )
    print(
        "pulled_allocation_endpoint_masks="
        f"descriptions={tuple(sorted(endpoint['allocation_descriptions'].items()))};"
        f"q0_q11_translation_deltas={endpoint['q0_q11_translation_deltas']};"
        f"individual_ranks={tuple(sorted((role, tuple(sorted(rows.items()))) for role, rows in endpoint['allocation_ranks'].items()))};"
        f"difference_rank_augmentation={tuple(sorted(endpoint['allocation_difference_ranks'].items()))};"
        "scope=full_endpoint_masks_on_pulled_intervals_not_a_q_action"
    )
    print(
        "exact_convolution_ranks="
        "prime:1000003,"
        f"endpoint_mask:{endpoint['ranks'][0]},"
        f"full_endpoint_plane_Z_derivative:{endpoint['ranks'][1]},"
        f"full_endpoint_plane_Z8_derivative:{endpoint['ranks'][2]},"
        f"vertical_projection:{endpoint['ranks'][3]}"
    )
    print(
        "rooted_path_endpoint_translates="
        f"{tuple(endpoint['translated_steps'].items())};"
        f"normalized_(a,1)_hits={endpoint['normalized_hits']};"
        "verdict=no canonically normalized endpoint Z vertex, but Z^8 occurs"
    )
    print(
        "endpoint_h8_h9_h4_leg_audit="
        f"{tuple((scope, data['paths'], tuple(sorted(data['leg_counts'].items())), tuple(sorted(data['leg_path_counts'].items())), tuple(sorted(data['leg_cells'].items())), tuple(sorted(data['leg_witnesses'].items())), tuple(sorted(data['triangle_counts'].items())), tuple(sorted(data['triangle_witnesses'].items()))) for scope, data in endpoint_legs.items())};"
        "typing=endpoint-address translations only, not allocation-q arrows"
    )
    print(
        "step2_to_step68_Z8_reference="
        f"cells={tuple(row[0] for row in step68['rows'])};"
        f"source_endpoint_deltas={tuple(sorted(step68['source_endpoint_deltas'].items()))};"
        f"target_endpoint_deltas={tuple(sorted(step68['target_endpoint_deltas'].items()))};"
        f"semantic={tuple(sorted(step68['semantic_pairs'].items()))};"
        f"ancestry_equal={tuple(sorted(step68['ancestry_pairs'].items()))};"
        f"ancestry_event_counts={step68['ancestry_event_counts']};"
        f"literal_ancestry_counts={step68['literal_ancestry_counts']};"
        f"literal_ancestry_digest={step68['literal_ancestry_digest']};"
        "literal_ancestry_set_differences="
        f"{step68['literal_ancestry_set_differences']};"
        "all_six_both_frames=1;source_target_carrier_delta0_preserved=1;"
        "labelled_rank=12;physical_interval_rank=1;"
        "Z8_generator_value=1_since_gcd(8,13)=1;"
        "physical_C13_action=0;endpoint_ancestry_intertwiner=0"
    )
    print(
        "step68_to_step90_Z8_plateau="
        f"steps={step68['plateau_steps']};"
        f"labelled_occurrences={step68['plateau_occurrences']};"
        "source_endpoint_deltas="
        f"{tuple(sorted(step68['plateau_source_endpoint_deltas'].items()))};"
        "target_endpoint_deltas="
        f"{tuple(sorted(step68['plateau_target_endpoint_deltas'].items()))};"
        f"semantic={tuple(sorted(step68['plateau_semantic_pairs'].items()))};"
        "even_steps_live_live=144;odd_steps_live_dead=132;"
        f"all_six_count={step68['plateau_factor_count']};"
        f"carrier_count={step68['plateau_carrier_count']};"
        f"ancestry_chamber_count={step68['plateau_ancestry_count']};"
        "source_frame_plateau=68..90;"
        "target_frame_plateau=68..89;"
        "step90_target_frame_boundary_failure=12;"
        "first_semantic_preserving_step=68;"
        "last_two_frame_semantic_preserving_step=88"
    )
    print(
        "THM2857_candidate_character_comparison="
        f"fixture={THM2857_CANDIDATE};"
        f"physical_delta={semilinear['delta']};"
        f"exact_X_Y_powers={semilinear['exact_powers']};"
        "physical_b_shift=8;torsor_reparameterization=r=10*b;"
        f"delta_r={semilinear['torsor_parameter_shift']};"
        f"candidate_centered_character_phase={semilinear['candidate_character']};"
        f"source_x_sweep_phase={semilinear['characters'][0]};"
        f"target_endpoint_sum_phase={semilinear['characters'][1]};"
        f"overlap={semilinear['overlap']};"
        f"two_embedding_rows={semilinear['specializations']};"
        "verdict=target_matches_character3_and_source_is_inverse;"
        "paired_phase_cancels;comparison_only_not_a_physical_clutch"
    )
    print(
        "TYPE_VERDICT="
        "positive labelled source section exists: each horn cell has a unique "
        "R with S=ad sending R to the fixed source I=M2; "
        "but R has empty source endpoint while M1=M2 have the same 81-point "
        "source endpoint mask, and the q7 horn needs source I retained while "
        "target E3 is absent.  The collar locks source birth to E3 repair; "
        "the horn needs the off-diagonal state (source present,E3 absent).  "
        "The endpoint address plane has the exact charged vertical C13 "
        "coordinate "
        "(rank 12 after vertical projection).  The +66h common-tangent move "
        "from step 2 to step 68 induces the generator-valued endpoint "
        "translation Z^8 on twelve horn cells while preserving both endpoint "
        "masks, all six factors, carrier, semantic value, and ancestry "
        "chamber.  It does not supply a C13 action or identify endpoint Z "
        "with ancestry carry T.  What remains missing is a lawful "
        "action/intertwiner and coupling from this endpoint-address arrow to "
        "the q3/q11/q7 allocation and E3/complement semantic-word horn."
    )
    print(f"endpoint_masks_evaluated={endpoint['mask_cache']}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
