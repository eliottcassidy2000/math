#!/usr/bin/env python3
"""Independent referee audit of the semilinear rectangle bridge scratch.

The audit uses a second canonical endpoint evaluator, enumerates affine
maps in one combined pass, exhausts the full coordinatewise PGL(2,13)
factor table, and checks the cyclotomic/Prony claims directly.

Scratch status only.  No executable Python assertion statement is used.
"""

from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge" / "results"
SOURCE = ROOT / ".scratch" / "lrc_semilinear_rectangle_bridge_20260728"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_digest(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINS = {
    SOURCE / "affine_rectangle_bridge_audit.py":
        "0b17358531fd69fc35160f0d6a8486786123d33fb2ec714bcedf2955340bd543",
    SOURCE / "REPORT.md":
        "cbc15b42dc66bbe8e5cefb3aaacd9bff639b715fb4ab49f8d5a7c97921a5d41d",
    SOURCE / "audit.out":
        "66c38cfd4c68f7c4ff3471d84ccad33c7329a142d5970c358d0d28ca07bc0ad3",
    COMP / "lrc14_horn_collar_endpoint_coboundary_thm2859.py":
        "3d7702641a2df258b829538d8fcf1d066cdf5f426cceef5781bbcfb37747bc15",
    COMP / "lrc14_horn_collar_endpoint_carry_thm2859.py":
        "6e062f3cc57c80fcff372c272bc138e280205bb953e484f1cc267340774260f0",
    COMP / "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    RESULTS / "lrc14_two_origin_endpoint_projective_kummer_thm2868.out":
        "ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9",
    COMP / "lrc14_endpoint_kummer_galois_bockstein_thm2874.py":
        "3f15c44dc5f66c660ac3605cc25814adc39594bf193aa130a0f5353d6a6178b0",
    RESULTS / "lrc14_endpoint_kummer_galois_bockstein_thm2874.out":
        "90b993b56508ef3603f94104596b899ed9ec7084a2b58ead1604882873ef5455",
}
for path, expected in PINS.items():
    require(lf_digest(path) == expected, f"pinned input changed: {path}")


import lrc14_horn_collar_endpoint_carry_thm2859 as endpoint
import lrc14_horn_collar_endpoint_coboundary_thm2859 as independent
import lrc14_q3_q11_transverse_endpoint_horn_thm2847 as horn


P = 13
KEYS = tuple((x, y) for x in range(P) for y in range(P))
ORIGINS = ((0, 0), (12, 0))
QS = (0, 3, 11)
STEPS = (2, 68)


def point_index(point):
    return P * point[0] + point[1]


def point_mask(points):
    mask = 0
    for point in points:
        mask |= 1 << point_index(point)
    return mask


def translate_mask(mask, translation):
    tx, ty = translation
    result = 0
    for x, y in KEYS:
        if (mask >> point_index((x, y))) & 1:
            result |= 1 << point_index(
                ((x + tx) % P, (y + ty) % P)
            )
    return result


def matrix_point(matrix, point):
    a, b, c, d = matrix
    x, y = point
    return (
        (a * x + b * y) % P,
        (c * x + d * y) % P,
    )


def affine_point(matrix, translation, point):
    linear = matrix_point(matrix, point)
    return (
        (linear[0] + translation[0]) % P,
        (linear[1] + translation[1]) % P,
    )


def linear_mask(points, matrix):
    return point_mask(matrix_point(matrix, point) for point in points)


def gl2():
    return tuple(
        (a, b, c, d)
        for a, b, c, d in product(range(P), repeat=4)
        if (a * d - b * c) % P
    )


def vertical_character_image(matrix):
    a, b, c, d = matrix
    determinant = (a * d - b * c) % P
    inverse_det = pow(determinant, -1, P)
    inverse_c = -c * inverse_det % P
    inverse_d = a * inverse_det % P
    return (3 * inverse_c % P, 3 * inverse_d % P)


def is_axis(matrix):
    return matrix[1] == 0 and matrix[2] == 0


def fixes_origins(matrix, translation):
    return tuple(
        affine_point(matrix, translation, point)
        for point in ORIGINS
    ) == ORIGINS


def preserves_origins_setwise(matrix, translation):
    return {
        affine_point(matrix, translation, point)
        for point in ORIGINS
    } == set(ORIGINS)


def new_stat(source_size, target_size):
    return {
        "source_size": source_size,
        "target_size": target_size,
        "best_overlap": -1,
        "witness": None,
        "exact": [],
    }


def update_stat(stat, overlap, row):
    if overlap > stat["best_overlap"]:
        stat["best_overlap"] = overlap
        stat["witness"] = row
    if (
        overlap == stat["source_size"]
        and overlap == stat["target_size"]
    ):
        stat["exact"].append(row)


def hamming(stat):
    return (
        stat["source_size"] + stat["target_size"]
        - 2 * stat["best_overlap"]
    )


def affine_defect(source, target, row):
    matrix, translation = row
    image = {
        affine_point(matrix, translation, point)
        for point in source
    }
    return (
        tuple(sorted(image - target)),
        tuple(sorted(target - image)),
    )


def sumset(first, second, left, right):
    return {
        (left * x + right * y) % P
        for x in first
        for y in second
    }


def pgl2():
    rows = []
    for matrix in product(range(P), repeat=4):
        a, b, c, d = matrix
        if (a * d - b * c) % P == 0:
            continue
        first = next(value for value in matrix if value)
        if first == 1:
            rows.append(matrix)
    return tuple(rows)


def mobius_maps(source, target, matrices):
    rows = []
    for a, b, c, d in matrices:
        if any((c * value + d) % P == 0 for value in source):
            continue
        image = {
            (a * value + b)
            * pow((c * value + d) % P, -1, P) % P
            for value in source
        }
        if image == target:
            rows.append((a, b, c, d))
    return tuple(rows)


def multiplicative_order(exponent, modulus):
    return modulus // gcd(exponent, modulus)


def main():
    allocation = horn.allocation
    _module, full, _details, _e3, _clocks, _q_pairs = (
        allocation.build_geometry()
    )
    unit = full.T // P
    atom = allocation.ATOM_INTERVAL
    step_bases = {
        2: atom,
        68: (
            atom[0] + 66 * endpoint.H,
            atom[1] + 66 * endpoint.H,
        ),
    }

    masks = {}
    points = {}
    factors = {}
    target_frame_equal = {}
    evaluator_equal = {}
    for q in QS:
        for step in STEPS:
            interval = horn.circular_shift_interval(
                step_bases[step], q * unit, full.T
            )
            primary_mask = endpoint.endpoint_mask(interval)
            independent_mask = independent.endpoint_mask(interval)
            target_mask = endpoint.endpoint_mask((
                interval[0] + endpoint.SHIFT,
                interval[1] + endpoint.SHIFT,
            ))
            evaluator_equal[q, step] = primary_mask == independent_mask
            target_frame_equal[q, step] = primary_mask == target_mask
            masks[q, step] = point_mask(
                point
                for point, value in zip(KEYS, primary_mask)
                if value
            )
            points[q, step] = frozenset(
                point
                for point, value in zip(KEYS, primary_mask)
                if value
            )
            first = frozenset(x for x, _y in points[q, step])
            second = frozenset(y for _x, y in points[q, step])
            require(
                points[q, step]
                == frozenset(product(first, second)),
                "one reconstructed mask ceased to be Cartesian",
            )
            factors[q, step] = (first, second)

    expected_factors = {
        (0, 2): (
            frozenset((0, 1, 2, 3, 4, 5, 6, 7, 12)),
            frozenset((0, 1, 2, 3, 4, 5, 8, 9, 10)),
        ),
        (0, 68): (
            frozenset((0, 1, 2, 3, 4, 5, 6, 7, 12)),
            frozenset((0, 1, 2, 5, 6, 7, 8, 9, 10)),
        ),
        (3, 2): (
            frozenset(range(10)),
            frozenset((0, 3, 4, 5, 8, 9, 10, 11, 12)),
        ),
        (3, 68): (
            frozenset(range(10)),
            frozenset((0, 5, 6, 7, 8, 9, 10, 11, 12)),
        ),
        (11, 2): (
            frozenset((0, 1, 2, 3, 4, 5, 8, 9, 12)),
            frozenset((0, 1, 2, 3, 4, 5, 8, 11, 12)),
        ),
        (11, 68): (
            frozenset((0, 1, 2, 3, 4, 5, 8, 9, 12)),
            frozenset((0, 1, 2, 5, 6, 7, 8, 11, 12)),
        ),
    }
    require(
        factors == expected_factors
        and all(evaluator_equal.values())
        and all(target_frame_equal.values()),
        "rectangle reconstruction or frame equality changed",
    )

    endpoint_base = endpoint.ENDPOINT
    right_origin_rep = endpoint_base.REPS[allocation.RIGHT_ORIGIN]

    def literal_factor_bits(interval):
        bits = []
        for index, speed in enumerate(endpoint_base.endpoint.W):
            if index == 0:
                dangerous = full.make_comb(
                    speed, 91,
                    -13 - 7 * right_origin_rep[index],
                    13 - 7 * right_origin_rep[index],
                )
            else:
                dangerous = full.make_comb(
                    speed, 182,
                    -13 - 14 * right_origin_rep[index],
                    13 - 14 * right_origin_rep[index],
                )
            hit = allocation.intersect_sorted((interval,), dangerous)
            if allocation.contained(interval, dangerous):
                bits.append("D")
            elif not hit:
                bits.append("S")
            else:
                raise RuntimeError(
                    f"literal factor {index} became partial"
                )
        return "".join(bits)

    expected_literal_orbit = (
        "SSSSSSSSD",
        "SSDSSSSSD",
        "SSDSSSSSD",
        "SSSSSSSSD",
        "SDSSSSSSD",
        "DDSSSSSSD",
        "DSSSSDSSD",
        "DSSSSDSSD",
        "DSSSSSSSD",
        "SSSSDSSSD",
        "SSSSDSSSD",
        "SSSSSSSSD",
        "SSSDSSSSD",
    )
    literal_orbits = {}
    for step in STEPS:
        target_base = (
            step_bases[step][0] + allocation.physical.SHIFT,
            step_bases[step][1] + allocation.physical.SHIFT,
        )
        literal_orbits[step] = tuple(
            literal_factor_bits(horn.circular_shift_interval(
                target_base, q * unit, full.T
            ))
            for q in range(P)
        )
    guard_q5_corners = {
        (bits[0], bits[5])
        for bits in literal_orbits[2]
    }
    require(
        literal_orbits
        == {2: expected_literal_orbit, 68: expected_literal_orbit}
        and guard_q5_corners == {("S", "S"), ("D", "S"), ("D", "D")}
        and ("S", "D") not in guard_q5_corners
        and all(
            representative[0] == representative[5] == 0
            for representative in endpoint_base.REPS.values()
        )
        and all(
            literal_orbits[step][q] == "SSSSSSSSD"
            for step in STEPS for q in QS
        ),
        "guard/q5 missing-corner comparison changed",
    )

    mixed_projection_sizes = {
        q: {
            len(sumset(
                factors[q, 2][0],
                factors[q, 2][1],
                left,
                right,
            ))
            for left in range(1, P)
            for right in range(1, P)
        }
        for q in QS
    }
    require(
        mixed_projection_sizes
        == {0: {13}, 3: {13}, 11: {13}},
        "mixed-row Cauchy-Davenport boundary changed",
    )

    matrices = gl2()
    require(len(matrices) == 26_208, "GL2 census changed")
    monomial_count = sum(
        (
            matrix[0] != 0 and matrix[3] != 0
            and matrix[1] == 0 and matrix[2] == 0
        ) or (
            matrix[1] != 0 and matrix[2] != 0
            and matrix[0] == 0 and matrix[3] == 0
        )
        for matrix in matrices
    )
    require(monomial_count == 288, "monomial GL2 census changed")

    comparison_pairs = {
        "static_0_3": ((0, 2), (3, 2)),
        "static_0_11": ((0, 2), (11, 2)),
        "edge_0": ((0, 2), (0, 68)),
        "edge_3": ((3, 2), (3, 68)),
        "edge_11": ((11, 2), (11, 68)),
    }
    stats = {
        name: new_stat(len(points[source]), len(points[target]))
        for name, (source, target) in comparison_pairs.items()
    }
    pair_stats = {
        q: new_stat(
            len(points[0, 2]) + len(points[0, 68]),
            len(points[q, 2]) + len(points[q, 68]),
        )
        for q in (3, 11)
    }
    constraint_names = (
        "axis", "chi3", "axis_chi3", "origins", "chi3_origins"
    )
    constrained_edges = {
        q: {
            name: new_stat(
                len(points[q, 2]), len(points[q, 68])
            )
            for name in constraint_names
        }
        for q in QS
    }
    constrained_static = {
        q: {
            name: new_stat(
                len(points[0, 2]), len(points[q, 2])
            )
            for name in constraint_names
        }
        for q in (3, 11)
    }

    inverse_target_banks = {
        key: tuple(
            translate_mask(
                masks[key],
                ((-translation[0]) % P, (-translation[1]) % P),
            )
            for translation in KEYS
        )
        for key in masks
    }
    chi3_linear_count = 0
    axis_linear_count = 0
    axis_chi3_linear_count = 0
    origin_affine_count = 0
    chi3_origin_affine_count = 0

    for matrix in matrices:
        linear = {
            key: linear_mask(points[key], matrix)
            for key in points
        }
        axis = is_axis(matrix)
        chi3 = vertical_character_image(matrix) == (0, 3)
        chi3_linear_count += int(chi3)
        axis_linear_count += int(axis)
        axis_chi3_linear_count += int(axis and chi3)
        for translation_index, translation in enumerate(KEYS):
            origin_fixed = fixes_origins(matrix, translation)
            origin_affine_count += int(origin_fixed)
            chi3_origin_affine_count += int(origin_fixed and chi3)
            overlaps = {}
            for name, (source, target) in comparison_pairs.items():
                overlap = (
                    linear[source]
                    & inverse_target_banks[target][translation_index]
                ).bit_count()
                overlaps[name] = overlap
                update_stat(stats[name], overlap, (matrix, translation))
            for q in (3, 11):
                pair_overlap = (
                    linear[0, 2]
                    & inverse_target_banks[q, 2][translation_index]
                ).bit_count() + (
                    linear[0, 68]
                    & inverse_target_banks[q, 68][translation_index]
                ).bit_count()
                update_stat(
                    pair_stats[q],
                    pair_overlap,
                    (matrix, translation),
                )
            predicates = {
                "axis": axis,
                "chi3": chi3,
                "axis_chi3": axis and chi3,
                "origins": origin_fixed,
                "chi3_origins": chi3 and origin_fixed,
            }
            for q in QS:
                for name, truth in predicates.items():
                    if truth:
                        update_stat(
                            constrained_edges[q][name],
                            overlaps[f"edge_{q}"],
                            (matrix, translation),
                        )
            for q in (3, 11):
                for name, truth in predicates.items():
                    if truth:
                        update_stat(
                            constrained_static[q][name],
                            overlaps[f"static_0_{q}"],
                            (matrix, translation),
                        )

    require(
        (
            chi3_linear_count,
            axis_linear_count,
            axis_chi3_linear_count,
            origin_affine_count,
            chi3_origin_affine_count,
        ) == (156, 144, 12, 156, 13),
        "typed affine subgroup counts changed",
    )
    expected_affine = {
        "static_0_3": (81, 9, 0),
        "static_0_11": (64, 34, 0),
        "edge_0": (81, 0, 8),
        "edge_3": (80, 20, 0),
        "edge_11": (72, 18, 0),
    }
    require(
        {
            name: (
                stat["best_overlap"],
                hamming(stat),
                len(stat["exact"]),
            )
            for name, stat in stats.items()
        } == expected_affine,
        "full AGL2 search summary changed",
    )
    require(
        {
            q: (
                stat["best_overlap"],
                hamming(stat),
                len(stat["exact"]),
            )
            for q, stat in pair_stats.items()
        } == {3: (144, 54, 0), 11: (128, 68, 0)},
        "pair-map AGL2 search summary changed",
    )
    expected_edge_hamming = {
        0: {
            "axis": 0, "chi3": 0, "axis_chi3": 0,
            "origins": 36, "chi3_origins": 36,
        },
        3: {
            "axis": 20, "chi3": 40, "axis_chi3": 40,
            "origins": 20, "chi3_origins": 40,
        },
        11: {
            "axis": 18, "chi3": 36, "axis_chi3": 36,
            "origins": 36, "chi3_origins": 36,
        },
    }
    expected_static_hamming = {
        3: {
            "axis": 9, "chi3": 9, "axis_chi3": 9,
            "origins": 27, "chi3_origins": 59,
        },
        11: {
            "axis": 34, "chi3": 50, "axis_chi3": 50,
            "origins": 64, "chi3_origins": 64,
        },
    }
    require(
        {
            q: {
                name: hamming(stat)
                for name, stat in rows.items()
            }
            for q, rows in constrained_edges.items()
        } == expected_edge_hamming
        and {
            q: {
                name: hamming(stat)
                for name, stat in rows.items()
            }
            for q, rows in constrained_static.items()
        } == expected_static_hamming,
        "constrained AGL2 minima changed",
    )

    q0_exact = tuple(stats["edge_0"]["exact"])
    q0_chi3 = tuple(
        row for row in q0_exact
        if vertical_character_image(row[0]) == (0, 3)
    )
    q0_pointwise = tuple(
        row for row in q0_exact if fixes_origins(*row)
    )
    q0_setwise = tuple(
        row for row in q0_exact if preserves_origins_setwise(*row)
    )
    require(
        len(q0_exact) == 8
        and len(q0_chi3) == 2
        and not q0_pointwise
        and not q0_setwise,
        "q0 exact map character/origin boundary changed",
    )
    q0_origin_images = tuple(
        tuple(affine_point(*row, origin) for origin in ORIGINS)
        for row in q0_chi3
    )
    require(
        q0_origin_images
        == (
            ((0, 5), (12, 5)),
            ((6, 5), (7, 5)),
        ),
        "q0 named-origin images changed",
    )

    selected_rows = {
        "static_0_3": ((1, 0, 0, 1), (1, 8)),
        "edge_3": ((1, 0, 0, 2), (0, 0)),
        "edge_11_axis": ((1, 0, 0, 6), (0, 1)),
    }
    selected_defects = {
        "static_0_3": affine_defect(
            points[0, 2], points[3, 2],
            selected_rows["static_0_3"],
        ),
        "edge_3": affine_defect(
            points[3, 2], points[3, 68],
            selected_rows["edge_3"],
        ),
        "edge_11_axis": affine_defect(
            points[11, 2], points[11, 68],
            selected_rows["edge_11_axis"],
        ),
    }
    require(
        selected_defects["static_0_3"]
        == (
            (),
            tuple((9, value) for value in sorted(factors[3, 2][1])),
        )
        and selected_defects["edge_3"]
        == (
            tuple((value, 3) for value in range(10)),
            tuple((value, 12) for value in range(10)),
        )
        and selected_defects["edge_11_axis"]
        == (
            tuple((value, 10) for value in sorted(factors[11, 2][0])),
            tuple((value, 11) for value in sorted(factors[11, 2][0])),
        ),
        "selected one-fibre defects changed",
    )
    require(
        vertical_character_image(selected_rows["edge_3"][0]) == (0, 8)
        and vertical_character_image(
            selected_rows["edge_11_axis"][0]
        ) == (0, 7)
        and (8 * 2) % P == 3
        and (7 * 6) % P == 3,
        "Fourier covector/coefficient rechart arithmetic changed",
    )

    factor_names = {
        "q0_A": factors[0, 2][0],
        "q0_B2": factors[0, 2][1],
        "q0_B68": factors[0, 68][1],
        "q3_A": factors[3, 2][0],
        "q3_B2": factors[3, 2][1],
        "q3_B68": factors[3, 68][1],
        "q11_A": factors[11, 2][0],
        "q11_B2": factors[11, 2][1],
        "q11_B68": factors[11, 68][1],
    }
    pgl_matrices = pgl2()
    require(len(pgl_matrices) == 2_184, "PGL2 line count changed")
    pgl_table = {
        (source, target): mobius_maps(
            factor_names[source], factor_names[target], pgl_matrices
        )
        for source in factor_names
        for target in factor_names
    }

    def coordinatewise_bridge(source_factors, target_factors):
        direct = (
            bool(pgl_table[source_factors[0], target_factors[0]])
            and bool(pgl_table[source_factors[1], target_factors[1]])
        )
        swapped = (
            bool(pgl_table[source_factors[0], target_factors[1]])
            and bool(pgl_table[source_factors[1], target_factors[0]])
        )
        return direct, swapped

    pgl_full_bridge_checks = {
        "q0_q11_step2": coordinatewise_bridge(
            ("q0_A", "q0_B2"), ("q11_A", "q11_B2")
        ),
        "q0_q11_step68": coordinatewise_bridge(
            ("q0_A", "q0_B68"), ("q11_A", "q11_B68")
        ),
        "q3_edge": coordinatewise_bridge(
            ("q3_A", "q3_B2"), ("q3_A", "q3_B68")
        ),
        "q11_edge": coordinatewise_bridge(
            ("q11_A", "q11_B2"), ("q11_A", "q11_B68")
        ),
    }
    require(
        all(
            not direct and not swapped
            for direct, swapped in pgl_full_bridge_checks.values()
        ),
        "a full coordinatewise PGL rectangle bridge appeared",
    )
    pgl_nonaffine_witness_counts = {
        pair: sum(matrix[2] != 0 for matrix in rows)
        for pair, rows in pgl_table.items()
        if any(matrix[2] != 0 for matrix in rows)
    }
    require(
        pgl_nonaffine_witness_counts
        == {
            ("q3_A", "q3_A"): 6,
            ("q11_A", "q11_A"): 4,
            ("q11_A", "q11_B2"): 4,
            ("q11_B2", "q11_A"): 4,
            ("q11_B2", "q11_B2"): 4,
        },
        "PGL non-affine factor-map correction census changed",
    )

    xi_order = 2_366
    node_exponents = (13, 169)
    require(
        tuple(
            multiplicative_order(exponent, xi_order)
            for exponent in node_exponents
        ) == (182, 14),
        "Prony node orders changed",
    )
    cyclotomic_rows = {}
    for slope in (2, 6):
        units = tuple(
            value for value in range(1, xi_order)
            if gcd(value, xi_order) == 1
            and value % P == slope
        )
        left_fixed = tuple(
            value for value in units
            if 13 * (value - 1) % xi_order == 0
        )
        right_fixed = tuple(
            value for value in units
            if 169 * (value - 1) % xi_order == 0
        )
        unordered_fixed = tuple(
            value for value in units
            if {
                13 * value % xi_order,
                169 * value % xi_order,
            } == {13, 169}
        )
        cyclotomic_rows[slope] = (
            len(units),
            len(left_fixed),
            len(right_fixed),
            len(unordered_fixed),
        )
    require(
        cyclotomic_rows == {2: (78, 0, 13, 0), 6: (78, 0, 13, 0)},
        "cyclotomic lift/Prony-node census changed",
    )

    theorem_2868 = (
        ROOT / "01-canon" / "theorems"
        / "THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent.md"
    ).read_text()
    theorem_2874 = (
        ROOT / "01-canon" / "theorems"
        / "THM-2874-endpoint-kummer-galois-clutch-and-bockstein-seam-transgression.md"
    ).read_text()
    require(
        "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
        in theorem_2868
        and "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
        in theorem_2874
        and "RESERVED / UNPROVED EMPTY STUB" not in theorem_2868,
        "promoted THM-2868/2874 status changed",
    )

    source_tree = ast.parse(Path(__file__).read_text())
    require(
        not any(
            isinstance(node, ast.Assert)
            for node in ast.walk(source_tree)
        ),
        "executable Python assertion statement found",
    )

    print("SEMILINEAR RECTANGLE BRIDGE INDEPENDENT REFEREE AUDIT")
    print(
        "pinned="
        f"{tuple((path.name, digest) for path, digest in PINS.items())}"
    )
    print(
        "rectangle_factors="
        f"{tuple(sorted((key, tuple(map(tuple, value))) for key, value in factors.items()))};"
        f"independent_evaluator_equal:{all(evaluator_equal.values())};"
        f"source_target_frame_equal:{all(target_frame_equal.values())}"
    )
    print(
        "group_counts="
        f"GL2:{len(matrices)},AGL2:{len(matrices) * P * P},"
        f"monomial:{monomial_count},chi3:{chi3_linear_count},"
        f"axis:{axis_linear_count},axis_chi3:{axis_chi3_linear_count},"
        f"origins:{origin_affine_count},"
        f"chi3_origins:{chi3_origin_affine_count}"
    )
    print(
        "mixed_projection_sizes="
        f"{tuple((q, tuple(sorted(values))) for q, values in mixed_projection_sizes.items())};"
        "Cauchy_Davenport_forces_monomial_exact_rectangle_maps"
    )
    print(
        "affine_summaries="
        f"{tuple((name, stat['best_overlap'], hamming(stat), len(stat['exact'])) for name, stat in stats.items())};"
        f"pair_summaries:{tuple((q, stat['best_overlap'], hamming(stat), len(stat['exact'])) for q, stat in pair_stats.items())}"
    )
    print(
        "constrained_edge_hamming="
        f"{expected_edge_hamming};"
        f"constrained_static_hamming:{expected_static_hamming}"
    )
    print(
        "selected_defects="
        f"{selected_defects};"
        "covectors=q0:(0,3),q3_near:(0,8),q11_near:(0,7);"
        "coefficient_recharts=2,6"
    )
    print(
        f"q0_exact_count:{len(q0_exact)};"
        f"q0_chi3_count:{len(q0_chi3)};"
        f"q0_named_origin_pointwise:{len(q0_pointwise)};"
        f"q0_named_origin_setwise:{len(q0_setwise)};"
        f"q0_chi3_origin_images:{q0_origin_images}"
    )
    print(
        f"PGL2_line:{len(pgl_matrices)};"
        f"full_coordinatewise_bridge_checks:{pgl_full_bridge_checks};"
        f"nonaffine_factor_map_correction:{pgl_nonaffine_witness_counts}"
    )
    print(
        f"Prony_node_orders:{tuple(multiplicative_order(exponent, xi_order) for exponent in node_exponents)};"
        f"cyclotomic_rows:{cyclotomic_rows}"
    )
    print(
        "guard_q5_comparison="
        f"corners:{tuple(sorted(guard_q5_corners))};"
        "missing:(S,D);step2_equals_step68:1;"
        "q0_q3_q11_all_map_to:(S,S);"
        "endpoint_address_projection_is_constant_on_all_169_addresses;"
        "therefore_nonzero_rectangle_edge_defects_map_to_zero under the "
        "natural guard/q5 quotient and do not explain the missing corner"
    )
    print(
        "promotion_verdict=finite affine support theorem survives; "
        "correct PGL claim to absence of full coordinatewise rectangle "
        "bridges, because non-affine individual factor maps do exist; "
        "replace stale RESERVED THM-2868 language by compatibility with "
        "promoted signed/projective THM-2868 and its THM-2874 clutch"
    )
    print("ALL INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    main()
