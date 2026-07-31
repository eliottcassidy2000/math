#!/usr/bin/env python3
"""Exact companion for THM-2877 endpoint-rectangle classification.

It reconstructs the q0, q3, and q11 step-2/step-68 masks from the
canonical endpoint evaluator, then searches:

* every affine map of one F_13 coordinate;
* every monomial affine map of F_13^2;
* every affine map of F_13^2, using exact bit-set correlations;
* the subgroup preserving the labelled vertical character-three line; and
* the subgroup fixing the two THM-2863 named origins.

No executable Python assert statement is used.
"""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from math import gcd
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_digest(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINS = {
    "lrc14_horn_collar_endpoint_carry_thm2859.py":
        "6e062f3cc57c80fcff372c272bc138e280205bb953e484f1cc267340774260f0",
    "lrc14_horn_collar_prony_typed_descent_gate_thm2859.py":
        "ff9a954e65209d0b96de7d9215ccc6a38dfdbb16245414564a63237924efad28",
    "lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
}
for name, expected in PINS.items():
    require(
        lf_digest(COMP / name) == expected,
        f"pinned canonical dependency changed: {name}",
    )

RESULT_PINS = {
    "lrc14_horn_collar_endpoint_carry_thm2859.out":
        "b33065698ee0ef4d3513ab51562244b331b01c1843c42c37949fb7406dbf239b",
    "lrc14_horn_collar_prony_typed_descent_gate_thm2859.out":
        "487363c8d8d34cf703dd83fa6d3867b5932796454792a6394da44354bd59278b",
    "lrc14_q3_q11_transverse_endpoint_horn_thm2847.out":
        "155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d",
}
for name, expected in RESULT_PINS.items():
    require(
        lf_digest(RESULTS / name) == expected,
        f"pinned canonical output changed: {name}",
    )


import lrc14_horn_collar_endpoint_carry_thm2859 as endpoint
import lrc14_q3_q11_transverse_endpoint_horn_thm2847 as horn


P = 13
QS = (0, 3, 11)
STEPS = (2, 68)
ORIGINS = ((0, 0), (12, 0))


def add_point(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def matrix_point(matrix, point):
    a, b, c, d = matrix
    x, y = point
    return ((a * x + b * y) % P, (c * x + d * y) % P)


def affine_point(matrix, translation, point):
    return add_point(matrix_point(matrix, point), translation)


def affine_1d_maps(source, target):
    source = frozenset(source)
    target = frozenset(target)
    return tuple(
        (multiplier, shift)
        for multiplier in range(1, P)
        for shift in range(P)
        if {
            (multiplier * value + shift) % P
            for value in source
        } == target
    )


def pgl1():
    rows = []
    for matrix in product(range(P), repeat=4):
        a, b, c, d = matrix
        if (a * d - b * c) % P == 0:
            continue
        first = next(value for value in matrix if value)
        inverse = pow(first, -1, P)
        normalized = tuple(value * inverse % P for value in matrix)
        if matrix == normalized:
            rows.append(matrix)
    return tuple(rows)


def projective_1d_maps(source, target, matrices):
    source = frozenset(source)
    target = frozenset(target)
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


def factor_description(mask):
    first, second, mass, cartesian = endpoint.cartesian_description(mask)
    require(cartesian, "canonical endpoint mask stopped being rectangular")
    return frozenset(first), frozenset(second), mass


def mask_points(mask):
    return frozenset(
        point
        for point, value in zip(endpoint.KEYS, mask)
        if value
    )


def point_mask(points):
    value = 0
    for x, y in points:
        value |= 1 << (P * x + y)
    return value


def translate_bit_mask(mask, translation):
    tx, ty = translation
    result = 0
    for x in range(P):
        for y in range(P):
            if (mask >> (P * x + y)) & 1:
                result |= 1 << (
                    P * ((x + tx) % P) + ((y + ty) % P)
                )
    return result


def linear_bit_mask(points, matrix):
    return point_mask(matrix_point(matrix, point) for point in points)


def gl2():
    return tuple(
        (a, b, c, d)
        for a in range(P)
        for b in range(P)
        for c in range(P)
        for d in range(P)
        if (a * d - b * c) % P
    )


def is_monomial(matrix):
    a, b, c, d = matrix
    return (
        a != 0 and d != 0 and b == 0 and c == 0
    ) or (
        b != 0 and c != 0 and a == 0 and d == 0
    )


def target_translate_bank(target_points):
    target = point_mask(target_points)
    return {
        translation: translate_bit_mask(
            target,
            ((-translation[0]) % P, (-translation[1]) % P),
        )
        for translation in endpoint.KEYS
    }


def best_affine_single(source_points, target_points, matrices):
    """Return exact maps, maximum overlap, and lexicographically first witness."""
    target_bank = target_translate_bank(target_points)
    target_size = len(target_points)
    source_size = len(source_points)
    best_overlap = -1
    best_witness = None
    exact = []
    for matrix in matrices:
        linear = linear_bit_mask(source_points, matrix)
        for translation in endpoint.KEYS:
            overlap = (
                linear & target_bank[translation]
            ).bit_count()
            if overlap > best_overlap:
                best_overlap = overlap
                best_witness = (matrix, translation)
            if overlap == source_size == target_size:
                exact.append((matrix, translation))
    return {
        "exact": tuple(exact),
        "best_overlap": best_overlap,
        "min_hamming": source_size + target_size - 2 * best_overlap,
        "witness": best_witness,
    }


def best_affine_pair(
    source_pair,
    target_pair,
    matrices,
):
    """Use one affine map on both masks and maximize total intersection."""
    source_left, source_right = source_pair
    target_left, target_right = target_pair
    left_bank = target_translate_bank(target_left)
    right_bank = target_translate_bank(target_right)
    source_total = len(source_left) + len(source_right)
    target_total = len(target_left) + len(target_right)
    best_overlap = -1
    best_witness = None
    exact = []
    for matrix in matrices:
        linear_left = linear_bit_mask(source_left, matrix)
        linear_right = linear_bit_mask(source_right, matrix)
        for translation in endpoint.KEYS:
            overlap = (
                linear_left & left_bank[translation]
            ).bit_count() + (
                linear_right & right_bank[translation]
            ).bit_count()
            if overlap > best_overlap:
                best_overlap = overlap
                best_witness = (matrix, translation)
            if overlap == source_total == target_total:
                exact.append((matrix, translation))
    return {
        "exact": tuple(exact),
        "best_overlap": best_overlap,
        "min_total_hamming":
            source_total + target_total - 2 * best_overlap,
        "witness": best_witness,
    }


def best_affine_candidates(source_points, target_points, candidates):
    target_mask = point_mask(target_points)
    source_size = len(source_points)
    target_size = len(target_points)
    best_overlap = -1
    best_witness = None
    exact = []
    for matrix, translation in candidates:
        image = translate_bit_mask(
            linear_bit_mask(source_points, matrix),
            translation,
        )
        overlap = (image & target_mask).bit_count()
        if overlap > best_overlap:
            best_overlap = overlap
            best_witness = (matrix, translation)
        if overlap == source_size == target_size:
            exact.append((matrix, translation))
    return {
        "exact": tuple(exact),
        "best_overlap": best_overlap,
        "min_hamming": source_size + target_size - 2 * best_overlap,
        "witness": best_witness,
    }


def affine_defect(source_points, target_points, witness):
    matrix, translation = witness
    image = frozenset(
        affine_point(matrix, translation, point)
        for point in source_points
    )
    return {
        "lost_from_image": tuple(sorted(image - target_points)),
        "gained_in_target": tuple(sorted(target_points - image)),
    }


def vertical_character_image(matrix):
    """Image of the vertical covector (0,3) under affine pushforward.

    If y=Mx+t and chi_k(x)=omega^(k.x), then pushforward carries k to
    k M^(-1); translation changes only its scalar normalization.
    """
    a, b, c, d = matrix
    determinant = (a * d - b * c) % P
    inverse_det = pow(determinant, -1, P)
    inverse = (
        d * inverse_det % P,
        -b * inverse_det % P,
        -c * inverse_det % P,
        a * inverse_det % P,
    )
    ia, ib, ic, id_ = inverse
    return (3 * ic % P, 3 * id_ % P)


def named_origin_image(matrix, translation):
    return tuple(
        affine_point(matrix, translation, origin)
        for origin in ORIGINS
    )


def complement(values):
    return tuple(sorted(set(range(P)) - set(values)))


def sumset(first, second, left_scalar, right_scalar):
    return {
        (left_scalar * x + right_scalar * y) % P
        for x in first
        for y in second
    }


def main():
    allocation = horn.allocation
    _module, full, _details, _e3, _clocks, _q_pairs = (
        allocation.build_geometry()
    )
    period = full.T
    unit = period // P
    atom = allocation.ATOM_INTERVAL
    translated_atom = (
        atom[0] + 66 * endpoint.H,
        atom[1] + 66 * endpoint.H,
    )

    masks = {}
    descriptions = {}
    points = {}
    for q in QS:
        for step, base in zip(STEPS, (atom, translated_atom)):
            interval = horn.circular_shift_interval(
                base, q * unit, period
            )
            mask = endpoint.endpoint_mask(interval)
            masks[(q, step)] = mask
            descriptions[(q, step)] = factor_description(mask)
            points[(q, step)] = mask_points(mask)

    expected = {
        (0, 2): (
            frozenset((0, 1, 2, 3, 4, 5, 6, 7, 12)),
            frozenset((0, 1, 2, 3, 4, 5, 8, 9, 10)),
            81,
        ),
        (0, 68): (
            frozenset((0, 1, 2, 3, 4, 5, 6, 7, 12)),
            frozenset((0, 1, 2, 5, 6, 7, 8, 9, 10)),
            81,
        ),
        (3, 2): (
            frozenset(range(10)),
            frozenset((0, 3, 4, 5, 8, 9, 10, 11, 12)),
            90,
        ),
        (3, 68): (
            frozenset(range(10)),
            frozenset((0, 5, 6, 7, 8, 9, 10, 11, 12)),
            90,
        ),
        (11, 2): (
            frozenset((0, 1, 2, 3, 4, 5, 8, 9, 12)),
            frozenset((0, 1, 2, 3, 4, 5, 8, 11, 12)),
            81,
        ),
        (11, 68): (
            frozenset((0, 1, 2, 3, 4, 5, 8, 9, 12)),
            frozenset((0, 1, 2, 5, 6, 7, 8, 11, 12)),
            81,
        ),
    }
    require(
        descriptions == expected,
        "canonical q0/q3/q11 rectangle bank changed",
    )

    factor_names = {
        "q0_A": expected[(0, 2)][0],
        "q0_B2": expected[(0, 2)][1],
        "q0_B68": expected[(0, 68)][1],
        "q3_A": expected[(3, 2)][0],
        "q3_B2": expected[(3, 2)][1],
        "q3_B68": expected[(3, 68)][1],
        "q11_A": expected[(11, 2)][0],
        "q11_B2": expected[(11, 2)][1],
        "q11_B68": expected[(11, 68)][1],
    }
    factor_map_table = {
        (left_name, right_name):
            affine_1d_maps(left, right)
        for left_name, left in factor_names.items()
        for right_name, right in factor_names.items()
    }
    pgl1_matrices = pgl1()
    require(len(pgl1_matrices) == 2_184, "PGL2(F13) census changed")
    projective_factor_table = {
        (left_name, right_name):
            projective_1d_maps(left, right, pgl1_matrices)
        for left_name, left in factor_names.items()
        for right_name, right in factor_names.items()
    }

    matrices = gl2()
    monomial_matrices = tuple(
        matrix for matrix in matrices if is_monomial(matrix)
    )
    require(
        len(matrices) == 26_208
        and len(monomial_matrices) == 288,
        "GL2 or monomial matrix census changed",
    )
    chi3_matrices = tuple(
        matrix for matrix in matrices
        if vertical_character_image(matrix) == (0, 3)
    )
    axis_chi3_matrices = tuple(
        matrix for matrix in chi3_matrices
        if matrix[1] == 0 and matrix[2] == 0
    )
    axis_matrices = tuple(
        matrix for matrix in matrices
        if matrix[1] == 0 and matrix[2] == 0
    )
    axis_candidates = tuple(
        (matrix, translation)
        for matrix in axis_matrices
        for translation in endpoint.KEYS
    )
    chi3_candidates = tuple(
        (matrix, translation)
        for matrix in chi3_matrices
        for translation in endpoint.KEYS
    )
    axis_chi3_candidates = tuple(
        (matrix, translation)
        for matrix in axis_chi3_matrices
        for translation in endpoint.KEYS
    )
    origin_fixed_candidates = tuple(
        (matrix, translation)
        for matrix in matrices
        for translation in ((0, 0),)
        if named_origin_image(matrix, translation) == ORIGINS
    )
    chi3_origin_candidates = tuple(
        row for row in origin_fixed_candidates
        if vertical_character_image(row[0]) == (0, 3)
    )
    require(
        len(chi3_matrices) == 156
        and len(axis_matrices) == 144
        and len(axis_chi3_matrices) == 12
        and len(origin_fixed_candidates) == 156
        and len(chi3_origin_candidates) == 13,
        "typed affine subgroup census changed",
    )

    # Cauchy-Davenport is sharp enough here: a row using both input
    # coordinates has full 13-point projection for every 9x9 source
    # rectangle, since 9+9-1>13.  Verify that exact finite statement on
    # every q0 factor pair and every nonzero mixed row.
    q0_first, q0_second, _mass = expected[(0, 2)]
    mixed_projection_sizes = {
        len(sumset(q0_first, q0_second, left, right))
        for left in range(1, P)
        for right in range(1, P)
    }
    require(
        mixed_projection_sizes == {13},
        "q0 mixed projection stopped being all of F13",
    )

    single_searches = {
        ((0, 2), (3, 2)): best_affine_single(
            points[(0, 2)], points[(3, 2)], matrices
        ),
        ((0, 2), (11, 2)): best_affine_single(
            points[(0, 2)], points[(11, 2)], matrices
        ),
    }

    direct_edge_searches = {
        q: best_affine_single(
            points[(q, 2)],
            points[(q, 68)],
            matrices,
        )
        for q in QS
    }
    q0_pair_searches = {
        q: best_affine_pair(
            (points[(0, 2)], points[(0, 68)]),
            (points[(q, 2)], points[(q, 68)]),
            matrices,
        )
        for q in (3, 11)
    }
    constrained_direct_edges = {
        q: {
            "axis": best_affine_candidates(
                points[(q, 2)],
                points[(q, 68)],
                axis_candidates,
            ),
            "chi3": best_affine_candidates(
                points[(q, 2)],
                points[(q, 68)],
                chi3_candidates,
            ),
            "axis_chi3": best_affine_candidates(
                points[(q, 2)],
                points[(q, 68)],
                axis_chi3_candidates,
            ),
            "origins": best_affine_candidates(
                points[(q, 2)],
                points[(q, 68)],
                origin_fixed_candidates,
            ),
            "chi3_and_origins": best_affine_candidates(
                points[(q, 2)],
                points[(q, 68)],
                chi3_origin_candidates,
            ),
        }
        for q in QS
    }
    constrained_static_bridges = {
        q: {
            "axis": best_affine_candidates(
                points[(0, 2)],
                points[(q, 2)],
                axis_candidates,
            ),
            "chi3": best_affine_candidates(
                points[(0, 2)],
                points[(q, 2)],
                chi3_candidates,
            ),
            "axis_chi3": best_affine_candidates(
                points[(0, 2)],
                points[(q, 2)],
                axis_chi3_candidates,
            ),
            "origins": best_affine_candidates(
                points[(0, 2)],
                points[(q, 2)],
                origin_fixed_candidates,
            ),
            "chi3_and_origins": best_affine_candidates(
                points[(0, 2)],
                points[(q, 2)],
                chi3_origin_candidates,
            ),
        }
        for q in (3, 11)
    }

    # Exact q0 edge maps can be typed more finely.
    q0_exact = direct_edge_searches[0]["exact"]
    require(
        len(q0_exact) == 8
        and all(is_monomial(matrix) for matrix, _translation in q0_exact),
        "q0 full-affine edge map census changed",
    )
    character_three_maps = tuple(
        row for row in q0_exact
        if vertical_character_image(row[0]) == (0, 3)
    )
    fixed_origin_maps = tuple(
        row for row in q0_exact
        if named_origin_image(*row) == ORIGINS
    )
    setwise_origin_maps = tuple(
        row for row in q0_exact
        if set(named_origin_image(*row)) == set(ORIGINS)
    )
    require(
        len(character_three_maps) == 2
        and len(fixed_origin_maps) == 0
        and fixed_origin_maps == setwise_origin_maps,
        "q0 character/origin compatible map census changed",
    )
    require(
        not set(character_three_maps) & set(fixed_origin_maps),
        "q0 acquired a map preserving both chi3 and named origins",
    )

    # The only semilinear Frobenius automorphism of the prime field F13 is
    # x -> x^13 = x.  Verify it pointwise rather than treating it as an
    # additional search family.
    require(
        all(pow(value, P, P) == value for value in range(P)),
        "prime-field Frobenius ceased to be the identity",
    )

    # Exact no-go expectations.
    require(
        not single_searches[((0, 2), (3, 2))]["exact"],
        "an 81-point q0 mask became affine-equivalent to a 90-point q3 mask",
    )
    require(
        not single_searches[((0, 2), (11, 2))]["exact"],
        "a q0 mask became affine-equivalent to a q11 mask",
    )
    require(
        not direct_edge_searches[3]["exact"]
        and not direct_edge_searches[11]["exact"]
        and not q0_pair_searches[3]["exact"]
        and not q0_pair_searches[11]["exact"],
        "q3/q11 acquired an affine endpoint edge or q0 pair bridge",
    )

    source_complement = complement(factor_names["q0_A"])
    target_complements = {
        name: complement(factor_names[name])
        for name in ("q11_A", "q11_B2", "q11_B68")
    }
    selected_defects = {
        "q0_step2_to_q3_step2": affine_defect(
            points[(0, 2)],
            points[(3, 2)],
            single_searches[((0, 2), (3, 2))]["witness"],
        ),
        "q3_direct_edge": affine_defect(
            points[(3, 2)],
            points[(3, 68)],
            direct_edge_searches[3]["witness"],
        ),
        "q11_direct_edge": affine_defect(
            points[(11, 2)],
            points[(11, 68)],
            direct_edge_searches[11]["witness"],
        ),
    }
    # A nontrivial vertical dilation can be made to preserve the chi3 label
    # only by enlarging "semilinear" from the prime-field state space to a
    # cyclotomic coefficient automorphism omega -> omega^s.  Such an
    # automorphism necessarily moves THM-2863's order-182 left Prony node.
    xi_order = 2_366
    coefficient_recharts = {}
    for slope in (2, 6):
        units = tuple(
            value for value in range(1, xi_order)
            if gcd(value, xi_order) == 1 and value % P == slope
        )
        fixes_left_node = tuple(
            value for value in units
            if 13 * (value - 1) % xi_order == 0
        )
        fixes_right_node = tuple(
            value for value in units
            if 169 * (value - 1) % xi_order == 0
        )
        coefficient_recharts[slope] = {
            "unit_count": len(units),
            "fixes_left_node": fixes_left_node,
            "fixes_right_node": fixes_right_node,
        }
    require(
        all(
            row["unit_count"] == 78
            and not row["fixes_left_node"]
            and len(row["fixes_right_node"]) == 13
            for row in coefficient_recharts.values()
        ),
        "cyclotomic character-rechart/Prony-node boundary changed",
    )

    source_tree = ast.parse(Path(__file__).read_text())
    require(
        not any(
            isinstance(node, ast.Assert)
            for node in ast.walk(source_tree)
        ),
        "executable Python assert statement found",
    )

    print("THM-2877 SEMILINEAR ENDPOINT RECTANGLE CLASSIFICATION")
    print(
        f"pinned_scripts={tuple(PINS.items())};"
        f"pinned_outputs={tuple(RESULT_PINS.items())}"
    )
    print(
        "rectangle_bank="
        + repr({
            key: (
                tuple(sorted(first)),
                tuple(sorted(second)),
                mass,
            )
            for key, (first, second, mass) in descriptions.items()
        })
    )
    print(
        f"group_census=GL2:{len(matrices)},"
        f"monomial_GL2:{len(monomial_matrices)},"
        f"AGL2:{len(matrices) * P * P};"
        f"PGL2_line:{len(pgl1_matrices)};"
        f"chi3_linear:{len(chi3_matrices)},"
        f"axis_linear:{len(axis_matrices)},"
        f"axis_chi3_linear:{len(axis_chi3_matrices)},"
        f"origin_fixed_affine:{len(origin_fixed_candidates)},"
        f"chi3_origin_affine:{len(chi3_origin_candidates)};"
        "F13_Frobenius=id"
    )
    print(
        "q0_mixed_row_projection_sizes="
        f"{tuple(sorted(mixed_projection_sizes))};"
        "therefore_exact_9x9_rectangle_images_are_monomial"
    )
    print(
        f"q0_horizontal_complement={source_complement};"
        f"q11_factor_complements={target_complements}"
    )
    print(
        "q0_to_q3_single="
        + repr(single_searches[((0, 2), (3, 2))])
    )
    print(
        "q0_to_q11_single="
        + repr(single_searches[((0, 2), (11, 2))])
    )
    print(f"direct_edge_searches={direct_edge_searches}")
    print(f"q0_pair_bridge_searches={q0_pair_searches}")
    print(f"constrained_direct_edges={constrained_direct_edges}")
    print(f"constrained_static_bridges={constrained_static_bridges}")
    print(f"selected_nearest_defects={selected_defects}")
    print(f"cyclotomic_character_recharts={coefficient_recharts}")
    print(
        "q0_exact_edge_maps="
        f"{q0_exact};"
        f"chi3_preserving={character_three_maps};"
        f"named_origins_pointwise={fixed_origin_maps}"
    )
    print(
        "q0_chi3_origin_intersection=empty;"
        "chi3 maps move named origins to="
        f"{tuple(named_origin_image(*row) for row in character_three_maps)};"
        "no exact affine q0 edge map fixes the named origins even setwise"
    )
    print(
        "selected_factor_map_table="
        + repr({
            key: factor_map_table[key]
            for key in (
                ("q0_A", "q11_A"),
                ("q0_A", "q11_B2"),
                ("q0_A", "q11_B68"),
                ("q0_B2", "q11_A"),
                ("q0_B2", "q11_B2"),
                ("q0_B2", "q11_B68"),
                ("q0_B68", "q11_A"),
                ("q0_B68", "q11_B2"),
                ("q0_B68", "q11_B68"),
                ("q0_B2", "q0_B68"),
                ("q3_B2", "q3_B68"),
                ("q11_B2", "q11_B68"),
            )
        })
    )
    print(
        "selected_projective_factor_map_table="
        + repr({
            key: projective_factor_table[key]
            for key in (
                ("q0_A", "q11_A"),
                ("q0_B2", "q11_B2"),
                ("q0_B68", "q11_B68"),
                ("q3_B2", "q3_B68"),
                ("q11_B2", "q11_B68"),
                ("q0_B2", "q0_B68"),
            )
        })
    )
    print(
        "TYPE_VERDICT=q3 fails cardinality before geometry; "
        "q11 fails full AGL2 because Cauchy-Davenport forces monomial "
        "rows and the required factor equivalences are absent; "
        "prime-field semilinearity adds nothing; coefficient phases do not "
        "change support"
    )
    print(
        "SIDECAR_VERDICT=q0,q3,q11 lie in the E3 block on the horn, but "
        "an endpoint-address affine map does not transport physical "
        "interval, U ancestry, owner/current provenance, or the q7 "
        "complement block"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
