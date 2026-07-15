#!/usr/bin/env python3
"""Exact B3/B2 address audit for THM-801 and HYP-6880.

The old seven-term staircase recursion is lifted from cardinalities to three
actual restriction maps.  In THM-553 pin coordinates these are the high end
face A, the gap-contraction face B, and the low end face C.  Their Cech nerve
reconstructs a tiling.  After quotienting by all-tile complement, the same
three faces give an exact line descent from n >= 6; the nonempty triple
overlap is the phase witness.

The script also:

* adds the missing B face to THM-796's Xi tensor, producing Omega;
* attaches the mirror crossing-layer B2 sidecar and the coarse B3 address;
* disintegrates those signatures over literal coloured node-pair fibres;
* audits the pure cubic upper/low/high blue interaction and the four-role
  upper/A/B/C colour law;
* treats information carriers, not runners, as the Tournament Analysis
  vertices.

The B face is a valid lower staircase tiling, but not induced deletion of one
tournament vertex.  It contracts the gap coordinate of every retained
interval root.  That preservation boundary is intentional.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from three_sorted_metagraph_recursion_codex_S9 import (
    SMALL_ATLAS,
    N7_ATLAS,
    apex_zero_endpoint,
    blue_fraction_exponent,
    carrier_tournament,
    complement,
    is_grid_symmetric,
    line_index,
    load_atlases,
    m_tiles,
    reflection_orbits,
    reflect,
    tile_index,
    tile_schema,
)
from tournament_tiling_metagraph_address_codex_S4 import c3_count, tiling_tournament


OUTPUT = Path("05-knowledge/results/mobius_cech_metagraph_codec_codex_S12.out")
JSON_OUTPUT = Path("05-knowledge/results/mobius_cech_metagraph_codec_codex_S12.json")
FACE_ORDER = ("A", "B", "C")
B3_WORDS = ("A", "B", "C", "AB", "AC", "BC", "ABC")


def face_target(n: int, side: str, tile: tuple[int, int]) -> tuple[int, int] | None:
    """THM-553's three size-(n-1) subtriangle coordinates."""
    a, b = tile
    if side == "A" and a < n:
        return a, b
    if side == "B" and a - b >= 3:
        return a - 1, b
    if side == "C" and b >= 2:
        return a - 1, b - 1
    return None


def b3_face_mask(mask: int, n: int, side: str) -> int:
    lower_index = tile_index(n - 1)
    value = 0
    for bit, tile in enumerate(tile_schema(n)[0]):
        target = face_target(n, side, tile)
        if target is not None:
            value |= ((mask >> bit) & 1) << lower_index[target]
    return value


def b3_faces(mask: int, n: int) -> tuple[int, int, int]:
    return tuple(b3_face_mask(mask, n, side) for side in FACE_ORDER)  # type: ignore[return-value]


def reconstruct_b3(faces: tuple[int, int, int], n: int) -> int:
    """Glue a compatible A/B/C face tuple; reject inconsistent overlaps."""
    lower_index = tile_index(n - 1)
    value = 0
    for bit, tile in enumerate(tile_schema(n)[0]):
        seen = []
        for face, side in zip(faces, FACE_ORDER):
            target = face_target(n, side, tile)
            if target is not None:
                seen.append((face >> lower_index[target]) & 1)
        if not seen or len(set(seen)) != 1:
            raise ValueError("incompatible B3 face tuple")
        value |= seen[0] << bit
    return value


def b3_line_face(line: int, n: int, side: str) -> int:
    return line_index(b3_face_mask(line, n, side), n - 1)


def pair_phase(e1: int, side1: str, e2: int, side2: str, n: int) -> int:
    """Relative canonical-endpoint phase on one pair overlap."""
    x = b3_face_mask(e1, n, side1)
    y = b3_face_mask(e2, n, side2)
    if x == y:
        return 0
    if x == complement(y, n - 1):
        return 1
    raise ValueError("bare overlap lines do not agree")


def compatible_line_triples(n: int) -> dict:
    """Enumerate the bare-line Cech matching object at upper size n."""
    lower_n = n - 1
    lower_lines = range(1 << max(0, m_tiles(lower_n) - 1))
    by_a_overlap: dict[int, list[int]] = defaultdict(list)
    for line in lower_lines:
        by_a_overlap[b3_line_face(line, lower_n, "A")].append(line)

    compatible = 0
    holonomy = Counter()
    for e_a in lower_lines:
        e_b_candidates = by_a_overlap[b3_line_face(e_a, lower_n, "B")]
        e_c_candidates = by_a_overlap[b3_line_face(e_a, lower_n, "C")]
        for e_b in e_b_candidates:
            phase_ab = pair_phase(e_a, "B", e_b, "A", lower_n)
            bc_target = b3_line_face(e_b, lower_n, "C")
            for e_c in e_c_candidates:
                if bc_target != b3_line_face(e_c, lower_n, "B"):
                    continue
                phase_ac = pair_phase(e_a, "C", e_c, "A", lower_n)
                phase_bc = pair_phase(e_b, "C", e_c, "B", lower_n)
                holonomy[phase_ab ^ phase_ac ^ phase_bc] += 1
                compatible += 1

    image = Counter(
        tuple(b3_line_face(line, n, side) for side in FACE_ORDER)
        for line in range(1 << (m_tiles(n) - 1))
    )
    return {
        "compatible_bare_triples": compatible,
        "holonomy_histogram": counter_json(holonomy),
        "upper_line_image_cells": len(image),
        "upper_line_image_fibre_histogram": counter_json(Counter(image.values())),
        "is_bijection_onto_matching_object": (
            compatible == len(image) and set(image.values()) == {1}
        ),
    }


def full_ie_word(n: int, tile: tuple[int, int]) -> str:
    return "".join(side for side in FACE_ORDER if face_target(n, side, tile) is not None)


def b3_signature(mask: int, n: int) -> tuple[int, ...]:
    counts = Counter()
    for bit, tile in enumerate(tile_schema(n)[0]):
        counts[full_ie_word(n, tile)] += (mask >> bit) & 1
    return tuple(counts[word] for word in B3_WORDS)


def b2_signature(mask: int, n: int) -> tuple[int, ...]:
    """Mirror-pair bit counts by the THM-553 crossing clock tau."""
    tiles, sigma = tile_schema(n)
    layers: dict[int, Counter[tuple[int, int]]] = defaultdict(Counter)
    fixed = Counter()
    for bit, (a, b) in enumerate(tiles):
        tau = a + b - 1
        if tau < n:
            layers[tau][((mask >> bit) & 1, (mask >> sigma[bit]) & 1)] += 1
        elif tau == n:
            assert sigma[bit] == bit
            fixed[(mask >> bit) & 1] += 1
    result = []
    for tau in range(3, n):
        result.extend(layers[tau][state] for state in ((0, 0), (0, 1), (1, 0), (1, 1)))
    result.extend((fixed[0], fixed[1]))
    return tuple(result)


def b2_skew_depth(signature: tuple[int, ...], n: int) -> int:
    return sum(signature[4 * (tau - 3) + 1] + signature[4 * (tau - 3) + 2] for tau in range(3, n))


def primitive(row: Counter[int]) -> tuple[tuple[int, int], ...]:
    divisor = 0
    for value in row.values():
        divisor = math.gcd(divisor, value)
    return tuple(sorted((key, value // divisor) for key, value in row.items()))


def compact_partition(values: dict[int, object]) -> dict[str, int]:
    cells = Counter(values.values())
    total = sum(cells.values())
    return {
        "cells": len(cells),
        "collision_cells": sum(size > 1 for size in cells.values()),
        "collision_excess": total - len(cells),
        "collision_pairs": sum(size * (size - 1) // 2 for size in cells.values()),
        "max_multiplicity": max(cells.values()),
        "separated_pairs": (total * total - sum(size * size for size in cells.values())) // 2,
    }


def counter_json(counter: Counter) -> dict[str, int]:
    return {str(key): value for key, value in sorted(counter.items(), key=lambda item: str(item[0]))}


def key_lines(values: dict[int, object], limit: int = 24) -> list[dict]:
    cells: dict[object, list[int]] = defaultdict(list)
    for line, value in values.items():
        cells[value].append(line)
    return [
        {"multiplicity": len(lines), "line_indices": lines}
        for lines in cells.values()
        if len(lines) > 1
    ][:limit]


def blue_word(mask: int, n: int, faces: tuple[int, int, int]) -> str:
    _, sigma = tile_schema(n)
    _, lower_sigma = tile_schema(n - 1)
    return ("B" if is_grid_symmetric(mask, sigma) else "K") + "".join(
        "B" if is_grid_symmetric(face, lower_sigma) else "K" for face in faces
    )


def predicted_four_role_atoms(n: int) -> Counter[str]:
    total = 1 << (m_tiles(n) - 1)
    upper = 1 << (reflection_orbits(n) - 1)
    face = 1 << (reflection_orbits(n - 1) + n - 3)
    face_pair = 1 << (n - 3 + (n - 2) // 2)
    upper_end = 1 << (n - 3)
    face_triple = 1 << (2 * ((n - 2) // 2))
    return Counter(
        {
            "BBBB": upper_end,
            "BKBK": upper - upper_end,
            "KBBB": face_triple - upper_end,
            "KBBK": face_pair - face_triple,
            "KBKB": face_pair - face_triple,
            "KBKK": face - 2 * face_pair + face_triple,
            "KKBB": face_pair - face_triple,
            "KKBK": face - 2 * face_pair + face_triple - upper + upper_end,
            "KKKB": face - 2 * face_pair + face_triple,
            "KKKK": total - 3 * face + 3 * face_pair - face_triple,
        }
    )


def cubic_colour_record(n: int) -> dict:
    total = 1 << (m_tiles(n) - 1)
    upper = 1 << (reflection_orbits(n) - 1)
    one_face = 1 << (reflection_orbits(n - 1) + n - 3)
    kappa = Fraction(upper * one_face * (total - one_face), total**3)
    return {
        "a_beta_n": f"{upper}/{total}",
        "b_beta_previous": f"{one_face}/{total}",
        "third_cumulant": f"{kappa.numerator}/{kappa.denominator}",
        "integral_numerator_U_L_TminusL": upper * one_face * (total - one_face),
        "product_plus_cubic_pgf": (
            "((1-a)+a*x)((1-b)+b*y)((1-b)+b*z)"
            "+a*b*(1-b)(x-1)(y-1)(z-1)"
        ),
        "exact_atom_defect": "(-1)^(3-|R|)*a*b*(1-b)",
        "legendre_mask_word_1_to_7": "++-+--+",
        "blue_conditional_face_covariance": f"{Fraction(one_face, total) * (1 - Fraction(one_face, total))}",
        "black_conditional_face_covariance": (
            "undefined" if upper == total else
            f"{-Fraction(upper * one_face * (total - one_face), total * total * (total - upper))}"
        ),
    }


def smith_horizontal(mask: int, n: int) -> Fraction:
    """BSST row current using THM-790 orientation bits b=1-t_atlas."""
    return sum(
        (Fraction(1 - ((mask >> bit) & 1), n - 1 - y)
         for bit, (_x, y) in enumerate(tile_schema(n)[0])),
        Fraction(),
    )


def interval_mobius_square(function, mask: int, n: int):
    """The full Möbius primitive of the legal marked-path interval faces."""
    face_a, _face_b, face_c = b3_faces(mask, n)
    core = b3_face_mask(face_c, n - 1, "A")
    return (
        function(mask, n)
        - function(face_c, n - 1)
        - function(face_a, n - 1)
        + function(core, n - 2)
    )


def endpoint_epsilon(mask: int, n: int) -> int:
    """THM-785 epsilon in THM-796 atlas bits: bottom ones minus top ones."""
    index = tile_index(n)
    bottom = sum((mask >> index[(x, 1)]) & 1 for x in range(3, n))
    top = sum((mask >> index[(n, y)]) & 1 for y in range(2, n - 1))
    return bottom - top


@lru_cache(maxsize=None)
def tiling_c3(mask: int, n: int) -> int:
    tiles, _sigma = tile_schema(n)
    return c3_count(tiling_tournament(mask, n, tiles), n)


def boundary_c3_curvature(mask: int, n: int) -> int:
    """Omega C3: cyclic triples meeting both fixed-path endpoints."""
    return interval_mobius_square(tiling_c3, mask, n)


def predicted_q_histogram(n: int, blue: bool | None = None) -> Counter[int]:
    """Closed apex-zero line polynomial in q=Omega C3."""
    degree = n - 4
    total_factor = 1 << (m_tiles(n) - 2 * n + 7)
    total = Counter(
        {k: total_factor * math.comb(degree, k) * 3 ** (degree - k) for k in range(degree + 1)}
    )
    r = reflection_orbits(n)
    blue_hist = Counter()
    if n % 2 == 0:
        factor = 1 << (r - n + 3)
        for k in range((n - 4) // 2 + 1):
            blue_hist[2 * k] += factor * math.comb((n - 4) // 2, k) * 3 ** ((n - 4) // 2 - k)
    else:
        factor = 1 << (r - n + 3)
        for k in range((n - 5) // 2 + 1):
            base = factor * math.comb((n - 5) // 2, k) * 3 ** ((n - 5) // 2 - k)
            blue_hist[2 * k] += base
            blue_hist[2 * k + 1] += base
    if blue is True:
        return blue_hist
    if blue is False:
        return total - blue_hist
    return total


def predicted_q_pair_histogram(n: int, blue: bool | None = None) -> Counter[tuple[int, int]]:
    """Closed edge polynomial for (Omega C3 at apex-zero, at its complement)."""
    degree = n - 4
    factor = 1 << (m_tiles(n) - 2 * n + 5)
    total = Counter()
    for k in range(degree + 1):
        base = factor * math.comb(degree, k) * 3 ** (degree - k)
        for boundary, coefficient in ((0, 1), (1, 2), (2, 1)):
            total[(k, k + boundary)] += base * coefficient

    r = reflection_orbits(n)
    blue_hist = Counter()
    factor = 1 << (r - n + 2)
    half_degree = (n - 4) // 2 if n % 2 == 0 else (n - 5) // 2
    for k in range(half_degree + 1):
        base = factor * math.comb(half_degree, k) * 3 ** (half_degree - k)
        middle_choices = ((0, 0),) if n % 2 == 0 else ((0, 0), (1, 1))
        for middle_q0, middle_q1 in middle_choices:
            for boundary in (0, 2):
                blue_hist[(2 * k + middle_q0, 2 * k + middle_q1 + boundary)] += base
    if blue is True:
        return blue_hist
    if blue is False:
        return total - blue_hist
    return total


def branch_census(n: int, node_by_mask: dict[int, list[int]]) -> dict:
    upper_nodes = node_by_mask[n]
    lower_nodes = node_by_mask[n - 1]
    fibres: dict[int, list[int]] = defaultdict(list)
    for mask, node in enumerate(upper_nodes):
        fibres[node].append(mask)
    rows = {side: {node: Counter() for node in fibres} for side in FACE_ORDER}
    for node, masks in fibres.items():
        for mask in masks:
            for side in FACE_ORDER:
                rows[side][node][lower_nodes[b3_face_mask(mask, n, side)]] += 1

    assert rows["A"] == rows["C"]
    result = {"A_equals_C_failures": 0, "faces": {}}
    for side in FACE_ORDER:
        result["faces"][side] = {
            "support": compact_partition({node: tuple(sorted(row)) for node, row in rows[side].items()}),
            "weighted": compact_partition({node: tuple(sorted(row.items())) for node, row in rows[side].items()}),
            "primitive": compact_partition({node: primitive(row) for node, row in rows[side].items()}),
        }
    result["joined_A_B"] = {
        "support": compact_partition(
            {
                node: (tuple(sorted(rows["A"][node])), tuple(sorted(rows["B"][node])))
                for node in fibres
            }
        ),
        "primitive": compact_partition(
            {node: (primitive(rows["A"][node]), primitive(rows["B"][node])) for node in fibres}
        ),
    }
    result["A_equals_B_weighted_rows"] = sum(rows["A"][node] == rows["B"][node] for node in fibres)
    return result


def refinement_summary(
    group_by_line: dict[int, tuple], signature_by_line: dict[int, object]
) -> dict:
    cells = Counter((group_by_line[line], signature) for line, signature in signature_by_line.items())
    group_cells: dict[tuple, list[int]] = defaultdict(list)
    for (group, _signature), count in cells.items():
        group_cells[group].append(count)
    stratified = {}
    for colour in ("B", "K"):
        for loop_type in ("loop", "cross"):
            groups = {
                group: counts
                for group, counts in group_cells.items()
                if group[0] == colour and ("loop" if group[1] == group[2] else "cross") == loop_type
            }
            stratified[colour + "_" + loop_type] = {
                "fibres": len(groups),
                "subcells": sum(len(counts) for counts in groups.values()),
                "collision_excess": sum(sum(counts) - len(counts) for counts in groups.values()),
                "fully_separated_fibres": sum(all(count == 1 for count in counts) for counts in groups.values()),
            }
    return {
        "subcells": len(cells),
        "collision_excess": sum(cells.values()) - len(cells),
        "multiplicity_histogram": counter_json(Counter(cells.values())),
        "max_multiplicity": max(cells.values()),
        "fully_separated_fibres": sum(all(count == 1 for count in counts) for counts in group_cells.values()),
        "total_fibres": len(group_cells),
        "stratified": stratified,
    }


def size_census(n: int, node_by_mask: dict[int, list[int]]) -> tuple[dict, dict[str, dict[int, object]]]:
    upper_nodes = node_by_mask[n]
    lower_nodes = node_by_mask[n - 1]
    line_count = 1 << (m_tiles(n) - 1)
    xi = {}
    xi_core = {}
    omega_plain = {}
    omega = {}
    omega_b2 = {}
    b2 = {}
    b3 = {}
    b23 = {}
    groups = {}
    q_pair = {}
    q_colour_pair = {}
    q_all = Counter()
    q_blue = Counter()
    q_black = Counter()
    q_pair_all = Counter()
    q_pair_blue = Counter()
    q_pair_black = Counter()
    node_curvature: dict[int, Counter[int]] = defaultdict(Counter)
    node_c3: dict[int, int] = {}
    epsilon_all_tilings = Counter()
    black_epsilon_by_group = Counter()
    colour_atoms = Counter()
    mirror_blue_failures = 0
    smith_antisymmetric_failures = 0
    smith_symmetric_failures = 0

    for line in range(line_count):
        mask = apex_zero_endpoint(line, n)
        other = complement(mask, n)
        faces = b3_faces(mask, n)
        word = blue_word(mask, n, faces)
        colour_atoms[word] += 1
        upper_pair = (upper_nodes[mask], upper_nodes[other])
        face_pairs = tuple(
            (lower_nodes[face], lower_nodes[complement(face, n - 1)]) for face in faces
        )
        xi_word = word[0] + word[3] + word[1]
        xi[line] = upper_pair + face_pairs[2] + face_pairs[0] + (xi_word,)
        if n >= 5:
            core = b3_face_mask(faces[2], n - 1, "A")
            core_nodes = node_by_mask[n - 2]
            xi_core[line] = xi[line] + (
                core_nodes[core], core_nodes[complement(core, n - 2)]
            )
        omega_plain[line] = upper_pair + sum(face_pairs, ())
        omega[line] = omega_plain[line] + (word,)
        b2[line] = b2_signature(mask, n)
        b3[line] = b3_signature(mask, n)
        b23[line] = (b2[line], b3[line])
        omega_b2[line] = (omega[line], b2[line])
        u, v = sorted(upper_pair)
        groups[line] = (word[0], u, v)
        mirror_blue_failures += (b2_skew_depth(b2[line], n) == 0) != (word[0] == "B")

        q0 = boundary_c3_curvature(mask, n)
        q1 = boundary_c3_curvature(other, n)
        q_pair[line] = (q0, q1)
        q_colour_pair[line] = (word[0], q0, q1)
        q_all[q0] += 1
        q_pair_all[(q0, q1)] += 1
        (q_blue if word[0] == "B" else q_black)[q0] += 1
        (q_pair_blue if word[0] == "B" else q_pair_black)[(q0, q1)] += 1
        node_curvature[upper_pair[0]][q0] += 1
        node_curvature[upper_pair[1]][q1] += 1
        for endpoint, node, q_value in ((mask, upper_pair[0], q0), (other, upper_pair[1], q1)):
            c3_value = tiling_c3(endpoint, n)
            if node in node_c3:
                assert node_c3[node] == c3_value
            node_c3[node] = c3_value
            epsilon = endpoint_epsilon(endpoint, n)
            epsilon_all_tilings[epsilon] += 1
            omega_h = interval_mobius_square(smith_horizontal, endpoint, n)
            omega_v = interval_mobius_square(smith_horizontal, reflect(endpoint, n), n)
            denominator = (n - 2) * (n - 3)
            smith_antisymmetric_failures += denominator * (omega_h - omega_v) != epsilon
            apex_orientation = 1 - ((endpoint >> tile_index(n)[(n, 1)]) & 1)
            index = tile_index(n)
            bottom = sum((endpoint >> index[(x, 1)]) & 1 for x in range(3, n))
            top = sum((endpoint >> index[(n, y)]) & 1 for y in range(2, n - 1))
            leg_current = bottom + top + 2 * ((endpoint >> index[(n, 1)]) & 1) - (n - 2)
            smith_symmetric_failures += (
                denominator * (omega_h + omega_v)
                != leg_current + (2 * apex_orientation - 1) * (n - 2)
            )
        if word[0] == "K":
            black_epsilon_by_group[(groups[line], endpoint_epsilon(mask, n))] += 1

    predicted_atoms = predicted_four_role_atoms(n)
    predicted_atoms += Counter()  # drop no keys yet; Counter equality ignores explicit zeroes
    assert colour_atoms == +predicted_atoms
    assert mirror_blue_failures == 0
    assert compact_partition(omega_b2)["cells"] == line_count
    assert q_all == predicted_q_histogram(n)
    assert q_blue == predicted_q_histogram(n, True)
    assert q_black == predicted_q_histogram(n, False)
    assert q_pair_all == predicted_q_pair_histogram(n)
    assert q_pair_blue == predicted_q_pair_histogram(n, True)
    assert q_pair_black == predicted_q_pair_histogram(n, False)
    assert smith_antisymmetric_failures == smith_symmetric_failures == 0

    black_epsilon_symmetry_failures = 0
    black_zero_epsilon_parity_failures = 0
    black_groups = {group for group, _epsilon in black_epsilon_by_group}
    for group in black_groups:
        epsilons = {epsilon for one_group, epsilon in black_epsilon_by_group if one_group == group}
        for epsilon in epsilons:
            black_epsilon_symmetry_failures += (
                black_epsilon_by_group[(group, epsilon)]
                != black_epsilon_by_group[(group, -epsilon)]
            )
        black_zero_epsilon_parity_failures += black_epsilon_by_group[(group, 0)] % 2
    assert black_epsilon_symmetry_failures == black_zero_epsilon_parity_failures == 0

    curvature_polynomial = {
        node: tuple(sorted(row.items())) for node, row in node_curvature.items()
    }
    curvature_with_c3 = {
        node: (node_c3[node], curvature_polynomial[node]) for node in node_curvature
    }

    refinements = {
        "B2": refinement_summary(groups, b2),
        "B3": refinement_summary(groups, b3),
        "B2_join_B3": refinement_summary(groups, b23),
        "C3_boundary_pair": refinement_summary(groups, q_pair),
    }
    result = {
        "n": n,
        "lines": line_count,
        "line_cech_descent": compatible_line_triples(n),
        "Xi": compact_partition(xi),
        "Xi_join_common_core_nodes": compact_partition(xi_core) if xi_core else None,
        "Omega_without_colour": compact_partition(omega_plain),
        "Omega": compact_partition(omega),
        "Omega_join_B2": compact_partition(omega_b2),
        "Omega_collision_witnesses": key_lines(omega),
        "B2_signature": compact_partition(b2),
        "B3_signature": compact_partition(b3),
        "B2_join_B3_signature": compact_partition(b23),
        "projected_node_pair_fibre_refinements": refinements,
        "four_role_colour_atoms_UABC": dict(sorted(colour_atoms.items())),
        "four_role_closed_formula_failures": int(colour_atoms != +predicted_atoms),
        "B2_zero_skew_iff_blue_failures": mirror_blue_failures,
        "pure_cubic_three_role_colour_law": cubic_colour_record(n),
        "interval_face_mobius_curvature": {
            "definition": "Omega f=f_n-f_(n-1)o d_L-f_(n-1)o d_H+f_(n-2)o core",
            "C3_positive_recursion": "C3_n=C3_low+C3_high-C3_core+q; q counts cyclic triples meeting both path endpoints",
            "E4_recursion": "Omega E4=2(n-1)-8q",
            "q_histogram": counter_json(q_all),
            "blue_q_histogram": counter_json(q_blue),
            "black_q_histogram": counter_json(q_black),
            "q_pair_histogram": counter_json(q_pair_all),
            "blue_q_pair_histogram": counter_json(q_pair_blue),
            "black_q_pair_histogram": counter_json(q_pair_black),
            "q_polynomial": "2^(M-2n+7)(3+z)^(n-4)",
            "q_pair_polynomial": "2^(M-2n+5)(1+w)^2(3+zw)^(n-4)",
            "node_curvature_polynomial_partition": compact_partition(curvature_polynomial),
            "node_curvature_plus_C3_partition": compact_partition(curvature_with_c3),
            "black_projected_fibre_q_coefficient_parity_failures": sum(
                count % 2
                for (group, q_value), count in Counter(
                    (groups[line], q_pair[line][0]) for line in range(line_count)
                    if groups[line][0] == "K"
                ).items()
            ),
        },
        "smith_interval_mobius_curvature": {
            "orientation_bit_convention": "b_xy=1-t_atlas_xy, so b=1 means x->y as in THM-790",
            "antisymmetric_identity": "den*(Omega J_h-Omega J_v)=epsilon=e1+e_n-(n-2)",
            "symmetric_identity": "den*(Omega J_h+Omega J_v)=lambda+(2a-1)(n-2), lambda=e1-e_n",
            "identity_failures": smith_antisymmetric_failures + smith_symmetric_failures,
            "epsilon_histogram_all_tilings": counter_json(epsilon_all_tilings),
            "black_projected_fibre_sign_symmetry_failures": black_epsilon_symmetry_failures,
            "black_zero_epsilon_evenness_failures": black_zero_epsilon_parity_failures,
            "warning": "blue implies antisymmetric curvature zero, but balanced black lines also exist",
        },
        "node_branch_correspondences": branch_census(n, node_by_mask),
    }
    carriers = {
        "node_pair_colour": groups,
        "B3_address": b3,
        "B2_address": b2,
        "B2_B3_join": b23,
        "C3_boundary_pair": q_colour_pair,
        "Xi": xi,
        "Omega": omega,
        "Omega_B2_join": omega_b2,
        "exact_line": {line: line for line in range(line_count)},
    }
    return result, carriers


def tournament_analysis(carriers: dict[str, dict[int, object]]) -> dict:
    stats = {name: compact_partition(values) for name, values in carriers.items()}
    retention = carrier_tournament(stats, "retention")
    economy = carrier_tournament(stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carriers)) for j in range(i + 1, len(carriers))
    )
    return {
        "vertices": list(carriers),
        "pairwise_observable": "number of unordered literal-line pairs separated by a carrier",
        "switches": ["retention", "retention_per_log2_cells"],
        "tie_hamiltonian_path": list(carriers),
        "carrier_stats": stats,
        "retention": retention,
        "economy": economy,
        "edge_flips_between_gauges": flips,
    }


def verify_regressions(result: dict) -> None:
    expected = {
        4: (4, 4, 4),
        5: (32, 32, 32),
        6: (509, 510, 512),
        7: (16031, 16308, 16384),
    }
    refinement_expected = {
        4: {"B2": (4, 0, 3), "B3": (4, 0, 3), "B2_join_B3": (4, 0, 3)},
        5: {"B2": (32, 0, 20), "B3": (32, 0, 20), "B2_join_B3": (32, 0, 20)},
        6: {"B2": (504, 8, 179), "B3": (468, 44, 156), "B2_join_B3": (512, 0, 187)},
        7: {"B2": (16212, 172, 5962), "B3": (15016, 1368, 5071), "B2_join_B3": (16368, 16, 6110)},
    }
    for size in result["sizes"]:
        n = size["n"]
        assert (size["Xi"]["cells"], size["Omega"]["cells"], size["Omega_join_B2"]["cells"]) == expected[n]
        for name, expected_row in refinement_expected[n].items():
            row = size["projected_node_pair_fibre_refinements"][name]
            assert (row["subcells"], row["collision_excess"], row["fully_separated_fibres"]) == expected_row
    n6, n7 = result["sizes"][2], result["sizes"][3]
    assert n6["Xi_join_common_core_nodes"]["cells"] == 510
    assert n7["Xi_join_common_core_nodes"]["cells"] == 16110
    assert n7["interval_face_mobius_curvature"]["node_curvature_polynomial_partition"]["cells"] == 238
    assert n7["interval_face_mobius_curvature"]["node_curvature_plus_C3_partition"]["cells"] == 249
    result["regression_failures"] = 0


def render(result: dict) -> str:
    lines = [
        "THM-801 MOBIUS/CECH METAGRAPH CODEC",
        "=" * 82,
        "A=high endpoint face; B=gap-contraction face; C=low endpoint face",
        "The three faces cover every tile and reconstruct compatible tilings exactly.",
        "Omega adds B to THM-796 Xi; B2 is the mirror crossing-layer sidecar.",
        "",
    ]
    for size in result["sizes"]:
        n = size["n"]
        descent = size["line_cech_descent"]
        refine = size["projected_node_pair_fibre_refinements"]
        branch = size["node_branch_correspondences"]
        curvature = size["interval_face_mobius_curvature"]
        smith = size["smith_interval_mobius_curvature"]
        xi_core = size["Xi_join_common_core_nodes"]
        lines.extend(
            [
                f"n={n}: lines={size['lines']}",
                f"  line Cech matching/image/holonomy/bijection="
                f"{descent['compatible_bare_triples']}/{descent['upper_line_image_cells']}/"
                f"{descent['holonomy_histogram']}/{descent['is_bijection_onto_matching_object']}",
                f"  Xi cells={size['Xi']['cells']} hist(max/excess)="
                f"{size['Xi']['max_multiplicity']}/{size['Xi']['collision_excess']}",
                f"  Omega cells={size['Omega']['cells']} hist(max/excess)="
                f"{size['Omega']['max_multiplicity']}/{size['Omega']['collision_excess']}",
                f"  Omega+B2 cells={size['Omega_join_B2']['cells']} "
                f"max/excess={size['Omega_join_B2']['max_multiplicity']}/"
                f"{size['Omega_join_B2']['collision_excess']}",
                f"  Xi+core cells/excess="
                f"{('degenerate' if xi_core is None else str(xi_core['cells']) + '/' + str(xi_core['collision_excess']))}",
                f"  projected fibres B2(subcells/excess/separated)="
                f"{refine['B2']['subcells']}/{refine['B2']['collision_excess']}/"
                f"{refine['B2']['fully_separated_fibres']}/{refine['B2']['total_fibres']}",
                f"  projected fibres B3(subcells/excess/separated)="
                f"{refine['B3']['subcells']}/{refine['B3']['collision_excess']}/"
                f"{refine['B3']['fully_separated_fibres']}/{refine['B3']['total_fibres']}",
                f"  projected fibres B2+B3(subcells/excess/separated)="
                f"{refine['B2_join_B3']['subcells']}/{refine['B2_join_B3']['collision_excess']}/"
                f"{refine['B2_join_B3']['fully_separated_fibres']}/"
                f"{refine['B2_join_B3']['total_fibres']}",
                f"  projected fibres q-pair(subcells/excess/separated)="
                f"{refine['C3_boundary_pair']['subcells']}/"
                f"{refine['C3_boundary_pair']['collision_excess']}/"
                f"{refine['C3_boundary_pair']['fully_separated_fibres']}/"
                f"{refine['C3_boundary_pair']['total_fibres']}",
                f"  UABC atoms={size['four_role_colour_atoms_UABC']}",
                f"  cubic cumulant={size['pure_cubic_three_role_colour_law']['third_cumulant']}",
                f"  q=OmegaC3 histogram={curvature['q_histogram']}; edge (q0,q1)="
                f"{curvature['q_pair_histogram']}",
                f"  node K_u/K_u+C3 cells="
                f"{curvature['node_curvature_polynomial_partition']['cells']}/"
                f"{curvature['node_curvature_plus_C3_partition']['cells']} of "
                f"{size['node_branch_correspondences']['faces']['A']['weighted']['cells']}",
                f"  Smith Omega identities/parity failures="
                f"{smith['identity_failures']}/"
                f"{smith['black_projected_fibre_sign_symmetry_failures']}/"
                f"{smith['black_zero_epsilon_evenness_failures']}",
                f"  A/B support cells="
                f"{branch['faces']['A']['support']['cells']}/"
                f"{branch['faces']['B']['support']['cells']}; joined="
                f"{branch['joined_A_B']['support']['cells']}; primitive="
                f"{branch['faces']['A']['primitive']['cells']}/"
                f"{branch['faces']['B']['primitive']['cells']}/"
                f"{branch['joined_A_B']['primitive']['cells']}",
                "",
            ]
        )
    ta = result["tournament_analysis_n7"]
    lines.extend(
        [
            "EXACT FRONTIER",
            "  n>=6: the full B3 line descent is bijective; G is the phase witness.",
            "  n=7: adding the gap face removes 277 of Xi's 353 collision excess.",
            "  n=7: Omega+B2 is an exact 16,384-line codec; B2+B3 alone leaves 16 collisions.",
            "  all proper upper/endpoint-face colour marginals are independent;",
            "  the only connected defect is a*b*(1-b)(x-1)(y-1)(z-1).",
            "  Omega is the full Mobius primitive of legal marked-path interval faces.",
            "  Omega C3 counts cycles meeting both endpoints and gives a positive C3 recursion.",
            "  primal/dual Smith Omega modes recover both THM-785 coordinates (lambda,epsilon).",
            "  signed epsilon cancels in every black node-pair fibre; only |epsilon| can bias drift.",
            "",
            "TOURNAMENT ANALYSIS (information carriers as vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  pairwise observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; edge flips={ta['edge_flips_between_gauges']}",
            f"  tie Hamiltonian path={tuple(ta['tie_hamiltonian_path'])}",
            f"  retention fingerprints: scores={ta['retention']['score_hist']} "
            f"C3={ta['retention']['directed_3cycles']} SCC={ta['retention']['scc_sizes']} "
            f"Hpaths={ta['retention']['hamiltonian_paths']}",
            f"  regression failures={result['regression_failures']}",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=SMALL_ATLAS)
    parser.add_argument("--n7-atlas", type=Path, default=N7_ATLAS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()

    _sizes, node_by_mask = load_atlases(args.small_atlas, args.n7_atlas)
    result = {
        "schema_version": 1,
        "theorem": "THM-801",
        "hypothesis": "HYP-6880",
        "face_semantics": {
            "A": "high endpoint deletion: a<n, (a,b)->(a,b)",
            "B": "gap contraction: a-b>=3, (a,b)->(a-1,b)",
            "C": "low endpoint deletion: b>=2, (a,b)->(a-1,b-1)",
        },
        "sizes": [],
    }
    carriers_n7 = None
    for n in range(4, 8):
        # Exhaust the tiling-level Cech reconstruction independently of lines.
        for mask in range(1 << m_tiles(n)):
            assert reconstruct_b3(b3_faces(mask, n), n) == mask
        size, carriers = size_census(n, node_by_mask)
        result["sizes"].append(size)
        if n == 7:
            carriers_n7 = carriers
    assert carriers_n7 is not None
    result["tournament_analysis_n7"] = tournament_analysis(carriers_n7)
    verify_regressions(result)
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
