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
    tile_index,
    tile_schema,
)


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
    omega_plain = {}
    omega = {}
    omega_b2 = {}
    b2 = {}
    b3 = {}
    b23 = {}
    groups = {}
    colour_atoms = Counter()
    mirror_blue_failures = 0

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
        omega_plain[line] = upper_pair + sum(face_pairs, ())
        omega[line] = omega_plain[line] + (word,)
        b2[line] = b2_signature(mask, n)
        b3[line] = b3_signature(mask, n)
        b23[line] = (b2[line], b3[line])
        omega_b2[line] = (omega[line], b2[line])
        u, v = sorted(upper_pair)
        groups[line] = (word[0], u, v)
        mirror_blue_failures += (b2_skew_depth(b2[line], n) == 0) != (word[0] == "B")

    predicted_atoms = predicted_four_role_atoms(n)
    predicted_atoms += Counter()  # drop no keys yet; Counter equality ignores explicit zeroes
    assert colour_atoms == +predicted_atoms
    assert mirror_blue_failures == 0
    assert compact_partition(omega_b2)["cells"] == line_count

    refinements = {
        "B2": refinement_summary(groups, b2),
        "B3": refinement_summary(groups, b3),
        "B2_join_B3": refinement_summary(groups, b23),
    }
    result = {
        "n": n,
        "lines": line_count,
        "line_cech_descent": compatible_line_triples(n),
        "Xi": compact_partition(xi),
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
        "node_branch_correspondences": branch_census(n, node_by_mask),
    }
    carriers = {
        "node_pair_colour": groups,
        "B3_address": b3,
        "B2_address": b2,
        "B2_B3_join": b23,
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
                f"  UABC atoms={size['four_role_colour_atoms_UABC']}",
                f"  cubic cumulant={size['pure_cubic_three_role_colour_law']['third_cumulant']}",
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
