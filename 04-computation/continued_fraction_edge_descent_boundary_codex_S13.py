#!/usr/bin/env python3
"""Exact X6->X7 centered-CF edge-descent boundary for THM-813."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from continued_fraction_metagraph_edge_transport_codex_S13 import (
    centered_increments,
    centered_height,
    coordinate_audit,
    coordinate_copy,
    mobius_coefficients,
    rho_5_6,
    truncated_signature,
)
from mobius_cech_metagraph_codec_codex_S12 import compact_partition
from three_sorted_metagraph_recursion_codex_S9 import (
    N7_ATLAS,
    SMALL_ATLAS,
    complement,
    is_grid_symmetric,
    line_index,
    load_atlases,
    m_tiles,
    reflect,
    tile_index,
    tile_schema,
)
from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament


OUTPUT = Path("05-knowledge/results/continued_fraction_edge_descent_boundary_codex_S13.out")
JSON_OUTPUT = Path("05-knowledge/results/continued_fraction_edge_descent_boundary_codex_S13.json")


RHO_6_7 = (0, 1, 2, 2, 3, 4, 5, 6, 6, 7, 8, 5, 7, 8, 9)


def edge_key(mask: int, n: int, node_map: list[int]) -> tuple[str, int, int]:
    colour = "B" if is_grid_symmetric(mask, tile_schema(n)[1]) else "K"
    nodes = sorted((node_map[mask], node_map[complement(mask, n)]))
    return colour, nodes[0], nodes[1]


def phi(mask: int, source_n: int, target_n: int, rho: tuple[int, ...]) -> int:
    return coordinate_copy(mask, source_n, target_n, rho)


def line_edge_audit(
    source_n: int,
    target_n: int,
    rho: tuple[int, ...],
    node_by_mask: dict[int, list[int]],
) -> dict:
    source_map, target_map = node_by_mask[source_n], node_by_mask[target_n]
    source_lines = 1 << (m_tiles(source_n) - 1)
    p_fibres: dict[tuple[str, int, int], list[int]] = defaultdict(list)
    p_images: dict[tuple[str, int, int], set[tuple[str, int, int]]] = defaultdict(set)
    q_fibres: dict[int, list[int]] = defaultdict(list)
    q_images: dict[int, set[tuple[str, int, int]]] = defaultdict(set)
    target_by_line = {}
    values = {}
    for line in range(source_lines):
        source_p = edge_key(line, source_n, source_map)
        q_rep = min(line, line_index(reflect(line, source_n), source_n))
        target_mask = phi(line, source_n, target_n, rho)
        target_p = edge_key(target_mask, target_n, target_map)
        p_fibres[source_p].append(line)
        p_images[source_p].add(target_p)
        q_fibres[q_rep].append(line)
        q_images[q_rep].add(target_p)
        target_by_line[line] = target_p
        values[line] = {
            "colour": source_p[0],
            "source_P": source_p,
            "source_Q": q_rep,
            "target_P": target_p,
            "source_P_target_P": (source_p, target_p),
            "literal_line": line,
        }

    p_failures = {key: images for key, images in p_images.items() if len(images) > 1}
    q_failures = {key: images for key, images in q_images.items() if len(images) > 1}
    multi_q_per_p = {}
    for key, lines in p_fibres.items():
        qs = {min(line, line_index(reflect(line, source_n), source_n)) for line in lines}
        if len(qs) > 1:
            multi_q_per_p[key] = sorted(qs)

    q_target_fibres: dict[tuple[str, int, int], list[int]] = defaultdict(list)
    for q_rep, images in q_images.items():
        assert len(images) == 1
        q_target_fibres[next(iter(images))].append(q_rep)
    q_target_collisions = {
        key: reps for key, reps in q_target_fibres.items() if len(reps) > 1
    }

    first_key = ("B", 15, 16)
    witness = {
        "source_P": list(first_key),
        "lines": p_fibres[first_key],
        "target_P": [list(target_by_line[line]) for line in p_fibres[first_key]],
        "source_one_tiles": [
            [list(tile) for bit, tile in enumerate(tile_schema(source_n)[0]) if line & (1 << bit)]
            for line in p_fibres[first_key]
        ],
        "target_masks": [phi(line, source_n, target_n, rho) for line in p_fibres[first_key]],
        "target_one_tiles": [
            [
                list(tile)
                for bit, tile in enumerate(tile_schema(target_n)[0])
                if phi(line, source_n, target_n, rho) & (1 << bit)
            ]
            for line in p_fibres[first_key]
        ],
    }

    carrier_order = (
        "colour", "source_P", "source_Q", "target_P", "source_P_target_P", "literal_line"
    )
    stats = {
        carrier: compact_partition({line: row[carrier] for line, row in values.items()})
        for carrier in carrier_order
    }
    retention = carrier_tournament(stats, "retention")
    economy = carrier_tournament(stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carrier_order)) for j in range(i + 1, len(carrier_order))
    )

    return {
        "source_lines": source_lines,
        "source_P_cells": len(p_fibres),
        "source_Q_cells": len(q_fibres),
        "target_P_cells_from_Q": len(q_target_fibres),
        "P_fibre_size_histogram": dict(sorted(Counter(map(len, p_fibres.values())).items())),
        "P_target_support_histogram": dict(sorted(Counter(map(len, p_images.values())).items())),
        "P_descent_failure_cells": len(p_failures),
        "P_descent_failure_blue": sum(key[0] == "B" for key in p_failures),
        "P_descent_failure_black": sum(key[0] == "K" for key in p_failures),
        "P_cells_with_multiple_Q_orbits": len(multi_q_per_p),
        "multi_Q_P_successes": [list(key) for key in multi_q_per_p if key not in p_failures],
        "Q_descent_failure_cells": len(q_failures),
        "Q_to_target_P_collision_cells": len(q_target_collisions),
        "Q_to_target_P_collision_excess": sum(len(reps) - 1 for reps in q_target_collisions.values()),
        "Q_to_target_P_collisions": [
            {"target_P": list(key), "Q_representatives": reps}
            for key, reps in sorted(q_target_collisions.items())
        ],
        "first_failure": witness,
        "tournament_analysis": {
            "vertices": list(carrier_order),
            "pairwise_observable": "number of unordered source-line pairs separated by a carrier",
            "switches": ["retention", "retention_per_log2_cells"],
            "tie_hamiltonian_path": list(carrier_order),
            "carrier_stats": stats,
            "retention": retention,
            "economy": economy,
            "edge_flips_between_gauges": flips,
        },
    }


def natural_core_lifts() -> list[tuple[int, ...]]:
    """All core-local surjective reflection-intertwining lifts with fixed legs."""
    fixed = {0: 0, 1: 1, 2: 2, 3: 2, 4: 3, 5: 4, 9: 7, 12: 7, 14: 9}
    pair_options = ((6, 8), (8, 6), (0, 0), (5, 5))
    lifts = []
    for first in pair_options:      # target pair (7,10)
        for second in pair_options: # target pair (8,13)
            for fixed6 in (0, 5):
                for fixed11 in (0, 5):
                    mapping = dict(fixed)
                    mapping[7], mapping[10] = first
                    mapping[8], mapping[13] = second
                    mapping[6], mapping[11] = fixed6, fixed11
                    rho = tuple(mapping[i] for i in range(m_tiles(7)))
                    if len(set(rho)) == m_tiles(6):
                        lifts.append(rho)
    return lifts


def robust_core_audit(node_by_mask: dict[int, list[int]]) -> dict:
    lifts = natural_core_lifts()
    failure_hist = Counter()
    witness_failures = 0
    for rho in lifts:
        coordinate = coordinate_audit(6, 7, rho)
        assert coordinate["rho_surjective"]
        assert coordinate["rho_sigma_failures"] == coordinate["reflection_failures"] == 0
        source_map, target_map = node_by_mask[6], node_by_mask[7]
        images: dict[tuple[str, int, int], set[tuple[str, int, int]]] = defaultdict(set)
        for line in range(1 << (m_tiles(6) - 1)):
            images[edge_key(line, 6, source_map)].add(
                edge_key(phi(line, 6, 7, rho), 7, target_map)
            )
        failures = sum(len(values) > 1 for values in images.values())
        failure_hist[failures] += 1
        witness_failures += len(images[("B", 15, 16)]) <= 1
    return {
        "candidate_lifts_before_surjectivity": 64,
        "surjective_equivariant_lifts": len(lifts),
        "P_failure_count_histogram": dict(sorted(failure_hist.items())),
        "B15_16_failure_missing_lifts": witness_failures,
    }


def compose_rho(outer: tuple[int, ...], inner: tuple[int, ...]) -> tuple[int, ...]:
    """outer: target->middle, inner: middle->source."""
    return tuple(inner[middle] for middle in outer)


def composed_audit(node_by_mask: dict[int, list[int]]) -> dict:
    rho = compose_rho(RHO_6_7, rho_5_6())
    coordinate = coordinate_audit(5, 7, rho)
    source_map, target_map = node_by_mask[5], node_by_mask[7]
    node_support: dict[int, set[int]] = defaultdict(set)
    edge_images: dict[tuple[str, int, int], set[tuple[str, int, int]]] = defaultdict(set)
    for mask in range(1 << m_tiles(5)):
        node_support[source_map[mask]].add(target_map[phi(mask, 5, 7, rho)])
    for line in range(1 << (m_tiles(5) - 1)):
        edge_images[edge_key(line, 5, source_map)].add(
            edge_key(phi(line, 5, 7, rho), 7, target_map)
        )

    source_index = tile_index(5)
    variable_bits = tuple(source_index[t] for t in ((5, 2), (5, 3), (3, 1), (4, 1), (5, 1)))
    core_bit = source_index[(4, 2)]
    core_rows = []
    for core in (0, 1):
        targets = []
        for assignment in range(32):
            mask = core << core_bit
            for variable, bit in enumerate(variable_bits):
                mask |= ((assignment >> variable) & 1) << bit
            targets.append(target_map[phi(mask, 5, 7, rho)])
        nodes = sorted(set(targets))
        coefficients = {
            node: mobius_coefficients([int(target == node) for target in targets])
            for node in nodes
        }
        cells_by_degree = []
        collisions = []
        for degree in range(6):
            cells: dict[tuple[int, ...], list[int]] = defaultdict(list)
            for node in nodes:
                cells[truncated_signature(coefficients[node], degree)].append(node)
            cells_by_degree.append(len(cells))
            collisions.append([cell for cell in cells.values() if len(cell) > 1])
        core_rows.append(
            {
                "core": core,
                "target_nodes": len(nodes),
                "cells_by_max_degree": cells_by_degree,
                "degree_three_collisions": collisions[3],
            }
        )
    return {
        "rho_5_7": list(rho),
        "coordinate_audit": coordinate,
        "source_node_count": len(node_support),
        "target_node_union": len(set().union(*node_support.values())),
        "node_support_size_histogram": dict(sorted(Counter(map(len, node_support.values())).items())),
        "source_P_cells": len(edge_images),
        "target_P_cells": len(set().union(*edge_images.values())),
        "P_descent_failures": sum(len(images) != 1 for images in edge_images.values()),
        "mobius": core_rows,
        "coefficient_composition_law": "(rho1)_* (rho2)_* mu = (rho1 o rho2)_* mu",
        "degree_saturation_warning": "target monomials with high |B| may have low |rho(B)|",
    }


def audit() -> dict:
    _, node_by_mask = load_atlases(SMALL_ATLAS, N7_ATLAS)
    schedule = {
        "p": 3,
        "q": 4,
        "phase": 0,
        "heights": [centered_height(4, 3, 0, i) for i in range(4)],
        "increments": list(centered_increments(4, 3, 0)),
        "partial_quotient": 1,
        "next_phase": 1,
        "has_midpoint_tie": False,
    }
    coordinate = coordinate_audit(6, 7, RHO_6_7)
    edges = line_edge_audit(6, 7, RHO_6_7, node_by_mask)
    robustness = robust_core_audit(node_by_mask)
    composed = composed_audit(node_by_mask)

    assert schedule["increments"] == [1, 2, 1]
    assert coordinate["image_size"] == 1024 and coordinate["rho_surjective"]
    assert coordinate["rho_sigma_failures"] == coordinate["complement_failures"] == 0
    assert coordinate["reflection_failures"] == coordinate["colour_failures"] == 0
    assert edges["source_P_cells"] == 187 and edges["source_Q_cells"] == 272
    assert edges["target_P_cells_from_Q"] == 268
    assert edges["P_fibre_size_histogram"] == {1: 24, 2: 115, 4: 31, 6: 8, 8: 6, 12: 2, 14: 1}
    assert edges["P_target_support_histogram"] == {1: 136, 2: 34, 3: 8, 4: 6, 6: 2, 7: 1}
    assert edges["P_descent_failure_cells"] == 51
    assert (edges["P_descent_failure_blue"], edges["P_descent_failure_black"]) == (4, 47)
    assert edges["P_cells_with_multiple_Q_orbits"] == 52
    assert edges["multi_Q_P_successes"] == [["K", 12, 29]]
    assert edges["Q_descent_failure_cells"] == 0
    assert edges["Q_to_target_P_collision_cells"] == edges["Q_to_target_P_collision_excess"] == 4
    assert edges["first_failure"]["lines"] == [132, 370]
    assert edges["first_failure"]["target_P"] == [["B", 143, 203], ["B", 27, 65]]
    assert robustness["surjective_equivariant_lifts"] == 40
    assert robustness["P_failure_count_histogram"] == {51: 8, 52: 32}
    assert robustness["B15_16_failure_missing_lifts"] == 0
    assert composed["rho_5_7"] == [0, 1, 2, 2, 2, 3, 0, 4, 4, 5, 4, 0, 5, 4, 5]
    assert composed["P_descent_failures"] == 0
    assert composed["source_P_cells"] == composed["target_P_cells"] == 20
    assert composed["mobius"][0]["cells_by_max_degree"] == [2, 5, 10, 16, 18, 18]
    assert composed["mobius"][1]["cells_by_max_degree"] == [2, 5, 11, 17, 18, 18]
    assert composed["mobius"][0]["degree_three_collisions"] == [[73, 143, 187]]
    assert composed["mobius"][1]["degree_three_collisions"] == [[35, 145]]

    return {
        "schema_version": 1,
        "theorem": "THM-813",
        "schedule": schedule,
        "coordinate_embedding": coordinate,
        "edge_descent": edges,
        "core_lift_robustness": robustness,
        "composed_X5_X7": composed,
        "general_Q_lemma": "any kappa/sigma-equivariant tiling map induces E_source/<sigma> -> P_target",
        "P_descent_criterion": "descent through P_source iff the induced Q-source function is constant inside every P-source fibre",
    }


def render(result: dict) -> str:
    s, c, e = result["schedule"], result["coordinate_embedding"], result["edge_descent"]
    lines = [
        "THM-813 CENTERED-CHRISTOFFEL EDGE-DESCENT BOUNDARY",
        "=" * 78,
        f"schedule (q,p,s)=({s['q']},{s['p']},{s['phase']}): F={tuple(s['heights'])}, "
        f"d={tuple(s['increments'])}, next phase={s['next_phase']}, tie={s['has_midpoint_tie']}",
        f"rho={tuple(c['rho'])}; image={c['image_size']}/{c['source_size']}; "
        f"intertwining/complement/reflection/colour failures={c['rho_sigma_failures']}/"
        f"{c['complement_failures']}/{c['reflection_failures']}/{c['colour_failures']}",
        "",
        "P=COLOUR+UNORDERED NODE PAIR; Q=COMPLEMENT LINE / STAIRCASE REFLECTION",
        f"  lines={e['source_lines']}; P cells={e['source_P_cells']}; Q cells={e['source_Q_cells']}",
        f"  P fibre histogram={e['P_fibre_size_histogram']}",
        f"  P target-support histogram={e['P_target_support_histogram']}",
        f"  P failures={e['P_descent_failure_cells']} "
        f"(blue={e['P_descent_failure_blue']}, black={e['P_descent_failure_black']})",
        f"  multi-Q P cells={e['P_cells_with_multiple_Q_orbits']}; "
        f"sole accidental success={e['multi_Q_P_successes']}",
        f"  Q descent failures={e['Q_descent_failure_cells']}; Q cells -> target P cells="
        f"{e['source_Q_cells']}->{e['target_P_cells_from_Q']} "
        f"({e['Q_to_target_P_collision_cells']} double collisions)",
        "",
        "FIRST EXACT FAILURE",
        f"  P={tuple(e['first_failure']['source_P'])}; lines={tuple(e['first_failure']['lines'])}",
        f"  target P={tuple(map(tuple, e['first_failure']['target_P']))}",
        f"  target masks={tuple(e['first_failure']['target_masks'])}",
    ]
    robust = result["core_lift_robustness"]
    lines.extend(
        [
            "",
            "CORE-LIFT ROBUSTNESS",
            f"  candidates/surjective equivariant={robust['candidate_lifts_before_surjectivity']}/"
            f"{robust['surjective_equivariant_lifts']}",
            f"  P-failure histogram={robust['P_failure_count_histogram']}; "
            f"B(15,16) missing failures={robust['B15_16_failure_missing_lifts']}",
        ]
    )
    comp = result["composed_X5_X7"]
    lines.extend(
        [
            "",
            "COMPOSED X5 -> X7 ACTION",
            f"  rho={tuple(comp['rho_5_7'])}",
            f"  source/target P cells/failures={comp['source_P_cells']}/"
            f"{comp['target_P_cells']}/{comp['P_descent_failures']}",
            f"  node support={comp['source_node_count']} source -> {comp['target_node_union']} target; "
            f"histogram={comp['node_support_size_histogram']}",
        ]
    )
    for row in comp["mobius"]:
        lines.append(
            f"  core={row['core']}: degree cells={tuple(row['cells_by_max_degree'])}; "
            f"degree-3 collisions={row['degree_three_collisions']}"
        )
    lines.append(f"  coefficient law: {comp['coefficient_composition_law']}")
    ta = e["tournament_analysis"]
    lines.extend(
        [
            "",
            "TOURNAMENT ANALYSIS (line-information carriers as vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; edge flips={ta['edge_flips_between_gauges']}",
            f"  retention scores={ta['retention']['score_hist']} C3="
            f"{ta['retention']['directed_3cycles']} SCC={ta['retention']['scc_sizes']} "
            f"Hpaths={ta['retention']['hamiltonian_paths']}",
            "",
            f"GENERAL LEMMA: {result['general_Q_lemma']}",
            f"DESCENT TEST: {result['P_descent_criterion']}",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()
    result = audit()
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
