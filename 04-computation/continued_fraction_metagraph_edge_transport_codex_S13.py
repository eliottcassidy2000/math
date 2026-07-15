#!/usr/bin/env python3
"""Exact centred-Christoffel transport from X_5 into X_6 for THM-812.

The (q,p,s)=(3,2,1) centred increment word is (1,2).  It replicates
the low leg in blocks (1,2), the reflected high leg in blocks (2,1), and
extends the literal core by (c,a,c).  The resulting coordinate-copy map is
audited on tilings, complement lines, converse-merged nodes, coloured
projected edge cells, and anchored Boolean Mobius coefficients.

Tournament Analysis deliberately uses line-information carriers as vertices:
the pairwise observable is separation of literal source complement lines,
the switches are retention and retention per log-cell, and the declared
carrier order is the tie Hamiltonian path.  Runners are payload, not vertices.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

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


OUTPUT = Path("05-knowledge/results/continued_fraction_metagraph_edge_transport_codex_S13.out")
JSON_OUTPUT = Path("05-knowledge/results/continued_fraction_metagraph_edge_transport_codex_S13.json")


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def centered_height(q: int, p: int, phase: int, i: int) -> int:
    return ceil_fraction(
        Fraction(q * i, p) + Fraction(q, 2 * p) - Fraction(phase, 2)
    )


def centered_increments(q: int, p: int, phase: int) -> tuple[int, ...]:
    return tuple(
        centered_height(q, p, phase, i + 1)
        - centered_height(q, p, phase, i)
        for i in range(p)
    )


def coordinate_copy(mask: int, source_n: int, target_n: int, rho: tuple[int, ...]) -> int:
    assert len(rho) == m_tiles(target_n)
    return sum(((mask >> source_bit) & 1) << target_bit for target_bit, source_bit in enumerate(rho))


def rho_5_6() -> tuple[int, ...]:
    """Target coordinate -> source coordinate for the (1,2) leg replication."""
    source = tile_index(5)
    mapping = {
        (6, 1): (5, 1),                         # apex a
        (3, 1): (3, 1), (4, 1): (3, 1),        # h0,h0
        (5, 1): (4, 1),                         # h1
        (6, 2): (5, 2),                         # l0
        (6, 3): (5, 3), (6, 4): (5, 3),        # l1,l1
        (4, 2): (4, 2), (5, 2): (5, 1), (5, 3): (4, 2),  # c,a,c
    }
    return tuple(source[mapping[tile]] for tile in tile_schema(6)[0])


def rho_4_6_tie() -> tuple[int, ...]:
    """Odd/odd tie sanity map: copy both one-bit legs thrice; core = apex."""
    source = tile_index(4)
    mapping = {
        (6, 1): (4, 1),
        (3, 1): (3, 1), (4, 1): (3, 1), (5, 1): (3, 1),
        (6, 2): (4, 2), (6, 3): (4, 2), (6, 4): (4, 2),
        (4, 2): (4, 1), (5, 2): (4, 1), (5, 3): (4, 1),
    }
    return tuple(source[mapping[tile]] for tile in tile_schema(6)[0])


def coordinate_audit(source_n: int, target_n: int, rho: tuple[int, ...]) -> dict:
    source_bits = m_tiles(source_n)
    target_bits = m_tiles(target_n)
    source_sigma = tile_schema(source_n)[1]
    target_sigma = tile_schema(target_n)[1]
    images = [coordinate_copy(mask, source_n, target_n, rho) for mask in range(1 << source_bits)]
    rho_sigma_failures = sum(
        rho[target_sigma[j]] != source_sigma[rho[j]] for j in range(target_bits)
    )
    complement_failures = sum(
        coordinate_copy(complement(mask, source_n), source_n, target_n, rho)
        != complement(images[mask], target_n)
        for mask in range(1 << source_bits)
    )
    reflection_failures = sum(
        coordinate_copy(reflect(mask, source_n), source_n, target_n, rho)
        != reflect(images[mask], target_n)
        for mask in range(1 << source_bits)
    )
    colour_failures = sum(
        is_grid_symmetric(mask, source_sigma)
        != is_grid_symmetric(images[mask], target_sigma)
        for mask in range(1 << source_bits)
    )
    return {
        "source_n": source_n,
        "target_n": target_n,
        "rho": list(rho),
        "rho_surjective": len(set(rho)) == source_bits,
        "image_size": len(set(images)),
        "source_size": 1 << source_bits,
        "rho_sigma_failures": rho_sigma_failures,
        "complement_failures": complement_failures,
        "reflection_failures": reflection_failures,
        "colour_failures": colour_failures,
    }


def edge_key(mask: int, n: int, node_map: list[int]) -> tuple[str, int, int]:
    _, sigma = tile_schema(n)
    colour = "B" if is_grid_symmetric(mask, sigma) else "K"
    nodes = sorted((node_map[mask], node_map[complement(mask, n)]))
    return colour, nodes[0], nodes[1]


def node_and_edge_transport(node_by_mask: dict[int, list[int]], rho: tuple[int, ...]) -> tuple[dict, dict[int, dict[str, object]]]:
    source_map = node_by_mask[5]
    target_map = node_by_mask[6]
    phi = lambda mask: coordinate_copy(mask, 5, 6, rho)

    node_support: dict[int, set[int]] = defaultdict(set)
    for mask in range(1 << m_tiles(5)):
        node_support[source_map[mask]].add(target_map[phi(mask)])
    support_records = [
        {"source_node": source, "target_nodes": sorted(targets), "support_size": len(targets)}
        for source, targets in sorted(node_support.items())
    ]

    source_fibres: dict[tuple[str, int, int], list[int]] = defaultdict(list)
    target_images: dict[tuple[str, int, int], set[tuple[str, int, int]]] = defaultdict(set)
    line_values: dict[int, dict[str, object]] = {}
    source_line_count = 1 << (m_tiles(5) - 1)
    for line in range(source_line_count):
        source_key = edge_key(line, 5, source_map)
        target_mask = phi(line)
        target_key = edge_key(target_mask, 6, target_map)
        source_fibres[source_key].append(line)
        target_images[source_key].add(target_key)
        line_values[line] = {
            "source_colour": source_key[0],
            "source_node_pair": source_key[1:],
            "source_edge_cell": source_key,
            "target_node_pair": target_key[1:],
            "target_edge_cell": target_key,
            "literal_line": line,
        }

    descent_failures = sum(len(images) != 1 for images in target_images.values())
    colour_descent_failures = sum(
        source[0] != next(iter(targets))[0]
        for source, targets in target_images.items()
    )
    reflection_fibre_failures = 0
    for source, lines in source_fibres.items():
        reflected = {line_index(reflect(line, 5), 5) for line in lines}
        reflection_fibre_failures += reflected != set(lines)

    mapping = []
    for source in sorted(target_images):
        targets = sorted(target_images[source])
        mapping.append(
            {
                "source": list(source),
                "target": list(targets[0]) if len(targets) == 1 else [list(x) for x in targets],
                "source_lines": source_fibres[source],
            }
        )

    carrier_order = (
        "source_colour",
        "source_node_pair",
        "source_edge_cell",
        "target_node_pair",
        "target_edge_cell",
        "literal_line",
    )
    values = {
        carrier: {line: record[carrier] for line, record in line_values.items()}
        for carrier in carrier_order
    }
    stats = {name: compact_partition(partition) for name, partition in values.items()}
    retention = carrier_tournament(stats, "retention")
    economy = carrier_tournament(stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carrier_order)) for j in range(i + 1, len(carrier_order))
    )

    audit = {
        "source_node_count": len(node_support),
        "target_node_union": len(set().union(*node_support.values())),
        "node_support_size_histogram": dict(sorted(Counter(map(len, node_support.values())).items())),
        "node_support": support_records,
        "bare_node_transport_is_function": all(len(targets) == 1 for targets in node_support.values()),
        "source_edge_cells": len(source_fibres),
        "target_edge_cells_in_image": len(set().union(*target_images.values())),
        "edge_cell_descent_failures": descent_failures,
        "colour_descent_failures": colour_descent_failures,
        "reflection_fibre_failures": reflection_fibre_failures,
        "source_fibre_size_histogram": dict(sorted(Counter(map(len, source_fibres.values())).items())),
        "blue_singleton_cells": sum(key[0] == "B" and len(lines) == 1 for key, lines in source_fibres.items()),
        "black_doubleton_cells": sum(key[0] == "K" and len(lines) == 2 for key, lines in source_fibres.items()),
        "mapping": mapping,
        "tournament_analysis": {
            "vertices": list(carrier_order),
            "pairwise_observable": "number of unordered source complement-line pairs separated by a carrier",
            "switches": ["retention", "retention_per_log2_cells"],
            "tie_hamiltonian_path": list(carrier_order),
            "carrier_stats": stats,
            "retention": retention,
            "economy": economy,
            "edge_flips_between_gauges": flips,
        },
    }
    return audit, line_values


def mobius_coefficients(truth: list[int]) -> list[int]:
    coefficients = truth.copy()
    dimension = (len(truth) - 1).bit_length()
    for bit in range(dimension):
        for mask in range(len(truth)):
            if mask & (1 << bit):
                coefficients[mask] -= coefficients[mask ^ (1 << bit)]
    return coefficients


def truncated_signature(coefficients: list[int], degree: int) -> tuple[int, ...]:
    return tuple(value for mask, value in enumerate(coefficients) if mask.bit_count() <= degree)


def mobius_transport(node_by_mask: dict[int, list[int]], rho: tuple[int, ...]) -> dict:
    target_map = node_by_mask[6]
    source_index = tile_index(5)
    variable_tiles = ((5, 2), (5, 3), (3, 1), (4, 1), (5, 1))
    variable_names = ("l0", "l1", "h0", "h1", "a")
    variable_bits = tuple(source_index[tile] for tile in variable_tiles)
    core_bit = source_index[(4, 2)]
    core_results = []
    coefficient_tables: dict[int, dict[int, list[int]]] = {}

    for core in (0, 1):
        target_at_assignment = []
        for assignment in range(1 << len(variable_bits)):
            source_mask = core << core_bit
            for variable, source_bit in enumerate(variable_bits):
                source_mask |= ((assignment >> variable) & 1) << source_bit
            target_mask = coordinate_copy(source_mask, 5, 6, rho)
            target_at_assignment.append(target_map[target_mask])
        target_nodes = sorted(set(target_at_assignment))
        coefficients = {}
        truth_sets = {}
        for node in target_nodes:
            truth = [int(value == node) for value in target_at_assignment]
            truth_sets[node] = [mask for mask, value in enumerate(truth) if value]
            coefficients[node] = mobius_coefficients(truth)
        coefficient_tables[core] = coefficients

        cells_by_degree = []
        collisions_by_degree = []
        for degree in range(len(variable_bits) + 1):
            cells: dict[tuple[int, ...], list[int]] = defaultdict(list)
            for node in target_nodes:
                cells[truncated_signature(coefficients[node], degree)].append(node)
            cells_by_degree.append(len(cells))
            collisions_by_degree.append([nodes for nodes in cells.values() if len(nodes) > 1])
        core_results.append(
            {
                "core": core,
                "target_nodes": target_nodes,
                "target_node_count": len(target_nodes),
                "cells_by_max_degree": cells_by_degree,
                "collisions_by_max_degree": collisions_by_degree,
                "truth_sets": {str(node): truth_sets[node] for node in target_nodes},
            }
        )

    collision_pair = coefficient_tables[1]
    first_separators = []
    for mask in range(1 << len(variable_bits)):
        left = collision_pair[12][mask]
        right = collision_pair[31][mask]
        if mask.bit_count() == 4 and left != right:
            first_separators.append(
                {
                    "subset_mask": mask,
                    "variables": [variable_names[i] for i in range(len(variable_names)) if mask & (1 << i)],
                    "bidegree": [
                        sum(mask & (1 << i) != 0 for i in (0, 1)),
                        sum(mask & (1 << i) != 0 for i in (2, 3)),
                        int(bool(mask & (1 << 4))),
                    ],
                    "node12": left,
                    "node31": right,
                }
            )

    return {
        "source_variable_order": list(variable_names),
        "source_variable_tiles": [list(tile) for tile in variable_tiles],
        "core_tile": [4, 2],
        "core_results": core_results,
        "unique_degree_three_collision": [12, 31],
        "collision_truth_sets": {
            "12": core_results[1]["truth_sets"]["12"],
            "31": core_results[1]["truth_sets"]["31"],
        },
        "first_degree_four_separators": first_separators,
        "general_pullback_law": "mu_(f o Phi)(A)=sum_{B subset S_target: rho(B)=A} mu_f(B)",
    }


def tie_sanity(node_by_mask: dict[int, list[int]], rho: tuple[int, ...]) -> dict:
    coordinate = coordinate_audit(4, 6, rho)
    source_map = node_by_mask[4]
    target_map = node_by_mask[6]
    support: dict[int, set[int]] = defaultdict(set)
    edge_images: dict[tuple[str, int, int], set[tuple[str, int, int]]] = defaultdict(set)
    for mask in range(1 << m_tiles(4)):
        target = coordinate_copy(mask, 4, 6, rho)
        support[source_map[mask]].add(target_map[target])
    for line in range(1 << (m_tiles(4) - 1)):
        target = coordinate_copy(line, 4, 6, rho)
        edge_images[edge_key(line, 4, source_map)].add(edge_key(target, 6, target_map))
    return {
        "schedule": {"p": 1, "q": 3, "phase": 1, "increments": [3], "has_midpoint_tie": True},
        "coordinate_audit": coordinate,
        "source_node_count": len(support),
        "node_support_size_histogram": dict(sorted(Counter(map(len, support.values())).items())),
        "source_lines": 1 << (m_tiles(4) - 1),
        "source_edge_cells": len(edge_images),
        "edge_cell_descent_failures": sum(len(images) != 1 for images in edge_images.values()),
    }


def audit() -> dict:
    _, node_by_mask = load_atlases(SMALL_ATLAS, N7_ATLAS)
    increments = centered_increments(3, 2, 1)
    heights = tuple(centered_height(3, 2, 1, i) for i in range(3))
    next_phase = (1 - 1) % 2
    rho = rho_5_6()
    coordinate = coordinate_audit(5, 6, rho)
    edge_transport, _ = node_and_edge_transport(node_by_mask, rho)
    mobius = mobius_transport(node_by_mask, rho)
    tie = tie_sanity(node_by_mask, rho_4_6_tie())

    assert increments == (1, 2) and heights == (1, 2, 4) and next_phase == 0
    assert coordinate["rho_surjective"] and coordinate["image_size"] == 64
    assert coordinate["rho_sigma_failures"] == 0
    assert coordinate["complement_failures"] == coordinate["reflection_failures"] == 0
    assert coordinate["colour_failures"] == 0
    assert edge_transport["source_node_count"] == 10
    assert edge_transport["target_node_union"] == 23
    assert edge_transport["node_support_size_histogram"] == {1: 3, 2: 1, 3: 2, 4: 1, 5: 1, 6: 2}
    assert not edge_transport["bare_node_transport_is_function"]
    assert edge_transport["source_edge_cells"] == edge_transport["target_edge_cells_in_image"] == 20
    assert edge_transport["edge_cell_descent_failures"] == 0
    assert edge_transport["colour_descent_failures"] == 0
    assert edge_transport["reflection_fibre_failures"] == 0
    assert edge_transport["blue_singleton_cells"] == 8
    assert edge_transport["black_doubleton_cells"] == 12
    assert mobius["core_results"][0]["cells_by_max_degree"] == [2, 5, 10, 13, 13, 13]
    assert mobius["core_results"][1]["cells_by_max_degree"] == [2, 5, 10, 14, 15, 15]
    assert mobius["core_results"][1]["collisions_by_max_degree"][3] == [[12, 31]]
    assert mobius["collision_truth_sets"] == {"12": [27, 29], "31": [15]}
    assert len(mobius["first_degree_four_separators"]) == 3
    assert tie["coordinate_audit"]["image_size"] == 8
    assert tie["coordinate_audit"]["reflection_failures"] == 0
    assert tie["edge_cell_descent_failures"] == 0 and tie["source_edge_cells"] == 3

    return {
        "schema_version": 1,
        "theorem": "THM-812",
        "schedule": {
            "p": 2,
            "q": 3,
            "phase": 1,
            "heights": list(heights),
            "increments": list(increments),
            "partial_quotient": 1,
            "next_phase": next_phase,
            "has_midpoint_tie": False,
        },
        "coordinate_embedding": coordinate,
        "node_and_edge_transport": edge_transport,
        "mobius_transport": mobius,
        "odd_odd_tie_sanity": tie,
    }


def render(result: dict) -> str:
    schedule = result["schedule"]
    coordinate = result["coordinate_embedding"]
    edge = result["node_and_edge_transport"]
    mobius = result["mobius_transport"]
    lines = [
        "THM-812 CENTERED-CHRISTOFFEL METAGRAPH EDGE TRANSPORT",
        "=" * 78,
        f"schedule (q,p,s)=({schedule['q']},{schedule['p']},{schedule['phase']}): "
        f"F={tuple(schedule['heights'])}, d={tuple(schedule['increments'])}, "
        f"next phase={schedule['next_phase']}, tie={schedule['has_midpoint_tie']}",
        "low replication=(l0,l1,l1); high=(h0,h0,h1); core=(c,a,c)",
        f"rho={tuple(coordinate['rho'])}; surjective={coordinate['rho_surjective']}; "
        f"image={coordinate['image_size']}/{coordinate['source_size']}",
        f"rho-sigma/complement/reflection/colour failures="
        f"{coordinate['rho_sigma_failures']}/{coordinate['complement_failures']}/"
        f"{coordinate['reflection_failures']}/{coordinate['colour_failures']}",
        "",
        "NODE AND COLOURED-EDGE TRANSPORT",
        f"  bare source nodes={edge['source_node_count']}; target union={edge['target_node_union']}; "
        f"support histogram={edge['node_support_size_histogram']}",
        f"  bare-node map is a function: {edge['bare_node_transport_is_function']}",
        f"  projected coloured edge cells source/target={edge['source_edge_cells']}/"
        f"{edge['target_edge_cells_in_image']}; descent failures={edge['edge_cell_descent_failures']}",
        f"  blue singleton cells={edge['blue_singleton_cells']}; "
        f"black reflection-doubletons={edge['black_doubleton_cells']}; "
        f"reflection failures={edge['reflection_fibre_failures']}",
    ]
    for record in edge["mapping"]:
        source = record["source"]
        target = record["target"]
        lines.append(
            f"    {source[0]}({source[1]},{source[2]}) -> "
            f"{target[0]}({target[1]},{target[2]}) lines={tuple(record['source_lines'])}"
        )
    lines.extend(["", "ANCHORED BOOLEAN MOBIUS PULLBACK"])
    for core in mobius["core_results"]:
        lines.append(
            f"  core={core['core']}: target nodes={core['target_node_count']}; "
            f"cells through degrees 0..5={tuple(core['cells_by_max_degree'])}"
        )
    lines.append(
        f"  unique degree<=3 collision={tuple(mobius['unique_degree_three_collision'])}; "
        f"truth sets={mobius['collision_truth_sets']}"
    )
    for separator in mobius["first_degree_four_separators"]:
        lines.append(
            f"    degree-4 {tuple(separator['variables'])} bidegree={tuple(separator['bidegree'])}: "
            f"mu_12={separator['node12']} mu_31={separator['node31']}"
        )
    lines.append(f"  law: {mobius['general_pullback_law']}")
    tie = result["odd_odd_tie_sanity"]
    lines.extend(
        [
            "",
            "ODD/ODD TIE SANITY",
            f"  schedule d={tuple(tie['schedule']['increments'])}; image="
            f"{tie['coordinate_audit']['image_size']}/8; source nodes={tie['source_node_count']}; "
            f"support histogram={tie['node_support_size_histogram']}",
            f"  lines/edge cells/descent failures={tie['source_lines']}/"
            f"{tie['source_edge_cells']}/{tie['edge_cell_descent_failures']}",
        ]
    )
    ta = edge["tournament_analysis"]
    lines.extend(
        [
            "",
            "TOURNAMENT ANALYSIS (line-information carriers as vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; edge flips={ta['edge_flips_between_gauges']}",
            f"  retention scores={ta['retention']['score_hist']} "
            f"C3={ta['retention']['directed_3cycles']} SCC={ta['retention']['scc_sizes']} "
            f"Hpaths={ta['retention']['hamiltonian_paths']}",
            "",
            "The CF replication is functorial on literal tilings, complement lines, and",
            "these 20 projected coloured edge fibres, but not on the 10 bare source",
            "nodes.  Degree <=3 also fails on one fixed-core node-indicator pair;",
            "this is a truncation obstruction, not a failure of the full Omega+B2 codec.",
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
