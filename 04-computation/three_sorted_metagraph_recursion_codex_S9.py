#!/usr/bin/env python3
"""Exact recursion among tilings, complement lines, and merged nodes.

The three sorts are kept separate:

  T_n  labelled fixed-Hamiltonian-path tilings,
  E_n  complement-line orbits {t,kappa(t)},
  N_n  converse-merged tournament isomorphism nodes.

At tiling level the script proves and audits the pullback

  T_n = (T_(n-1) x_[T_(n-2)] T_(n-1)) x {0,1}_apex.

Choosing the unique apex-zero endpoint of a line removes the last factor and
gives E_n the compatible-pair model.  Projection T_n -> N_n produces exact
weighted branching and blue/black endpoint kernels; it does not produce a
Markov recursion on bare nodes, and the script records the first failures of
strong lumpability.

Tournament Analysis uses information carriers as vertices.  Its pairwise
observable is the number of node pairs separated by each carrier; retention
and retention-per-cell are the two gauges, with the declared carrier list as
the tie Hamiltonian path.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

from tournament_tiling_metagraph_address_codex_S4 import (
    carrier_tournament,
    is_grid_symmetric,
    partition_stats,
    tile_schema,
)


SMALL_ATLAS = Path("05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json")
N7_ATLAS = Path("05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json")


def m_tiles(n: int) -> int:
    return math.comb(n - 1, 2)


def reflection_orbits(n: int) -> int:
    return (m_tiles(n) + (n - 1) // 2) // 2


@lru_cache(maxsize=None)
def tile_index(n: int) -> dict[tuple[int, int], int]:
    tiles, _ = tile_schema(n)
    return {tile: k for k, tile in enumerate(tiles)}


@lru_cache(maxsize=None)
def face_bit_map(n: int, side: str) -> tuple[tuple[int, int], ...]:
    """Map upper tile positions to endpoint-deleted face positions."""
    tiles, _ = tile_schema(n)
    lower_index = tile_index(n - 1)
    result = []
    for upper_bit, (x, y) in enumerate(tiles):
        if side == "low":
            if y == 1:
                continue
            lower_tile = (x - 1, y - 1)
        elif side == "high":
            if x == n:
                continue
            lower_tile = (x, y)
        else:
            raise ValueError(side)
        result.append((upper_bit, lower_index[lower_tile]))
    return tuple(result)


def face_mask(mask: int, n: int, side: str) -> int:
    return sum(((mask >> upper) & 1) << lower for upper, lower in face_bit_map(n, side))


def apex_bit(mask: int, n: int) -> int:
    return (mask >> tile_index(n)[(n, 1)]) & 1


def reconstruct_tiling(low: int, high: int, apex: int, n: int) -> int:
    """Inverse of (low face, high face, apex), requiring core compatibility."""
    low_core = face_mask(low, n - 1, "high")
    high_core = face_mask(high, n - 1, "low")
    if low_core != high_core:
        raise ValueError("incompatible endpoint faces")
    low_index = tile_index(n - 1)
    value = 0
    for bit, (x, y) in enumerate(tile_schema(n)[0]):
        if (x, y) == (n, 1):
            tile_value = apex
        elif y > 1:
            tile_value = (low >> low_index[(x - 1, y - 1)]) & 1
        else:
            tile_value = (high >> low_index[(x, y)]) & 1
        value |= tile_value << bit
    return value


def complement(mask: int, n: int) -> int:
    return mask ^ ((1 << m_tiles(n)) - 1)


def reflect(mask: int, n: int) -> int:
    _, sigma = tile_schema(n)
    return sum(((mask >> bit) & 1) << sigma[bit] for bit in range(len(sigma)))


def line_index(mask: int, n: int) -> int:
    """Canonical line endpoint; it is also the dense index 0..|E_n|-1."""
    return min(mask, complement(mask, n))


def line_to_tilings(index: int, n: int) -> tuple[int, int]:
    return index, complement(index, n)


def line_to_nodes(index: int, n: int, node_by_mask: list[int]) -> tuple[int, int]:
    a, b = line_to_tilings(index, n)
    return node_by_mask[a], node_by_mask[b]


def apex_zero_endpoint(index: int, n: int) -> int:
    """The unique endpoint of a complement line whose apex coordinate is zero."""
    a, b = line_to_tilings(index, n)
    return a if apex_bit(a, n) == 0 else b


def line_face(index: int, n: int, side: str) -> int:
    """Endpoint-independent face line of an upper complement line."""
    return line_index(face_mask(index, n, side), n - 1)


def line_phase_mate(index: int, n: int) -> int:
    """Other upper line with the same two bare face lines, for n >= 5."""
    endpoint = apex_zero_endpoint(index, n)
    low = face_mask(endpoint, n, "low")
    high = face_mask(endpoint, n, "high")
    mate = reconstruct_tiling(complement(low, n - 1), complement(high, n - 1), 0, n)
    return line_index(mate, n)


def normalized_counter(counter: Counter[int]) -> tuple[tuple[int, int], ...]:
    divisor = 0
    for value in counter.values():
        divisor = math.gcd(divisor, value)
    return tuple(sorted((key, value // divisor) for key, value in counter.items()))


def partition_summary(values: dict[int, object]) -> dict:
    cells: dict[object, list[int]] = defaultdict(list)
    for node, value in values.items():
        cells[value].append(node)
    return {
        "cells": len(cells),
        "largest_cell": max(map(len, cells.values())),
        "colliding_pairs": sum(len(cell) * (len(cell) - 1) // 2 for cell in cells.values()),
        "collision_cells": [cell for cell in cells.values() if len(cell) > 1],
    }


def counter_json(counter: Counter) -> dict[str, int]:
    def text(key: object) -> str:
        if isinstance(key, tuple):
            return "|".join(map(str, key))
        return str(key)

    return {text(key): value for key, value in sorted(counter.items(), key=lambda item: repr(item[0]))}


def predicted_colour_descent(n: int) -> Counter[str]:
    """Closed line counts by (upper, low-face, high-face) colours."""
    total = 1 << (m_tiles(n) - 1)
    upper_blue = 1 << (reflection_orbits(n) - 1)
    one_face_blue = 1 << (reflection_orbits(n - 1) + n - 3)
    both_faces_blue = 1 << (n - 3 + (n - 2) // 2)
    all_three_blue = 1 << (n - 3)
    return Counter(
        {
            "BBB": all_three_blue,
            "BKK": upper_blue - all_three_blue,
            "KBB": both_faces_blue - all_three_blue,
            "KBK": one_face_blue - both_faces_blue,
            "KKB": one_face_blue - both_faces_blue,
            "KKK": total
            - upper_blue
            - 2 * one_face_blue
            + both_faces_blue
            + all_three_blue,
        }
    )


def blue_fraction_exponent(n: int) -> int:
    """Return e with beta_n = 2^e, the blue fraction among E_n lines."""
    return reflection_orbits(n) - m_tiles(n)


def load_atlases(small_path: Path, n7_path: Path) -> tuple[dict[int, dict], dict[int, list[int]]]:
    small = json.loads(small_path.read_text())
    sizes = {int(size["n"]): size for size in small["sizes"]}
    sizes[7] = json.loads(n7_path.read_text())
    node_by_mask = {}
    for n, size in sizes.items():
        if n <= 6:
            ordered = sorted(size["tiling_map"], key=lambda record: int(record["mask"]))
            node_by_mask[n] = [int(record["node_rank"]) for record in ordered]
        else:
            node_by_mask[n] = list(map(int, size["node_rank_by_mask"]))
    return sizes, node_by_mask


def node_records(size: dict) -> dict[int, dict]:
    return {int(record["rank"]): record for record in size["nodes"]}


def basic_three_sorted_census(n: int, size: dict, node_map: list[int]) -> dict:
    tiles, sigma = tile_schema(n)
    tiling_count = 1 << len(tiles)
    full = tiling_count - 1
    records = node_records(size)
    masks_by_node: dict[int, list[int]] = defaultdict(list)
    for mask, node in enumerate(node_map):
        masks_by_node[node].append(mask)

    directed_blue: Counter[tuple[int, int]] = Counter()
    directed_black: Counter[tuple[int, int]] = Counter()
    projected_fibres: dict[tuple[str, int, int], list[int]] = defaultdict(list)
    line_colour_counts: Counter[str] = Counter()
    loop_line_counts: Counter[str] = Counter()
    loop_node_support: dict[str, set[int]] = {"blue": set(), "black": set()}
    for mask in range(tiling_count):
        other = mask ^ full
        u, v = node_map[mask], node_map[other]
        blue = is_grid_symmetric(mask, sigma)
        (directed_blue if blue else directed_black)[(u, v)] += 1
        if mask < other:
            colour = "blue" if blue else "black"
            line_colour_counts[colour] += 1
            a, b = sorted((u, v))
            projected_fibres[(colour, a, b)].append(mask)
            if u == v:
                loop_line_counts[colour] += 1
                loop_node_support[colour].add(u)

    fibre_size = {node: len(masks) for node, masks in masks_by_node.items()}
    incidence_failures = symmetry_failures = diagonal_failures = 0
    for kernel in (directed_blue, directed_black):
        symmetry_failures += sum(value != kernel[(v, u)] for (u, v), value in kernel.items())
    for node in masks_by_node:
        row_sum = sum(directed_blue[(node, v)] + directed_black[(node, v)] for v in masks_by_node)
        incidence_failures += row_sum != fibre_size[node]
        diagonal_failures += directed_blue[(node, node)] % 2
        diagonal_failures += directed_black[(node, node)] % 2
    assert incidence_failures == symmetry_failures == diagonal_failures == 0

    multiplicity_histogram = Counter(len(lines) for lines in projected_fibres.values())
    edge_support_signature = {}
    edge_weighted_signature = {}
    colour_degree_signature = {}
    for node in masks_by_node:
        support = []
        weighted = []
        colour_degrees = Counter()
        for colour, kernel in (("B", directed_blue), ("K", directed_black)):
            for target in masks_by_node:
                count = kernel[(node, target)]
                if count:
                    support.append((colour, target))
                    weighted.append((colour, target, count))
                    colour_degrees[(colour, "loop" if node == target else "cross")] += count
        edge_support_signature[node] = tuple(support)
        edge_weighted_signature[node] = tuple(weighted)
        colour_degree_signature[node] = tuple(sorted(colour_degrees.items()))

    projected_json = []
    for (colour, u, v), lines in sorted(projected_fibres.items()):
        projected_json.append(
            {
                "colour": colour,
                "node_u": records[u]["id"],
                "node_v": records[v]["id"],
                "line_count": len(lines),
                "line_indices": lines,
            }
        )

    return {
        "n": n,
        "tile_count": len(tiles),
        "tilings": tiling_count,
        "lines": tiling_count // 2,
        "nodes": len(masks_by_node),
        "node_fibre_size_histogram": counter_json(Counter(fibre_size.values())),
        "line_colour_counts": dict(line_colour_counts),
        "projected_coloured_supports": len(projected_fibres),
        "projected_edge_fibre_size_histogram": counter_json(multiplicity_histogram),
        "projected_loop_line_counts": dict(loop_line_counts),
        "projected_loop_node_support_counts": {
            colour: len(nodes) for colour, nodes in loop_node_support.items()
        },
        "incidence_failures": incidence_failures,
        "kernel_symmetry_failures": symmetry_failures,
        "kernel_diagonal_parity_failures": diagonal_failures,
        "carriers_base": {
            "fibre_size": fibre_size,
            "colour_degree": colour_degree_signature,
            "line_support_row": edge_support_signature,
            "line_weighted_row": edge_weighted_signature,
        },
        "projected_edge_fibres": projected_json,
    }


def recursive_census(
    n: int,
    sizes: dict[int, dict],
    node_by_mask: dict[int, list[int]],
) -> dict:
    upper_map = node_by_mask[n]
    lower_map = node_by_mask[n - 1]
    upper_records = node_records(sizes[n])
    lower_records = node_records(sizes[n - 1])
    tiles, sigma = tile_schema(n)
    _, lower_sigma = tile_schema(n - 1)
    full = (1 << len(tiles)) - 1
    upper_fibres: dict[int, list[int]] = defaultdict(list)
    lower_fibres: dict[int, list[int]] = defaultdict(list)
    for mask, node in enumerate(upper_map):
        upper_fibres[node].append(mask)
    for mask, node in enumerate(lower_map):
        lower_fibres[node].append(mask)

    low_branch: dict[int, Counter[int]] = {node: Counter() for node in upper_fibres}
    high_branch: dict[int, Counter[int]] = {node: Counter() for node in upper_fibres}
    lift_blocks: dict[tuple[int, int], Counter[int]] = defaultdict(Counter)
    pullback_failures = complement_failures = reflection_failures = core_failures = 0
    colour_descent: Counter[str] = Counter()
    line_face_pair_fibres: Counter[tuple[int, int]] = Counter()
    line_torsor_failures = 0
    apex_zero_lines: set[int] = set()
    apex_zero_bijection_failures = 0

    for mask, upper_node in enumerate(upper_map):
        low = face_mask(mask, n, "low")
        high = face_mask(mask, n, "high")
        apex = apex_bit(mask, n)
        low_core = face_mask(low, n - 1, "high")
        high_core = face_mask(high, n - 1, "low")
        core_failures += low_core != high_core
        pullback_failures += reconstruct_tiling(low, high, apex, n) != mask
        complement_failures += face_mask(mask ^ full, n, "low") != complement(low, n - 1)
        complement_failures += face_mask(mask ^ full, n, "high") != complement(high, n - 1)
        reflected = reflect(mask, n)
        reflection_failures += face_mask(reflected, n, "low") != reflect(high, n - 1)
        reflection_failures += face_mask(reflected, n, "high") != reflect(low, n - 1)
        reflection_failures += apex_bit(reflected, n) != apex
        low_node, high_node = lower_map[low], lower_map[high]
        low_branch[upper_node][low_node] += 1
        high_branch[upper_node][high_node] += 1
        lift_blocks[(upper_node, low_node)][low] += 1

        if apex == 0:
            line = line_index(mask, n)
            apex_zero_bijection_failures += line in apex_zero_lines
            apex_zero_lines.add(line)
            upper_colour = "B" if is_grid_symmetric(mask, sigma) else "K"
            low_colour = "B" if is_grid_symmetric(low, lower_sigma) else "K"
            high_colour = "B" if is_grid_symmetric(high, lower_sigma) else "K"
            colour_descent[upper_colour + low_colour + high_colour] += 1

    assert not (
        pullback_failures
        or complement_failures
        or reflection_failures
        or core_failures
        or apex_zero_bijection_failures
        or len(apex_zero_lines) != 1 << (m_tiles(n) - 1)
    )
    assert low_branch == high_branch

    row_failures = column_failures = 0
    for node, row in low_branch.items():
        row_failures += sum(row.values()) != len(upper_fibres[node])
    for lower_node, masks in lower_fibres.items():
        column_sum = sum(row[lower_node] for row in low_branch.values())
        column_failures += column_sum != (1 << (n - 2)) * len(masks)
    assert row_failures == column_failures == 0

    if n >= 5:
        for line in range(1 << (m_tiles(n) - 1)):
            low_line = line_face(line, n, "low")
            high_line = line_face(line, n, "high")
            low_core_line = line_face(low_line, n - 1, "high")
            high_core_line = line_face(high_line, n - 1, "low")
            line_torsor_failures += low_core_line != high_core_line
            line_face_pair_fibres[(low_line, high_line)] += 1
            mate = line_phase_mate(line, n)
            line_torsor_failures += mate == line
            line_torsor_failures += line_phase_mate(mate, n) != line
            line_torsor_failures += line_face(mate, n, "low") != low_line
            line_torsor_failures += line_face(mate, n, "high") != high_line
            line_torsor_failures += is_grid_symmetric(
                apex_zero_endpoint(mate, n), sigma
            ) != is_grid_symmetric(apex_zero_endpoint(line, n), sigma)
        line_torsor_failures += any(
            multiplicity != 2 for multiplicity in line_face_pair_fibres.values()
        )
        line_torsor_failures += len(line_face_pair_fibres) != 1 << (m_tiles(n) - 2)
        assert line_torsor_failures == 0

    support_signature = {node: tuple(sorted(row)) for node, row in low_branch.items()}
    weighted_signature = {
        node: tuple(sorted(row.items())) for node, row in low_branch.items()
    }
    normalized_signature = {
        node: normalized_counter(row) for node, row in low_branch.items()
    }
    support_partition = partition_summary(support_signature)
    weighted_partition = partition_summary(weighted_signature)
    normalized_partition = partition_summary(normalized_signature)

    nonuniform_blocks = []
    for (upper_node, lower_node), lifts in lift_blocks.items():
        values = [lifts[mask] for mask in lower_fibres[lower_node]]
        if len(set(values)) > 1:
            minimum, maximum = min(values), max(values)
            nonuniform_blocks.append(
                {
                    "upper_node": upper_records[upper_node]["id"],
                    "lower_node": lower_records[lower_node]["id"],
                    "minimum_lifts_per_lower_tiling": minimum,
                    "maximum_lifts_per_lower_tiling": maximum,
                    "distinct_lift_counts": sorted(set(values)),
                    "minimum_witness_masks": [
                        mask for mask in lower_fibres[lower_node] if lifts[mask] == minimum
                    ][:4],
                    "maximum_witness_masks": [
                        mask for mask in lower_fibres[lower_node] if lifts[mask] == maximum
                    ][:4],
                }
            )

    predicted = predicted_colour_descent(n)
    colour_failures = colour_descent != predicted
    assert not colour_failures

    branch_rows_json = []
    for upper_node, row in sorted(low_branch.items()):
        branch_rows_json.append(
            {
                "upper_node": upper_records[upper_node]["id"],
                "lower_nodes": [
                    {"node": lower_records[lower]["id"], "tiling_count": count}
                    for lower, count in sorted(row.items())
                ],
            }
        )

    return {
        "n": n,
        "pullback_formula": "T_n=(T_(n-1) x_[T_(n-2)] T_(n-1)) x {0,1}_apex",
        "line_formula": "E_n=T_(n-1) x_[T_(n-2)] T_(n-1) via apex-zero endpoint",
        "bare_line_tower_formula": (
            "E_n is a C2-torsor over "
            "E_(n-1) x_[E_(n-2)] E_(n-1), n>=5"
        ),
        "compatible_face_pairs": 1 << (m_tiles(n) - 1),
        "pullback_failures": pullback_failures,
        "core_compatibility_failures": core_failures,
        "complement_naturality_failures": complement_failures,
        "reflection_face_swap_failures": reflection_failures,
        "apex_zero_line_bijection_failures": apex_zero_bijection_failures,
        "branch_row_conservation_failures": row_failures,
        "branch_column_conservation_failures": column_failures,
        "low_high_branch_matrix_mismatches": int(low_branch != high_branch),
        "line_face_pair_support": len(line_face_pair_fibres),
        "line_face_pair_fibre_size_histogram": counter_json(
            Counter(line_face_pair_fibres.values())
        ),
        "line_torsor_failures": line_torsor_failures,
        "branch_support_partition": support_partition,
        "branch_weighted_partition": weighted_partition,
        "branch_normalized_partition": normalized_partition,
        "lift_blocks": len(lift_blocks),
        "nonuniform_lift_blocks": len(nonuniform_blocks),
        "strong_lumpability_failures": len(nonuniform_blocks),
        "maximum_lift_range": max(
            (
                witness["maximum_lifts_per_lower_tiling"]
                - witness["minimum_lifts_per_lower_tiling"]
                for witness in nonuniform_blocks
            ),
            default=0,
        ),
        "nonuniform_lift_witnesses": nonuniform_blocks[:24],
        "colour_descent_histogram": dict(colour_descent),
        "closed_form_colour_descent_histogram": dict(predicted),
        "blue_fraction_power_of_two_exponent": blue_fraction_exponent(n),
        "previous_blue_fraction_power_of_two_exponent": blue_fraction_exponent(n - 1),
        "blue_event_structure": (
            "upper/low/high blue are pairwise independent; upper blue forces "
            "low and high colours equal; their three-way interaction is nonzero"
        ),
        "colour_descent_formula_failures": int(colour_failures),
        "branch_rows": branch_rows_json,
        "carriers_recursive": {
            "lower_face_support": support_signature,
            "lower_face_weighted": weighted_signature,
            "lower_face_normalized": normalized_signature,
        },
    }


def attach_tournament_analysis(size_result: dict, recursive: dict | None) -> None:
    base = size_result.pop("carriers_base")
    carriers = {
        "fibre_size": base["fibre_size"],
        "colour_degree": base["colour_degree"],
        "line_support_row": base["line_support_row"],
        "line_weighted_row": base["line_weighted_row"],
    }
    if recursive is not None:
        carriers.update(recursive.pop("carriers_recursive"))
    carriers["exact_node"] = {node: node for node in base["fibre_size"]}
    stats = {name: partition_stats(values) for name, values in carriers.items()}
    retention = carrier_tournament(stats, "retention")
    economy = carrier_tournament(stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carriers))
        for j in range(i + 1, len(carriers))
    )
    size_result["carrier_stats"] = stats
    size_result["tournament_analysis"] = {
        "vertices": list(carriers),
        "pairwise_observable": "number of unordered merged-node pairs separated by the carrier",
        "switches": ["retention", "retention_per_cell"],
        "tie_hamiltonian_path": list(carriers),
        "retention": retention,
        "economy": economy,
        "edge_flips_between_gauges": flips,
    }


def render(result: dict) -> str:
    lines = [
        "THM-796 THREE-SORTED RECURSIVE METAGRAPH INCIDENCE",
        "=" * 82,
        "T_n = (T_(n-1) x_[T_(n-2)] T_(n-1)) x {0,1}_apex",
        "E_n = T_(n-1) x_[T_(n-2)] T_(n-1) (choose apex-zero endpoint)",
        "A^colour_uv = #{t in fibre(u): complement(t) in fibre(v), colour(t)=colour}",
        "D_n(u,v) = #{t in fibre(u): low_face(t) in fibre(v)}",
        "",
    ]
    recursive_by_n = {row["n"]: row for row in result["recursive_sizes"]}
    for size in result["sizes"]:
        n = size["n"]
        lines.extend(
            [
                f"n={n}: tilings={size['tilings']} lines={size['lines']} nodes={size['nodes']}",
                f"  line colours={size['line_colour_counts']} projected supports="
                f"{size['projected_coloured_supports']} loops(lines/nodes)="
                f"{size['projected_loop_line_counts']}/"
                f"{size['projected_loop_node_support_counts']}",
                f"  edge-fibre multiplicities={size['projected_edge_fibre_size_histogram']}",
                f"  incidence/symmetry/diagonal failures="
                f"{size['incidence_failures']}/{size['kernel_symmetry_failures']}/"
                f"{size['kernel_diagonal_parity_failures']}",
            ]
        )
        if n in recursive_by_n:
            rec = recursive_by_n[n]
            lines.extend(
                [
                    f"  branch partitions support/weighted/normalized="
                    f"{rec['branch_support_partition']['cells']}/"
                    f"{rec['branch_weighted_partition']['cells']}/"
                    f"{rec['branch_normalized_partition']['cells']} of {size['nodes']}",
                    f"  lift blocks nonuniform/total={rec['nonuniform_lift_blocks']}/"
                    f"{rec['lift_blocks']} max_range={rec['maximum_lift_range']}",
                    f"  bare-line face-pair support/fibres/torsor-failures="
                    f"{rec['line_face_pair_support']}/"
                    f"{rec['line_face_pair_fibre_size_histogram']}/"
                    f"{rec['line_torsor_failures']}",
                    f"  colour descent={rec['colour_descent_histogram']}",
                    f"  recursion failures pullback/core/complement/reflection/apex/row/column/colour="
                    f"{rec['pullback_failures']}/{rec['core_compatibility_failures']}/"
                    f"{rec['complement_naturality_failures']}/"
                    f"{rec['reflection_face_swap_failures']}/"
                    f"{rec['apex_zero_line_bijection_failures']}/"
                    f"{rec['branch_row_conservation_failures']}/"
                    f"{rec['branch_column_conservation_failures']}/"
                    f"{rec['colour_descent_formula_failures']}",
                ]
            )
        lines.append("")

    n7 = result["sizes"][-1]
    r7 = recursive_by_n[7]
    lines.extend(
        [
            "N=7 INFORMATION BOUNDARY",
            f"  lower-face support alone: {r7['branch_support_partition']['cells']}/272 cells, "
            f"{r7['branch_support_partition']['colliding_pairs']} twin pairs",
            f"  weighted row: {r7['branch_weighted_partition']['cells']}/272; "
            f"normalized row: {r7['branch_normalized_partition']['cells']}/272",
            f"  node-only strong lumpability fails on {r7['nonuniform_lift_blocks']}/"
            f"{r7['lift_blocks']} nonzero node blocks",
            "  conclusion: relative lower-face multiplicities identify every audited node,",
            "  but a lower node does not make its tilings exchangeable for the next lift",
            "",
            "TOURNAMENT ANALYSIS (n=7)",
            f"  carriers={tuple(n7['tournament_analysis']['vertices'])}",
            f"  retention/economy edge flips="
            f"{n7['tournament_analysis']['edge_flips_between_gauges']}",
            f"  tie Hamiltonian path="
            f"{tuple(n7['tournament_analysis']['tie_hamiltonian_path'])}",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=SMALL_ATLAS)
    parser.add_argument("--n7-atlas", type=Path, default=N7_ATLAS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    sizes, node_by_mask = load_atlases(args.small_atlas, args.n7_atlas)
    result = {
        "schema_version": 1,
        "theorem": "THM-796",
        "sorts": {
            "tilings": "T_n={0,1}^binom(n-1,2)",
            "lines": "E_n=T_n/complement",
            "nodes": "N_n=tournament isomorphism classes modulo converse",
        },
        "sizes": [],
        "recursive_sizes": [],
    }
    recursive_by_n = {}
    for n in range(4, 8):
        recursive_by_n[n] = recursive_census(n, sizes, node_by_mask)
        result["recursive_sizes"].append(recursive_by_n[n])
    for n in range(3, 8):
        size_result = basic_three_sorted_census(n, sizes[n], node_by_mask[n])
        attach_tournament_analysis(size_result, recursive_by_n.get(n))
        result["sizes"].append(size_result)

    output = render(result)
    print(output, end="")
    if args.output:
        args.output.write_text(output)
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
