#!/usr/bin/env python3
"""Exact factorization refinement of THM-830's node-coloured defect groupoid.

For a fixed trace component, a tiling is an arrow (x, delta) in the elementary
abelian pair groupoid.  Its base colour retains the *full* reflection defect
and its converse-merged tournament node.  We refine an arrow colour by the
multiset of colour pairs occurring in all two-arrow factorizations.

Two versions are audited:

* ordered: (colour(first), colour(second)), the ordinary composition tensor;
* symmetric: {colour(first), colour(second)}, the reversal-compatible Jordan
  tensor appropriate before quotienting tilings by grid reflection.

The exact committed tiling-to-node atlases make this finite and exhaustive for
3 <= n <= 7.  Tournament Analysis uses information carriers, rather than
runners or original arcs, as vertices.  The pairwise observable is separated
arrow pairs; the switch divides this retention by the number of address bits.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Hashable, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SMALL = ROOT / "05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json"
DEFAULT_N7 = ROOT / "05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json"
DEFAULT_OUTPUT = ROOT / "05-knowledge/results/node_coloured_defect_factorization_codex_S15.out"
DEFAULT_JSON = ROOT / "05-knowledge/results/node_coloured_defect_factorization_codex_S15.json"


@dataclass(frozen=True)
class Atlas:
    n: int
    tiles: tuple[tuple[int, int], ...]
    node_by_mask: tuple[int, ...]


@dataclass(frozen=True)
class Chart:
    n: int
    tiles: tuple[tuple[int, int], ...]
    reflection: tuple[int, ...]
    fixed: tuple[int, ...]
    pairs: tuple[tuple[int, int], ...]
    embedded_defect: tuple[int, ...]

    @property
    def h(self) -> int:
        return len(self.fixed) + len(self.pairs)

    @property
    def r(self) -> int:
        return len(self.pairs)

    @property
    def objects(self) -> int:
        return 1 << self.h

    @property
    def defects(self) -> int:
        return 1 << self.r

    def decode(self, obj: int, defect: int) -> int:
        mask = 0
        for coordinate, tile_bit in enumerate(self.fixed):
            mask |= ((obj >> coordinate) & 1) << tile_bit
        offset = len(self.fixed)
        for coordinate, (left, right) in enumerate(self.pairs):
            value = (obj >> (offset + coordinate)) & 1
            mask |= value << left
            mask |= (value ^ ((defect >> coordinate) & 1)) << right
        return mask

    def reflect_mask(self, mask: int) -> int:
        return sum(
            ((mask >> source) & 1) << target
            for source, target in enumerate(self.reflection)
        )


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_atlases(small_path: Path, n7_path: Path) -> dict[int, Atlas]:
    small = json.loads(small_path.read_text())
    answer: dict[int, Atlas] = {}
    for size in small["sizes"]:
        n = int(size["n"])
        nodes = [-1] * (1 << int(size["tile_count"]))
        for row in size["tiling_map"]:
            nodes[int(row["mask"])] = int(row["node_rank"])
        assert min(nodes) >= 0
        answer[n] = Atlas(
            n=n,
            tiles=tuple(tuple(map(int, tile)) for tile in size["tile_order"]),
            node_by_mask=tuple(nodes),
        )

    n7 = json.loads(n7_path.read_text())
    assert int(n7["n"]) == 7
    answer[7] = Atlas(
        n=7,
        tiles=tuple(tuple(map(int, tile)) for tile in n7["tile_order"]),
        node_by_mask=tuple(map(int, n7["node_rank_by_mask"])),
    )
    assert set(answer) == set(range(3, 8))
    return answer


def make_chart(atlas: Atlas) -> Chart:
    position = {tile: bit for bit, tile in enumerate(atlas.tiles)}
    reflection = tuple(
        position[(atlas.n - y + 1, atlas.n - x + 1)]
        for x, y in atlas.tiles
    )
    assert all(reflection[reflection[bit]] == bit for bit in range(len(reflection)))
    fixed = tuple(bit for bit, target in enumerate(reflection) if bit == target)
    pairs = tuple((bit, target) for bit, target in enumerate(reflection) if bit < target)
    embedded = tuple(
        sum(
            ((defect >> coordinate) & 1) << (len(fixed) + coordinate)
            for coordinate in range(len(pairs))
        )
        for defect in range(1 << len(pairs))
    )
    chart = Chart(atlas.n, atlas.tiles, reflection, fixed, pairs, embedded)
    assert chart.h == ((atlas.n - 1) ** 2) // 4
    assert chart.r == ((atlas.n - 2) ** 2) // 4
    assert chart.h + chart.r == len(atlas.tiles)
    return chart


def dense_rank(values: Iterable[Hashable]) -> list[int]:
    materialized = list(values)
    keys = sorted(set(materialized))
    rank = {key: index for index, key in enumerate(keys)}
    return [rank[value] for value in materialized]


def fibre_profile(values: Iterable[Hashable]) -> dict[int, int]:
    fibres = Counter(Counter(values).values())
    return dict(sorted(fibres.items()))


def separated_pairs(values: Iterable[Hashable]) -> int:
    counts = Counter(values)
    total = sum(counts.values())
    return total * (total - 1) // 2 - sum(size * (size - 1) // 2 for size in counts.values())


def refinement_step(
    colours: list[int], chart: Chart, symmetric: bool
) -> tuple[list[int], list[tuple[int, tuple[tuple[tuple[int, int], int], ...]]]]:
    """One exact two-factor refinement step.

    The returned semantic signature includes the old colour, so the operation
    can only split cells.  In symmetric mode each factor-colour pair is sorted;
    this forces invariance under arrow reversal.
    """
    dcount = chart.defects
    signatures = []
    for obj in range(chart.objects):
        for result_defect in range(dcount):
            factors: Counter[tuple[int, int]] = Counter()
            for left_defect in range(dcount):
                left = colours[obj * dcount + left_defect]
                middle = obj ^ chart.embedded_defect[left_defect]
                right_defect = left_defect ^ result_defect
                right = colours[middle * dcount + right_defect]
                pair = (min(left, right), max(left, right)) if symmetric else (left, right)
                factors[pair] += 1
            signatures.append((colours[obj * dcount + result_defect], tuple(sorted(factors.items()))))
    return dense_rank(signatures), signatures


def unlabelled_node_factor_step(
    base: list[int],
    base_semantic: list[tuple[int, int]],
    chart: Chart,
    symmetric: bool,
) -> list[int]:
    """Forget which exact factor defect carries each intermediate node pair."""
    dcount = chart.defects
    signatures = []
    for obj in range(chart.objects):
        for result_defect in range(dcount):
            factors: Counter[tuple[int, int]] = Counter()
            for left_defect in range(dcount):
                left_node = base_semantic[obj * dcount + left_defect][1]
                middle = obj ^ chart.embedded_defect[left_defect]
                right_node = base_semantic[
                    middle * dcount + (left_defect ^ result_defect)
                ][1]
                pair = (
                    (min(left_node, right_node), max(left_node, right_node))
                    if symmetric
                    else (left_node, right_node)
                )
                factors[pair] += 1
            signatures.append((base[obj * dcount + result_defect], tuple(sorted(factors.items()))))
    return dense_rank(signatures)


def first_coarse_collision(
    coarse: list[int], fine: list[int], masks: list[int], chart: Chart
) -> dict | None:
    fibres: dict[int, dict[int, int]] = defaultdict(dict)
    for index, (coarse_colour, fine_colour) in enumerate(zip(coarse, fine)):
        fibres[coarse_colour].setdefault(fine_colour, index)
    for images in fibres.values():
        if len(images) < 2:
            continue
        left, right = sorted(images.values(), key=lambda index: masks[index])[:2]
        left_obj, left_defect = divmod(left, chart.defects)
        right_obj, right_defect = divmod(right, chart.defects)
        assert left_defect == right_defect
        return {
            "left_mask": masks[left],
            "right_mask": masks[right],
            "defect": left_defect,
            "left_object": left_obj,
            "right_object": right_obj,
        }
    return None


def split_cell_count(old: list[int], new: list[int]) -> int:
    images: dict[int, set[int]] = defaultdict(set)
    for before, after in zip(old, new):
        images[before].add(after)
    return sum(len(targets) > 1 for targets in images.values())


def arrow_factor_semantics(
    index: int,
    chart: Chart,
    base_semantic: list[tuple[int, int]],
    symmetric: bool,
) -> Counter[tuple[tuple[int, int], tuple[int, int]]]:
    dcount = chart.defects
    obj, result_defect = divmod(index, dcount)
    answer: Counter[tuple[tuple[int, int], tuple[int, int]]] = Counter()
    for left_defect in range(dcount):
        left_index = obj * dcount + left_defect
        middle = obj ^ chart.embedded_defect[left_defect]
        right_index = middle * dcount + (left_defect ^ result_defect)
        pair = (base_semantic[left_index], base_semantic[right_index])
        if symmetric and pair[1] < pair[0]:
            pair = (pair[1], pair[0])
        answer[pair] += 1
    return answer


def first_witness(
    atlas: Atlas,
    chart: Chart,
    masks: list[int],
    base: list[int],
    base_semantic: list[tuple[int, int]],
    refined: list[int],
    symmetric: bool,
) -> dict | None:
    fibres: dict[int, list[int]] = defaultdict(list)
    for index, colour in enumerate(base):
        fibres[colour].append(index)
    for indices in fibres.values():
        by_new: dict[int, int] = {}
        for index in indices:
            by_new.setdefault(refined[index], index)
        if len(by_new) < 2:
            continue
        left, right = sorted(by_new.values(), key=lambda index: masks[index])[:2]
        left_factors = arrow_factor_semantics(left, chart, base_semantic, symmetric)
        right_factors = arrow_factor_semantics(right, chart, base_semantic, symmetric)
        separator = next(
            key
            for key in sorted(set(left_factors) | set(right_factors))
            if left_factors[key] != right_factors[key]
        )

        def arrow_record(index: int) -> dict:
            obj, defect = divmod(index, chart.defects)
            mask = masks[index]
            return {
                "mask": mask,
                "bit_word": format(mask, f"0{len(atlas.tiles)}b"),
                "selected_tiles": [list(atlas.tiles[bit]) for bit in range(len(atlas.tiles)) if (mask >> bit) & 1],
                "object": obj,
                "defect": defect,
                "defect_word": format(defect, f"0{chart.r}b"),
                "node": atlas.node_by_mask[mask],
            }

        return {
            "left": arrow_record(left),
            "right": arrow_record(right),
            "separating_factor_colour_pair": [
                {"defect": separator[0][0], "node": separator[0][1]},
                {"defect": separator[1][0], "node": separator[1][1]},
            ],
            "left_multiplicity": left_factors[separator],
            "right_multiplicity": right_factors[separator],
        }
    return None


def strongly_connected_components(vertices: list[str], adjacency: dict[str, set[str]]) -> list[int]:
    reach = {v: {v} | set(adjacency[v]) for v in vertices}
    for pivot in vertices:
        for left in vertices:
            if pivot in reach[left]:
                reach[left] |= reach[pivot]
    unseen = set(vertices)
    sizes = []
    while unseen:
        seed = min(unseen)
        component = {v for v in unseen if v in reach[seed] and seed in reach[v]}
        sizes.append(len(component))
        unseen -= component
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(vertices: list[str], adjacency: dict[str, set[str]]) -> int:
    count = len(vertices)
    dp = [[0] * count for _ in range(1 << count)]
    for vertex in range(count):
        dp[1 << vertex][vertex] = 1
    for subset in range(1 << count):
        for last in range(count):
            ways = dp[subset][last]
            if not ways:
                continue
            for nxt in range(count):
                if (subset >> nxt) & 1 or vertices[nxt] not in adjacency[vertices[last]]:
                    continue
                dp[subset | (1 << nxt)][nxt] += ways
    return sum(dp[-1])


def tournament_fingerprint(vertices: list[str], metric: dict[str, Fraction]) -> tuple[dict, dict[frozenset[str], str]]:
    adjacency = {vertex: set() for vertex in vertices}
    winners: dict[frozenset[str], str] = {}
    for left, right in itertools.combinations(vertices, 2):
        winner, loser = (
            (left, right)
            if (metric[left], left) > (metric[right], right)
            else (right, left)
        )
        adjacency[winner].add(loser)
        winners[frozenset((left, right))] = winner
    scores = Counter(len(adjacency[vertex]) for vertex in vertices)
    triangles = sum(
        b in adjacency[a] and c in adjacency[b] and a in adjacency[c]
        for a, b, c in itertools.permutations(vertices, 3)
    ) // 3
    return (
        {
            "score_histogram": dict(sorted(scores.items())),
            "directed_3cycles": triangles,
            "scc_sizes": strongly_connected_components(vertices, adjacency),
            "hamiltonian_path_count": hamiltonian_path_count(vertices, adjacency),
            "tie_hamiltonian_path": sorted(vertices, key=lambda name: (metric[name], name), reverse=True),
        },
        winners,
    )


def carrier_tournament(carriers: dict[str, list[Hashable]]) -> dict:
    vertices = list(carriers)
    retention = {name: separated_pairs(values) for name, values in carriers.items()}
    cells = {name: len(set(values)) for name, values in carriers.items()}
    raw_metric = {name: Fraction(retention[name]) for name in vertices}
    economy_metric = {
        name: Fraction(retention[name], max(1, math.ceil(math.log2(max(2, cells[name])))))
        for name in vertices
    }
    raw, raw_winners = tournament_fingerprint(vertices, raw_metric)
    economy, economy_winners = tournament_fingerprint(vertices, economy_metric)
    return {
        "vertices": vertices,
        "pairwise_observable": "number of unordered arrow pairs separated",
        "switch": "retention divided by ceiling(log2(number of cells))",
        "tie_rule": "carrier name",
        "cells": cells,
        "separated_pairs": retention,
        "raw_retention": raw,
        "retention_per_address_bit": economy,
        "edge_flips_between_gauges": sum(
            raw_winners[pair] != economy_winners[pair] for pair in raw_winners
        ),
    }


def analyze_size(atlas: Atlas) -> tuple[dict, dict[str, list[Hashable]]]:
    chart = make_chart(atlas)
    masks = [
        chart.decode(obj, defect)
        for obj in range(chart.objects)
        for defect in range(chart.defects)
    ]
    assert len(masks) == len(set(masks)) == len(atlas.node_by_mask)
    index_by_mask = {mask: index for index, mask in enumerate(masks)}
    assert all(
        atlas.node_by_mask[mask] == atlas.node_by_mask[chart.reflect_mask(mask)]
        for mask in masks
    )

    base_semantic = [
        (defect, atlas.node_by_mask[masks[obj * chart.defects + defect]])
        for obj in range(chart.objects)
        for defect in range(chart.defects)
    ]
    base = dense_rank(base_semantic)
    ordered, _ordered_signatures = refinement_step(base, chart, symmetric=False)
    symmetric, _symmetric_signatures = refinement_step(base, chart, symmetric=True)
    unlabelled_ordered = unlabelled_node_factor_step(base, base_semantic, chart, symmetric=False)
    unlabelled_symmetric = unlabelled_node_factor_step(base, base_semantic, chart, symmetric=True)

    endpoint_ordered = []
    endpoint_symmetric = []
    for obj in range(chart.objects):
        source_node = base_semantic[obj * chart.defects][1]
        for defect in range(chart.defects):
            target = obj ^ chart.embedded_defect[defect]
            target_node = base_semantic[target * chart.defects][1]
            endpoint_ordered.append((base[obj * chart.defects + defect], source_node, target_node))
            endpoint_symmetric.append(
                (
                    base[obj * chart.defects + defect],
                    min(source_node, target_node),
                    max(source_node, target_node),
                )
            )
    endpoint_ordered = dense_rank(endpoint_ordered)
    endpoint_symmetric = dense_rank(endpoint_symmetric)

    assert all(
        symmetric[obj * chart.defects + defect]
        == symmetric[(obj ^ chart.embedded_defect[defect]) * chart.defects + defect]
        for obj in range(chart.objects)
        for defect in range(chart.defects)
    )
    arrow_count = len(masks)
    reflection_orbits = chart.objects + (chart.defects - 1) * chart.objects // 2
    if atlas.n >= 5:
        assert len(set(ordered)) == arrow_count
        assert len(set(symmetric)) == reflection_orbits
    else:
        ordered_next, _ = refinement_step(ordered, chart, symmetric=False)
        symmetric_next, _ = refinement_step(symmetric, chart, symmetric=True)
        assert len(set(ordered_next)) == len(set(ordered))
        assert len(set(symmetric_next)) == len(set(symmetric))

    full = arrow_count - 1
    line_keys = []
    for mask in range(arrow_count):
        mate = mask ^ full
        if mask > mate:
            continue
        left = symmetric[index_by_mask[mask]]
        right = symmetric[index_by_mask[mate]]
        line_keys.append((min(left, right), max(left, right)))
    line_reflection_orbits = len(set(line_keys))
    literal_line_representatives = [
        mask for mask in range(arrow_count) if mask < (mask ^ full)
    ]
    direct_line_orbits = {
        min(
            mask,
            min(chart.reflect_mask(mask), chart.reflect_mask(mask) ^ full),
        )
        for mask in literal_line_representatives
    }
    expected_line_reflection_orbits = len(direct_line_orbits)
    # Burnside form: every blue line is fixed and black lines form pairs.
    blue_lines = 1 << (chart.h - 1)
    all_lines = 1 << (len(atlas.tiles) - 1)
    burnside_line_orbits = blue_lines + (all_lines - blue_lines) // 2
    assert expected_line_reflection_orbits == burnside_line_orbits
    assert line_reflection_orbits == burnside_line_orbits

    node_values = [semantic[1] for semantic in base_semantic]
    weight_values = [(semantic[1], semantic[0].bit_count()) for semantic in base_semantic]
    line_colour = [semantic[0] == 0 for semantic in base_semantic]
    carriers = {
        "blue/black": line_colour,
        "merged node": node_values,
        "node + defect weight": weight_values,
        "node + exact defect": base,
        "node + defect + endpoints": endpoint_symmetric,
        "unlabelled midpoint-node deck": unlabelled_symmetric,
        "symmetric factor colour": symmetric,
        "ordered factor colour": ordered,
    }

    initial_orbit_sizes: dict[int, int] = {}
    for colour in set(base):
        members = [index for index, value in enumerate(base) if value == colour]
        defect = base_semantic[members[0]][0]
        initial_orbit_sizes[colour] = len(members) if defect == 0 else len(members) // 2

    result = {
        "n": atlas.n,
        "tile_count": len(atlas.tiles),
        "fixed_trace_bits": len(chart.fixed),
        "defect_rank": chart.r,
        "objects": chart.objects,
        "arrows_tilings": arrow_count,
        "reflection_orbit_tilings": reflection_orbits,
        "base_cells_node_plus_exact_defect": len(set(base)),
        "base_reflection_orbit_fibre_profile": dict(sorted(Counter(initial_orbit_sizes.values()).items())),
        "base_cells_split_by_ordered_factorization": split_cell_count(base, ordered),
        "base_cells_split_by_symmetric_factorization": split_cell_count(base, symmetric),
        "ordered_one_step_cells": len(set(ordered)),
        "symmetric_one_step_cells": len(set(symmetric)),
        "ordered_endpoint_node_cells": len(set(endpoint_ordered)),
        "symmetric_endpoint_node_cells": len(set(endpoint_symmetric)),
        "unlabelled_ordered_factor_cells": len(set(unlabelled_ordered)),
        "unlabelled_symmetric_factor_cells": len(set(unlabelled_symmetric)),
        "unlabelled_symmetric_residual_reflection_orbits": (
            reflection_orbits - len(set(unlabelled_symmetric))
        ),
        "first_unlabelled_factor_collision": first_coarse_collision(
            unlabelled_symmetric, symmetric, masks, chart
        ),
        "ordered_one_step_is_discrete": len(set(ordered)) == arrow_count,
        "symmetric_one_step_is_exactly_reflection_orbits": len(set(symmetric)) == reflection_orbits,
        "complement_lines": all_lines,
        "blue_lines": blue_lines,
        "line_reflection_orbits_burnside": burnside_line_orbits,
        "line_keys_from_symmetric_factor_colour": line_reflection_orbits,
        "ordered_witness": first_witness(atlas, chart, masks, base, base_semantic, ordered, False),
        "symmetric_witness": first_witness(atlas, chart, masks, base, base_semantic, symmetric, True),
    }
    return result, carriers


def run(small_path: Path, n7_path: Path) -> dict:
    atlases = load_atlases(small_path, n7_path)
    sizes = []
    carrier_by_n = {}
    for n in range(3, 8):
        size, carriers = analyze_size(atlases[n])
        sizes.append(size)
        carrier_by_n[n] = carriers

    result = {
        "schema_version": 1,
        "theorem": "THM-851",
        "inputs": {
            str(small_path.relative_to(ROOT)): sha256(small_path),
            str(n7_path.relative_to(ROOT)): sha256(n7_path),
        },
        "sizes": sizes,
        "tournament_analysis_n7": carrier_tournament(carrier_by_n[7]),
        "preservation": {
            "preserves": [
                "fixed-path tiling", "full reflection defect", "merged tournament node",
                "two-arrow factorization multiplicity", "grid reflection", "complement line incidence",
                "blue/black line colour",
            ],
            "destroys": [
                "runner speeds", "metric scale", "LRC gaps", "owners", "wall chronology",
                "continued-fraction carry", "loneliness predicate",
            ],
            "challenged_assumption": (
                "the useful vertices are not runners or tournament vertices: they are coloured "
                "groupoid arrows and factorization obligations"
            ),
        },
    }

    expected = {
        3: (2, 2, 2),
        4: (4, 6, 5),
        5: (23, 64, 40),
        6: (255, 1024, 544),
        7: (7926, 32768, 16640),
    }
    for size in sizes:
        assert (
            size["base_cells_node_plus_exact_defect"],
            size["ordered_one_step_cells"],
            size["symmetric_one_step_cells"],
        ) == expected[size["n"]]
    expected_unlabelled = {
        3: (2, 2),
        4: (6, 5),
        5: (62, 39),
        6: (1022, 543),
        7: (32688, 16604),
    }
    for size in sizes:
        assert (
            size["unlabelled_ordered_factor_cells"],
            size["unlabelled_symmetric_factor_cells"],
        ) == expected_unlabelled[size["n"]]
    return result


def render(result: dict) -> str:
    lines = [
        "THM-851 NODE-COLOURED DEFECT FACTORIZATION CLOSURE",
        "=" * 72,
        "n  (f,r) objects arrows  C0  ordered-C1  symmetric-C1  sigma-orbits  Q-lines",
    ]
    for row in result["sizes"]:
        lines.append(
            f"{row['n']}  ({row['fixed_trace_bits']},{row['defect_rank']})"
            f" {row['objects']:7d} {row['arrows_tilings']:6d}"
            f" {row['base_cells_node_plus_exact_defect']:5d}"
            f" {row['ordered_one_step_cells']:11d}"
            f" {row['symmetric_one_step_cells']:13d}"
            f" {row['reflection_orbit_tilings']:13d}"
            f" {row['line_keys_from_symmetric_factor_colour']:7d}"
        )
        lines.append(
            "   C0 orbit-fibres=" + str(row["base_reflection_orbit_fibre_profile"])
            + f"; split cells ordered/symmetric="
            f"{row['base_cells_split_by_ordered_factorization']}/"
            f"{row['base_cells_split_by_symmetric_factorization']}"
        )
        lines.append(
            f"   endpoint O/S={row['ordered_endpoint_node_cells']}/"
            f"{row['symmetric_endpoint_node_cells']}; unlabelled-factor O/S="
            f"{row['unlabelled_ordered_factor_cells']}/"
            f"{row['unlabelled_symmetric_factor_cells']}; residual sigma-orbits="
            f"{row['unlabelled_symmetric_residual_reflection_orbits']}"
        )
        if row["symmetric_witness"]:
            witness = row["symmetric_witness"]
            lines.append(
                f"   first symmetric witness masks="
                f"{witness['left']['mask']},{witness['right']['mask']}"
                f" delta={witness['left']['defect']} node={witness['left']['node']}"
                f" separating multiplicity={witness['left_multiplicity']}/"
                f"{witness['right_multiplicity']}"
            )

    ta = result["tournament_analysis_n7"]
    lines.extend(
        [
            "",
            "n=7 TOURNAMENT ANALYSIS ON INFORMATION CARRIERS",
            f"vertices={ta['vertices']}",
            f"cells={ta['cells']}",
            f"raw scores/cycles/SCC/Hpaths="
            f"{ta['raw_retention']['score_histogram']}/"
            f"{ta['raw_retention']['directed_3cycles']}/"
            f"{ta['raw_retention']['scc_sizes']}/"
            f"{ta['raw_retention']['hamiltonian_path_count']}",
            f"economy scores/cycles/SCC/Hpaths="
            f"{ta['retention_per_address_bit']['score_histogram']}/"
            f"{ta['retention_per_address_bit']['directed_3cycles']}/"
            f"{ta['retention_per_address_bit']['scc_sizes']}/"
            f"{ta['retention_per_address_bit']['hamiltonian_path_count']}",
            f"edge flips between gauges={ta['edge_flips_between_gauges']}",
            "",
            "PRESERVES: tilings, exact defect, merged nodes, factor counts, reflection, line incidence, colour",
            "DESTROYS: all runner/metric/owner/wall/carry/loneliness data",
            "CHALLENGED ASSUMPTION: factorization obligations, not runners or original arcs, are the vertices",
            "ALL ASSERTIONS PASSED",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    if not __debug__:
        raise RuntimeError("THM-851 verification requires assertions; do not run with python -O")
    parser = argparse.ArgumentParser()
    parser.add_argument("--small-atlas", type=Path, default=DEFAULT_SMALL)
    parser.add_argument("--n7-atlas", type=Path, default=DEFAULT_N7)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    args = parser.parse_args()
    result = run(args.small_atlas, args.n7_atlas)
    args.output.write_text(render(result))
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(render(result), end="")


if __name__ == "__main__":
    main()
