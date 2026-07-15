#!/usr/bin/env python3
"""Exact R7 -> n8 regression for THM-818's ordered Cech join.

This checks the join construction one size below the n=9 frontier, where
THM-809 independently gives the answer 418.  The overlap keys are ordered
pairs of literal masks and every role has its own coordinate projection.

Assumption challenge: the proof-bearing vertices are equality witnesses
(x,y), not tournament vertices.  The regression preserves literal ordered
overlaps and upper apex/color, but destroys LRC metric gaps, owners, and all
continuation outside the audited face map.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament


ATLAS = Path("05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json")


def tiles(n: int) -> list[tuple[int, int]]:
    return [(a, b) for b in range(1, n - 1) for a in range(n, b + 1, -1)]


def reflection_map(n: int, ts: list[tuple[int, int]]) -> list[int]:
    position = {tile: i for i, tile in enumerate(ts)}
    return [position[(n - b + 1, n - a + 1)] for a, b in ts]


def symmetric(mask: int, reflection: list[int]) -> bool:
    return all(((mask >> i) & 1) == ((mask >> j) & 1) for i, j in enumerate(reflection))


def face_maps(n: int) -> list[dict[int, int]]:
    upper, lower = tiles(n), tiles(n - 1)
    position = {tile: i for i, tile in enumerate(lower)}
    maps = [dict(), dict(), dict()]
    for upper_bit, (a, b) in enumerate(upper):
        if a < n:
            maps[0][upper_bit] = position[(a, b)]
        if a - b >= 3:
            maps[1][upper_bit] = position[(a - 1, b)]
        if b >= 2:
            maps[2][upper_bit] = position[(a - 1, b - 1)]
    return maps


def pair_projections(maps: list[dict[int, int]], lower_size: int) -> dict[tuple[int, int], list[list[int]]]:
    result: dict[tuple[int, int], list[list[int]]] = {}
    for s, t in ((0, 1), (0, 2), (1, 2)):
        shared = sorted(set(maps[s]) & set(maps[t]))
        assert len(shared) == 10  # the literal size-six staircase overlap
        tables = []
        for role in (s, t):
            table = [0] * (1 << lower_size)
            for mask in range(1 << lower_size):
                projected = 0
                for target_bit, upper_bit in enumerate(shared):
                    projected |= ((mask >> maps[role][upper_bit]) & 1) << target_bit
                table[mask] = projected
            tables.append(table)
        result[(s, t)] = tables
    return result


def glue(faces: tuple[int, int, int], maps: list[dict[int, int]]) -> int:
    assigned: dict[int, int] = {}
    for role, mask in enumerate(faces):
        for upper_bit, lower_bit in maps[role].items():
            value = (mask >> lower_bit) & 1
            previous = assigned.setdefault(upper_bit, value)
            assert previous == value
    assert len(assigned) == len(tiles(8))
    return sum(value << bit for bit, value in assigned.items())


def audit(atlas_path: Path) -> dict:
    atlas = json.loads(atlas_path.read_text())
    node = atlas["node_rank_by_mask"]
    assert len(node) == 1 << 15

    lower_n, lower_bits = 7, 15
    full = (1 << lower_bits) - 1
    sigma7 = reflection_map(lower_n, tiles(lower_n))
    fibres: dict[tuple[int, int, bool], list[int]] = defaultdict(list)
    for mask in range(1 << lower_bits):
        fibres[(node[mask], node[mask ^ full], symmetric(mask, sigma7))].append(mask)

    rows: list[tuple[int, int, bool]] = []
    for fibre in fibres.values():
        rows.extend((x, y, x == y) for x in fibre for y in fibre)
    relation_rows = len(rows)
    off_diagonal = relation_rows - (1 << lower_bits)
    assert relation_rows == 113_632 and off_diagonal == 80_864

    maps = face_maps(8)
    projections = pair_projections(maps, lower_bits)
    ab_a, ab_b = projections[(0, 1)]
    ac_a, ac_c = projections[(0, 2)]
    bc_b, bc_c = projections[(1, 2)]

    # The two coordinates of every relation row remain ordered in all keys.
    def ordered_key(table: list[int], row: tuple[int, int, bool]) -> int:
        return table[row[0]] | (table[row[1]] << 10)

    b_index: dict[int, dict[int, list[int]]] = defaultdict(lambda: defaultdict(list))
    c_index: dict[int, dict[int, list[int]]] = defaultdict(lambda: defaultdict(list))
    for index, row in enumerate(rows):
        b_index[ordered_key(ab_b, row)][ordered_key(bc_b, row)].append(index)
        c_index[ordered_key(ac_c, row)][ordered_key(bc_c, row)].append(index)

    sigma8 = reflection_map(8, tiles(8))
    nontrivial = apex_zero = upper_colour = 0
    first_off_diagonal = {"A": 0, "B": 0, "C": 0}
    first_off_diagonal_apex = {"A": 0, "B": 0, "C": 0}
    first_off_diagonal_colour = {"A": 0, "B": 0, "C": 0}
    canonical_by_first = {"A": set(), "B": set(), "C": set()}
    canonical_pairs: set[tuple[int, int]] = set()
    for row_a in rows:
        b_groups = b_index.get(ordered_key(ab_a, row_a))
        c_groups = c_index.get(ordered_key(ac_a, row_a))
        if not b_groups or not c_groups:
            continue
        for bc_key in b_groups.keys() & c_groups.keys():
            for b_index_value in b_groups[bc_key]:
                row_b = rows[b_index_value]
                for c_index_value in c_groups[bc_key]:
                    row_c = rows[c_index_value]
                    if row_a[2] and row_b[2] and row_c[2]:
                        continue
                    first = "A" if not row_a[2] else ("B" if not row_b[2] else "C")
                    nontrivial += 1
                    first_off_diagonal[first] += 1
                    upper_x = glue((row_a[0], row_b[0], row_c[0]), maps)
                    upper_y = glue((row_a[1], row_b[1], row_c[1]), maps)
                    assert upper_x != upper_y
                    if (upper_x & 1) or (upper_y & 1):
                        continue
                    apex_zero += 1
                    first_off_diagonal_apex[first] += 1
                    if symmetric(upper_x, sigma8) != symmetric(upper_y, sigma8):
                        continue
                    upper_colour += 1
                    first_off_diagonal_colour[first] += 1
                    canonical = (min(upper_x, upper_y), max(upper_x, upper_y))
                    canonical_pairs.add(canonical)
                    canonical_by_first[first].add(canonical)

    assert nontrivial == 1_672
    assert apex_zero == 836 and upper_colour == 836
    assert len(canonical_pairs) == 418
    assert first_off_diagonal == {"A": 1_656, "B": 16, "C": 0}
    assert first_off_diagonal_apex == {"A": 828, "B": 8, "C": 0}
    assert first_off_diagonal_colour == {"A": 828, "B": 8, "C": 0}
    assert {key: len(value) for key, value in canonical_by_first.items()} == {
        "A": 414, "B": 4, "C": 0
    }

    stages = {
        "nontrivial_triangles": {"separated_pairs": 0, "cells": nontrivial},
        "both_apex_zero": {"separated_pairs": nontrivial - apex_zero, "cells": apex_zero},
        "upper_colour": {"separated_pairs": nontrivial - upper_colour, "cells": upper_colour},
        "unordered_pairs": {"separated_pairs": nontrivial - len(canonical_pairs), "cells": len(canonical_pairs)},
    }
    retention = carrier_tournament(stages, "retention")
    economy = carrier_tournament(stages, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(stages)) for j in range(i + 1, len(stages))
    )
    return {
        "R7_rows": relation_rows,
        "R7_off_diagonal": off_diagonal,
        "nontrivial_ordered_triangles": nontrivial,
        "both_apex_zero": apex_zero,
        "upper_colour_equal": upper_colour,
        "canonical_unordered_pairs": len(canonical_pairs),
        "first_off_diagonal": first_off_diagonal,
        "first_off_diagonal_apex": first_off_diagonal_apex,
        "first_off_diagonal_colour": first_off_diagonal_colour,
        "canonical_by_first": {key: len(value) for key, value in canonical_by_first.items()},
        "tournament_analysis": {
            "vertices": list(stages),
            "pairwise_observable": "number of raw nontrivial triangles rejected by the filter carrier",
            "switches": ["rejection retention", "rejection per surviving cell"],
            "retention": retention,
            "economy": economy,
            "edge_flips": flips,
        },
    }


def render(result: dict) -> str:
    ta = result["tournament_analysis"]
    return "\n".join(
        [
            "THM-818 ORDERED CECH-JOIN REGRESSION R7 -> N8",
            "=" * 72,
            f"R7 rows={result['R7_rows']} off-diagonal={result['R7_off_diagonal']}",
            f"nontrivial ordered triangles={result['nontrivial_ordered_triangles']}",
            f"both-apex-zero={result['both_apex_zero']}",
            f"upper-colour-equal={result['upper_colour_equal']}",
            f"canonical unordered pairs={result['canonical_unordered_pairs']}",
            f"first-off-diagonal A/B/C={result['first_off_diagonal']}",
            f"  after apex={result['first_off_diagonal_apex']}",
            f"  after colour={result['first_off_diagonal_colour']}",
            f"  canonical={result['canonical_by_first']}",
            "",
            "TOURNAMENT ANALYSIS (filter carriers as planning vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; flips={ta['edge_flips']}",
            f"  retention scores={ta['retention']['score_hist']} C3={ta['retention']['directed_3cycles']} "
            f"SCC={ta['retention']['scc_sizes']} Hpaths={ta['retention']['hamiltonian_paths']}",
            "",
            "Verdict: exact replay recovers THM-809's 418 base collisions.",
        ]
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, default=ATLAS)
    args = parser.parse_args()
    print(render(audit(args.atlas)), end="")


if __name__ == "__main__":
    main()
