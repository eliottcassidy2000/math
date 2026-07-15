#!/usr/bin/env python3
"""Exact R8 census for THM-818's collision-driven n=9 Cech join.

This intentionally reuses the certified 6,880-class n=8 classifier from
metagraph_flow_n8_check_opus_S305.py.  It persists only the exact equality-
relation census and resource preflight, not a 2^21 node-rank atlas.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import runpy
from collections import Counter
from pathlib import Path

import numpy as np

from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament


LEGACY = Path("04-computation/metagraph_flow_n8_check_opus_S305.py")
OUTPUT = Path("05-knowledge/results/mobius_cech_n9_relation_preflight_codex_S13.out")
JSON_OUTPUT = Path("05-knowledge/results/mobius_cech_n9_relation_preflight_codex_S13.json")


def audit(legacy: Path) -> dict:
    with contextlib.redirect_stdout(io.StringIO()):
        namespace = runpy.run_path(str(legacy))
    bits = namespace["bits_all"]
    cls = namespace["cls_of"]
    reflection = namespace["refl"]
    m = int(namespace["m"])
    tilings = int(namespace["N"])
    grid_symmetric = namespace["gs"]
    full = tilings - 1

    reflected = np.zeros(tilings, dtype=np.int64)
    for source_bit in range(m):
        reflected |= ((bits >> source_bit) & 1) << int(reflection[source_bit])

    converse_pairs = np.stack(
        (np.minimum(cls, cls[reflected]), np.maximum(cls, cls[reflected])), axis=1
    )
    unique_nodes, node_rank = np.unique(converse_pairs, axis=0, return_inverse=True)
    merged_nodes = len(unique_nodes)
    other = bits ^ full
    packed_h = (
        ((node_rank.astype(np.int64) * merged_nodes + node_rank[other]) << 1)
        | grid_symmetric.astype(np.int64)
    )
    keys, counts = np.unique(packed_h, return_counts=True)
    counts_object = counts.astype(object)
    relation_rows = int(np.dot(counts_object, counts_object))
    off_diagonal = relation_rows - tilings
    blue = (keys & 1) == 1
    blue_relation_rows = int(np.dot(counts[blue].astype(object), counts[blue].astype(object)))
    black_relation_rows = relation_rows - blue_relation_rows

    histogram = Counter(map(int, counts))
    memory = {
        "packed_relation_rows_bytes": 8 * relation_rows,
        "rows_plus_two_16byte_indexes": 40 * relation_rows,
        "rows_plus_three_16byte_indexes": 56 * relation_rows,
    }
    assert len(set(map(int, cls))) == 6880
    assert merged_nodes == 3528
    assert len(keys) == 876512 and int(counts.max()) == 26
    assert relation_rows == 5997416 and off_diagonal == 3900264
    assert blue_relation_rows == 5216 and black_relation_rows == 5992200
    assert sum(multiplicity * cells for multiplicity, cells in histogram.items()) == tilings
    assert sum(multiplicity * multiplicity * cells for multiplicity, cells in histogram.items()) == relation_rows

    carrier_stats = {
        "fibre_count_table": {"separated_pairs": 1, "cells": len(keys)},
        "packed_relation_rows": {"separated_pairs": 2, "cells": relation_rows},
        "two_projection_indexes": {"separated_pairs": 3, "cells": memory["rows_plus_two_16byte_indexes"]},
        "three_projection_indexes": {"separated_pairs": 4, "cells": memory["rows_plus_three_16byte_indexes"]},
    }
    capability = carrier_tournament(carrier_stats, "retention")
    economy = carrier_tournament(carrier_stats, "economy")
    flips = sum(
        capability["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(carrier_stats)) for j in range(i + 1, len(carrier_stats))
    )

    return {
        "schema_version": 1,
        "theorem": "THM-818",
        "n8_tilings": tilings,
        "n8_unmerged_classes": 6880,
        "n8_merged_nodes": merged_nodes,
        "H8_support": len(keys),
        "H8_max_fibre": int(counts.max()),
        "H8_fibre_histogram": dict(sorted(histogram.items())),
        "R8_rows": relation_rows,
        "R8_off_diagonal_rows": off_diagonal,
        "R8_blue_rows": blue_relation_rows,
        "R8_black_rows": black_relation_rows,
        "R8_rows_per_tiling": relation_rows / tilings,
        "memory": memory,
        "single_machine_preflight": relation_rows <= 30_000_000,
        "n9_key_bits": {
            "six_n8_node_ranks": 72,
            "UABC": 4,
            "S2": 22,
            "total": 98,
            "S2_radix_product": 3_200_000,
        },
        "exact_reduction": "ordered equal-face-H8 pairs are role-specific overlap-compatible triangles in R8^3; lower-codec collisions then require the non-diagonal, both-apex-zero, upper-UABC, and full-S2 filters",
        "tournament_analysis": {
            "vertices": list(carrier_stats),
            "pairwise_observable": "number of declared Cech-join obligations executable by the representation",
            "switches": ["capability", "capability_per_log2_storage_cells"],
            "tie_hamiltonian_path": list(carrier_stats),
            "carrier_stats": carrier_stats,
            "capability": capability,
            "economy": economy,
            "edge_flips_between_gauges": flips,
        },
    }


def render(result: dict) -> str:
    lines = [
        "THM-818 N=9 LOWER-COLLISION RELATION PREFLIGHT",
        "=" * 78,
        f"n8 tilings/classes/merged nodes={result['n8_tilings']}/"
        f"{result['n8_unmerged_classes']}/{result['n8_merged_nodes']}",
        f"H8 support={result['H8_support']} max fibre={result['H8_max_fibre']}",
        f"H8 fibre histogram={result['H8_fibre_histogram']}",
        "",
        "EXACT EQUALITY RELATION",
        f"  R8 rows={result['R8_rows']} off-diagonal={result['R8_off_diagonal_rows']} "
        f"rows/tiling={result['R8_rows_per_tiling']:.12g}",
        f"  blue/black rows={result['R8_blue_rows']}/{result['R8_black_rows']}",
        f"  reduction: {result['exact_reduction']}",
        "",
        "MEMORY PREFLIGHT",
        f"  packed uint64 rows={result['memory']['packed_relation_rows_bytes']} bytes",
        f"  + two 16-byte indexes={result['memory']['rows_plus_two_16byte_indexes']} bytes",
        f"  + three 16-byte indexes={result['memory']['rows_plus_three_16byte_indexes']} bytes",
        f"  single-machine gate R8<=30m: {result['single_machine_preflight']}",
        "",
        "N=9 LOWER KEY",
        f"  bits node/UABC/S2/total={result['n9_key_bits']['six_n8_node_ranks']}/"
        f"{result['n9_key_bits']['UABC']}/{result['n9_key_bits']['S2']}/"
        f"{result['n9_key_bits']['total']}",
        f"  S2 radix product={result['n9_key_bits']['S2_radix_product']} (<2^22)",
    ]
    ta = result["tournament_analysis"]
    lines.extend(
        [
            "",
            "TOURNAMENT ANALYSIS (join representations as vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; flips={ta['edge_flips_between_gauges']}",
            f"  capability scores={ta['capability']['score_hist']} C3="
            f"{ta['capability']['directed_3cycles']} SCC={ta['capability']['scc_sizes']} "
            f"Hpaths={ta['capability']['hamiltonian_paths']}",
            "",
            "Verdict: execute the exact three-way Cech join; do not enumerate 2^27 lines.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy", type=Path, default=LEGACY)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()
    result = audit(args.legacy)
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
