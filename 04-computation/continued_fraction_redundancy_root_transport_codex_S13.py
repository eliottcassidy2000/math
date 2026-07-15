#!/usr/bin/env python3
"""Audit THM-808's centered-Christoffel redundancy-root cocycle.

The input is THM-778's exact eight-owner wall movie.  A gap between covered
walls is represented by its owner-count vector c.  The token update

    k_a -> k_a - c_a w_a^{-1}  (mod p)

implies an affine update of the unique duplicated residue.  The script checks
that law, identifies the five exact block types, and tests which quotients are
continuation-deterministic on the ten-wall skeleton.

Tournament Analysis uses preservation carriers as vertices, literal gap pairs
as the pairwise observable, retention/economy as switches, and increasing
carrier richness as the tie Hamiltonian path.  The runner labels are payload,
not the tournament vertices.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from tournament_tiling_metagraph_address_codex_S4 import carrier_tournament
from mobius_cech_metagraph_codec_codex_S12 import compact_partition


INPUT = Path("05-knowledge/results/lrc14_centered_christoffel_endpoint_skew_product_codex_S7.json")
OUTPUT = Path("05-knowledge/results/continued_fraction_redundancy_root_transport_codex_S13.out")
JSON_OUTPUT = Path("05-knowledge/results/continued_fraction_redundancy_root_transport_codex_S13.json")


def deterministic_transition(records: list[dict], state_fields: tuple[str, ...], target: str) -> dict:
    images: dict[tuple, set] = defaultdict(set)
    indices: dict[tuple, list[int]] = defaultdict(list)
    for record in records:
        key = tuple(record[field] for field in state_fields)
        images[key].add(record[target])
        indices[key].append(record["index"])
    conflicts = [
        {
            "state": repr(key),
            "indices": indices[key],
            "targets": sorted(values),
        }
        for key, values in images.items()
        if len(values) > 1
    ]
    return {
        "state_fields": list(state_fields),
        "target": target,
        "cells": len(images),
        "conflict_cells": len(conflicts),
        "conflicts": conflicts,
        "descends_to_function": not conflicts,
    }


def audit(data: dict, prime: int = 7) -> dict:
    movie = next(movie for movie in data["movies"] if len(movie["speeds"]) == prime + 1)
    speeds = tuple(movie["speeds"])
    walls = movie["covered_walls"]
    gaps = movie["covered_wall_gaps"]
    assert len(walls) == 10 and len(gaps) == 9

    chambers = movie["exact_chambers"]
    left_of = {chamber["right"]: chamber["position"] for chamber in chambers}
    right_of = {chamber["left"]: chamber["position"] for chamber in chambers}
    root_sum_failures = 0
    for wall in walls:
        x = wall["x"]
        root_sum_failures += sum(left_of[x]) % prime != wall["duplicate_sheet_before"]
        root_sum_failures += sum(right_of[x]) % prime != wall["duplicate_sheet_after"]
    assert root_sum_failures == 0

    records = []
    cocycle_failures = 0
    for index, gap in enumerate(gaps):
        source = walls[index]
        target = walls[index + 1]
        counts = tuple(gap["owner_counts"])
        coordinate_shifts = tuple(
            (-count * pow(speed, -1, prime)) % prime
            for count, speed in zip(counts, speeds)
        )
        cocycle_delta = sum(coordinate_shifts) % prime
        observed_delta = (
            target["duplicate_sheet_before"] - source["duplicate_sheet_after"]
        ) % prime
        cocycle_failures += cocycle_delta != observed_delta
        records.append(
            {
                "index": index,
                "wall_blocks": gap["wall_blocks"],
                "individual_events": gap["individual_events"],
                "owner_counts": counts,
                "owner_counts_mod_prime": tuple(c % prime for c in counts),
                "coordinate_shifts": coordinate_shifts,
                "cocycle_delta": cocycle_delta,
                "observed_delta": observed_delta,
                "source_mask": source["a267_mask"],
                "target_mask": target["a267_mask"],
                "source_root": source["duplicate_sheet_after"],
                "target_root": target["duplicate_sheet_before"],
            }
        )
    assert cocycle_failures == 0

    by_block: dict[tuple[int, ...], list[dict]] = defaultdict(list)
    for record in records:
        by_block[record["owner_counts"]].append(record)
    block_types = []
    block_action_failures = 0
    for rank, (counts, members) in enumerate(
        sorted(by_block.items(), key=lambda item: min(row["index"] for row in item[1]))
    ):
        deltas = sorted({row["cocycle_delta"] for row in members})
        block_action_failures += len(deltas) != 1
        block_types.append(
            {
                "rank": rank,
                "wall_blocks": members[0]["wall_blocks"],
                "individual_events": members[0]["individual_events"],
                "owner_counts": list(counts),
                "occurrences": [row["index"] for row in members],
                "root_translation_mod_prime": deltas,
                "coordinate_shifts": list(members[0]["coordinate_shifts"]),
            }
        )
    assert len(block_types) == 5 and block_action_failures == 0

    transition_audits = {
        "block_to_root": deterministic_transition(records, ("owner_counts",), "target_root"),
        "root_block_to_root": deterministic_transition(
            records, ("source_root", "owner_counts"), "target_root"
        ),
        "mask_block_to_mask": deterministic_transition(
            records, ("source_mask", "owner_counts"), "target_mask"
        ),
        "mask_root_block_to_mask": deterministic_transition(
            records, ("source_mask", "source_root", "owner_counts"), "target_mask"
        ),
    }
    assert transition_audits["root_block_to_root"]["descends_to_function"]
    assert not transition_audits["mask_block_to_mask"]["descends_to_function"]
    assert transition_audits["mask_root_block_to_mask"]["descends_to_function"]

    values = {
        "block_length": {row["index"]: row["wall_blocks"] for row in records},
        "owner_count_block": {row["index"]: row["owner_counts"] for row in records},
        "source_root": {row["index"]: row["source_root"] for row in records},
        "source_mask": {row["index"]: row["source_mask"] for row in records},
        "mask_plus_block": {
            row["index"]: (row["source_mask"], row["owner_counts"]) for row in records
        },
        "root_plus_block": {
            row["index"]: (row["source_root"], row["owner_counts"]) for row in records
        },
        "mask_root_block": {
            row["index"]: (row["source_mask"], row["source_root"], row["owner_counts"])
            for row in records
        },
        "exact_gap": {row["index"]: row["index"] for row in records},
    }
    stats = {name: compact_partition(partition) for name, partition in values.items()}
    retention = carrier_tournament(stats, "retention")
    economy = carrier_tournament(stats, "economy")
    flips = sum(
        retention["adjacency"][i][j] != economy["adjacency"][i][j]
        for i in range(len(values)) for j in range(i + 1, len(values))
    )

    return {
        "schema_version": 1,
        "theorem": "THM-808",
        "prime": prime,
        "speeds": list(speeds),
        "general_identity": "d'=d-sum_a c_a*w_a^(-1) (mod p)",
        "root_equals_token_sum_failures": root_sum_failures,
        "block_cocycle_failures": cocycle_failures,
        "block_types": block_types,
        "records": records,
        "transition_audits": transition_audits,
        "tournament_analysis": {
            "vertices": list(values),
            "pairwise_observable": "number of unordered covered-gap pairs separated by a carrier",
            "switches": ["retention", "retention_per_log2_cells"],
            "tie_hamiltonian_path": list(values),
            "carrier_stats": stats,
            "retention": retention,
            "economy": economy,
            "edge_flips_between_gauges": flips,
        },
    }


def render(result: dict) -> str:
    lines = [
        "THM-808 CENTERED-CHRISTOFFEL REDUNDANCY-ROOT TRANSPORT",
        "=" * 82,
        f"prime={result['prime']} speeds={tuple(result['speeds'])}",
        f"law: {result['general_identity']}",
        f"root-sum failures={result['root_equals_token_sum_failures']}; "
        f"block-cocycle failures={result['block_cocycle_failures']}",
        "",
        "FIVE EXACT FIRST-RETURN BLOCK TYPES",
    ]
    for block in result["block_types"]:
        lines.append(
            f"  B{block['rank']}: walls/events={block['wall_blocks']}/{block['individual_events']} "
            f"occurrences={tuple(block['occurrences'])} root translation="
            f"{tuple(block['root_translation_mod_prime'])} counts={tuple(block['owner_counts'])}"
        )
    lines.append("")
    lines.append("CONTINUATION DESCENT")
    for name, audit in result["transition_audits"].items():
        lines.append(
            f"  {name}: cells={audit['cells']} conflicts={audit['conflict_cells']} "
            f"function={audit['descends_to_function']}"
        )
        for conflict in audit["conflicts"]:
            lines.append(
                f"    witness indices={tuple(conflict['indices'])} targets={tuple(conflict['targets'])} "
                f"state={conflict['state']}"
            )
    ta = result["tournament_analysis"]
    lines.extend(
        [
            "",
            "TOURNAMENT ANALYSIS (preservation carriers as vertices)",
            f"  vertices={tuple(ta['vertices'])}",
            f"  observable={ta['pairwise_observable']}",
            f"  switches={tuple(ta['switches'])}; edge flips={ta['edge_flips_between_gauges']}",
            f"  retention scores={ta['retention']['score_hist']} "
            f"C3={ta['retention']['directed_3cycles']} SCC={ta['retention']['scc_sizes']} "
            f"Hpaths={ta['retention']['hamiltonian_paths']}",
            "",
            "The repeated three-wall block starts twice from mask 31115 but reaches",
            "different masks.  The source redundancy roots differ and restore the",
            "observed transition.  Thus the CF block acts on the labelled/root stalk,",
            "not on the metagraph mask quotient alone.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--json", type=Path, default=JSON_OUTPUT)
    args = parser.parse_args()
    result = audit(json.loads(args.input.read_text()))
    text = render(result)
    print(text, end="")
    args.output.write_text(text)
    args.json.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
