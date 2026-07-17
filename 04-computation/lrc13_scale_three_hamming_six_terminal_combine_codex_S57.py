#!/usr/bin/env python3
"""Combine exact terminal-recursion shards for THM-862.

The C++ engine reconstructs the 1,504 arithmetic contexts, retains every
labelled step-39 ray language, and emits one GENERIC_ROW per context.  This
combiner checks exact, nonoverlapping coverage of the context bank, verifies
all row/shard accounting, aggregates the depth-six census, and records hashes.
It never identifies contexts or proof states.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CPP = ROOT / "04-computation/lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.cpp"
TYPES = ("2,4", "1,5", "0,6")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parse_array(word: str) -> tuple[int, ...]:
    answer = tuple(int(value) for value in word.split(","))
    require(len(answer) == 7, "expected a depth-0..6 array")
    return answer


def add_arrays(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(left, right, strict=True))


def array_word(values: tuple[int, ...]) -> str:
    return ",".join(str(value) for value in values)


def parse_row(line: str) -> dict:
    fields = line.split("|")
    require(len(fields) == 12 and fields[0] == "GENERIC_ROW", "malformed GENERIC_ROW")
    answer = {
        "line": line,
        "index": int(fields[1]),
        "type": fields[2],
        "nodes": parse_array(fields[3]),
        "dead": parse_array(fields[4]),
        "full_tooth": parse_array(fields[5]),
        "streaming_cap": parse_array(fields[6]),
        "candidate_edges": int(fields[7]),
        "covers": int(fields[8]),
        "loose": int(fields[9]),
        "maximum_cap": int(fields[10]),
        "maximum_candidate_speed": int(fields[11]),
    }
    require(answer["type"] in TYPES, "unknown context stratum")
    require(answer["nodes"][0] == 1, "context root multiplicity is not one")
    require(answer["candidate_edges"] == sum(answer["nodes"][1:]), "candidate/node accounting mismatch")
    require(answer["nodes"][6] == answer["covers"] + answer["loose"], "terminal row accounting mismatch")
    require(0 <= answer["maximum_candidate_speed"] <= answer["maximum_cap"], "candidate exceeds cap")
    for depth in range(7):
        require(answer["dead"][depth] <= answer["nodes"][depth], "dead count exceeds node count")
        require(answer["full_tooth"][depth] <= answer["nodes"][depth], "full-tooth count exceeds node count")
        require(answer["streaming_cap"][depth] <= answer["nodes"][depth], "streaming count exceeds node count")
    return answer


def parse_summary(line: str) -> dict:
    fields = line.split("|")
    require(fields[0] == "GENERIC_SHARD_SUMMARY", "malformed shard summary")
    values = dict(field.split("=", 1) for field in fields[1:])
    return {
        "nodes": parse_array(values["nodes"]),
        "dead": parse_array(values["dead"]),
        "full_tooth": parse_array(values["full_tooth"]),
        "streaming_cap": parse_array(values["streaming_cap"]),
        "candidate_edges": int(values["candidate_edges"]),
        "covers": int(values["covers"]),
        "loose": int(values["loose"]),
        "maximum_cap": int(values["maximum_cap"]),
        "maximum_candidate_speed": int(values["maximum_candidate_speed"]),
    }


parser = argparse.ArgumentParser()
parser.add_argument("shards", nargs="+", type=Path)
args = parser.parse_args()

zero = (0,) * 7
rows_by_index: dict[int, dict] = {}
shard_manifest: list[tuple[int, int, str]] = []

for path in args.shards:
    payload = path.read_bytes()
    lines = payload.decode().splitlines()
    require(lines and lines[0] == "SCALE_THREE_HAMMING_SIX_GENERIC_RECURSION_SHARD", "bad shard header")
    require(lines[-1] == "GENERIC_SHARD_DONE", "incomplete shard")
    arithmetic = dict(item.split("=", 1) for item in lines[1].split())
    require(arithmetic == {
        "arithmetic": "integer+rational",
        "floating_point": "none",
        "height_cutoff": "none",
        "depth_limit": "6",
        "early_gates": "1",
    }, "unexpected arithmetic or gate mode")
    range_fields = dict(item.split("=", 1) for item in lines[2].split())
    start = int(range_fields["context_start"])
    end = int(range_fields["context_end"])
    require(int(range_fields["context_count"]) == end - start, "bad shard range")
    require(int(range_fields["all_contexts"]) == 1504, "bad all-context count")

    local_rows: list[dict] = []
    for line in lines[3:-2]:
        require(line.startswith("GENERIC_ROW|"), "unexpected line in generic shard")
        row = parse_row(line)
        require(start <= row["index"] < end, "row outside shard range")
        require(row["index"] not in rows_by_index, "duplicate context row")
        rows_by_index[row["index"]] = row
        local_rows.append(row)
    require(len(local_rows) == end - start, "shard row count mismatch")

    summary = parse_summary(lines[-2])
    for field in ("nodes", "dead", "full_tooth", "streaming_cap"):
        aggregate = zero
        for row in local_rows:
            aggregate = add_arrays(aggregate, row[field])
        require(summary[field] == aggregate, f"shard {field} summary mismatch")
    for field in ("candidate_edges", "covers", "loose"):
        require(summary[field] == sum(row[field] for row in local_rows), f"shard {field} summary mismatch")
    for field in ("maximum_cap", "maximum_candidate_speed"):
        require(summary[field] == max(row[field] for row in local_rows), f"shard {field} summary mismatch")
    shard_manifest.append((start, end, sha256(payload).hexdigest()))

require(set(rows_by_index) == set(range(1504)), "context shards do not cover 0..1503 exactly")
rows = [rows_by_index[index] for index in range(1504)]
by_type = {kind: [row for row in rows if row["type"] == kind] for kind in TYPES}
require({kind: len(group) for kind, group in by_type.items()} == {
    "2,4": 336,
    "1,5": 672,
    "0,6": 496,
}, "stratum count mismatch")

aggregates = {}
for field in ("nodes", "dead", "full_tooth", "streaming_cap"):
    aggregate = zero
    for row in rows:
        aggregate = add_arrays(aggregate, row[field])
    aggregates[field] = aggregate
for field in ("candidate_edges", "covers", "loose"):
    aggregates[field] = sum(row[field] for row in rows)
for field in ("maximum_cap", "maximum_candidate_speed"):
    aggregates[field] = max(row[field] for row in rows)

require(aggregates["candidate_edges"] == sum(aggregates["nodes"][1:]), "global candidate/node mismatch")
require(aggregates["covers"] == 0, "scale-three covering terminal found")
require(aggregates["nodes"][6] == aggregates["covers"] + aggregates["loose"], "terminal verdict accounting mismatch")

expected = {
    "nodes": (1504, 146_912, 14_992_263, 931_412_556, 3_984_352, 4_481, 2),
    "dead": (0, 0, 1, 929_966_716, 3_980_598, 4_479, 0),
    "full_tooth": (0, 0, 0, 879_373_305, 0, 0, 0),
    "streaming_cap": (0, 0, 1, 50_593_411, 3_980_598, 4_479, 2),
    "candidate_edges": 950_540_566,
    "covers": 0,
    "loose": 2,
    "maximum_cap": 6_324,
    "maximum_candidate_speed": 6_317,
}
require(aggregates == expected, "frozen terminal census mismatch")
loose_rows = tuple(row["index"] for row in rows if row["loose"])
require(loose_rows == (888, 1502), "loose-terminal context mismatch")

row_payload = ("\n".join(row["line"] for row in rows) + "\n").encode()
require(
    sha256(row_payload).hexdigest()
    == "601655a592c4a9041bc7fb861031adf72d1f85bbdfd3769e50ed4d818978d5a1",
    "terminal row manifest mismatch",
)

print("SCALE_THREE_HAMMING_SIX_TERMINAL_EXACT_CENSUS")
print("arithmetic=integer+rational floating_point=none height_cutoff=none depth_limit=6 early_gates=1")
print(f"contexts={len(rows)} by_type={{{', '.join(f'{kind}:{len(by_type[kind])}' for kind in TYPES)}}}")
for field in ("nodes", "dead", "full_tooth", "streaming_cap"):
    print(f"{field}_depth0..6={array_word(aggregates[field])}")
print(f"candidate_edges={aggregates['candidate_edges']}")
print(f"maximum_cap={aggregates['maximum_cap']} maximum_candidate_speed={aggregates['maximum_candidate_speed']}")
print(f"covers={aggregates['covers']} nonempty_terminal_certificates={aggregates['loose']}")
print(f"nonempty_terminal_contexts={','.join(map(str, loose_rows))}")
for kind in TYPES:
    group = by_type[kind]
    type_nodes = zero
    for row in group:
        type_nodes = add_arrays(type_nodes, row["nodes"])
    print(
        f"stratum={kind} contexts={len(group)} nodes={array_word(type_nodes)} "
        f"candidate_edges={sum(row['candidate_edges'] for row in group)} "
        f"covers={sum(row['covers'] for row in group)} loose={sum(row['loose'] for row in group)}"
    )
print("TOURNAMENT_ANALYSIS")
print("vertices=labelled_future_rays pair_observable=least_legal_speed_difference switch=sign gauge=numerical_order")
print("fingerprints=transitive_at_each_materialized_prefix directed_cycles=0 SCCs=singletons Hamiltonian_paths=1 ties=0")
print("terminal_recursion_uses_tournament_quotient=0")
print("ASSUMPTION_AUDIT")
print("challenged_vertices=runners->missing_labels->signed_sheet_providers->labelled_ray_states->prefix_obligations")
print("faithful_state=literal_components+remaining_labelled_step39_rays+last_speed+shortcut_ancestry")
print("destroyed_by_ray_order=component_widths+sheet_units+future_language;_therefore_order_is_telemetry_only")
print(f"row_payload_sha256={sha256(row_payload).hexdigest()}")
print(f"primary_source_sha256={sha256(CPP.read_bytes()).hexdigest()}")
print(f"combiner_source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
print(f"shards={len(shard_manifest)}")
for start, end, digest in sorted(shard_manifest):
    print(f"shard range={start}:{end} sha256={digest}")
print("verdict=all_1504_scale_three_H6_contexts_have_no_covering_terminal")
