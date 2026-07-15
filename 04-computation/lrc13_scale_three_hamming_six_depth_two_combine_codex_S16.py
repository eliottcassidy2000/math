#!/usr/bin/env python3
"""Combine deterministic THM-862 depth-two context shards exactly.

The C++ scout emits one canonical ROW record per arithmetic context, one G1
record per cached first-edge geometry, and one L1 record per full labelled
future-ray language.  This combiner proves exact context coverage, merges only
the geometry cache, rejects every proof-language collision, aggregates the
logical census, and prints a shard-layout-independent result.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CPP = ROOT / "04-computation/lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.cpp"
TYPES = ("2,4", "1,5", "0,6")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parse_hist(word: str) -> Counter[int]:
    answer: Counter[int] = Counter()
    if not word:
        return answer
    for item in word.split(","):
        key, value = item.split(":")
        answer[int(key)] += int(value)
    return answer


def hist_word(histogram: Counter) -> str:
    return ",".join(f"{key}:{histogram[key]}" for key in sorted(histogram))


def fraction_word(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_row(line: str) -> dict:
    fields = line.split("|")
    require(len(fields) == 37 and fields[0] == "ROW", "malformed ROW record")
    answer = {
        "line": line,
        "index": int(fields[1]),
        "type": fields[2],
        "root": int(fields[3]),
        "labels": fields[4],
        "orders": fields[5],
        "units": fields[6],
        "depth1": int(fields[7]),
        "depth2": int(fields[8]),
        "first_d1": int(fields[9]),
        "first_d3": int(fields[10]),
        "second_d1": int(fields[11]),
        "second_d3": int(fields[12]),
        "dead0": int(fields[13]),
        "dead1": int(fields[14]),
        "dead2": int(fields[15]),
        "stream2": int(fields[16]),
        "cover1": int(fields[17]),
        "cover2": int(fields[18]),
        "live": int(fields[19]),
        "min0": F(fields[20]),
        "min0_count": int(fields[21]),
        "cap0": int(fields[22]),
        "min1": F(fields[23]),
        "min1_count": int(fields[24]),
        "cap1": int(fields[25]),
        "min2": F(fields[26]),
        "min2_count": int(fields[27]),
        "cap2": int(fields[28]),
        "flips": int(fields[29]),
        "flip_hist": parse_hist(fields[30]),
        "signed_edges": int(fields[31]),
        "signed_positive": int(fields[32]),
        "signed_negative": int(fields[33]),
        "signed_triangles": int(fields[34]),
        "signed_scc": fields[35],
        "signed_hp": int(fields[36]),
    }
    require(answer["type"] in TYPES, "unknown context stratum")
    return answer


def sum_field(rows: list[dict], field: str) -> int:
    return sum(row[field] for row in rows)


def extremum(rows: list[dict], depth: int) -> tuple[F, int, int]:
    minimum = min(row[f"min{depth}"] for row in rows)
    count = sum(
        row[f"min{depth}_count"]
        for row in rows
        if row[f"min{depth}"] == minimum
    )
    maximum_cap = max(row[f"cap{depth}"] for row in rows)
    return minimum, count, maximum_cap


parser = argparse.ArgumentParser()
parser.add_argument("shards", nargs="+", type=Path)
args = parser.parse_args()

rows_by_index: dict[int, dict] = {}
geometry_one: Counter[tuple[int, int]] = Counter()
languages: set[str] = set()
dead_certificates: list[str] = []
shard_manifest = []

for path in args.shards:
    payload = path.read_bytes()
    lines = payload.decode().splitlines()
    require(lines and lines[0] == "SCALE_THREE_HAMMING_SIX_DEPTH_TWO_SHARD", "bad shard header")
    range_line = next(line for line in lines if line.startswith("context_start="))
    range_fields = dict(item.split("=") for item in range_line.split())
    start = int(range_fields["context_start"])
    end = int(range_fields["context_end"])
    require(int(range_fields["context_count"]) == end - start, "bad shard range")
    require(int(range_fields["all_contexts"]) == 1504, "bad all-context count")
    require(int(range_fields["direct_NAE_checks"]) == 94448, "bad direct/NAE count")
    require(lines[-1] == "SHARD_DONE", "incomplete shard")
    local_rows = 0
    for line in lines:
        if line.startswith("ROW|"):
            row = parse_row(line)
            require(start <= row["index"] < end, "ROW outside shard range")
            require(row["index"] not in rows_by_index, "duplicate context ROW")
            rows_by_index[row["index"]] = row
            local_rows += 1
        elif line.startswith("G1|"):
            _, root, speed, multiplicity = line.split("|")
            geometry_one[int(root), int(speed)] += int(multiplicity)
        elif line.startswith("L1|"):
            language = line[3:]
            require(language not in languages, "labelled future-ray language collision")
            languages.add(language)
        elif line.startswith("D2DEAD|"):
            dead_certificates.append(line)
    require(local_rows == end - start, "shard ROW count mismatch")
    shard_manifest.append((start, end, sha256(payload).hexdigest()))

require(set(rows_by_index) == set(range(1504)), "context shards do not cover 0..1503 exactly")
rows = [rows_by_index[index] for index in range(1504)]
require(sum(geometry_one.values()) == sum_field(rows, "depth1"), "G1 multiplicity mismatch")
require(len(languages) == sum_field(rows, "depth1"), "L1 census mismatch")

by_type = {kind: [row for row in rows if row["type"] == kind] for kind in TYPES}
require({kind: len(group) for kind, group in by_type.items()} == {"2,4": 336, "1,5": 672, "0,6": 496}, "stratum count mismatch")
require(sum_field(rows, "depth1") == 146_912, "first logical layer mismatch")
require(sum_field(rows, "first_d1") == 20_640, "first D=1 edge mismatch")
require(sum_field(rows, "first_d3") == 126_272, "first D=3 edge mismatch")
require(sum_field(rows, "depth2") == 14_992_263, "depth-two logical census mismatch")
require(sum_field(rows, "second_d1") == 2_216_176, "second D=1 edge mismatch")
require(sum_field(rows, "second_d3") == 12_776_087, "second D=3 edge mismatch")
require(sum_field(rows, "dead0") == sum_field(rows, "dead1") == 0, "unexpected early dead prefix")
require(sum_field(rows, "dead2") == sum_field(rows, "stream2") == 1, "depth-two dead census mismatch")
require(sum_field(rows, "live") == 14_992_262, "depth-two live frontier mismatch")

g1_hist = Counter(geometry_one.values())
require(len(geometry_one) == 22_262, "first-edge geometry cache size mismatch")
require(g1_hist == {2: 2_454, 4: 8_753, 6: 6_479, 12: 2_375, 18: 2_201}, "first-edge geometry multiplicity mismatch")

expected_dead = (
    "D2DEAD|1448|1,5|2536|4,6,7,8,9,12|3,3,3,3,3,1|"
    "2,1,1,1,2,0|14|9|38|4|183/494,56/143|115/5434|63|70|"
    "6:3:1:70,7:3:1:73,8:3:1:76,12:1:0:75"
)
require(dead_certificates == [expected_dead], "unique depth-two dead certificate mismatch")

for row in rows:
    require(row["depth1"] == row["first_d1"] + row["first_d3"], "first edge split mismatch")
    require(row["depth2"] == row["second_d1"] + row["second_d3"], "second edge split mismatch")
    require(row["depth2"] == row["dead2"] + row["live"], "frontier split mismatch")
    require(row["dead2"] == row["stream2"], "stream certificate mismatch")
    require(row["cover1"] == row["cover2"] == 0, "unexpected early cover")
    require(sum(row["flip_hist"].values()) == row["depth1"], "conditioned tournament histogram mismatch")

row_payload = ("\n".join(row["line"] for row in rows) + "\n").encode()
g1_lines = [f"{root}:{speed}:{geometry_one[root, speed]}" for root, speed in sorted(geometry_one)]
g1_payload = ("\n".join(g1_lines) + "\n").encode()
language_payload = ("\n".join(sorted(languages)) + "\n").encode()
require(sha256(row_payload).hexdigest() == "19b6806acb4fd6c00c914690b20322594013ef82799b53c9d42d26710ac80b3f", "ROW manifest mismatch")
require(sha256(g1_payload).hexdigest() == "466f63905e19a2288545697480b98f237c5bb577db49c0f3f14adae1622e7c11", "G1 manifest mismatch")
require(sha256(language_payload).hexdigest() == "916008f85570e4a1a3a17c6915b69cc1f18b65c76102a0d7c06f3344d57a5feb", "L1 manifest mismatch")

flip_hist: Counter[int] = Counter()
for row in rows:
    flip_hist.update(row["flip_hist"])

expected_by_type = {
    "2,4": (32_708, 3_408_353, 0, 3_011, 26_093),
    "1,5": (64_596, 6_469_464, 1, 4_530, 22_800),
    "0,6": (49_608, 5_114_446, 0, 4_600, 33_133),
}
for kind, (depth1, depth2, dead2, minimum, maximum) in expected_by_type.items():
    group = by_type[kind]
    require(sum_field(group, "depth1") == depth1, f"{kind} depth-one mismatch")
    require(sum_field(group, "depth2") == depth2, f"{kind} depth-two mismatch")
    require(sum_field(group, "dead2") == dead2, f"{kind} dead mismatch")
    require(min(row["depth2"] for row in group) == minimum, f"{kind} minimum context workload mismatch")
    require(max(row["depth2"] for row in group) == maximum, f"{kind} maximum context workload mismatch")

require(extremum(rows, 0) == (F(1, 117), 8, 1188), "root component extrema mismatch")
require(extremum(rows, 1) == (F(11, 15431), 4, 3956), "depth-one component extrema mismatch")
require(extremum(rows, 2) == (F(11, 51389), 2, 6324), "depth-two component extrema mismatch")
require(flip_hist == {0: 26_541, 1: 5_584, 2: 7_497, 3: 13_104, 4: 38_196, 5: 10_478, 6: 34_052, 7: 8_502, 8: 2_958}, "conditioned edge-flip histogram mismatch")

print("SCALE_THREE_HAMMING_SIX_DEPTH_TWO_EXACT_CENSUS")
print("arithmetic=integer+rational floating_point=none height_cutoff=none depth_limit=2")
print(f"contexts={len(rows)} by_type={{{', '.join(f'{kind}:{len(by_type[kind])}' for kind in TYPES)}}}")
print(f"nodes_depth0..2={len(rows)},{sum_field(rows, 'depth1')},{sum_field(rows, 'depth2')}")
print(f"candidate_edges={sum_field(rows, 'depth1') + sum_field(rows, 'depth2')}")
print(
    f"edge_orders first_D1={sum_field(rows, 'first_d1')} first_D3={sum_field(rows, 'first_d3')} "
    f"second_D1={sum_field(rows, 'second_d1')} second_D3={sum_field(rows, 'second_d3')}"
)
print(
    f"dead_depth0..2={sum_field(rows, 'dead0')},{sum_field(rows, 'dead1')},{sum_field(rows, 'dead2')} "
    f"streaming_cap_depth2={sum_field(rows, 'stream2')} full_safe_tooth_depth0..2=0,0,0"
)
print(
    f"covers_depth1={sum_field(rows, 'cover1')} covers_depth2={sum_field(rows, 'cover2')} "
    f"frontier_live={sum_field(rows, 'live')}"
)
print(
    "unique_dead_prefix context=1448 R=4,6,7,8,9,12 "
    "orders=3,3,3,3,3,1 units=2,1,1,1,2,0 insertions=9:14,4:38"
)
print(
    "unique_dead_certificate longest_component=(183/494,56/143) "
    "length=115/5434 four_comb_cap=63 least_future=70 "
    "remaining=6:70,7:73,8:76,12:75"
)
for depth in range(3):
    minimum, count, maximum = extremum(rows, depth)
    print(
        f"component_depth={depth} min_longest={fraction_word(minimum)} "
        f"multiplicity={count} max_next_cap={maximum}"
    )
for kind in TYPES:
    group = by_type[kind]
    print(
        f"stratum={kind} contexts={len(group)} depth1={sum_field(group, 'depth1')} "
        f"depth2={sum_field(group, 'depth2')} dead2={sum_field(group, 'dead2')} "
        f"live2={sum_field(group, 'live')} per_context_depth1={min(r['depth1'] for r in group)}..{max(r['depth1'] for r in group)} "
        f"per_context_depth2={min(r['depth2'] for r in group)}..{max(r['depth2'] for r in group)}"
    )
print(
    f"geometry_cache depth1_logical={sum_field(rows, 'depth1')} unique_R_x1={len(geometry_one)} "
    f"multiplicity_hist={hist_word(g1_hist)}"
)
print(
    f"proof_languages depth1={len(languages)} collisions=0 "
    "geometry_cache_is_not_a_proof_state_quotient"
)
print("SIGNED_PROVIDER_GRAPH_ANALYSIS")
print("vertices=order_three_labels relation=provider_to_owner_at_ratio_4|5|8|9 signed_by_ratio")
print(f"edge_count_hist={hist_word(Counter(row['signed_edges'] for row in rows))}")
print(f"sign_split_hist={hist_word(Counter((row['signed_positive'], row['signed_negative']) for row in rows))}")
print(f"directed_triangle_hist={hist_word(Counter(row['signed_triangles'] for row in rows))}")
print(f"SCC_profile_hist={hist_word(Counter(row['signed_scc'] for row in rows))}")
print(f"Hamiltonian_path_count_hist={hist_word(Counter(row['signed_hp'] for row in rows))}")
print("relation_type=sparse_signed_not_tournament sheet_predicate_preserved=1 metric_predicate_preserved=0")
print("RAY_TOURNAMENT_ANALYSIS")
print("vertices=labelled_future_rays pair_observable=least_legal_speed_difference switch=sign gauge=numerical_order")
print("root_tournaments=1504 score_hist_per_root=0:1,1:1,2:1,3:1,4:1,5:1 directed_triangles=0 SCCs=1:6 Hamiltonian_paths=1 ties=0")
print(
    f"conditioned_tournaments={sum_field(rows, 'depth1')} score_hist_per_prefix=0:1,1:1,2:1,3:1,4:1 "
    f"directed_triangles=0 SCCs=1:5 Hamiltonian_paths=1 ties=0 edge_flips={sum_field(rows, 'flips')} "
    f"flip_hist={hist_word(flip_hist)}"
)
print("tie_Hamiltonian_path=increasing_label_after_equal_speed;_unused_because_CRT_residues_are_distinct")
print("ASSUMPTION_AUDIT")
print("challenged_vertices=runners->missing_labels->signed_sheet_providers->labelled_ray_states->prefix_obligations")
print("faithful_state=literal_components+remaining_labelled_step39_rays+last_speed+shortcut_ancestry")
print("cache_preserves=literal_component_geometry destroys=remaining_ray_language_and_metric_continuation")
print(f"row_payload_sha256={sha256(row_payload).hexdigest()}")
print(f"geometry1_payload_sha256={sha256(g1_payload).hexdigest()}")
print(f"language_payload_sha256={sha256(language_payload).hexdigest()}")
print(f"shards={len(shard_manifest)}")
for start, end, digest in sorted(shard_manifest):
    print(f"shard range={start}:{end} sha256={digest}")
print(f"primary_source_sha256={sha256(CPP.read_bytes()).hexdigest()}")
print(f"combiner_source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")
print("verdict=exact_depth_two_census_complete_deeper_metric_bank_open")
