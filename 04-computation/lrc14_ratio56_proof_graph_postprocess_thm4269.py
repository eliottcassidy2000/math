#!/usr/bin/env python3
"""Semantic proof-graph novelty for the 4:5 then 5:6 outsider rays.

Input is the independently frozen THM-4266 candidate residual.  No edge list
for either adjacent ratio is copied: both are selected by their exact ratio,
primitive scale, and theorem bridge start.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


THM4266_RESIDUAL = Path(
    "05-knowledge/results/lrc14_three_round_learned_carrier_thm4266/"
    "thm4266_final_residual.csv"
)
THM4267_AUDIT = Path(
    "05-knowledge/results/"
    "lrc14_ratio45_common_deck_proof_graph_independent_audit_thm4267.out"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_add(value: int, word: int) -> int:
    for shift in range(0, 64, 8):
        value ^= (word >> shift) & 0xFF
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


def edge_fnv(edges: list[tuple[int, int]]) -> int:
    value = 0xCBF29CE484222325
    for left, right in edges:
        value = fnv_add(value, left)
        value = fnv_add(value, right)
    return value


def edge_bytes(edges: list[tuple[int, int]]) -> bytes:
    return b"".join(f"{left},{right}\n".encode("ascii") for left, right in edges)


def edge_sha(edges: list[tuple[int, int]]) -> str:
    return hashlib.sha256(edge_bytes(edges)).hexdigest()


def load_edges(path: Path) -> list[tuple[int, int]]:
    edges = []
    for line in path.read_text(encoding="ascii").splitlines():
        left_text, right_text = line.split(",")
        edges.append((int(left_text), int(right_text)))
    require(edges == sorted(edges), "input residual order changed")
    require(len(edges) == len(set(edges)), "input residual has duplicates")
    return edges


parser = argparse.ArgumentParser()
parser.add_argument("--repo", type=Path, default=Path("."))
args = parser.parse_args()
repo = args.repo.resolve()

post4266 = load_edges(repo / THM4266_RESIDUAL)
require(len(post4266) == 177_585, "THM-4266 residual count changed")
require(edge_fnv(post4266) == 0x6CE05D05EB01DAED,
        "THM-4266 residual FNV changed")
require(edge_sha(post4266) ==
        "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9",
        "THM-4266 residual SHA changed")

canonical_4267 = (repo / THM4267_AUDIT).read_bytes()
require(hashlib.sha256(canonical_4267).hexdigest() ==
        "211740c11a47d2a5f4da757a266b31523794a0570b362cb686810fcf08dd5001",
        "canonical THM-4267 proof-graph audit bytes changed")
require(
    b"NOVEL_AFTER_THM4266 COUNT 63 FNV 2280174481469c0d "
    b"SHA256 a63b3309dc30775c82ab7533472f4523ae4bda69525bc262ac06e749011dde23\n"
    in canonical_4267,
    "canonical THM-4267 novel ledger changed",
)
require(
    b"FINAL_AFTER_THM4267 COUNT 177522 FNV 33142f955cc93379 "
    b"SHA256 d277aebe296153ead14a77207ea1499c961c8b06796b7e62f324e34f7a9ef087\n"
    in canonical_4267,
    "canonical THM-4267 final ledger changed",
)

ratio45 = []
scales45 = []
for left, right in post4266:
    if left % 4 != 0:
        continue
    scale = left // 4
    if scale >= 73 and right == 5 * scale:
        ratio45.append((left, right))
        scales45.append(scale)

require(len(ratio45) == 63, "post-THM-4266 4:5 novelty changed")
require(scales45 == [*range(73, 132), 133, 134, 135, 137],
        "post-THM-4266 4:5 scale ledger changed")

ratio45_set = set(ratio45)
post4267 = [edge for edge in post4266 if edge not in ratio45_set]
require(len(post4267) == 177_522, "post-THM-4267 count changed")

ratio56 = []
scales56 = []
for left, right in post4267:
    if left % 5 != 0:
        continue
    scale = left // 5
    if scale >= 59 and right == 6 * scale:
        ratio56.append((left, right))
        scales56.append(scale)

require(len(ratio56) == 53, "post-THM-4267 5:6 novelty changed")
ratio56_set = set(ratio56)
require(not (ratio45_set & ratio56_set), "4:5/5:6 novelty overlap")
post4269 = [edge for edge in post4267 if edge not in ratio56_set]
require(len(post4269) == 177_469, "post-THM-4269 count changed")

maximum = max(right for _, right in post4269)
top = [edge for edge in post4269 if edge[1] == maximum]
require(maximum == 688 and top == [(520, 688)],
        "post-THM-4269 top layer changed")

print("LRC14_RATIO_5_6_PROOF_GRAPH_POSTPROCESS")
print(
    f"POST_THM4266 COUNT {len(post4266)} FNV {edge_fnv(post4266):016x} "
    f"SHA256 {edge_sha(post4266)}"
)
print(
    f"THM4267_RATIO_4_5_NOVEL COUNT {len(ratio45)} "
    f"FNV {edge_fnv(ratio45):016x} SHA256 {edge_sha(ratio45)}"
)
print("THM4267_SCALES " + ",".join(map(str, scales45)))
print(
    f"POST_THM4267 COUNT {len(post4267)} FNV {edge_fnv(post4267):016x} "
    f"SHA256 {edge_sha(post4267)}"
)
print(
    f"THM4269_RATIO_5_6_NOVEL COUNT {len(ratio56)} "
    f"FNV {edge_fnv(ratio56):016x} SHA256 {edge_sha(ratio56)}"
)
print("THM4269_SCALES " + ",".join(map(str, scales56)))
print(
    f"POST_THM4269 COUNT {len(post4269)} FNV {edge_fnv(post4269):016x} "
    f"SHA256 {edge_sha(post4269)}"
)
print(
    f"TOP_LAYER ENDPOINT {maximum} COUNT {len(top)} EDGES "
    + " ".join(f"{left},{right}" for left, right in top)
)
print("VERDICT EXACT_SEMANTIC_NOVELTY_CONDITIONAL_ON_THM4266_4267 LRC14_OPEN")
