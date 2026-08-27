"""Exact current proof-graph postprocessor for THM-4270's four rays."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RATIOS = ((3, 5), (7, 8), (8, 9), (11, 12))
EXPECTED = {
    (3, 5): (32, "e1d17d394cb9d341",
             "30619f841239e9ed8b96da32a5ca554abde89197689d413585abd4cfbdc88b8f"),
    (7, 8): (44, "3b0afe245a171fcd",
             "f87090f6b399a0192cbfefaab70b0272292f0dd2786a329f5fab8caae19a0c67"),
    (8, 9): (40, "56718b95ea9dbd7c",
             "fd8815dccd0a7653646120c1879ba1c90cf4ecb82b082a626b87709dd12b7ca6"),
    (11, 12): (30, "ca69d142e1f493b0",
               "25a31b32ccf9c92bec53b61ebd4623b0c93688785b1f148500fa9a28e445167f"),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_edges(path: Path) -> list[tuple[int, int]]:
    rows = []
    for line in path.read_text(encoding="ascii").splitlines():
        left, right = line.split(",")
        rows.append((int(left), int(right)))
    require(rows == sorted(set(rows)), f"unordered or duplicate edge file: {path}")
    return rows


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


def edge_sha(edges: list[tuple[int, int]]) -> str:
    raw = b"".join(f"{left},{right}\n".encode("ascii") for left, right in edges)
    return hashlib.sha256(raw).hexdigest()


source = (
    ROOT / "05-knowledge/results/lrc14_three_round_learned_carrier_thm4266"
    / "thm4266_final_residual.csv"
)
post4266 = read_edges(source)
require(len(post4266) == 177_585, "THM-4266 residual count changed")
require(edge_fnv(post4266) == 0x6CE05D05EB01DAED,
        "THM-4266 residual FNV changed")
require(edge_sha(post4266) ==
        "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9",
        "THM-4266 residual SHA256 changed")

# Reconstruct THM-4267 semantically rather than importing its postprocessor.
ray45 = []
for left, right in post4266:
    if left % 4:
        continue
    scale = left // 4
    if 73 <= scale < 154 and right == 5 * scale:
        ray45.append((left, right))
require(len(ray45) == 63, "THM-4267 novel ray count changed")
require(edge_fnv(ray45) == 0x2280174481469C0D,
        "THM-4267 novel ray FNV changed")
require(edge_sha(ray45) ==
        "a63b3309dc30775c82ab7533472f4523ae4bda69525bc262ac06e749011dde23",
        "THM-4267 novel ray SHA256 changed")
ray45_set = set(ray45)
post4267 = [edge for edge in post4266 if edge not in ray45_set]
require(len(post4267) == 177_522, "THM-4267 residual count changed")
require(edge_fnv(post4267) == 0x33142F955CC93379,
        "THM-4267 residual FNV changed")
require(edge_sha(post4267) ==
        "d277aebe296153ead14a77207ea1499c961c8b06796b7e62f324e34f7a9ef087",
        "THM-4267 residual SHA256 changed")

# Reconstruct the subsequently promoted THM-4269 5:6 component as well.
ray56 = []
for left, right in post4267:
    if left % 5:
        continue
    scale = left // 5
    if 59 <= scale < 129 and right == 6 * scale:
        ray56.append((left, right))
require(len(ray56) == 53, "THM-4269 novel ray count changed")
require(edge_fnv(ray56) == 0x09A794C92E5AB5C9,
        "THM-4269 novel ray FNV changed")
require(edge_sha(ray56) ==
        "94a21039c9e379433890e370a3e6985cf5d85e2f61054bc6897d778bffb34823",
        "THM-4269 novel ray SHA256 changed")
ray56_set = set(ray56)
post4269 = [edge for edge in post4267 if edge not in ray56_set]
require(len(post4269) == 177_469, "THM-4269 residual count changed")
require(edge_fnv(post4269) == 0x4D1FEAE0C1E653D5,
        "THM-4269 residual FNV changed")
require(edge_sha(post4269) ==
        "289cede32347b364123827e7dea02d728b71e8c87d079a9892d3e0492b4a08ae",
        "THM-4269 residual SHA256 changed")

selected: dict[tuple[int, int], list[tuple[int, int]]] = {}
for ratio in RATIOS:
    u, v = ratio
    first = 290 // u + 1
    tail = (770 + v - 1) // v
    rows = []
    for left, right in post4269:
        if left % u:
            continue
        scale = left // u
        if first <= scale < tail and right == v * scale:
            rows.append((left, right))
    selected[ratio] = rows
    expected_count, expected_fnv, expected_sha = EXPECTED[ratio]
    require(len(rows) == expected_count, f"{ratio} selected count changed")
    require(f"{edge_fnv(rows):016x}" == expected_fnv,
            f"{ratio} selected FNV changed")
    require(edge_sha(rows) == expected_sha, f"{ratio} selected SHA256 changed")

all_rows = [edge for ratio in RATIOS for edge in selected[ratio]]
require(len(all_rows) == len(set(all_rows)), "new ray contributions overlap")
require(len(all_rows) == 146, "new ray union count changed")
require(edge_fnv(all_rows) == 0x3115735824BDB7F5,
        "new ray union FNV changed")
require(edge_sha(all_rows) ==
        "b2fcda2e0e602a4284d243d486ae426f7376de0be429c30d1bb40ee6c455f750",
        "new ray union SHA256 changed")
row_set = set(all_rows)
updated = [edge for edge in post4269 if edge not in row_set]
require(len(updated) == 177_323, "updated residual count changed")
require(edge_fnv(updated) == 0xF1DCC8033FA727D9,
        "updated residual FNV changed")
require(edge_sha(updated) ==
        "8c0b1fac01d00bd54784178034f4e5f21c2a29ea95a9cb0ed5a63b06fbc20872",
        "updated residual SHA256 changed")
require(max(right for _, right in updated) == 688,
        "updated residual maximum endpoint changed")
require([edge for edge in updated if edge[1] == 688] == [(520, 688)],
        "updated residual top layer changed")

print("LRC14_FOUR_PRIMITIVE_RAYS_CURRENT_PROOF_GRAPH_THM4270")
print(
    f"POST_THM4266 COUNT {len(post4266)} FNV {edge_fnv(post4266):016x} "
    f"SHA256 {edge_sha(post4266)}"
)
print(
    f"POST_THM4267 COUNT {len(post4267)} FNV {edge_fnv(post4267):016x} "
    f"SHA256 {edge_sha(post4267)}"
)
print(
    f"POST_THM4269 COUNT {len(post4269)} FNV {edge_fnv(post4269):016x} "
    f"SHA256 {edge_sha(post4269)}"
)
for ratio in RATIOS:
    rows = selected[ratio]
    scales = [left // ratio[0] for left, _ in rows]
    print(
        f"RATIO {ratio[0]}:{ratio[1]} NOVEL_COUNT {len(rows)} "
        f"FNV {edge_fnv(rows):016x} SHA256 {edge_sha(rows)} "
        f"SCALES {','.join(map(str, scales))}"
    )
print(
    f"NOVEL_UNION COUNT {len(all_rows)} FNV {edge_fnv(all_rows):016x} "
    f"SHA256 {edge_sha(all_rows)}"
)
print(
    f"UPDATED COUNT {len(updated)} FNV {edge_fnv(updated):016x} "
    f"SHA256 {edge_sha(updated)} MAX {max(right for _, right in updated)} "
    "TOP 520,688"
)
print("VERDICT EXACT_CURRENT_PROOF_GRAPH LRC14_OPEN")
