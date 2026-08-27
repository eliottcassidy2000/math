"""Exact THM-4270 proof-graph postprocessor for four fixed-pool rays."""

from __future__ import annotations

import contextlib
import hashlib
import io
import os
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[1]
RATIOS = ((3, 5), (7, 8), (8, 9), (11, 12))
EXPECTED = {
    (3, 5): (32, "e1d17d394cb9d341",
             "30619f841239e9ed8b96da32a5ca554abde89197689d413585abd4cfbdc88b8f"),
    (7, 8): (46, "fd0e23c7ec1ddd2f",
             "d34d828d3bcd929a2ea7ad478aa2c01112562a0e1442f95dbc63da13bb6d4880"),
    (8, 9): (43, "f193ad0e40f85f36",
             "11e449a85df94fb65833c440727f7ef74d378d6e564ac6d4e155294936d8df38"),
    (11, 12): (31, "34a2ac9fdaeac516",
               "fbc6dcff82543ee76ad91be79e9ca52f010072eec7a95b557a84019ccfa4ee67"),
}


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


def edge_sha(edges: list[tuple[int, int]]) -> str:
    raw = b"".join(f"{left},{right}\n".encode("ascii") for left, right in edges)
    return hashlib.sha256(raw).hexdigest()


def load(relative: str):
    prior = Path.cwd()
    try:
        os.chdir(ROOT)
        with contextlib.redirect_stdout(io.StringIO()):
            return runpy.run_path(str(ROOT / relative))
    finally:
        os.chdir(prior)


band_namespace = load(
    "04-computation/"
    "lrc14_prefix_union_semantic_band_residual_postprocess_thm4261.py"
)
ratio34_namespace = load(
    "04-computation/lrc14_ratio34_residual_postprocess_thm4262.py"
)

post4261 = band_namespace["updated"]
post4262_only = ratio34_namespace["updated"]
ray34 = set(ratio34_namespace["ray"])
current = [edge for edge in post4261 if edge not in ray34]
require(len(current) == 180_622, "combined THM-4261/4262 residual changed")
require(set(current) == set(post4262_only) - set(band_namespace["band"]),
        "combined residual order-independent check failed")
require(edge_fnv(current) == 0x0CEF4E2887C8F24E,
        "inherited residual FNV changed")
require(edge_sha(current) ==
        "fa1c5672b0f2cd2490413e9b69a4720bf1dc4eef8aee694c1c73d390aba58e11",
        "inherited residual SHA256 changed")

selected: dict[tuple[int, int], list[tuple[int, int]]] = {}
for ratio in RATIOS:
    u, v = ratio
    first = 290 // u + 1
    tail = (770 + v - 1) // v
    rows = []
    for left, right in current:
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
row_set = set(all_rows)
updated = [edge for edge in current if edge not in row_set]
require(len(all_rows) == 152, "new ray union count changed")
require(edge_fnv(all_rows) == 0x80ED45E7E179D8FF,
        "new ray union FNV changed")
require(edge_sha(all_rows) ==
        "46b1c987fe337d879a721efab99da8265d126d12798788295b0829f6bd5741fd",
        "new ray union SHA256 changed")
require(len(updated) == 180_470, "updated residual count changed")
require(edge_fnv(updated) == 0xD9AFDC10F8E2AA88,
        "updated residual FNV changed")
require(edge_sha(updated) ==
        "e2ba42307e4c628ea9ef517a858456e4dd64e0b7a034fd6dd74ff707a8838f3f",
        "updated residual SHA256 changed")
require(max(right for _, right in updated) == 732,
        "updated residual maximum endpoint changed")

print("LRC14_FOUR_PRIMITIVE_RAYS_RESIDUAL_POSTPROCESS_THM4270")
print(
    f"INHERITED COUNT {len(current)} FNV {edge_fnv(current):016x} "
    f"SHA256 {edge_sha(current)}"
)
for ratio in RATIOS:
    rows = selected[ratio]
    scales = [left // ratio[0] for left, _ in rows]
    print(
        f"RATIO {ratio[0]}:{ratio[1]} COUNT {len(rows)} "
        f"FNV {edge_fnv(rows):016x} SHA256 {edge_sha(rows)} "
        f"SCALES {','.join(map(str, scales))}"
    )
print(
    f"UNION COUNT {len(all_rows)} FNV {edge_fnv(all_rows):016x} "
    f"SHA256 {edge_sha(all_rows)}"
)
print(
    f"UPDATED COUNT {len(updated)} FNV {edge_fnv(updated):016x} "
    f"SHA256 {edge_sha(updated)} MAX {max(right for _, right in updated)}"
)
print("VERDICT EXACT_POSTPROCESS LRC14_OPEN")
