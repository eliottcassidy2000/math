#!/usr/bin/env python3
"""Exact current proof-graph accounting for the 2D rectangle certificate."""

from __future__ import annotations

import argparse
import hashlib
import math
from collections import Counter
from pathlib import Path

OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
BOUNDARY_4271 = {(256, 671), (384, 671)}
BOUNDARY_4276 = {(256, 670), (384, 670)}
RATIOS_4270 = ((3, 5), (7, 8), (8, 9), (11, 12))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_edges(path: Path) -> list[tuple[int, int]]:
    edges: list[tuple[int, int]] = []
    for number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        fields = line.split(",")
        require(len(fields) == 2, f"bad edge line {number}")
        edges.append((int(fields[0]), int(fields[1])))
    require(edges == sorted(set(edges)), "edge source is not a sorted set")
    return edges


def fnv(edges: list[tuple[int, int]]) -> str:
    state = OFFSET
    for q, r in edges:
        for word in (q, r):
            for shift in range(0, 64, 8):
                state ^= (word >> shift) & 0xFF
                state = state * PRIME & MASK64
    return f"{state:016x}"


def sha(edges: list[tuple[int, int]]) -> str:
    raw = b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)
    return hashlib.sha256(raw).hexdigest()


def fp(label: str, edges: list[tuple[int, int]]) -> None:
    print(f"{label} COUNT {len(edges)} FNV {fnv(edges)} SHA256 {sha(edges)}")


def remove(before: list[tuple[int, int]], deletion: list[tuple[int, int]]) -> list[tuple[int, int]]:
    deleted = set(deletion)
    require(len(deleted) == len(deletion), "deletion contains duplicates")
    require(deleted <= set(before), "deletion is not contained in source")
    return [edge for edge in before if edge not in deleted]


def ray(edges: list[tuple[int, int]], u: int, v: int) -> list[tuple[int, int]]:
    return [(q, r) for q, r in edges if q > 290 and v * q == u * r]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    args = parser.parse_args()
    source = (
        args.repo.resolve()
        / "05-knowledge/results/lrc14_three_round_learned_carrier_thm4266"
        / "thm4266_final_residual.csv"
    )
    post4266 = read_edges(source)
    require((len(post4266), fnv(post4266), sha(post4266)) == (
        177_585,
        "6ce05d05eb01daed",
        "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9",
    ), "post-THM-4266 ledger changed")

    post4267 = remove(post4266, ray(post4266, 4, 5))
    require((len(post4267), fnv(post4267), sha(post4267)) == (
        177_522,
        "33142f955cc93379",
        "d277aebe296153ead14a77207ea1499c961c8b06796b7e62f324e34f7a9ef087",
    ), "post-THM-4267 ledger changed")

    post4269 = remove(post4267, ray(post4267, 5, 6))
    require((len(post4269), fnv(post4269), sha(post4269)) == (
        177_469,
        "4d1feae0c1e653d5",
        "289cede32347b364123827e7dea02d728b71e8c87d079a9892d3e0492b4a08ae",
    ), "post-THM-4269 ledger changed")

    components4270 = [ray(post4269, u, v) for u, v in RATIOS_4270]
    deletion4270 = [edge for component in components4270 for edge in component]
    require(len(deletion4270) == 146 and len(set(deletion4270)) == 146,
            "THM-4270 component changed")
    post4270 = remove(post4269, deletion4270)
    require((len(post4270), fnv(post4270), sha(post4270)) == (
        177_323,
        "f1dcc8033fa727d9",
        "8c0b1fac01d00bd54784178034f4e5f21c2a29ea95a9cb0ed5a63b06fbc20872",
    ), "post-THM-4270 ledger changed")

    deletion4271 = [edge for edge in post4270
                    if edge[1] >= 671 and edge not in BOUNDARY_4271]
    post4271 = remove(post4270, deletion4271)
    require((len(deletion4271), fnv(deletion4271), sha(deletion4271)) == (
        2_419,
        "ff169664a6750abe",
        "7980354ba5d9dde9ce994eda992bf7030d79bd0562cf9c8f2742d4bb53653e89",
    ), "THM-4271 deletion changed")
    require((len(post4271), fnv(post4271), sha(post4271)) == (
        174_904,
        "b3db855040bcf19e",
        "07d84c1572baeb89d9f88e095788e52e3916dab6074f89cf3b164e5ebea3a5a6",
    ), "post-THM-4271 ledger changed")

    deletion4276 = [edge for edge in post4271
                    if edge[1] >= 670 and edge not in BOUNDARY_4276]
    post4276 = remove(post4271, deletion4276)
    require((len(deletion4276), fnv(deletion4276), sha(deletion4276)) == (
        163,
        "edc377490e3f58bb",
        "49960a608af498bf385808f3fe234b75923ae203ea9a6e95875e61df1c43de26",
    ), "THM-4276 deletion changed")
    require((len(post4276), fnv(post4276), sha(post4276)) == (
        174_741,
        "f13745b05320f83c",
        "51d5723b146cb108a2e11627924a2fd6af46435564e2460ab78af936bfb12dd0",
    ), "post-THM-4276 ledger changed")

    full_rectangle = [(q, r) for q in range(450, 500) for r in range(600, 651)]
    ratio_counts = Counter((q // math.gcd(q, r), r // math.gcd(q, r))
                           for q, r in full_rectangle)
    multiplicities = Counter(ratio_counts.values())
    gcd_values = {math.gcd(q, r) for q, r in full_rectangle}
    require(len(full_rectangle) == 2_550, "full rectangle census changed")
    require(len(ratio_counts) == 2_500, "primitive-ratio census changed")
    require(multiplicities == {1: 2_474, 2: 16, 3: 7, 5: 2, 13: 1},
            "primitive-ratio multiplicity profile changed")
    require((len(gcd_values), min(gcd_values), max(gcd_values)) == (77, 1, 162),
            "rectangle gcd profile changed")

    rectangle = [(q, r) for q, r in post4276
                 if 450 <= q <= 499 and 600 <= r <= 650]
    require((len(rectangle), fnv(rectangle), sha(rectangle)) == (
        2_419,
        "67b373dc22ac870d",
        "e9ee49675d5345b06157e64a060506be3f6dd6a94835cc92f6eb5138f346cffb",
    ), "rectangle novelty ledger changed")
    require(not any(v * q == u * r for q, r in rectangle
                    for u, v in ((4, 5), (5, 6), *RATIOS_4270)),
            "rectangle reclaims a proved ray")
    require(not set(rectangle) & set(deletion4276),
            "rectangle overlaps THM-4276 deletion")
    final = remove(post4276, rectangle)
    require((len(final), fnv(final), sha(final)) == (
        172_322,
        "30b2a7e597ac548c",
        "7228658eae4067db4bbcceb6c9b1ccebf1bd3e6f128e202ea184854acc53f309",
    ), "final residual ledger changed")
    require(max(r for _, r in final) == 670, "final maximum endpoint changed")
    require([edge for edge in final if edge[1] == 670] == [(256, 670), (384, 670)],
            "final top layer changed")

    print("LRC14_2D_RECTANGLE_CURRENT_PROOF_GRAPH_V1")
    print("FULL_RECTANGLE PAIRS 2550 PRIMITIVE_RATIOS 2500 "
          "MULTIPLICITIES 1:2474,2:16,3:7,5:2,13:1 "
          "DISTINCT_GCDS 77 MIN_GCD 1 MAX_GCD 162")
    fp("POST_THM4270", post4270)
    fp("THM4271_NOVEL", deletion4271)
    fp("POST_THM4271", post4271)
    fp("THM4276_NOVEL", deletion4276)
    fp("POST_THM4276", post4276)
    fp("RECTANGLE_NOVEL", rectangle)
    fp("FINAL_AFTER_RECTANGLE", final)
    print("RECTANGLE_DISJOINT_PROVED_RAYS 4:5,5:6,3:5,7:8,8:9,11:12")
    print("RECTANGLE_DISJOINT_THM4276 YES")
    print("FINAL_TOP ENDPOINT 670 COUNT 2 EDGES 256,670;384,670")
    print("VERDICT PASS EXACT_CURRENT_PROOF_GRAPH LRC14_OPEN")


if __name__ == "__main__":
    main()
