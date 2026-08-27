#!/usr/bin/env python3
"""Portable semantic proof-graph audit for the compact round-five carrier."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
RATIOS_4270 = ((3, 5), (7, 8), (8, 9), (11, 12))
HOSTILES_4271 = {(256, 671), (384, 671)}
HOSTILES_ROUND5 = {(256, 670), (384, 670)}


def read_edges(path: Path) -> list[tuple[int, int]]:
    edges = [tuple(map(int, line.split(",")))
             for line in path.read_text(encoding="ascii").splitlines()]
    assert edges == sorted(set(edges))
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
    return hashlib.sha256(b"".join(
        f"{q},{r}\n".encode("ascii") for q, r in edges
    )).hexdigest()


def fp(label: str, edges: list[tuple[int, int]]) -> None:
    print(f"{label} COUNT {len(edges)} FNV {fnv(edges)} SHA256 {sha(edges)}")


def remove(edges: list[tuple[int, int]],
           deletion: set[tuple[int, int]]) -> list[tuple[int, int]]:
    assert deletion <= set(edges)
    return [edge for edge in edges if edge not in deletion]


def ray(edges: list[tuple[int, int]],
        u: int, v: int) -> list[tuple[int, int]]:
    return [(q, r) for q, r in edges if v * q == u * r and q > 290]


def audit_compact_transcript(path: Path) -> None:
    text = path.read_text(encoding="ascii")
    assert (
        "COMPACT_ROUND5 NOVEL 6 UNION 8524 UNION_FNV 5ddb84a44f5d2ad7 "
        "SEED256_ACTIVE 4296 SEED256_FAILURES 0 "
        "SEED384_ACTIVE 5994 SEED384_FAILURES 0"
    ) in text
    layer = re.search(
        r"^LAYER 670 ROWS (\d+) RESISTANT (\d+) ACTIVE_SUM (\d+) "
        r"CHECKS (\d+) MAX_CHECKS (\d+) ROW_LEDGER ([0-9a-f]+)$",
        text,
        re.MULTILINE,
    )
    assert layer is not None
    assert layer.groups() == (
        "163", "2", "1331051", "69626501969", "3851",
        "90e990cdced5a71b",
    )
    hostiles = [tuple(map(int, match)) for match in re.findall(
        r"^RESISTANT PAIR (\d+),(670) ", text, re.MULTILINE
    )]
    assert hostiles == [(256, 670), (384, 670)]
    assert text.rstrip().endswith(
        "STOP FIRST_RESISTANT_LAYER 670 VERDICT BOUNDARY_FOUND"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--compact-transcript", type=Path, required=True)
    parser.add_argument("--emit-post4271", type=Path)
    args = parser.parse_args()
    audit_compact_transcript(args.compact_transcript)

    source = (
        args.repo.resolve()
        / "05-knowledge/results/lrc14_three_round_learned_carrier_thm4266"
        / "thm4266_final_residual.csv"
    )
    post4266 = read_edges(source)
    assert (len(post4266), fnv(post4266), sha(post4266)) == (
        177_585,
        "6ce05d05eb01daed",
        "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9",
    )
    post4267 = remove(post4266, set(ray(post4266, 4, 5)))
    assert (len(post4267), fnv(post4267), sha(post4267)) == (
        177_522,
        "33142f955cc93379",
        "d277aebe296153ead14a77207ea1499c961c8b06796b7e62f324e34f7a9ef087",
    )
    post4269 = remove(post4267, set(ray(post4267, 5, 6)))
    assert (len(post4269), fnv(post4269), sha(post4269)) == (
        177_469,
        "4d1feae0c1e653d5",
        "289cede32347b364123827e7dea02d728b71e8c87d079a9892d3e0492b4a08ae",
    )
    deletion4270: set[tuple[int, int]] = set()
    for u, v in RATIOS_4270:
        component = set(ray(post4269, u, v))
        assert not deletion4270 & component
        deletion4270 |= component
    assert len(deletion4270) == 146
    post4270 = remove(post4269, deletion4270)
    assert (len(post4270), fnv(post4270), sha(post4270)) == (
        177_323,
        "f1dcc8033fa727d9",
        "8c0b1fac01d00bd54784178034f4e5f21c2a29ea95a9cb0ed5a63b06fbc20872",
    )

    raw4271 = {
        edge for edge in post4267
        if edge[1] >= 671 and edge not in HOSTILES_4271
    }
    current4271 = raw4271 & set(post4270)
    assert len(current4271) == 2_419
    post4271 = remove(post4270, current4271)
    assert (len(post4271), fnv(post4271), sha(post4271)) == (
        174_904,
        "b3db855040bcf19e",
        "07d84c1572baeb89d9f88e095788e52e3916dab6074f89cf3b164e5ebea3a5a6",
    )
    assert [edge for edge in post4271 if edge[1] == 671] == [
        (256, 671), (384, 671)
    ]
    if args.emit_post4271 is not None:
        args.emit_post4271.write_text(
            "".join(f"{q},{r}\n" for q, r in post4271), encoding="ascii"
        )

    compact_closure = [
        edge for edge in post4271
        if edge[1] >= 670 and edge not in HOSTILES_ROUND5
    ]
    assert len(compact_closure) == 163
    final = remove(post4271, set(compact_closure))
    assert (len(compact_closure), fnv(compact_closure),
            sha(compact_closure)) == (
        163,
        "edc377490e3f58bb",
        "49960a608af498bf385808f3fe234b75923ae203ea9a6e95875e61df1c43de26",
    )
    assert (len(final), fnv(final), sha(final)) == (
        174_741,
        "f13745b05320f83c",
        "51d5723b146cb108a2e11627924a2fd6af46435564e2460ab78af936bfb12dd0",
    )
    assert max(r for _, r in final) == 670
    assert [edge for edge in final if edge[1] == 670] == [
        (256, 670), (384, 670)
    ]

    print("ROUND5_COMPACT_CURRENT_PROOF_GRAPH_V1")
    fp("POST_THM4270", post4270)
    fp("THM4271_CURRENT_DELETION", sorted(current4271))
    fp("POST_THM4271", post4271)
    fp("ROUND5_COMPACT_DELETION", compact_closure)
    fp("FINAL_AFTER_ROUND5", final)
    print("ROUND5_DELETION_LAYERS R671=2 R670=161")
    print("FINAL_TOP ENDPOINT 670 COUNT 2 EDGES 256,670;384,670")
    print("VERDICT PASS EXACT_CURRENT_PROOF_GRAPH LRC14_OPEN")


if __name__ == "__main__":
    main()
