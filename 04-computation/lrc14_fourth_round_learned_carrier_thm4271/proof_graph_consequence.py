#!/usr/bin/env python3
"""Portable exact proof-graph postprocessor for THM-4271.

The script starts only from the canonical THM-4266 residual.  It rebuilds
the THM-4267, THM-4269, and THM-4270 deletions by their semantic ray
relations, then accounts for the round-four carrier component in every
requested proof order.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
BOUNDARY_HOSTILES = {(256, 671), (384, 671)}
RATIOS_4270 = ((3, 5), (7, 8), (8, 9), (11, 12))
EXPECTED_LAYERS = [
    (688, 1, 0), (687, 126, 0), (686, 126, 0), (685, 127, 0),
    (684, 130, 0), (683, 133, 0), (682, 135, 0), (681, 136, 0),
    (680, 139, 0), (679, 143, 0), (678, 144, 0), (677, 146, 0),
    (676, 149, 0), (675, 153, 0), (674, 156, 0), (673, 159, 0),
    (672, 161, 0), (671, 162, 2),
]


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
    data = b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)
    return hashlib.sha256(data).hexdigest()


def fingerprint(label: str, edges: list[tuple[int, int]]) -> None:
    print(f"{label} COUNT {len(edges)} FNV {fnv(edges)} SHA256 {sha(edges)}")


def remove(before: list[tuple[int, int]],
           deletion: list[tuple[int, int]]) -> list[tuple[int, int]]:
    deletion_set = set(deletion)
    assert len(deletion_set) == len(deletion)
    assert deletion_set <= set(before)
    return [edge for edge in before if edge not in deletion_set]


def semantic_ray(edges: list[tuple[int, int]],
                 u: int, v: int) -> list[tuple[int, int]]:
    return [(q, r) for q, r in edges if v * q == u * r and q > 290]


def audit_descent(path: Path) -> None:
    text = path.read_text(encoding="ascii")
    layers = [tuple(map(int, match)) for match in re.findall(
        r"^LAYER (\d+) ROWS (\d+) RESISTANT (\d+) ", text, re.MULTILINE
    )]
    hostiles = [tuple(map(int, match)) for match in re.findall(
        r"^RESISTANT PAIR (\d+),(\d+) ", text, re.MULTILINE
    )]
    assert layers == EXPECTED_LAYERS
    assert hostiles == [(256, 671), (384, 671)]
    assert sum(rows - resistant for _, rows, resistant in layers) == 2_424
    assert text.rstrip().endswith(
        "STOP FIRST_RESISTANT_LAYER 671 VERDICT BOUNDARY_FOUND"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--emit-post4267", type=Path)
    parser.add_argument("--descent-output", type=Path)
    args = parser.parse_args()
    source = (
        args.repo.resolve()
        / "05-knowledge/results/lrc14_three_round_learned_carrier_thm4266"
        / "thm4266_final_residual.csv"
    )
    descent_output = args.descent_output or (
        args.repo.resolve()
        / "05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271"
        / "round4_descent_688_671_O3.out"
    )
    audit_descent(descent_output)
    post4266 = read_edges(source)
    assert (len(post4266), fnv(post4266), sha(post4266)) == (
        177_585,
        "6ce05d05eb01daed",
        "009614651bb81e9763b2a9ff4b580497bfb6978a6c69d18cf986346e369374d9",
    )

    ray45 = semantic_ray(post4266, 4, 5)
    post4267 = remove(post4266, ray45)
    assert (len(ray45), fnv(ray45), sha(ray45)) == (
        63,
        "2280174481469c0d",
        "a63b3309dc30775c82ab7533472f4523ae4bda69525bc262ac06e749011dde23",
    )
    assert (len(post4267), fnv(post4267), sha(post4267)) == (
        177_522,
        "33142f955cc93379",
        "d277aebe296153ead14a77207ea1499c961c8b06796b7e62f324e34f7a9ef087",
    )
    if args.emit_post4267 is not None:
        args.emit_post4267.write_text(
            "".join(f"{q},{r}\n" for q, r in post4267), encoding="ascii"
        )

    raw4271 = [edge for edge in post4267
               if edge[1] >= 671 and edge not in BOUNDARY_HOSTILES]
    final_after4267 = remove(post4267, raw4271)
    assert (len(raw4271), fnv(raw4271), sha(raw4271)) == (
        2_424,
        "3ec97f3ae7e5d142",
        "44f00e2b22adbf071eaebe3edf474337abf789ac9dad867cc48d0e157cfebc94",
    )
    assert (len(final_after4267), fnv(final_after4267),
            sha(final_after4267)) == (
        175_098,
        "917c15bec89ad602",
        "e7f5c682295979ee3d15b41d5e3b53de0e9557e7494d6d5355a180c2643be39e",
    )

    ray56 = semantic_ray(post4267, 5, 6)
    post4269 = remove(post4267, ray56)
    assert (len(ray56), fnv(ray56), sha(ray56)) == (
        53,
        "09a794c92e5ab5c9",
        "94a21039c9e379433890e370a3e6985cf5d85e2f61054bc6897d778bffb34823",
    )
    assert (len(post4269), fnv(post4269), sha(post4269)) == (
        177_469,
        "4d1feae0c1e653d5",
        "289cede32347b364123827e7dea02d728b71e8c87d079a9892d3e0492b4a08ae",
    )
    ray56_set = set(ray56)
    overlap4269 = [edge for edge in raw4271 if edge in ray56_set]
    assert overlap4269 == []

    rays4270 = [semantic_ray(post4269, u, v) for u, v in RATIOS_4270]
    expected = [
        (32, "e1d17d394cb9d341",
         "30619f841239e9ed8b96da32a5ca554abde89197689d413585abd4cfbdc88b8f"),
        (44, "3b0afe245a171fcd",
         "f87090f6b399a0192cbfefaab70b0272292f0dd2786a329f5fab8caae19a0c67"),
        (40, "56718b95ea9dbd7c",
         "fd8815dccd0a7653646120c1879ba1c90cf4ecb82b082a626b87709dd12b7ca6"),
        (30, "ca69d142e1f493b0",
         "25a31b32ccf9c92bec53b61ebd4623b0c93688785b1f148500fa9a28e445167f"),
    ]
    for ray, wanted in zip(rays4270, expected):
        assert (len(ray), fnv(ray), sha(ray)) == wanted
    union4270_component_order = [edge for ray in rays4270 for edge in ray]
    assert len(set(union4270_component_order)) == 146
    assert (fnv(union4270_component_order), sha(union4270_component_order)) == (
        "3115735824bdb7f5",
        "b2fcda2e0e602a4284d243d486ae426f7376de0be429c30d1bb40ee6c455f750",
    )
    post4270 = remove(post4269, union4270_component_order)
    assert (len(post4270), fnv(post4270), sha(post4270)) == (
        177_323,
        "f1dcc8033fa727d9",
        "8c0b1fac01d00bd54784178034f4e5f21c2a29ea95a9cb0ed5a63b06fbc20872",
    )

    union4270_set = set(union4270_component_order)
    overlap4270 = [edge for edge in raw4271 if edge in union4270_set]
    assert overlap4270 == [
        (588, 672),
        (595, 680),
        (600, 675),
        (608, 684),
        (616, 672),
    ]
    post4270_set = set(post4270)
    novel_after4270 = [edge for edge in raw4271 if edge in post4270_set]
    final_after4271 = remove(post4270, novel_after4270)
    assert (len(overlap4270), fnv(overlap4270), sha(overlap4270)) == (
        5,
        "7484ad77ffddd129",
        "7eafa6774730c49213c08eb7b595952ca64c4e224db30b348933bfa422f82060",
    )
    assert (len(novel_after4270), fnv(novel_after4270),
            sha(novel_after4270)) == (
        2_419,
        "ff169664a6750abe",
        "7980354ba5d9dde9ce994eda992bf7030d79bd0562cf9c8f2742d4bb53653e89",
    )
    assert (len(final_after4271), fnv(final_after4271),
            sha(final_after4271)) == (
        174_904,
        "b3db855040bcf19e",
        "07d84c1572baeb89d9f88e095788e52e3916dab6074f89cf3b164e5ebea3a5a6",
    )
    assert max(r for _, r in final_after4271) == 671
    assert [edge for edge in final_after4271 if edge[1] == 671] == [
        (256, 671), (384, 671)
    ]

    print("THM4271_ROUND4_CURRENT_PROOF_GRAPH_V1")
    fingerprint("POST_THM4266", post4266)
    fingerprint("POST_THM4267", post4267)
    fingerprint("RAW_THM4271_AFTER4267", raw4271)
    fingerprint("FINAL_IF_AFTER4267", final_after4267)
    fingerprint("THM4269_5_6", ray56)
    print("OVERLAP_THM4269 COUNT 0")
    fingerprint("POST_THM4269", post4269)
    fingerprint("THM4270_COMPONENT_ORDER", union4270_component_order)
    fingerprint("POST_THM4270", post4270)
    fingerprint("OVERLAP_THM4270", overlap4270)
    print("OVERLAP_THM4270_EDGES " + ";".join(
        f"{q},{r}" for q, r in overlap4270
    ))
    fingerprint("THM4271_NOVEL_AFTER4270", novel_after4270)
    fingerprint("FINAL_AFTER_THM4271", final_after4271)
    print("FINAL_TOP ENDPOINT 671 COUNT 2 EDGES 256,671;384,671")
    print("VERDICT PASS EXACT_CURRENT_PROOF_GRAPH LRC14_OPEN")


if __name__ == "__main__":
    main()
