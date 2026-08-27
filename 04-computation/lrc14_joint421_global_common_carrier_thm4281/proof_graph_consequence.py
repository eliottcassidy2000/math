#!/usr/bin/env python3
"""Exact proof-graph consequence of the joint421 carrier descent."""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from pathlib import Path

OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
BOUNDARY_4276 = {(256, 670), (384, 670)}
BOUNDARY_4281: set[tuple[int, int]] = set()


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--post4271", type=Path, required=True)
    parser.add_argument("--common", type=Path, required=True)
    parser.add_argument("--emit-final", type=Path)
    args = parser.parse_args()

    post4271 = read_edges(args.post4271)
    require((len(post4271), fnv(post4271), sha(post4271)) == (
        174_904,
        "b3db855040bcf19e",
        "07d84c1572baeb89d9f88e095788e52e3916dab6074f89cf3b164e5ebea3a5a6",
    ), "post-THM-4271 ledger changed")

    deletion4276 = [edge for edge in post4271
                    if edge[1] >= 670 and edge not in BOUNDARY_4276]
    post4276 = remove(post4271, deletion4276)
    require((len(post4276), fnv(post4276), sha(post4276)) == (
        174_741,
        "f13745b05320f83c",
        "51d5723b146cb108a2e11627924a2fd6af46435564e2460ab78af936bfb12dd0",
    ), "post-THM-4276 ledger changed")

    rectangle = [(q, r) for q, r in post4276
                 if 450 <= q <= 499 and 600 <= r <= 650]
    post4277 = remove(post4276, rectangle)
    require((len(rectangle), fnv(rectangle), sha(rectangle)) == (
        2_419,
        "67b373dc22ac870d",
        "e9ee49675d5345b06157e64a060506be3f6dd6a94835cc92f6eb5138f346cffb",
    ), "THM-4277 rectangle ledger changed")
    require((len(post4277), fnv(post4277), sha(post4277)) == (
        172_322,
        "30b2a7e597ac548c",
        "7228658eae4067db4bbcceb6c9b1ccebf1bd3e6f128e202ea184854acc53f309",
    ), "post-THM-4277 ledger changed")

    common = read_edges(args.common)
    require(set(common) <= set(post4277),
            "common-active deletion is not contained in post-THM4277")
    require((len(common), fnv(common), sha(common)) == (
        148_063,
        "465fb756a183167e",
        "412b942759ed7afde1bffaaeabc9e6ae31b8b8bc25f26d73c71712523a057aaa",
    ), "global common-active ledger changed")
    complement = remove(post4277, common)
    require((len(complement), fnv(complement), sha(complement)) == (
        24_259,
        "78b212469c336f37",
        "77f3d21d127bb5e21f583556314f74032271e3a4903696415885308b363624ef",
    ), "post-common complement ledger changed")

    carrier_raw = [edge for edge in post4277
                   if edge[1] >= 664 and edge not in BOUNDARY_4281]
    common_set = set(common)
    carrier_overlap = [edge for edge in carrier_raw if edge in common_set]
    carrier_novel = [edge for edge in carrier_raw if edge not in common_set]
    carrier_novel_set = set(carrier_novel)
    deletion4281 = [edge for edge in post4277
                    if edge in common_set or edge in carrier_novel_set]
    final = remove(post4277, deletion4281)
    require((len(carrier_raw), fnv(carrier_raw), sha(carrier_raw)) == (
        1_068,
        "c61906c5ef460bb8",
        "ffbb4e8030dae146267a623009603df4cd631468d80a5fa0a8965ec7a5232b93",
    ), "raw endpoint-carrier deletion ledger changed")
    require((len(carrier_overlap), fnv(carrier_overlap), sha(carrier_overlap)) == (
        1_032,
        "edb538ac3fb4c990",
        "c3357b39a4b356b07c51654114f3b6c2e8e2455fe5ded630a0526686982f7a69",
    ), "carrier/common overlap ledger changed")
    require((len(carrier_novel), fnv(carrier_novel), sha(carrier_novel)) == (
        36,
        "a8027022c6436919",
        "b14f2cb2ff92ea4f14676e3b6cf6b63fcff1b3378d2ece0461a9a4c89bccfd0e",
    ), "carrier-novel ledger changed")
    require(len(common) + len(carrier_novel) == len(deletion4281),
            "THM4281 union count changed")
    require((len(deletion4281), fnv(deletion4281), sha(deletion4281)) == (
        148_099,
        "0308f3f6b8d7a23e",
        "ed993db4e2b53a416e6031696c4a0b01dfd3f9b9e47fb2607f5f191059752836",
    ), "THM4281 total deletion ledger changed")
    require((len(final), fnv(final), sha(final)) == (
        24_223,
        "80ec0687d8c7dba7",
        "75a3c7616c982538363c7801ed2dab3fe9aa775ab601f7a7119dd9fb5d301552",
    ), "final residual ledger changed")
    require(max(r for _, r in final) == 663, "final maximum endpoint changed")
    require([edge for edge in final if edge[1] == 663] ==
                [(256, 663), (366, 663), (520, 663)],
            "final top layer changed")
    carrier_layers = Counter(r for _, r in carrier_raw)
    require(carrier_layers == {670: 2, 669: 168, 668: 170, 667: 176,
                               666: 179, 665: 183, 664: 190},
            "THM-4281 raw carrier layer census changed")

    print("JOINT421_CARRIER_CURRENT_PROOF_GRAPH_V1")
    fp("POST_THM4271", post4271)
    fp("POST_THM4276", post4276)
    fp("THM4277_RECTANGLE", rectangle)
    fp("POST_THM4277", post4277)
    fp("THM4281_COMMON_ACTIVE", common)
    fp("POST_COMMON_COMPLEMENT", complement)
    fp("THM4281_CARRIER_RAW", carrier_raw)
    fp("THM4281_CARRIER_COMMON_OVERLAP", carrier_overlap)
    fp("THM4281_CARRIER_NOVEL", carrier_novel)
    fp("THM4281_TOTAL_NOVEL", deletion4281)
    fp("FINAL_AFTER_THM4281", final)
    print("THM4281_CARRIER_RAW_LAYERS R670=2 R669=168 R668=170 R667=176 "
          "R666=179 R665=183 R664=190")
    print("THM4281_DISJOINT_THM4277 YES")
    print("FINAL_TOP ENDPOINT 663 COUNT 3 EDGES 256,663;366,663;520,663")
    print("VERDICT PASS EXACT_CURRENT_PROOF_GRAPH LRC14_OPEN")

    if args.emit_final is not None:
        args.emit_final.write_text(
            "".join(f"{q},{r}\n" for q, r in final), encoding="ascii")


if __name__ == "__main__":
    main()
