#!/usr/bin/env python3
"""Exact set ledger for the joint421 common-active/carrier-closure union.

This verifier independently reconstructs the post-THM-4277 universe from the
post-THM-4271 CSV.  It treats the primary common-active CSV as an exact input,
and reconstructs the raw carrier-closed set from the endpoint descent certified
by the supplied literal/scan transcripts.  All edge ledgers use sorted ``q,r``
ASCII rows and 64-bit little-endian FNV-1a over the two integer coordinates.
Layer-table FNV uses the same encoding over every numeric field in column order.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1

BOUNDARY_4276 = {(256, 670), (384, 670)}
BOUNDARY_4281: set[tuple[int, int]] = set()

EXPECTED_POST4271 = (
    174_904,
    "b3db855040bcf19e",
    "07d84c1572baeb89d9f88e095788e52e3916dab6074f89cf3b164e5ebea3a5a6",
)
EXPECTED_POST4276 = (
    174_741,
    "f13745b05320f83c",
    "51d5723b146cb108a2e11627924a2fd6af46435564e2460ab78af936bfb12dd0",
)
EXPECTED_RECTANGLE = (
    2_419,
    "67b373dc22ac870d",
    "e9ee49675d5345b06157e64a060506be3f6dd6a94835cc92f6eb5138f346cffb",
)
EXPECTED_POST4277 = (
    172_322,
    "30b2a7e597ac548c",
    "7228658eae4067db4bbcceb6c9b1ccebf1bd3e6f128e202ea184854acc53f309",
)
EXPECTED_COMMON = (
    148_063,
    "465fb756a183167e",
    "412b942759ed7afde1bffaaeabc9e6ae31b8b8bc25f26d73c71712523a057aaa",
)
EXPECTED_PARTITIONS = {
    "complement": (
        24_259,
        "78b212469c336f37",
        "77f3d21d127bb5e21f583556314f74032271e3a4903696415885308b363624ef",
    ),
    "carrier": (
        1_068,
        "c61906c5ef460bb8",
        "ffbb4e8030dae146267a623009603df4cd631468d80a5fa0a8965ec7a5232b93",
    ),
    "overlap": (
        1_032,
        "edb538ac3fb4c990",
        "c3357b39a4b356b07c51654114f3b6c2e8e2455fe5ded630a0526686982f7a69",
    ),
    "common_not_carrier": (
        147_031,
        "cacf2bef3f0b32bb",
        "922a4d172d00e5d3eb3424cb0ca0356f2cf1c32986c308b110b3179b41e92f11",
    ),
    "carrier_not_common": (
        36,
        "a8027022c6436919",
        "b14f2cb2ff92ea4f14676e3b6cf6b63fcff1b3378d2ece0461a9a4c89bccfd0e",
    ),
    "total_deletion": (
        148_099,
        "0308f3f6b8d7a23e",
        "ed993db4e2b53a416e6031696c4a0b01dfd3f9b9e47fb2607f5f191059752836",
    ),
    "post_carrier": (
        171_254,
        "197cdaf61aff1315",
        "25b9939804be1707d6ebfbd4841181b0ac9252a2a4a3252c467794022e7df4ff",
    ),
    "final": (
        24_223,
        "80ec0687d8c7dba7",
        "75a3c7616c982538363c7801ed2dab3fe9aa775ab601f7a7119dd9fb5d301552",
    ),
}
EXPECTED_COMPLEMENT_LAYERS = (
    639,
    "244707d5e67a299e",
    "0704be8ec0e23cf9f9f90c373476d666e329f2b4ba1b019182c754a72673ad87",
)
EXPECTED_OVERLAP_LAYERS = (
    639,
    "d19efeaf22f776f9",
    "b52c468020501e8dec4622dbecfd9b1afdcdd543d0e265c4fc5d4525e742d3e3",
)
EXPECTED_CARRIER_LAYERS = {
    670: 2,
    669: 168,
    668: 170,
    667: 176,
    666: 179,
    665: 183,
    664: 190,
}


Edge = tuple[int, int]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_words(rows: Iterable[Sequence[int]]) -> str:
    state = OFFSET
    for row in rows:
        for word in row:
            require(word >= 0, "FNV input must be nonnegative")
            for shift in range(0, 64, 8):
                state ^= (word >> shift) & 0xFF
                state = state * PRIME & MASK64
    return f"{state:016x}"


def edge_bytes(edges: Sequence[Edge]) -> bytes:
    return b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)


def fingerprint(edges: Sequence[Edge]) -> tuple[int, str, str]:
    raw = edge_bytes(edges)
    return len(edges), fnv_words(edges), hashlib.sha256(raw).hexdigest()


def read_edges(path: Path) -> list[Edge]:
    edges: list[Edge] = []
    for number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        fields = line.split(",")
        require(len(fields) == 2, f"bad edge row {number} in {path}")
        q, r = map(int, fields)
        require(0 < q < r, f"bad edge coordinates at row {number} in {path}")
        edges.append((q, r))
    require(edges == sorted(set(edges)), f"edge order/distinctness changed in {path}")
    return edges


def ordered_subset(universe: Sequence[Edge], members: set[Edge]) -> list[Edge]:
    require(members <= set(universe), "requested set is not contained in universe")
    return [edge for edge in universe if edge in members]


def show(label: str, edges: Sequence[Edge]) -> None:
    count, fnv_value, sha = fingerprint(edges)
    print(f"{label} COUNT {count} FNV {fnv_value} SHA256 {sha}")


def show_top(label: str, edges: Sequence[Edge]) -> None:
    require(edges, f"cannot report top of empty set {label}")
    endpoint = max(r for _, r in edges)
    top = [edge for edge in edges if edge[1] == endpoint]
    count, fnv_value, sha = fingerprint(top)
    text = ";".join(f"{q},{r}" for q, r in top)
    if len(top) <= 20:
        detail = f"EDGES {text}"
    else:
        detail = f"FIRST {top[0][0]},{endpoint} LAST {top[-1][0]},{endpoint}"
    print(
        f"{label}_TOP ENDPOINT {endpoint} COUNT {count} "
        f"FNV {fnv_value} SHA256 {sha} {detail}"
    )


def require_text(path: Path, patterns: Sequence[str]) -> str:
    text = path.read_text(encoding="ascii")
    for pattern in patterns:
        require(re.search(pattern, text, flags=re.MULTILINE) is not None,
                f"missing transcript evidence {pattern!r} in {path}")
    return text


def verify_carrier_transcripts(args: argparse.Namespace) -> None:
    # The detached literal replay closes the two retained endpoint-670 rows.
    require_text(args.literal670, [
        r"^PAIR 256,670 AUGMENTED_ACTIVE 2569 BODIES 14307150 FAILURES 0 CHECKS 613731713 MAX_CHECKS 2497$",
        r"^PAIR 384,670 AUGMENTED_ACTIVE 4218 BODIES 14307150 FAILURES 0 CHECKS 541846700 MAX_CHECKS 4205$",
        r"^ENDPOINT664_WITNESS c0c125 ACTIVE 1 DISJOINT 2 FINAL_CARRIER 8951 FINAL_FNV 188f82ab9dd1695a FINAL_ACTIVE 4528 BODIES 14307150 FAILURES 0 CHECKS 536018179 MAX_CHECKS 4528$",
        r"^VERDICT PASS DETACHED_LITERAL_WALL_DESCENT_THROUGH_664$",
    ])

    # At odd boundary layers the sole resistant row is closed by the next
    # append; the intervening even layer is already closed.  The final 664
    # scan leaves exactly (256,664).
    require_text(args.scan669, [
        r"^AUGMENTED_CARRIER 8945 FNV 3212efa05dd18c00$",
        r"^LAYER 669 ROWS 168 RESISTANT 1 ACTIVE_SUM 1455317 CHECKS 71524673001 MAX_CHECKS 3359 ROW_LEDGER 1b12eb2fb263db7e$",
        r"^RESISTANT PAIR 256,669 ACTIVE 3359 FAILURES 4 FIRST_FAILURE f047001 ",
    ])
    require_text(args.scan668, [
        r"^AUGMENTED_CARRIER 8947 FNV b246d6377503ae7$",
        r"^REPAIRED PAIR 256,669 ACTIVE 3361 FAILURES 0 CHECKS 569666955 MAX_CHECKS 3361$",
        r"^LAYER 668 ROWS 170 RESISTANT 0 ACTIVE_SUM 1502443 CHECKS 71368851159 MAX_CHECKS 6511 ROW_LEDGER 5b5c447f482c7d7$",
    ])
    require_text(args.scan667, [
        r"^AUGMENTED_CARRIER 8947 FNV b246d6377503ae7$",
        r"^LAYER 667 ROWS 176 RESISTANT 1 ACTIVE_SUM 1536773 CHECKS 74521089350 MAX_CHECKS 4199 ROW_LEDGER 94237267aec84f12$",
        r"^RESISTANT PAIR 256,667 ACTIVE 3982 FAILURES 18 FIRST_FAILURE 514a488 ",
    ])
    require_text(args.scan666, [
        r"^AUGMENTED_CARRIER 8949 FNV 1496ebecca72684b$",
        r"^REPAIRED PAIR 256,667 ACTIVE 3984 FAILURES 0 CHECKS 572210490 MAX_CHECKS 3984$",
        r"^LAYER 666 ROWS 179 RESISTANT 0 ACTIVE_SUM 1583997 CHECKS 75099164361 MAX_CHECKS 5344 ROW_LEDGER 6b834428d86c0231$",
    ])
    require_text(args.scan665, [
        r"^AUGMENTED_CARRIER 8949 FNV 1496ebecca72684b$",
        r"^LAYER 665 ROWS 183 RESISTANT 1 ACTIVE_SUM 1565443 CHECKS 78347641341 MAX_CHECKS 5907 ROW_LEDGER 22942b499df769d0$",
        r"^RESISTANT PAIR 520,665 ACTIVE 5156 FAILURES 1 FIRST_FAILURE 151a3400 ",
    ])
    require_text(args.scan664, [
        r"^AUGMENTED_CARRIER 8950 FNV f07022300e266930$",
        r"^REPAIRED PAIR 520,665 ACTIVE 5157 FAILURES 0 CHECKS 531031657 MAX_CHECKS 5157$",
        r"^LAYER 664 ROWS 190 RESISTANT 1 ACTIVE_SUM 1660345 CHECKS 80180889034 MAX_CHECKS 6307 ROW_LEDGER 172cd995dd68cdd6$",
        r"^RESISTANT PAIR 256,664 ACTIVE 4527 FAILURES 2 FIRST_FAILURE d143408 ",
    ])
    require_text(args.scan664_final, [
        r"^AUGMENTED_CARRIER 8951 FNV 188f82ab9dd1695a$",
        r"^REPAIRED PAIR 256,664 ACTIVE 4528 FAILURES 0 CHECKS 536018179 MAX_CHECKS 4528$",
        r"^LAYER 664 ROWS 190 RESISTANT 0 ACTIVE_SUM 1660533 CHECKS 80180889036 MAX_CHECKS 6307 ROW_LEDGER 2973e6c6d40d5af8$",
        r"^VERDICT RANGE_CLOSED EXACT_PRIMARY_ENDPOINT_DESCENT$",
    ])


def parse_claimed_fingerprints(path: Path) -> dict[str, tuple[int, str, str]]:
    claims: dict[str, tuple[int, str, str]] = {}
    pattern = re.compile(
        r"^([A-Z0-9_]+) COUNT (\d+) FNV ([0-9a-f]{16}) SHA256 ([0-9a-f]{64})$"
    )
    for line in path.read_text(encoding="ascii").splitlines():
        match = pattern.fullmatch(line)
        if match:
            claims[match.group(1)] = (
                int(match.group(2)), match.group(3), match.group(4)
            )
    return claims


def csv_bytes(header: Sequence[str], rows: Sequence[Sequence[int]]) -> bytes:
    lines = [",".join(header)]
    lines.extend(",".join(str(value) for value in row) for row in rows)
    return ("\n".join(lines) + "\n").encode("ascii")


def write_exact(path: Path | None, raw: bytes) -> None:
    if path is not None:
        path.write_bytes(raw)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--post4271", type=Path, required=True)
    parser.add_argument("--common", type=Path, required=True)
    parser.add_argument("--proof-graph", type=Path, required=True)
    parser.add_argument("--literal670", type=Path, required=True)
    for endpoint in range(664, 670):
        parser.add_argument(f"--scan{endpoint}", type=Path, required=True)
    parser.add_argument("--scan664-final", dest="scan664_final",
                        type=Path, required=True)
    parser.add_argument("--emit-complement", type=Path)
    parser.add_argument("--emit-complement-layers", type=Path)
    parser.add_argument("--emit-overlap-layers", type=Path)
    parser.add_argument("--emit-carrier", type=Path)
    parser.add_argument("--emit-overlap", type=Path)
    parser.add_argument("--emit-common-not-carrier", type=Path)
    parser.add_argument("--emit-carrier-not-common", type=Path)
    parser.add_argument("--emit-total-deletion", type=Path)
    parser.add_argument("--emit-final", type=Path)
    args = parser.parse_args()

    post4271 = read_edges(args.post4271)
    require(fingerprint(post4271) == EXPECTED_POST4271,
            "post-THM-4271 ledger changed")
    post4276 = [edge for edge in post4271
                if edge[1] < 670 or edge in BOUNDARY_4276]
    require(fingerprint(post4276) == EXPECTED_POST4276,
            "post-THM-4276 ledger changed")
    rectangle = [edge for edge in post4276
                 if 450 <= edge[0] <= 499 and 600 <= edge[1] <= 650]
    require(fingerprint(rectangle) == EXPECTED_RECTANGLE,
            "THM-4277 rectangle ledger changed")
    rectangle_set = set(rectangle)
    post4277 = [edge for edge in post4276 if edge not in rectangle_set]
    require(fingerprint(post4277) == EXPECTED_POST4277,
            "post-THM-4277 ledger changed")

    common = read_edges(args.common)
    require(fingerprint(common) == EXPECTED_COMMON,
            "primary common-active ledger changed")
    universe_set = set(post4277)
    common_set = set(common)
    require(common_set <= universe_set,
            "primary common-active set escapes post-THM-4277 universe")

    verify_carrier_transcripts(args)
    carrier = [edge for edge in post4277
               if edge[1] >= 664 and edge not in BOUNDARY_4281]
    require(Counter(r for _, r in carrier) == EXPECTED_CARRIER_LAYERS,
            "reconstructed carrier-closed layer census changed")
    carrier_set = set(carrier)

    complement = ordered_subset(post4277, universe_set - common_set)
    overlap = ordered_subset(post4277, common_set & carrier_set)
    common_not_carrier = ordered_subset(post4277, common_set - carrier_set)
    carrier_not_common = ordered_subset(post4277, carrier_set - common_set)
    total_deletion = ordered_subset(post4277, common_set | carrier_set)
    final = ordered_subset(post4277, universe_set - (common_set | carrier_set))
    post_carrier = ordered_subset(post4277, universe_set - carrier_set)

    require(len(overlap) + len(common_not_carrier) == len(common),
            "common partition count failed")
    require(len(overlap) + len(carrier_not_common) == len(carrier),
            "carrier partition count failed")
    require(len(total_deletion) + len(final) == len(post4277),
            "deletion/final partition count failed")
    require(len(complement) == len(carrier_not_common) + len(final),
            "common complement partition count failed")

    frozen_partitions = {
        "complement": complement,
        "carrier": carrier,
        "overlap": overlap,
        "common_not_carrier": common_not_carrier,
        "carrier_not_common": carrier_not_common,
        "total_deletion": total_deletion,
        "post_carrier": post_carrier,
        "final": final,
    }
    for label, edges in frozen_partitions.items():
        require(fingerprint(edges) == EXPECTED_PARTITIONS[label],
                f"frozen set fingerprint changed: {label}")

    claimed = parse_claimed_fingerprints(args.proof_graph)
    expected_claims = {
        "POST_THM4271": fingerprint(post4271),
        "POST_THM4276": fingerprint(post4276),
        "THM4277_RECTANGLE": fingerprint(rectangle),
        "POST_THM4277": fingerprint(post4277),
        "THM4281_COMMON_ACTIVE": fingerprint(common),
        "POST_COMMON_COMPLEMENT": fingerprint(complement),
        "THM4281_CARRIER_RAW": fingerprint(carrier),
        "THM4281_CARRIER_COMMON_OVERLAP": fingerprint(overlap),
        "THM4281_CARRIER_NOVEL": fingerprint(carrier_not_common),
        "THM4281_TOTAL_NOVEL": fingerprint(total_deletion),
        "FINAL_AFTER_THM4281": fingerprint(final),
    }
    for label, expected in expected_claims.items():
        require(claimed.get(label) == expected,
                f"proof-graph claim mismatch or missing: {label}")

    residual_layers = Counter(r for _, r in post4277)
    common_layers = Counter(r for _, r in common)
    complement_layers = Counter(r for _, r in complement)
    carrier_layers = Counter(r for _, r in carrier)
    overlap_layers = Counter(r for _, r in overlap)
    common_not_carrier_layers = Counter(r for _, r in common_not_carrier)
    carrier_not_common_layers = Counter(r for _, r in carrier_not_common)
    deletion_layers = Counter(r for _, r in total_deletion)
    final_layers = Counter(r for _, r in final)
    endpoints = sorted(residual_layers)

    complement_rows = [
        (r, residual_layers[r], common_layers[r], complement_layers[r])
        for r in endpoints
    ]
    overlap_rows = [
        (
            r,
            residual_layers[r],
            common_layers[r],
            complement_layers[r],
            carrier_layers[r],
            overlap_layers[r],
            common_not_carrier_layers[r],
            carrier_not_common_layers[r],
            deletion_layers[r],
            final_layers[r],
        )
        for r in endpoints
    ]
    complement_raw = csv_bytes(
        ("endpoint", "post_thm4277_pairs", "common_active_pairs",
         "not_common_active_pairs"),
        complement_rows,
    )
    overlap_raw = csv_bytes(
        ("endpoint", "post_thm4277_pairs", "common_active_pairs",
         "not_common_active_pairs", "carrier_closed_pairs",
         "common_and_carrier_pairs", "common_not_carrier_pairs",
         "carrier_not_common_pairs", "total_deletion_pairs",
         "final_residual_pairs"),
        overlap_rows,
    )
    require(
        (len(complement_rows), fnv_words(complement_rows),
         hashlib.sha256(complement_raw).hexdigest()) == EXPECTED_COMPLEMENT_LAYERS,
        "complement layer ledger changed",
    )
    require(
        (len(overlap_rows), fnv_words(overlap_rows),
         hashlib.sha256(overlap_raw).hexdigest()) == EXPECTED_OVERLAP_LAYERS,
        "overlap layer ledger changed",
    )

    write_exact(args.emit_complement, edge_bytes(complement))
    write_exact(args.emit_complement_layers, complement_raw)
    write_exact(args.emit_overlap_layers, overlap_raw)
    write_exact(args.emit_carrier, edge_bytes(carrier))
    write_exact(args.emit_overlap, edge_bytes(overlap))
    write_exact(args.emit_common_not_carrier, edge_bytes(common_not_carrier))
    write_exact(args.emit_carrier_not_common, edge_bytes(carrier_not_common))
    write_exact(args.emit_total_deletion, edge_bytes(total_deletion))
    write_exact(args.emit_final, edge_bytes(final))

    print("JOINT421_COMMON_COMPLEMENT_CARRIER_OVERLAP_VERIFY_V1")
    show("POST_THM4271", post4271)
    show("POST_THM4276", post4276)
    show("THM4277_RECTANGLE", rectangle)
    show("POST_THM4277", post4277)
    show("COMMON_ACTIVE_A", common)
    show("NOT_COMMON_ACTIVE_U_MINUS_A", complement)
    show("CARRIER_CLOSED_C", carrier)
    show("COMMON_AND_CARRIER_A_INTERSECT_C", overlap)
    show("COMMON_NOT_CARRIER_A_MINUS_C", common_not_carrier)
    show("CARRIER_NOT_COMMON_C_MINUS_A", carrier_not_common)
    show("TOTAL_DELETION_D_EQ_A_UNION_C", total_deletion)
    show("POST_CARRIER_U_MINUS_C", post_carrier)
    show("FINAL_U_MINUS_D", final)
    show("COMMON_INTERSECT_DELETION_A_INTERSECT_D", common)
    show("COMMON_MINUS_DELETION_A_MINUS_D", [])
    show("DELETION_MINUS_COMMON_D_MINUS_A", carrier_not_common)
    show_top("NOT_COMMON_ACTIVE_U_MINUS_A", complement)
    show_top("COMMON_ACTIVE_A", common)
    show_top("CARRIER_CLOSED_C", carrier)
    show_top("COMMON_AND_CARRIER_A_INTERSECT_C", overlap)
    show_top("COMMON_NOT_CARRIER_A_MINUS_C", common_not_carrier)
    show_top("CARRIER_NOT_COMMON_C_MINUS_A", carrier_not_common)
    show_top("TOTAL_DELETION_D_EQ_A_UNION_C", total_deletion)
    show_top("FINAL_U_MINUS_D", final)
    print(
        f"COMPLEMENT_LAYER_CENSUS ROWS {len(complement_rows)} "
        f"FNV {fnv_words(complement_rows)} "
        f"SHA256 {hashlib.sha256(complement_raw).hexdigest()}"
    )
    print(
        f"OVERLAP_LAYER_CENSUS ROWS {len(overlap_rows)} "
        f"FNV {fnv_words(overlap_rows)} "
        f"SHA256 {hashlib.sha256(overlap_raw).hexdigest()}"
    )
    for endpoint in (670, 669, 668, 667, 666, 665, 664, 663, 650,
                     600, 500, 400, 300, 200, 100, 25, 2):
        if endpoint in residual_layers:
            print(
                f"LAYER {endpoint} U {residual_layers[endpoint]} "
                f"A {common_layers[endpoint]} U_MINUS_A {complement_layers[endpoint]} "
                f"C {carrier_layers[endpoint]} A_INTERSECT_C {overlap_layers[endpoint]} "
                f"A_MINUS_C {common_not_carrier_layers[endpoint]} "
                f"C_MINUS_A {carrier_not_common_layers[endpoint]} "
                f"D {deletion_layers[endpoint]} FINAL {final_layers[endpoint]}"
            )
    print("DEFINITIONS U=POST_THM4277 A=COMMON_ACTIVE C=CARRIER_CLOSED "
          "D=A_UNION_C FINAL=U_MINUS_D")
    print("CARRIER_EVIDENCE R670_LITERAL R669_TO_R664_EXACT_ENDPOINT_SCANS")
    print("VERDICT PASS EXACT_SET_PARTITION_AND_LAYER_REPLAY")


if __name__ == "__main__":
    main()
