#!/usr/bin/env python3
"""Exact THM-4283 carrier/surgery proof-graph consequence.

The script treats the common-family union and the descending carrier scan as
separate certificate nodes, audits their exact overlap, and writes the union
and complementary residual in canonical lexicographic serialization.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1

EXPECTED_RESIDUAL_COUNT = 23_373
EXPECTED_RESIDUAL_FNV = 0xC6AB0AE49EE32273
EXPECTED_RESIDUAL_SHA = "c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3"
EXPECTED_COMMON_COUNT = 640
EXPECTED_COMMON_FNV = 0x45E9ECDF240E6417
EXPECTED_COMMON_SHA = "3246d76e82e9e19d07e5851810da3107ad8fe98a1dfbd087edb2b9c5d8b27fa1"

EXPECTED_LAYERS = {
    644: (7, 0x195E5D7D703B7D4C, 2, 0),
    643: (14, 0xA14EB5B1EE96EDB4, 0, 0),
    642: (9, 0x32318BDCCA33A0F4, 0, 0),
    641: (7, 0x399BD7D3C8D4B81D, 0, 0),
    640: (14, 0xB548BEB345C96F2, 0, 0),
    639: (4, 0xFC85AA6B3AB0D13F, 0, 0),
    638: (9, 0xD03393624CA9D75F, 40, 40),
}


def fnv_words(words: list[int] | tuple[int, ...]) -> int:
    state = FNV_OFFSET
    for word in words:
        if word < 0:
            raise ValueError("negative FNV word")
        for shift in range(0, 64, 8):
            state ^= (word >> shift) & 0xFF
            state = (state * FNV_PRIME) & MASK64
    return state


def pair_fnv(rows: list[tuple[int, int]]) -> int:
    return fnv_words([coordinate for row in rows for coordinate in row])


def raw_rows(rows: list[tuple[int, int]]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")


def sha_rows(rows: list[tuple[int, int]]) -> str:
    return hashlib.sha256(raw_rows(rows)).hexdigest()


def read_pairs(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        if not line:
            raise ValueError(f"{path}:{number}: blank row")
        fields = line.split(",")
        if len(fields) != 2:
            raise ValueError(f"{path}:{number}: malformed pair")
        pair = (int(fields[0]), int(fields[1]))
        if not 0 < pair[0] < pair[1]:
            raise ValueError(f"{path}:{number}: invalid ordered pair")
        if rows and not rows[-1] < pair:
            raise ValueError(f"{path}:{number}: ledger not strictly ordered")
        rows.append(pair)
    return rows


def write_pairs(path: Path, rows: list[tuple[int, int]]) -> None:
    path.write_bytes(raw_rows(rows))


def tagged(tokens: list[str], tag: str) -> str:
    try:
        return tokens[tokens.index(tag) + 1]
    except (ValueError, IndexError) as error:
        raise ValueError(f"missing transcript tag {tag}") from error


def parse_scan(path: Path, residual: list[tuple[int, int]]) -> tuple[
        list[tuple[int, int]], list[tuple[int, int]]]:
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0] != "THM4283_ENDPOINT_TOP_BAND_SCAN_V1":
        raise ValueError("carrier scan banner changed")
    header = lines[1].split()
    if (
        int(tagged(header, "RESIDUAL")) != EXPECTED_RESIDUAL_COUNT
        or int(tagged(header, "FNV"), 16) != EXPECTED_RESIDUAL_FNV
        or int(tagged(header, "REQUESTED_LOWER")) != 630
        or int(tagged(header, "SELECTED")) != 108
        or int(tagged(header, "PAIR_FNV"), 16) != 0x734EB66EAD03881A
        or int(tagged(header, "BASE8996_FNV"), 16) != 0xFD899660F14B311C
        or int(tagged(header, "NESTED8997_FNV"), 16) != 0x8E1860A25D0FCF87
    ):
        raise ValueError("carrier scan header identity changed")

    scanned: list[tuple[int, int]] = []
    layers: dict[int, tuple[int, int, int, int]] = {}
    hostile_638: list[tuple[int, int]] = []
    saw_exact_644 = False
    saw_exact_638 = False
    for line in lines[2:]:
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "PAIR":
            q_text, r_text = tokens[1].split(",")
            pair = (int(q_text), int(r_text))
            scanned.append(pair)
            nested_failures = int(tagged(tokens, "NESTED_FAILURES"))
            if pair[1] >= 639 and nested_failures != 0:
                raise ValueError(f"nested carrier failure inside closed band: {pair}")
            if pair[1] == 638:
                hostile_638.append(pair)
            if pair == (256, 644):
                saw_exact_644 = (
                    int(tagged(tokens, "BASE_FAILURES")) == 2
                    and nested_failures == 0
                    and int(tagged(tokens, "RESPONSE_CLASSES")) == 4
                    and int(tagged(tokens, "FULL_RESPONDERS")) == 367
                    and tagged(tokens, "LEAST_FULL") == "014c9084"
                    and int(tagged(tokens, "EXACT_MINIMUM")) == 1
                )
            if pair == (256, 638):
                saw_exact_638 = (
                    int(tagged(tokens, "BASE_FAILURES")) == 40
                    and nested_failures == 40
                    and int(tagged(tokens, "RESPONSE_CLASSES")) == 315
                    and int(tagged(tokens, "FULL_RESPONDERS")) == 0
                    and tagged(tokens, "LEAST_FULL") == "00000000"
                    and int(tagged(tokens, "EXACT_MINIMUM")) == 9
                )
        elif tokens[0] == "LAYER":
            endpoint = int(tokens[1])
            layers[endpoint] = (
                int(tagged(tokens, "ROWS")),
                int(tagged(tokens, "FNV"), 16),
                int(tagged(tokens, "BASE_FAILURES")),
                int(tagged(tokens, "NESTED_FAILURES")),
            )
        elif tokens[0] == "PAIR_LEDGER_FNV":
            if (
                int(tokens[1], 16) != 0xAB0A8C09B7AAFF5
                or int(tagged(tokens, "TOTAL_BASE_FAILURES")) != 42
                or int(tagged(tokens, "TOTAL_NESTED_FAILURES")) != 40
                or int(tagged(tokens, "STOPPED_AT_FAILURE")) != 1
            ):
                raise ValueError("carrier scan final ledger changed")
        elif tokens[0] == "VERDICT" and line != "VERDICT PASS EXACT_DESCENDING_SCAN":
            raise ValueError("carrier scan verdict changed")

    if layers != EXPECTED_LAYERS or not saw_exact_644 or not saw_exact_638:
        raise ValueError("carrier layer or exact-response controls changed")
    expected_scanned = sorted(
        (pair for pair in residual if pair[1] >= 638),
        key=lambda pair: (-pair[1], pair[0]),
    )
    if scanned != expected_scanned:
        raise ValueError("carrier scan did not exhaust the exact 638+ band")
    carrier = sorted(pair for pair in scanned if pair[1] >= 639)
    hostile = sorted(hostile_638)
    if len(carrier) != 55 or len(hostile) != 9:
        raise ValueError("carrier/hostile layer counts changed")
    return carrier, hostile


def audit_boundary_witness(path: Path) -> None:
    lines = path.read_text(encoding="ascii").splitlines()
    if len(lines) != 6 or lines[0] != "THM4283_ENDPOINT638_EXACT_RESPONSE_WITNESS_V1":
        raise ValueError("endpoint-638 witness banner or row count changed")
    pair = lines[1].split()
    if (
        pair[1] != "256,638"
        or int(tagged(pair, "BASE8997")) != 8997
        or int(tagged(pair, "FNV"), 16) != 0x8E1860A25D0FCF87
        or int(tagged(pair, "FAILURES")) != 40
        or int(tagged(pair, "FAILURE_FNV"), 16) != 0x917D107C4536EFC9
        or int(tagged(pair, "RESPONSE_CLASSES")) != 315
        or int(tagged(pair, "FULL_RESPONDERS")) != 0
        or int(tagged(pair, "EXACT_MINIMUM")) != 9
    ):
        raise ValueError("endpoint-638 witness boundary changed")
    witnesses = lines[2].split()
    if (
        int(tagged(witnesses, "WITNESSES")) != 9
        or int(tagged(witnesses, "FNV"), 16) != 0x02B936529030E4BC
        or witnesses[witnesses.index("MASKS") + 1:] != [
            "02203226", "081e1084", "08a89440", "180a8281", "18261042",
            "18a0d040", "1a82a200", "202a9440", "280a0a88",
        ]
    ):
        raise ValueError("endpoint-638 explicit witness changed")
    replay = lines[3].split()
    if (
        int(tagged(replay, "REPAIRED9006_FNV"), 16) != 0xFDC1C57AE4DC1BB6
        or int(tagged(replay, "REPLAY_FAILURES")) != 0
        or lines[4] != "SCOPE LOCAL_PAIR_RESPONSE_ONLY_NO_LOWER_LAYER_CLAIM"
        or lines[5] != "VERDICT PASS EXACT_MINIMUM_NINE_AND_EXPLICIT_WITNESS"
    ):
        raise ValueError("endpoint-638 witness replay changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("residual", type=Path)
    parser.add_argument("common_union", type=Path)
    parser.add_argument("carrier_scan", type=Path)
    parser.add_argument("boundary_witness", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    residual = read_pairs(args.residual)
    if (
        len(residual) != EXPECTED_RESIDUAL_COUNT
        or pair_fnv(residual) != EXPECTED_RESIDUAL_FNV
        or sha_rows(residual) != EXPECTED_RESIDUAL_SHA
    ):
        raise ValueError("post-THM-4282 residual identity changed")
    common = read_pairs(args.common_union)
    if (
        len(common) != EXPECTED_COMMON_COUNT
        or pair_fnv(common) != EXPECTED_COMMON_FNV
        or sha_rows(common) != EXPECTED_COMMON_SHA
        or not set(common) <= set(residual)
    ):
        raise ValueError("all-127 common union identity or typing changed")

    carrier55, boundary638 = parse_scan(args.carrier_scan, residual)
    audit_boundary_witness(args.boundary_witness)
    carrier = sorted(carrier55 + boundary638)
    if len(carrier) != 64 or carrier != [pair for pair in residual if pair[1] >= 638]:
        raise ValueError("repaired carrier band is not the exact 638+ layer union")
    residual_set = set(residual)
    common_set = set(common)
    carrier_set = set(carrier)
    overlap = sorted(common_set & carrier_set)
    proof_union = sorted(common_set | carrier_set)
    final_residual = sorted(residual_set - set(proof_union))
    top_endpoint = max(r for _, r in final_residual)
    top = [(q, r) for q, r in final_residual if r == top_endpoint]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_pairs(args.output_dir / "carrier_band64.csv", carrier)
    write_pairs(args.output_dir / "carrier_common_overlap.csv", overlap)
    write_pairs(args.output_dir / "proof_union.csv", proof_union)
    write_pairs(args.output_dir / "final_residual.csv", final_residual)
    write_pairs(args.output_dir / "boundary_endpoint638.csv", boundary638)

    print("THM4283_PROOF_GRAPH_CONSEQUENCE_V1")
    print(
        f"INPUT residual={len(residual)} fnv={pair_fnv(residual):016x} "
        f"sha256={sha_rows(residual)}"
    )
    for label, rows in (
        ("SIGNATURE_FIBRE", common),
        ("CARRIER638_644", carrier),
        ("OVERLAP", overlap),
        ("UNION", proof_union),
        ("BOUNDARY638", boundary638),
        ("FINAL_RESIDUAL", final_residual),
    ):
        print(
            f"{label} count={len(rows)} fnv={pair_fnv(rows):016x} "
            f"sha256={sha_rows(rows)}"
        )
    print(
        "FINAL_BOUNDARY max_r=" + str(top_endpoint)
        + " top=" + ";".join(f"{q},{r}" for q, r in top)
    )
    print("SCOPE TYPED_SET_UNION_NO_SINGLE_GLOBAL_DECK_NO_LRC14_CONSEQUENCE")
    print("VERDICT PASS EXACT_CARRIER_SURGERY_PROOF_GRAPH")


if __name__ == "__main__":
    main()
