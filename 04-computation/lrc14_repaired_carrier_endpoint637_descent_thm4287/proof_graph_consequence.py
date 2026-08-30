#!/usr/bin/env python3
"""Exact typed proof-graph subtraction for THM-4287.

The two closure ledgers are independent certificate nodes.  This consumer
checks their types, their exact overlap, and the endpoint-637 boundary before
subtracting their set union from the pinned THM-4283 residual.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1

EXPECTED_INPUT_COUNT = 22_682
EXPECTED_INPUT_FNV = 0xF7563445F15EFEBF
EXPECTED_INPUT_SHA = "7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102"

EXPECTED_CARRIER = [(100, 637), (294, 637), (520, 637)]
EXPECTED_CARRIER_FNV = 0xF4B48B16A28826AD
EXPECTED_CARRIER_SHA = "e31697893a8a4dbeb4b22ea62d7ef2bfed7fd55bbe35138a1bfc6f460638601b"

EXPECTED_SIGNATURE_COUNT = 33
EXPECTED_SIGNATURE_FNV = 0xF11425C1E1B17094
EXPECTED_SIGNATURE_SHA = "5b5de4e502802287ff2c504d34cfc4de21c10adffd1694c4bfa1c9a1d735c85b"
EXPECTED_OVERLAP = [(520, 637)]
EXPECTED_OVERLAP_FNV = 0x3CEFEA79822890F4
EXPECTED_OVERLAP_SHA = "0d8ccfa43b9f309187ba867180e667eb9aac32341991ede0ea34a4f4b3b2fc12"
EXPECTED_UNION_COUNT = 35
EXPECTED_UNION_FNV = 0x0B5D8D9C28D7325D
EXPECTED_UNION_SHA = "1d36b9e6899496f93d44785d389523d4f3a85efe4eaadaf2e1297904116bbd78"
EXPECTED_FINAL_COUNT = 22_647
EXPECTED_FINAL_FNV = 0xDF5374D4ACA67677
EXPECTED_FINAL_SHA = "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317"
EXPECTED_FINAL_MAX = 636
EXPECTED_FINAL_TOP = [
    (100, 636), (256, 636), (294, 636),
    (338, 636), (372, 636), (384, 636),
]


def raw_rows(rows: list[tuple[int, int]]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")


def fnv_words(words: list[int]) -> int:
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


def pair_sha(rows: list[tuple[int, int]]) -> str:
    return hashlib.sha256(raw_rows(rows)).hexdigest()


def read_pairs(path: Path) -> list[tuple[int, int]]:
    try:
        data = path.read_bytes()
        text = data.decode("ascii")
    except UnicodeDecodeError as error:
        raise ValueError(f"{path}: ledger is not ASCII") from error

    rows: list[tuple[int, int]] = []
    for number, line in enumerate(text.splitlines(), 1):
        if not line:
            raise ValueError(f"{path}:{number}: blank row")
        fields = line.split(",")
        if len(fields) != 2 or any(not field.isdigit() for field in fields):
            raise ValueError(f"{path}:{number}: malformed pair")
        pair = (int(fields[0]), int(fields[1]))
        if not 0 < pair[0] < pair[1]:
            raise ValueError(f"{path}:{number}: invalid ordered pair")
        if rows and not rows[-1] < pair:
            raise ValueError(f"{path}:{number}: ledger not strictly ordered")
        rows.append(pair)

    if data != raw_rows(rows):
        raise ValueError(f"{path}: ledger serialization is not canonical")
    return rows


def identity(label: str, rows: list[tuple[int, int]]) -> str:
    return (
        f"{label} count={len(rows)} fnv={pair_fnv(rows):016x} "
        f"sha256={pair_sha(rows)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("residual_u4283", type=Path)
    parser.add_argument("carrier_endpoint637", type=Path)
    parser.add_argument("signature_fibre_closure", type=Path)
    parser.add_argument("output_residual", type=Path)
    args = parser.parse_args()

    residual = read_pairs(args.residual_u4283)
    carrier = read_pairs(args.carrier_endpoint637)
    signature = read_pairs(args.signature_fibre_closure)

    if (
        len(residual) != EXPECTED_INPUT_COUNT
        or pair_fnv(residual) != EXPECTED_INPUT_FNV
        or pair_sha(residual) != EXPECTED_INPUT_SHA
    ):
        raise ValueError("THM-4283 final residual identity changed")
    if (
        carrier != EXPECTED_CARRIER
        or pair_fnv(carrier) != EXPECTED_CARRIER_FNV
        or pair_sha(carrier) != EXPECTED_CARRIER_SHA
    ):
        raise ValueError("endpoint-637 carrier closure identity changed")
    if (
        len(signature) != EXPECTED_SIGNATURE_COUNT
        or pair_fnv(signature) != EXPECTED_SIGNATURE_FNV
        or pair_sha(signature) != EXPECTED_SIGNATURE_SHA
    ):
        raise ValueError("signature-fibre closure identity changed")

    residual_set = set(residual)
    carrier_set = set(carrier)
    signature_set = set(signature)
    if not carrier_set <= residual_set or not signature_set <= residual_set:
        raise ValueError("a closure ledger is not typed in the THM-4283 residual")

    input_max = max(r for _, r in residual)
    input_top = [pair for pair in residual if pair[1] == input_max]
    if input_max != 637 or input_top != carrier:
        raise ValueError("carrier closure is not exactly the maximal endpoint layer")

    overlap = sorted(carrier_set & signature_set)
    closure_union = sorted(carrier_set | signature_set)
    removed = sorted(residual_set & set(closure_union))
    final_residual = sorted(residual_set - set(closure_union))
    if (
        overlap != EXPECTED_OVERLAP
        or pair_fnv(overlap) != EXPECTED_OVERLAP_FNV
        or pair_sha(overlap) != EXPECTED_OVERLAP_SHA
    ):
        raise ValueError("closure-ledger overlap changed")
    if (
        len(closure_union) != EXPECTED_UNION_COUNT
        or pair_fnv(closure_union) != EXPECTED_UNION_FNV
        or pair_sha(closure_union) != EXPECTED_UNION_SHA
        or removed != closure_union
    ):
        raise ValueError("typed closure union/removal changed")
    if (
        len(final_residual) != EXPECTED_FINAL_COUNT
        or pair_fnv(final_residual) != EXPECTED_FINAL_FNV
        or pair_sha(final_residual) != EXPECTED_FINAL_SHA
    ):
        raise ValueError("post-THM-4287 residual count changed")

    final_max = max(r for _, r in final_residual)
    final_top = [pair for pair in final_residual if pair[1] == final_max]
    if final_max != EXPECTED_FINAL_MAX or final_top != EXPECTED_FINAL_TOP:
        raise ValueError("post-THM-4287 residual boundary changed")

    args.output_residual.parent.mkdir(parents=True, exist_ok=True)
    args.output_residual.write_bytes(raw_rows(final_residual))

    print("THM4287_PROOF_GRAPH_CONSEQUENCE_V1")
    for label, rows in (
        ("INPUT_U4283", residual),
        ("CARRIER_ENDPOINT637", carrier),
        ("SIGNATURE_FIBRE_CLOSURE", signature),
        ("OVERLAP", overlap),
        ("CLOSURE_UNION", closure_union),
        ("REMOVED", removed),
        ("FINAL_RESIDUAL", final_residual),
    ):
        print(identity(label, rows))
    print(
        "INPUT_BOUNDARY max_r=" + str(input_max)
        + " top=" + ";".join(f"{q},{r}" for q, r in input_top)
    )
    print(
        "FINAL_BOUNDARY max_r=" + str(final_max)
        + " top=" + ";".join(f"{q},{r}" for q, r in final_top)
    )
    print("SCOPE TYPED_RESIDUAL_SUBTRACTION_NO_PHYSICAL_ENTRY_NO_LRC14_CONSEQUENCE")
    print("VERDICT PASS EXACT_THM4287_PROOF_GRAPH_CONSEQUENCE")


if __name__ == "__main__":
    main()
