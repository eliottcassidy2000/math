#!/usr/bin/env python3
"""Exact proof-graph consequence of the index-367 singleton-fibre surgery."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from typing import Iterable, Sequence


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
TARGET_INDEX = 367
TARGET_MASK = 0x02188125
REPLACEMENT_MASK = 0x0A188803

EXPECTED_POST4281 = (
    24_223,
    "80ec0687d8c7dba7",
    "75a3c7616c982538363c7801ed2dab3fe9aa775ab601f7a7119dd9fb5d301552",
)
EXPECTED_SIGNATURE_SHA = "4f0e8da3fdab8bd5a0e14f3b4fa30050602025f63486aa35e0cf03374e6e3832"
EXPECTED_FIBRE = (
    26,
    "9758d0de94c51f75",
    "67f57d4ee960adfcd365c85f4baeb4a8ff2c393b211ea6ac561a5b0ad233300a",
)
EXPECTED_POST_FIBRE = (
    24_197,
    "c9746d48cbf3fc37",
    "8d59477a996503c338434dda77ae6ad897637d9d090101be118d179631c73e09",
)
EXPECTED_ORIGINAL_DECK = (421, "20d63dd42fe8150e")
EXPECTED_REBUILT_DECK = (421, "87b42cf8a2069177")


Edge = tuple[int, int]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv_words(rows: Iterable[Sequence[int]]) -> str:
    state = OFFSET
    for row in rows:
        for word in row:
            require(word >= 0, "FNV word must be nonnegative")
            for shift in range(0, 64, 8):
                state ^= (word >> shift) & 0xFF
                state = state * PRIME & MASK64
    return f"{state:016x}"


def edge_bytes(edges: Sequence[Edge]) -> bytes:
    return b"".join(f"{q},{r}\n".encode("ascii") for q, r in edges)


def edge_fingerprint(edges: Sequence[Edge]) -> tuple[int, str, str]:
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


def read_deck(path: Path) -> list[int]:
    deck = [int(token, 16) for token in path.read_text(encoding="ascii").split()]
    require(len(deck) == len(set(deck)), f"deck duplicates in {path}")
    require(all(mask.bit_count() == 8 and mask < (1 << 30) for mask in deck),
            f"invalid rank-eight mask in {path}")
    return deck


def read_signature_edges(path: Path) -> tuple[list[Edge], list[int]]:
    edges: list[Edge] = []
    signatures: list[int] = []
    with path.open(newline="", encoding="ascii") as handle:
        rows = csv.reader(handle)
        require(next(rows) == ["q", "r", "inactive_count", "w0", "w1",
                               "w2", "w3", "w4", "w5", "w6"],
                "signature header changed")
        for number, fields in enumerate(rows, 2):
            require(len(fields) == 10, f"bad signature row {number}")
            q, r, weight = map(int, fields[:3])
            words = [int(value, 16) for value in fields[3:]]
            signature = sum(word << (64 * index)
                            for index, word in enumerate(words))
            require(signature >> 421 == 0 and signature.bit_count() == weight,
                    f"signature weight/padding changed at row {number}")
            edges.append((q, r))
            signatures.append(signature)
    require(edges == sorted(set(edges)), "signature edge order changed")
    return edges, signatures


def show(label: str, edges: Sequence[Edge]) -> None:
    count, fnv_value, sha = edge_fingerprint(edges)
    print(f"{label} COUNT {count} FNV {fnv_value} SHA256 {sha}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("post4281", type=Path)
    parser.add_argument("signatures", type=Path)
    parser.add_argument("fibre", type=Path)
    parser.add_argument("original_deck", type=Path)
    parser.add_argument("rebuilt_deck", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    post4281 = read_edges(args.post4281)
    require(edge_fingerprint(post4281) == EXPECTED_POST4281,
            "post-THM-4281 residual changed")
    require(hashlib.sha256(args.signatures.read_bytes()).hexdigest() ==
            EXPECTED_SIGNATURE_SHA, "full signature CSV changed")
    signature_edges, signatures = read_signature_edges(args.signatures)
    require(signature_edges == post4281,
            "signature universe differs from post-THM-4281 residual")

    fibre = read_edges(args.fibre)
    require(edge_fingerprint(fibre) == EXPECTED_FIBRE,
            "singleton fibre ledger changed")
    derived_fibre = [edge for edge, signature in zip(signature_edges, signatures)
                     if signature == (1 << TARGET_INDEX)]
    require(derived_fibre == fibre,
            "supplied fibre is not the exact singleton index-367 fibre")

    original_deck = read_deck(args.original_deck)
    rebuilt_deck = read_deck(args.rebuilt_deck)
    require((len(original_deck), fnv_words((mask,) for mask in original_deck)) ==
            EXPECTED_ORIGINAL_DECK, "original deck changed")
    require((len(rebuilt_deck), fnv_words((mask,) for mask in rebuilt_deck)) ==
            EXPECTED_REBUILT_DECK, "rebuilt deck changed")
    require(original_deck[TARGET_INDEX] == TARGET_MASK,
            "target mask/index changed")
    require(rebuilt_deck == original_deck[:TARGET_INDEX] +
            original_deck[TARGET_INDEX + 1:] + [REPLACEMENT_MASK],
            "rebuilt deck is not the stated one-mask surgery")

    fibre_set = set(fibre)
    post_fibre = [edge for edge in post4281 if edge not in fibre_set]
    require(edge_fingerprint(post_fibre) == EXPECTED_POST_FIBRE,
            "post-fibre residual changed")
    args.output.write_bytes(edge_bytes(post_fibre))

    print("INDEX367_SINGLETON_FIBRE_PROOF_GRAPH_V1")
    show("POST_THM4281", post4281)
    print(f"SIGNATURE_CSV SHA256 {EXPECTED_SIGNATURE_SHA}")
    show("INDEX367_SINGLETON_FIBRE", fibre)
    print("DECK_SURGERY INDEX 367 REMOVE 02188125 APPEND 0a188803 "
          "ORIGINAL_FNV 20d63dd42fe8150e REBUILT_FNV 87b42cf8a2069177")
    show("POST_INDEX367_FIBRE", post_fibre)
    endpoint = max(r for _, r in post_fibre)
    top = [edge for edge in post_fibre if edge[1] == endpoint]
    show("POST_INDEX367_FIBRE_TOP", top)
    print("TOP_EDGES " + ";".join(f"{q},{r}" for q, r in top))
    print("VERDICT PASS EXACT_PROOF_GRAPH_CONSEQUENCE_ONLY")


if __name__ == "__main__":
    main()
