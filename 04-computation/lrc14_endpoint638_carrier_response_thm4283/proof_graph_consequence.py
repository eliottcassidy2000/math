#!/usr/bin/env python3
"""Exact set/carrier consequence for THM-4283."""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
REPAIR = 0x014C9084
TOP638_Q = {100, 282, 294, 366, 416, 420, 512, 520}


def require(condition: bool, message: str) -> None:
    """Optimization-safe audit gate."""
    if not condition:
        raise RuntimeError(message)


def fnv_words(words: list[int]) -> int:
    state = OFFSET
    for word in words:
        for byte in int(word).to_bytes(8, "little", signed=False):
            state ^= byte
            state = (state * PRIME) & MASK64
    return state


def fnv_rows(rows: list[tuple[int, int]]) -> int:
    return fnv_words([coordinate for row in rows for coordinate in row])


def raw_rows(rows: list[tuple[int, int]]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")


def sha_rows(rows: list[tuple[int, int]]) -> str:
    return hashlib.sha256(raw_rows(rows)).hexdigest()


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for number, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        pieces = line.split(",")
        require(len(pieces) == 2, f"bad row at {path}:{number}")
        rows.append((int(pieces[0]), int(pieces[1])))
    require(len(rows) == len(set(rows)), f"duplicate row in {path}")
    return rows


def read_masks(path: Path) -> list[int]:
    return [int(token, 16) for token in path.read_text(encoding="ascii").split()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--post4282", type=Path, required=True)
    parser.add_argument("--closed", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--additions", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.post4282)
    closed = read_rows(args.closed)
    require(universe == sorted(universe), "post-THM4282 rows lost lexicographic order")
    require(closed == sorted(closed), "closed rows lost lexicographic order")
    require(
        len(universe) == 23_373
        and fnv_rows(universe) == 0xC6AB0AE49EE32273
        and sha_rows(universe)
        == "c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3",
        "post-THM4282 universe identity changed",
    )

    expected_closed = [
        (q, r)
        for q, r in universe
        if r >= 639 or (r == 638 and q in TOP638_Q)
    ]
    require(closed == expected_closed, "closed63 definition changed")
    require(
        len(closed) == 63
        and fnv_rows(closed) == 0x3ED7002531EFFA14
        and sha_rows(closed)
        == "2c48a1b46f3ac2fce4b9fdaf426373288bd2ee15649c228de2a8e4d9fc6f46fd",
        "closed63 identity changed",
    )

    closed_set = set(closed)
    residual = [row for row in universe if row not in closed_set]
    require(
        len(residual) == 23_310
        and fnv_rows(residual) == 0x89EA6A588DA4BA0A
        and sha_rows(residual)
        == "75c5881e9de9622c1627de4da07dca1df0be82f3366ef1e7eb78e36ff0f71d14",
        "post-THM4283 residual identity changed",
    )
    maximum = max(r for _, r in residual)
    top = [row for row in residual if row[1] == maximum]
    require(maximum == 638 and top == [(256, 638)], "top singleton changed")

    base = read_masks(args.base)
    additions = read_masks(args.additions)
    require(
        len(base) == 8_951
        and fnv_words(base) == 0x188F82AB9DD1695A,
        "base carrier identity changed",
    )
    require(
        len(additions) == 45
        and fnv_words(additions) == 0xEC083B65CC8C34E3,
        "addition ledger identity changed",
    )
    prior = base + additions
    require(
        len(set(prior)) == 8_996
        and fnv_words(prior) == 0xFD899660F14B311C,
        "inherited carrier identity changed",
    )
    require(REPAIR not in set(prior) and REPAIR.bit_count() == 8, "repair changed")
    repaired = prior + [REPAIR]
    require(
        fnv_words(repaired) == 0x8E1860A25D0FCF87,
        "repaired carrier identity changed",
    )

    pool = [
        8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
        80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
        168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
    ]
    support = [pool[index] for index in range(30) if REPAIR >> index & 1]
    require(
        support == [15, 42, 85, 120, 143, 145, 176, 193],
        "repair support changed",
    )

    layers = Counter(r for _, r in closed)
    print("THM4283_PROOF_GRAPH_CONSEQUENCE_V1")
    print("PRIOR_CARRIER", len(prior), f"{fnv_words(prior):016x}")
    print("REPAIR", f"{REPAIR:08x}", support)
    print("REPAIRED_CARRIER", len(repaired), f"{fnv_words(repaired):016x}")
    print(
        "CLOSED63",
        len(closed),
        f"FNV={fnv_rows(closed):016x}",
        f"SHA256={sha_rows(closed)}",
    )
    print("LAYERS", " ".join(f"{r}:{layers[r]}" for r in sorted(layers, reverse=True)))
    print(
        "RESIDUAL",
        len(residual),
        f"FNV={fnv_rows(residual):016x}",
        f"SHA256={sha_rows(residual)}",
    )
    print("TOP", maximum, top)
    print("VERDICT PASS EXACT SET DIFFERENCE; LRC14 OPEN")


if __name__ == "__main__":
    main()
