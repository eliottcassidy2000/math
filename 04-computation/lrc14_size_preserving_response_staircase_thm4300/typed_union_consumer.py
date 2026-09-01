#!/usr/bin/env python3
"""Typed consequence join for THM-4300.

This consumer unions row consequences only.  The endpoint carrier, the
index-297 deck, and every inherited THM-4295/4296 consumer remain distinct.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_pairs(path: Path) -> tuple[tuple[int, int], ...]:
    rows = tuple(
        tuple(map(int, line.split(",")))
        for line in path.read_text(encoding="ascii").splitlines()
        if line and line != "q,r"
    )
    require(rows == tuple(sorted(set(rows))), f"unordered/duplicate: {path}")
    require(all(0 < q < r for q, r in rows), f"invalid pair: {path}")
    return rows


def payload(rows: list[tuple[int, int]] | tuple[tuple[int, int], ...]) -> bytes:
    return b"".join(f"{q},{r}\n".encode("ascii") for q, r in rows)


def fnv(rows: list[tuple[int, int]] | tuple[tuple[int, int], ...]) -> str:
    state = FNV_OFFSET
    for q, r in rows:
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state ^= byte
                state = state * FNV_PRIME & MASK64
    return f"{state:016x}"


def identity(label: str, values: set[tuple[int, int]]) -> None:
    rows = sorted(values)
    maximum = max(r for _, r in rows) if rows else 0
    print(
        f"{label} rows={len(rows)} fnv={fnv(rows)} "
        f"sha256={hashlib.sha256(payload(rows)).hexdigest()} "
        f"max_endpoint={maximum}"
    )


def show(values: set[tuple[int, int]]) -> str:
    return " ".join(f"({q},{r})" for q, r in sorted(values)) or "EMPTY"


def write_pairs(path: Path, values: set[tuple[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload(sorted(values)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("live", type=Path)
    parser.add_argument("old_union", type=Path)
    parser.add_argument("h297", type=Path)
    parser.add_argument("--prefix-out", type=Path, required=True)
    parser.add_argument("--union-out", type=Path, required=True)
    parser.add_argument("--residual-out", type=Path, required=True)
    args = parser.parse_args()

    live_rows = read_pairs(args.live)
    old_rows = read_pairs(args.old_union)
    h_rows = read_pairs(args.h297)
    live = set(live_rows)
    old = set(old_rows)
    h297 = set(h_rows)
    prefix = {(q, r) for q, r in live if r >= 597}

    require(len(live) == 22647 and fnv(live_rows) == "df5374d4aca67677",
            "inherited live residual changed")
    require(len(old) == 1324 and fnv(old_rows) == "f55ee025df29bb65",
            "inherited typed union changed")
    require(len(h297) == 42 and fnv(h_rows) == "211843ee21a19170",
            "index-297 ideal changed")
    require(old <= live and h297 <= live and prefix <= live,
            "typed node escaped inherited residual")

    print("THM4300_TYPED_UNION_CONSUMER_V1")
    identity("LIVE", live)
    identity("OLD_UNION", old)
    identity("K597", prefix)
    identity("H297", h297)
    require(len(prefix) == 354 and fnv(sorted(prefix)) == "33b0069ca994b786",
            "complete endpoint prefix changed")
    require(hashlib.sha256(payload(sorted(prefix))).hexdigest() ==
            "e653857c9bfe1ef50e7724cfad05232b3695f88534ca289844a46d914ca52df5",
            "endpoint-prefix SHA changed")

    old_k = old & prefix
    old_h = old & h297
    k_h = prefix & h297
    triple = old & prefix & h297
    print(f"OVERLAPS old_k597={len(old_k)} old_h297={len(old_h)} "
          f"k597_h297={len(k_h)} triple={len(triple)}")
    require(len(old_k) == 96 and not old_h and not k_h and not triple,
            "typed overlap profile changed")

    union = old | prefix | h297
    residual = live - union
    identity("TYPED_UNION", union)
    identity("FINAL_RESIDUAL", residual)
    require(len(union) == 1624 and fnv(sorted(union)) == "11414a33ab91fef6",
            "final typed union changed")
    require(hashlib.sha256(payload(sorted(union))).hexdigest() ==
            "ef9102553cd030f67ab1bdb7d6965c3efaf4b0d8aa85daa1092354c9703caf26",
            "final typed-union SHA changed")
    require(len(residual) == 21023 and fnv(sorted(residual)) ==
            "e93e8089e9dc58c0", "final residual changed")
    require(hashlib.sha256(payload(sorted(residual))).hexdigest() ==
            "ce215cb53a742e5e0d0d4f16e344687da9583ce5dc1417c7ea70399fb0bf70ba",
            "final residual SHA changed")
    require(not union & residual and union | residual == live,
            "typed union does not partition live residual")
    maximum = max(r for _, r in residual)
    top = {(q, r) for q, r in residual if r == maximum}
    expected_top = {(96,596), (100,596), (192,596), (210,596),
                    (256,596), (260,596), (294,596), (306,596),
                    (384,596)}
    require(maximum == 596 and top == expected_top,
            "residual boundary changed")
    print(f"FINAL_MAX {maximum} TOP {show(top)}")
    print("INCREMENTS K597_NEW 258 H297_NEW 42 TOTAL_NEW 300")
    print("TYPE_GUARD separate_certificates_row_consequences_only")
    print("SCOPE FINITE_EXACT_TYPED_ROW_SET_UNION_NO_COMMON_DECK_NO_LRC14")

    write_pairs(args.prefix_out, prefix)
    write_pairs(args.union_out, union)
    write_pairs(args.residual_out, residual)


if __name__ == "__main__":
    main()
