#!/usr/bin/env python3
"""Independent typed-set arithmetic for the THM-4287 LRC residual.

Inputs are four plain ``q,r`` ledgers: current residual, singleton consumer,
J19 consumer, and endpoint-carrier consumer.  This script does no geometry
and intentionally makes no claim that the three consumer types share a deck.
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


def read_pairs(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    with path.open("r", encoding="ascii", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line or line == "q,r":
                continue
            fields = line.split(",")
            require(len(fields) == 2, f"malformed pair row in {path}: {line}")
            q, r = map(int, fields)
            require(0 < q < r, f"invalid pair in {path}: {line}")
            rows.append((q, r))
    require(rows == sorted(rows), f"ledger is not lexicographically sorted: {path}")
    require(len(rows) == len(set(rows)), f"ledger contains duplicates: {path}")
    return rows


def fnv_add(state: int, word: int) -> int:
    for shift in range(0, 64, 8):
        state ^= (word >> shift) & 0xFF
        state = (state * FNV_PRIME) & MASK64
    return state


def fnv(rows: list[tuple[int, int]]) -> str:
    state = FNV_OFFSET
    for q, r in rows:
        state = fnv_add(state, q)
        state = fnv_add(state, r)
    return f"{state:016x}"


def sha(rows: list[tuple[int, int]]) -> str:
    payload = b"".join(f"{q},{r}\n".encode("ascii") for q, r in rows)
    return hashlib.sha256(payload).hexdigest()


def max_top(rows: list[tuple[int, int]]) -> tuple[int, list[tuple[int, int]]]:
    require(rows, "cannot find maximum endpoint of empty ledger")
    endpoint = max(r for _, r in rows)
    return endpoint, [(q, r) for q, r in rows if r == endpoint]


def describe(name: str, values: set[tuple[int, int]]) -> None:
    rows = sorted(values)
    endpoint, top = max_top(rows)
    print(f"{name} rows={len(rows)} fnv={fnv(rows)} sha256={sha(rows)} "
          f"max_endpoint={endpoint} top=" +
          " ".join(f"({q},{r})" for q, r in top))


def write_pairs(path: Path, values: set[tuple[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for q, r in sorted(values):
            handle.write(f"{q},{r}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("current", type=Path)
    parser.add_argument("singleton", type=Path)
    parser.add_argument("j19", type=Path)
    parser.add_argument("endpoint", type=Path)
    parser.add_argument("--union-out", type=Path)
    parser.add_argument("--residual-out", type=Path)
    args = parser.parse_args()

    current = set(read_pairs(args.current))
    singleton = set(read_pairs(args.singleton))
    j19 = set(read_pairs(args.j19))
    endpoint = set(read_pairs(args.endpoint))
    for name, component in (("singleton", singleton), ("j19", j19),
                            ("endpoint", endpoint)):
        require(component <= current, f"{name} is not a subset of current")

    sj = singleton & j19
    se = singleton & endpoint
    je = j19 & endpoint
    triple = singleton & j19 & endpoint
    union = singleton | j19 | endpoint
    residual = current - union

    print("INDEPENDENT_TYPED_PROOF_GRAPH_CONSUMER_V1")
    describe("CURRENT", current)
    describe("SINGLETON", singleton)
    describe("J19", j19)
    describe("ENDPOINT", endpoint)
    print(f"OVERLAPS singleton_j19={len(sj)} singleton_endpoint={len(se)} "
          f"j19_endpoint={len(je)} triple={len(triple)}")
    for name, values in (("singleton_j19", sj), ("singleton_endpoint", se),
                         ("j19_endpoint", je), ("triple", triple)):
        print(name.upper() + "_ROWS " +
              (" ".join(f"({q},{r})" for q, r in sorted(values))
               if values else "EMPTY"))
    print("PROFILE "
          f"singleton_only={len(singleton - j19 - endpoint)} "
          f"j19_only={len(j19 - singleton - endpoint)} "
          f"endpoint_only={len(endpoint - singleton - j19)} "
          f"singleton_j19_only={len(sj - endpoint)} "
          f"singleton_endpoint_only={len(se - j19)} "
          f"j19_endpoint_only={len(je - singleton)} triple={len(triple)}")
    describe("TYPED_UNION", union)
    describe("FINAL_RESIDUAL", residual)
    require(not (union & residual), "union/residual are not disjoint")
    require(union | residual == current, "union/residual do not partition current")
    require(len(union) == (len(singleton) + len(j19) + len(endpoint)
                           - len(sj) - len(se) - len(je) + len(triple)),
            "inclusion-exclusion count failed")
    print("PARTITION PASS inclusion_exclusion=PASS subset_checks=PASS")
    print("TYPE_GUARD separate_consumers_no_common_deck_claim")

    if args.union_out is not None:
        write_pairs(args.union_out, union)
    if args.residual_out is not None:
        write_pairs(args.residual_out, residual)


if __name__ == "__main__":
    main()
