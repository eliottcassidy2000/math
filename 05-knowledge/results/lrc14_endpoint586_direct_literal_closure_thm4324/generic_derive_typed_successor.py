#!/usr/bin/env python3
"""Exact typed-partition successor for a frozen complete endpoint layer."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)), f"noncanonical rows: {path}")
    return rows


def fnv(rows: list[tuple[int, int]]) -> int:
    state = OFFSET
    for row in rows:
        for word in row:
            for shift in range(0, 64, 8):
                state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def write(path: Path, rows: list[tuple[int, int]]) -> str:
    data = "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--prior-union", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--top", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    universe = read_rows(args.universe)
    prior_union = read_rows(args.prior_union)
    prior_residual = read_rows(args.prior_residual)
    top = read_rows(args.top)
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe),
            "prior partition invalid")
    endpoint = max(r for _, r in prior_residual)
    require(top == [row for row in prior_residual if row[1] == endpoint],
            "top is not the complete current endpoint")
    union = sorted(set(prior_union) | set(top))
    residual = sorted(set(prior_residual) - set(top))
    require(set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "successor partition invalid")
    next_endpoint = max(r for _, r in residual)
    next_top = [row for row in residual if row[1] == next_endpoint]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_name = f"typed_union{len(union)}.csv"
    residual_name = f"final_residual{len(residual)}.csv"
    top_name = f"residual_top{next_endpoint}.csv"
    union_sha = write(args.output_dir / union_name, union)
    residual_sha = write(args.output_dir / residual_name, residual)
    top_sha = write(args.output_dir / top_name, next_top)
    print("LRC14_GENERIC_TYPED_SUCCESSOR_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv(prior_residual):016x}")
    print(f"CONSUMED_TOP {endpoint} ROWS {len(top)} FNV {fnv(top):016x}")
    print(f"NEW_UNION {len(union)} FNV {fnv(union):016x} SHA256 {union_sha} FILE {union_name}")
    print(f"NEW_RESIDUAL {len(residual)} FNV {fnv(residual):016x} SHA256 {residual_sha} FILE {residual_name}")
    print(f"NEW_TOP {next_endpoint} ROWS {len(next_top)} FNV {fnv(next_top):016x} SHA256 {top_sha} FILE {top_name}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
