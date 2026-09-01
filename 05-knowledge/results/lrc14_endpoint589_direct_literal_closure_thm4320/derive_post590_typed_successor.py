#!/usr/bin/env python3
"""Consume THM-4314's endpoint-590 layer into the frozen typed partition."""

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


def fnv_rows(rows: list[tuple[int, int]]) -> int:
    state = OFFSET
    for row in rows:
        for word in row:
            for shift in range(0, 64, 8):
                state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def read_rows(path: Path, expected: int) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)), f"noncanonical row ledger: {path}")
    require(len(rows) == expected, f"row count changed: {path}")
    return rows


def write_rows(path: Path, rows: list[tuple[int, int]]) -> str:
    data = "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--prior-union", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--top590", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe, 22647)
    prior_union = read_rows(args.prior_union, 2100)
    prior_residual = read_rows(args.prior_residual, 20547)
    top590 = read_rows(args.top590, 13)
    require(fnv_rows(universe) == 0xDF5374D4ACA67677,
            "universe identity changed")
    require(fnv_rows(prior_union) == 0x3B2D991DA091A7DF,
            "prior union identity changed")
    require(fnv_rows(prior_residual) == 0x59CA49A11D140EC5,
            "prior residual identity changed")
    require(fnv_rows(top590) == 0x44AA8A793D162CF9,
            "top590 identity changed")
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe),
            "prior partition invalid")
    require(set(top590).issubset(prior_residual) and
            all(r == 590 for _, r in top590),
            "top590 escaped prior residual/endpoint")

    union = sorted(set(prior_union) | set(top590))
    residual = sorted(set(prior_residual) - set(top590))
    require(len(union) == 2113 and len(residual) == 20534,
            "successor partition counts changed")
    require(set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "successor partition invalid")
    endpoint = max(r for _, r in residual)
    top = [row for row in residual if row[1] == endpoint]
    require(fnv_rows(union) == 0xC806CCE6B836FDFF,
            "successor union identity changed")
    require(fnv_rows(residual) == 0x11285B5A49F4150D,
            "successor residual identity changed")
    require(endpoint == 589 and len(top) == 28 and
            fnv_rows(top) == 0x5D9429C9F9971322,
            "successor top frontier changed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_sha = write_rows(args.output_dir / "typed_union2113.csv", union)
    residual_sha = write_rows(args.output_dir / "final_residual20534.csv", residual)
    top_sha = write_rows(args.output_dir / f"residual_top{endpoint}.csv", top)
    require(union_sha ==
            "0006156b27e794783a9ddc65a932ef005064b6bbd168193ddb5d092511881aa3" and
            residual_sha ==
            "664ffe4f24d281ccb6cbd0de250f50b008929f836696319e78fa345b591799bb" and
            top_sha ==
            "005306a28c756a93862fd1745414fc058032f4d99e32086c079c29607bfed0c6",
            "successor typed SHA changed")
    print("LRC14_POST590_TYPED_SUCCESSOR_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv_rows(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv_rows(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv_rows(prior_residual):016x}")
    print(f"CONSUMED_TOP590 {len(top590)} FNV {fnv_rows(top590):016x}")
    print(f"NEW_UNION {len(union)} FNV {fnv_rows(union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(residual)} FNV {fnv_rows(residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {endpoint} ROWS {len(top)} FNV {fnv_rows(top):016x} SHA256 {top_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
