#!/usr/bin/env python3
"""Consume the finite-exact endpoint-588 layer into the typed partition."""

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
    parser.add_argument("--top588", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe, 22647)
    prior_union = read_rows(args.prior_union, 2141)
    prior_residual = read_rows(args.prior_residual, 20506)
    top588 = read_rows(args.top588, 66)
    require(fnv_rows(universe) == 0xDF5374D4ACA67677,
            "universe identity changed")
    require(fnv_rows(prior_union) == 0xC84BB7F7EAA0F230,
            "prior union identity changed")
    require(fnv_rows(prior_residual) == 0x3CD0863A93C7602E,
            "prior residual identity changed")
    require(fnv_rows(top588) == 0x18CF9A572CF9A5BE,
            "top588 identity changed")
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe),
            "prior partition invalid")
    require(set(top588).issubset(prior_residual) and
            all(r == 588 for _, r in top588),
            "top588 escaped prior residual/endpoint")

    union = sorted(set(prior_union) | set(top588))
    residual = sorted(set(prior_residual) - set(top588))
    require(len(union) == 2207 and len(residual) == 20440,
            "successor partition counts changed")
    require(set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "successor partition invalid")
    endpoint = max(r for _, r in residual)
    top = [row for row in residual if row[1] == endpoint]
    require(fnv_rows(union) == 0x18D067B5614CF47F,
            "successor union identity changed")
    require(fnv_rows(residual) == 0x794BD808E92E27CD,
            "successor residual identity changed")
    require(endpoint == 587 and len(top) == 10 and
            fnv_rows(top) == 0xF48CA5F1904D6F52,
            "successor top frontier changed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_sha = write_rows(args.output_dir / "typed_union2207.csv", union)
    residual_sha = write_rows(
        args.output_dir / "final_residual20440.csv", residual)
    top_sha = write_rows(args.output_dir / "residual_top587.csv", top)
    require(union_sha ==
            "f03c84f15d9a149b0a083b50e922118e814d5644f5fa21e7011ae1c414ff3675" and
            residual_sha ==
            "be2d63e98beefb062e9ae4436d490ee2f630352989bf309cf85b5a1ffc44278c" and
            top_sha ==
            "e2e841bccef0773513cc71d6b60ed03aa227cc701e19bc8a4673b4b7971d2a63",
            "successor typed SHA changed")

    print("LRC14_POST588_TYPED_SUCCESSOR_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv_rows(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv_rows(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv_rows(prior_residual):016x}")
    print(f"CONSUMED_TOP588 {len(top588)} FNV {fnv_rows(top588):016x}")
    print(f"NEW_UNION {len(union)} FNV {fnv_rows(union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(residual)} FNV {fnv_rows(residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {endpoint} ROWS {len(top)} FNV {fnv_rows(top):016x} SHA256 {top_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
