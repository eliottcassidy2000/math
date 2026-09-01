#!/usr/bin/env python3
"""Consume the finite-exact endpoint-587 layer into the typed partition."""

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
    parser.add_argument("--top587", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe, 22647)
    prior_union = read_rows(args.prior_union, 2207)
    prior_residual = read_rows(args.prior_residual, 20440)
    top587 = read_rows(args.top587, 10)
    require(fnv_rows(universe) == 0xDF5374D4ACA67677,
            "universe identity changed")
    require(fnv_rows(prior_union) == 0x18D067B5614CF47F,
            "prior union identity changed")
    require(fnv_rows(prior_residual) == 0x794BD808E92E27CD,
            "prior residual identity changed")
    require(fnv_rows(top587) == 0xF48CA5F1904D6F52,
            "top587 identity changed")
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe),
            "prior partition invalid")
    require(set(top587).issubset(prior_residual) and
            all(r == 587 for _, r in top587),
            "top587 escaped prior residual/endpoint")

    union = sorted(set(prior_union) | set(top587))
    residual = sorted(set(prior_residual) - set(top587))
    require(len(union) == 2217 and len(residual) == 20430,
            "successor partition counts changed")
    require(set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "successor partition invalid")
    endpoint = max(r for _, r in residual)
    top = [row for row in residual if row[1] == endpoint]
    require(fnv_rows(union) == 0xE6592CBEF9B616D8,
            "successor union identity changed")
    require(fnv_rows(residual) == 0x4710F750DFCF91EA,
            "successor residual identity changed")
    require(endpoint == 586 and len(top) == 12 and
            fnv_rows(top) == 0xA1B617FAA2E7F63F,
            "successor top frontier changed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_sha = write_rows(args.output_dir / "typed_union2217.csv", union)
    residual_sha = write_rows(
        args.output_dir / "final_residual20430.csv", residual)
    top_sha = write_rows(args.output_dir / "residual_top586.csv", top)
    require(union_sha ==
            "d465a2f62c77ddaf921e7f3d7f32c674ea45e46fcdc348c90c711d7ba8a7a6b6" and
            residual_sha ==
            "8203dcfd6cc26b67bfdb648c7d8b50f94d7af1ab7ddbd6af2ff68acb15941f0b" and
            top_sha ==
            "d38b7fd9ea2aa9afdd12446d646cc8a9466cc5d4429612f03c8dff3165edf7ea",
            "successor typed SHA changed")

    print("LRC14_POST587_TYPED_SUCCESSOR_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv_rows(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv_rows(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv_rows(prior_residual):016x}")
    print(f"CONSUMED_TOP587 {len(top587)} FNV {fnv_rows(top587):016x}")
    print(f"NEW_UNION {len(union)} FNV {fnv_rows(union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(residual)} FNV {fnv_rows(residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {endpoint} ROWS {len(top)} FNV {fnv_rows(top):016x} SHA256 {top_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
