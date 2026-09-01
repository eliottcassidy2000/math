#!/usr/bin/env python3
"""Consume the finite-exact endpoint-589 layer into the typed partition."""

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
    parser.add_argument("--top589", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe, 22647)
    prior_union = read_rows(args.prior_union, 2113)
    prior_residual = read_rows(args.prior_residual, 20534)
    top589 = read_rows(args.top589, 28)
    require(fnv_rows(universe) == 0xDF5374D4ACA67677,
            "universe identity changed")
    require(fnv_rows(prior_union) == 0xC806CCE6B836FDFF,
            "prior union identity changed")
    require(fnv_rows(prior_residual) == 0x11285B5A49F4150D,
            "prior residual identity changed")
    require(fnv_rows(top589) == 0x5D9429C9F9971322,
            "top589 identity changed")
    require(set(prior_union).isdisjoint(prior_residual) and
            set(prior_union) | set(prior_residual) == set(universe),
            "prior partition invalid")
    require(set(top589).issubset(prior_residual) and
            all(r == 589 for _, r in top589),
            "top589 escaped prior residual/endpoint")

    union = sorted(set(prior_union) | set(top589))
    residual = sorted(set(prior_residual) - set(top589))
    require(len(union) == 2141 and len(residual) == 20506,
            "successor partition counts changed")
    require(set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "successor partition invalid")
    endpoint = max(r for _, r in residual)
    top = [row for row in residual if row[1] == endpoint]
    require(fnv_rows(union) == 0xC84BB7F7EAA0F230,
            "successor union identity changed")
    require(fnv_rows(residual) == 0x3CD0863A93C7602E,
            "successor residual identity changed")
    require(endpoint == 588 and len(top) == 66 and
            fnv_rows(top) == 0x18CF9A572CF9A5BE,
            "successor top frontier changed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_sha = write_rows(args.output_dir / "typed_union2141.csv", union)
    residual_sha = write_rows(
        args.output_dir / "final_residual20506.csv", residual)
    top_sha = write_rows(args.output_dir / "residual_top588.csv", top)
    require(union_sha ==
            "09f6d765afd78c03b60e7c4d047cae7df883ef8b2782f4256ee0e7c7be38ee09" and
            residual_sha ==
            "c782f56439163ac6e9f4b5f230cde1397c048e7ae3308bcc08d449aeb4ede2d9" and
            top_sha ==
            "a490630678d4ba088b98e24f2ca3098fe2b9407e197ace1a4696df823cfeda69",
            "successor typed SHA changed")

    print("LRC14_POST589_TYPED_SUCCESSOR_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv_rows(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv_rows(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv_rows(prior_residual):016x}")
    print(f"CONSUMED_TOP589 {len(top589)} FNV {fnv_rows(top589):016x}")
    print(f"NEW_UNION {len(union)} FNV {fnv_rows(union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(residual)} FNV {fnv_rows(residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {endpoint} ROWS {len(top)} FNV {fnv_rows(top):016x} SHA256 {top_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
