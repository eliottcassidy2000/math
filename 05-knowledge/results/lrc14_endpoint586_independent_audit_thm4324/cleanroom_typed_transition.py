#!/usr/bin/env python3
"""Import-free typed-partition reconstruction at endpoint 586."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split(",")
        require(len(fields) == 2, f"malformed row in {path.name}")
        q, r = map(int, fields)
        require(0 < q < r, f"invalid typed row in {path.name}")
        rows.append((q, r))
    require(rows == sorted(set(rows)), f"noncanonical rows in {path.name}")
    return rows


def fnv(rows: list[tuple[int, int]]) -> int:
    state = OFFSET
    for row in rows:
        for word in row:
            for shift in range(0, 64, 8):
                state = ((state ^ ((word >> shift) & 255)) * PRIME) & MASK64
    return state


def emit(path: Path, rows: list[tuple[int, int]]) -> str:
    path.write_bytes("".join(f"{q},{r}\n" for q, r in rows).encode("ascii"))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--prior-union", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe)
    prior_union = read_rows(args.prior_union)
    prior_residual = read_rows(args.prior_residual)
    require(len(universe) == 22647 and len(prior_union) == 2217 and
            len(prior_residual) == 20430, "typed input counts changed")
    require(set(prior_union).isdisjoint(prior_residual),
            "prior partition overlaps")
    require(set(prior_union) | set(prior_residual) == set(universe),
            "prior partition does not equal universe")

    endpoint = max(r for _, r in prior_residual)
    top = [row for row in prior_residual if row[1] == endpoint]
    require(endpoint == 586 and len(top) == 12,
            "reconstructed endpoint/frontier changed")
    union = sorted(set(prior_union) | set(top))
    residual = sorted(set(prior_residual) - set(top))
    next_endpoint = max(r for _, r in residual)
    next_top = [row for row in residual if row[1] == next_endpoint]
    require(next_endpoint == 585 and len(next_top) == 23,
            "next endpoint/frontier changed")
    require(set(union).isdisjoint(residual) and
            set(union) | set(residual) == set(universe),
            "successor partition failed")

    args.output.mkdir(parents=True, exist_ok=True)
    top_sha = emit(args.output / "reconstructed_top586.csv", top)
    union_sha = emit(args.output / "typed_union2229.csv", union)
    residual_sha = emit(args.output / "final_residual20418.csv", residual)
    next_sha = emit(args.output / "residual_top585.csv", next_top)

    print("LRC14_ENDPOINT586_CLEANROOM_TYPED_TRANSITION_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv(prior_residual):016x}")
    print(f"RECONSTRUCTED_TOP {endpoint} ROWS {len(top)} FNV {fnv(top):016x} SHA256 {top_sha}")
    print(f"NEW_UNION {len(union)} FNV {fnv(union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(residual)} FNV {fnv(residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {next_endpoint} ROWS {len(next_top)} FNV {fnv(next_top):016x} SHA256 {next_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
