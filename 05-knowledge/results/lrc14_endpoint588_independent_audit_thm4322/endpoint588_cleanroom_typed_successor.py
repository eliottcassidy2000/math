#!/usr/bin/env python3
"""Independent set-theoretic consumer for the endpoint-588 typed layer."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

FNV_BASIS = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_rows(path: Path) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw:
            continue
        fields = raw.split(",")
        need(len(fields) == 2, f"malformed typed row in {path}")
        row = (int(fields[0]), int(fields[1]))
        need(0 < row[0] < row[1], f"invalid typed row in {path}")
        rows.append(row)
    need(rows == sorted(set(rows)), f"typed rows are not strict lexicographic set: {path}")
    return rows


def fnv(rows: list[tuple[int, int]]) -> int:
    state = FNV_BASIS
    for row in rows:
        for word in row:
            for byte in word.to_bytes(8, "little", signed=False):
                state ^= byte
                state = (state * FNV_PRIME) & ((1 << 64) - 1)
    return state


def write_rows(path: Path, rows: list[tuple[int, int]]) -> str:
    payload = "".join(f"{q},{r}\n" for q, r in rows).encode()
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--prior-union", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--top588", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe)
    prior_union = read_rows(args.prior_union)
    prior_residual = read_rows(args.prior_residual)
    top = read_rows(args.top588)
    U, C, R, T = map(set, (universe, prior_union, prior_residual, top))
    need(len(universe) == 22647 and fnv(universe) == 0xDF5374D4ACA67677,
         "universe identity changed")
    need(len(prior_union) == 2141 and fnv(prior_union) == 0xC84BB7F7EAA0F230,
         "prior union identity changed")
    need(len(prior_residual) == 20506 and fnv(prior_residual) == 0x3CD0863A93C7602E,
         "prior residual identity changed")
    need(len(top) == 66 and fnv(top) == 0x18CF9A572CF9A5BE,
         "endpoint-588 layer identity changed")
    need(C.isdisjoint(R) and C | R == U, "prior typed partition failed")
    need(T <= R and all(r == 588 for _, r in T), "top layer is not residual endpoint 588")
    need(T == {row for row in R if row[1] == max(r for _, r in R)},
         "top588 is not the complete maximum residual endpoint")

    new_union = sorted(C | T)
    new_residual = sorted(R - T)
    next_endpoint = max(r for _, r in new_residual)
    next_top = sorted(row for row in new_residual if row[1] == next_endpoint)
    need(next_endpoint == 587 and len(next_top) == 10,
         "unexpected typed successor frontier")
    need(set(new_union).isdisjoint(new_residual) and
         set(new_union) | set(new_residual) == U,
         "successor typed partition failed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_sha = write_rows(args.output_dir / "typed_union2207.csv", new_union)
    residual_sha = write_rows(args.output_dir / "final_residual20440.csv", new_residual)
    top_sha = write_rows(args.output_dir / "residual_top587.csv", next_top)
    print("LRC14_ENDPOINT588_CLEANROOM_TYPED_SUCCESSOR_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv(prior_residual):016x}")
    print(f"CONSUMED_TOP588 {len(top)} FNV {fnv(top):016x}")
    print(f"NEW_UNION {len(new_union)} FNV {fnv(new_union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(new_residual)} FNV {fnv(new_residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {next_endpoint} ROWS {len(next_top)} FNV {fnv(next_top):016x} SHA256 {top_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
