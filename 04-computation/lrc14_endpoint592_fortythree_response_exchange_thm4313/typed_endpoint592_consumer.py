#!/usr/bin/env python3
"""Consume the repaired endpoint-592 layer into the inherited typed partition."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def fnv_rows(rows: list[tuple[int, int]]) -> int:
    state = OFFSET
    for q, r in rows:
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state = ((state ^ byte) * PRIME) & MASK64
    return state


def read_rows(path: Path, expected: int | None = None) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        q_text, r_text = line.split(",")
        row = (int(q_text), int(r_text))
        if not 0 < row[0] < row[1]:
            raise RuntimeError(f"invalid row in {path}")
        rows.append(row)
    if rows != sorted(set(rows)):
        raise RuntimeError(f"row identity/order changed in {path}")
    if expected is not None and len(rows) != expected:
        raise RuntimeError(f"row count changed in {path}")
    return rows


def write_rows(path: Path, rows: list[tuple[int, int]]) -> str:
    data = "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def read_audit(path: Path, expected: int) -> dict[tuple[int, int], int]:
    audited: dict[tuple[int, int], int] = {}
    with path.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["q"]), int(row["r"]))
            if key in audited:
                raise RuntimeError(f"duplicate audit row in {path}")
            audited[key] = int(row["failures"])
    if len(audited) != expected:
        raise RuntimeError(f"audit row count changed in {path}")
    return audited


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--prior-union", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--top592", type=Path, required=True)
    parser.add_argument("--prefix391-audit", type=Path, required=True)
    parser.add_argument("--top594-audit", type=Path, required=True)
    parser.add_argument("--top593", type=Path, required=True)
    parser.add_argument("--pair-audit", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe, 22647)
    prior_union = read_rows(args.prior_union, 2052)
    prior_residual = read_rows(args.prior_residual, 20595)
    top592 = read_rows(args.top592, 35)
    if not set(prior_union).isdisjoint(prior_residual):
        raise RuntimeError("prior partition overlaps")
    if set(prior_union) | set(prior_residual) != set(universe):
        raise RuntimeError("prior partition does not equal universe")
    if not set(top592).issubset(prior_residual):
        raise RuntimeError("top592 escaped prior residual")
    if any(r != 592 for _, r in top592):
        raise RuntimeError("top592 endpoint changed")

    prefix391 = read_audit(args.prefix391_audit, 391)
    top594 = read_audit(args.top594_audit, 25)
    top593 = set(read_rows(args.top593, 16))
    expected_audited = set(prefix391) | set(top594) | top593 | set(top592)
    if len(expected_audited) != 467:
        raise RuntimeError("inherited/new audit targets overlap unexpectedly")
    audited = read_audit(args.pair_audit, 467)
    if set(audited) != expected_audited or any(audited.values()):
        raise RuntimeError("exchange pair audit did not close exact full467 target")
    with args.failures.open(newline="", encoding="ascii") as handle:
        if list(csv.DictReader(handle)):
            raise RuntimeError("exchange failure ledger is nonempty")

    new_union = sorted(set(prior_union) | set(top592))
    new_residual = sorted(set(prior_residual) - set(top592))
    if len(new_union) != 2087 or len(new_residual) != 20560:
        raise RuntimeError("new partition counts changed")
    max_endpoint = max(r for _, r in new_residual)
    new_top = [row for row in new_residual if row[1] == max_endpoint]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_sha = write_rows(args.output_dir / "typed_union2087.csv", new_union)
    residual_sha = write_rows(args.output_dir / "final_residual20560.csv", new_residual)
    top_sha = write_rows(args.output_dir / f"residual_top{max_endpoint}.csv", new_top)
    print("LRC14_ENDPOINT592_TYPED_CONSUMER_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv_rows(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv_rows(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv_rows(prior_residual):016x}")
    print(f"CONSUMED_TOP592 {len(top592)} FNV {fnv_rows(top592):016x}")
    print(f"NEW_UNION {len(new_union)} FNV {fnv_rows(new_union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(new_residual)} FNV {fnv_rows(new_residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {max_endpoint} ROWS {len(new_top)} FNV {fnv_rows(new_top):016x} SHA256 {top_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_EXCHANGE_CERTIFICATE_NOT_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
