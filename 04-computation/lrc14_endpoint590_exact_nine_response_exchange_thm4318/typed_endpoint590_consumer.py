#!/usr/bin/env python3
"""Consume the endpoint-590 exchange-closed layer into the typed partition."""

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
    parser.add_argument("--top590", type=Path, required=True)
    parser.add_argument("--pair-audit", type=Path, required=True)
    parser.add_argument("--failures", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_rows(args.universe, 22647)
    prior_union = read_rows(args.prior_union, 2100)
    prior_residual = read_rows(args.prior_residual, 20547)
    top590 = read_rows(args.top590, 13)
    if fnv_rows(universe) != 0xDF5374D4ACA67677:
        raise RuntimeError("universe identity changed")
    if fnv_rows(prior_union) != 0x3B2D991DA091A7DF:
        raise RuntimeError("prior union identity changed")
    if fnv_rows(prior_residual) != 0x59CA49A11D140EC5:
        raise RuntimeError("prior residual identity changed")
    if fnv_rows(top590) != 0x44AA8A793D162CF9:
        raise RuntimeError("top590 identity changed")
    if not set(prior_union).isdisjoint(prior_residual):
        raise RuntimeError("prior partition overlaps")
    if set(prior_union) | set(prior_residual) != set(universe):
        raise RuntimeError("prior partition does not equal universe")
    if not set(top590).issubset(prior_residual):
        raise RuntimeError("top590 escaped prior residual")
    if any(r != 590 for _, r in top590):
        raise RuntimeError("top590 endpoint changed")

    audited = read_audit(args.pair_audit, 493)
    if not set(top590).issubset(audited) or any(audited.values()):
        raise RuntimeError("exchange audit did not close the full target493")
    with args.failures.open(newline="", encoding="ascii") as handle:
        if list(csv.DictReader(handle)):
            raise RuntimeError("exchange failure ledger is nonempty")

    new_union = sorted(set(prior_union) | set(top590))
    new_residual = sorted(set(prior_residual) - set(top590))
    if len(new_union) != 2113 or len(new_residual) != 20534:
        raise RuntimeError("successor partition counts changed")
    if not set(new_union).isdisjoint(new_residual):
        raise RuntimeError("successor partition overlaps")
    if set(new_union) | set(new_residual) != set(universe):
        raise RuntimeError("successor partition does not equal universe")
    max_endpoint = max(r for _, r in new_residual)
    new_top = [row for row in new_residual if row[1] == max_endpoint]
    if fnv_rows(new_union) != 0xC806CCE6B836FDFF:
        raise RuntimeError("successor union identity changed")
    if fnv_rows(new_residual) != 0x11285B5A49F4150D:
        raise RuntimeError("successor residual identity changed")
    if max_endpoint != 589 or len(new_top) != 28:
        raise RuntimeError("successor top endpoint/count changed")
    if fnv_rows(new_top) != 0x5D9429C9F9971322:
        raise RuntimeError("successor top identity changed")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    union_sha = write_rows(args.output_dir / "typed_union2113.csv", new_union)
    residual_sha = write_rows(args.output_dir / "final_residual20534.csv", new_residual)
    top_sha = write_rows(args.output_dir / "residual_top589.csv", new_top)
    if union_sha != "0006156b27e794783a9ddc65a932ef005064b6bbd168193ddb5d092511881aa3":
        raise RuntimeError("successor union SHA changed")
    if residual_sha != "664ffe4f24d281ccb6cbd0de250f50b008929f836696319e78fa345b591799bb":
        raise RuntimeError("successor residual SHA changed")
    if top_sha != "005306a28c756a93862fd1745414fc058032f4d99e32086c079c29607bfed0c6":
        raise RuntimeError("successor top SHA changed")

    print("LRC14_ENDPOINT590_TYPED_CONSUMER_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv_rows(universe):016x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv_rows(prior_union):016x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv_rows(prior_residual):016x}")
    print(f"CONSUMED_TOP590 {len(top590)} FNV {fnv_rows(top590):016x}")
    print(f"NEW_UNION {len(new_union)} FNV {fnv_rows(new_union):016x} SHA256 {union_sha}")
    print(f"NEW_RESIDUAL {len(new_residual)} FNV {fnv_rows(new_residual):016x} SHA256 {residual_sha}")
    print(f"NEW_TOP {max_endpoint} ROWS {len(new_top)} FNV {fnv_rows(new_top):016x} SHA256 {top_sha}")
    print("SCOPE TYPED_ROW_SET_CONSUMER_EXCHANGE_CERTIFICATE_NOT_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
