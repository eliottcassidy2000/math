#!/usr/bin/env python3
"""Set-theoretic endpoint-587 frontier and successor reconstruction."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

FNV_BASIS = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
WORD_MASK = (1 << 64) - 1


def read_pairs(path: Path) -> tuple[tuple[int, int], ...]:
    pairs: list[tuple[int, int]] = []
    with path.open("r", encoding="ascii", newline="") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line:
                raise AssertionError(f"blank row in {path}")
            fields = line.split(",")
            if len(fields) != 2:
                raise AssertionError(f"malformed pair in {path}")
            pair = (int(fields[0]), int(fields[1]))
            if not (0 < pair[0] < pair[1]):
                raise AssertionError(f"invalid pair {pair}")
            if pairs and pairs[-1] >= pair:
                raise AssertionError(f"unsorted/duplicate pair in {path}")
            pairs.append(pair)
    return tuple(pairs)


def fnv(pairs: tuple[tuple[int, int], ...]) -> int:
    value = FNV_BASIS
    for pair in pairs:
        for word in pair:
            for shift in range(0, 64, 8):
                value ^= (word >> shift) & 0xFF
                value = (value * FNV_PRIME) & WORD_MASK
    return value


def write_pairs(path: Path, pairs: tuple[tuple[int, int], ...]) -> str:
    payload = "".join(f"{q},{r}\n" for q, r in pairs).encode("ascii")
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--universe", type=Path, required=True)
    parser.add_argument("--prior-union", type=Path, required=True)
    parser.add_argument("--prior-residual", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    universe = read_pairs(args.universe)
    prior_union = read_pairs(args.prior_union)
    prior_residual = read_pairs(args.prior_residual)
    universe_set = set(universe)
    union_set = set(prior_union)
    residual_set = set(prior_residual)
    if union_set & residual_set:
        raise AssertionError("prior union and residual overlap")
    if union_set | residual_set != universe_set:
        raise AssertionError("prior partition does not equal universe")
    if len(universe) != 22647 or fnv(universe) != 0xDF5374D4ACA67677:
        raise AssertionError("universe identity changed")
    if len(prior_union) != 2207 or fnv(prior_union) != 0x18D067B5614CF47F:
        raise AssertionError("prior union identity changed")
    if len(prior_residual) != 20440 or fnv(prior_residual) != 0x794BD808E92E27CD:
        raise AssertionError("prior residual identity changed")

    endpoint = max(r for _, r in prior_residual)
    top = tuple(pair for pair in prior_residual if pair[1] == endpoint)
    if endpoint != 587 or len(top) != 10 or fnv(top) != 0xF48CA5F1904D6F52:
        raise AssertionError("derived endpoint-587 frontier identity changed")

    new_union = tuple(sorted(union_set | set(top)))
    new_residual = tuple(pair for pair in prior_residual if pair not in set(top))
    next_endpoint = max(r for _, r in new_residual)
    next_top = tuple(pair for pair in new_residual if pair[1] == next_endpoint)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    top_sha = write_pairs(args.output_dir / "residual_top587.csv", top)
    union_sha = write_pairs(args.output_dir / "typed_union2217.csv", new_union)
    residual_sha = write_pairs(args.output_dir / "final_residual20430.csv", new_residual)
    next_sha = write_pairs(args.output_dir / "residual_top586.csv", next_top)

    print("LRC14_ENDPOINT587_CLEANROOM_TYPED_SUCCESSOR_V1")
    print(f"UNIVERSE {len(universe)} FNV {fnv(universe):x}")
    print(f"PRIOR_UNION {len(prior_union)} FNV {fnv(prior_union):x}")
    print(f"PRIOR_RESIDUAL {len(prior_residual)} FNV {fnv(prior_residual):x}")
    print(f"DERIVED_TOP587 {len(top)} FNV {fnv(top):x} SHA256 {top_sha}")
    print(
        f"NEW_UNION {len(new_union)} FNV {fnv(new_union):x} SHA256 {union_sha}"
    )
    print(
        f"NEW_RESIDUAL {len(new_residual)} FNV {fnv(new_residual):x} "
        f"SHA256 {residual_sha}"
    )
    print(
        f"NEW_TOP {next_endpoint} ROWS {len(next_top)} FNV {fnv(next_top):x} "
        f"SHA256 {next_sha}"
    )
    print("SCOPE TYPED_ROW_SET_CONSUMER_NO_PHYSICAL_ENTRY_NO_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
