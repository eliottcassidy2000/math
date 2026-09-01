#!/usr/bin/env python3
"""Independent row-set/hash consumer for H297, old union1324, and K597."""

from __future__ import annotations

import csv
import hashlib
import pathlib
import sys

FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1
TARGET_INDEX = 297


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized(rows: list[tuple[int, int]]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in rows).encode("ascii")


def sha256_rows(rows: list[tuple[int, int]]) -> str:
    return hashlib.sha256(normalized(rows)).hexdigest()


def fnv_add(state: int, word: int) -> int:
    for byte in word.to_bytes(8, "little", signed=False):
        state ^= byte
        state = (state * FNV_PRIME) & MASK64
    return state


def fnv_rows(rows: list[tuple[int, int]]) -> str:
    state = FNV_OFFSET
    for q, r in rows:
        state = fnv_add(state, q)
        state = fnv_add(state, r)
    return f"{state:016x}"


def read_pairs(path: pathlib.Path) -> list[tuple[int, int]]:
    with path.open(newline="", encoding="ascii") as source:
        rows = [(int(fields[0]), int(fields[1])) for fields in csv.reader(source)]
    if len(rows) != len(set(rows)):
        raise ValueError(f"duplicate pair in {path}")
    return rows


def derive_h297(
    signatures_path: pathlib.Path, live: set[tuple[int, int]]
) -> tuple[list[tuple[int, int]], int]:
    rows: list[tuple[int, int]] = []
    atlas_count = 0
    with signatures_path.open(newline="", encoding="ascii") as source:
        reader = csv.reader(source)
        header = next(reader)
        expected = ["q", "r", "inactive_count"] + [f"w{i}" for i in range(7)]
        if header != expected:
            raise ValueError("signature header changed")
        for fields in reader:
            atlas_count += 1
            if len(fields) != 10:
                raise ValueError("malformed signature row")
            pair = (int(fields[0]), int(fields[1]))
            words = [int(word, 16) for word in fields[3:]]
            weight = sum(word.bit_count() for word in words)
            if weight != int(fields[2]):
                raise ValueError(f"signature weight mismatch at {pair}")
            if pair not in live or weight != 1:
                continue
            target_word, target_bit = divmod(TARGET_INDEX, 64)
            expected_words = [0] * 7
            expected_words[target_word] = 1 << target_bit
            if words == expected_words:
                rows.append(pair)
    return rows, atlas_count


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: audit SIGNATURES LIVE TYPED_UNION1324")
    signatures_path, live_path, union_path = map(pathlib.Path, sys.argv[1:])
    live_rows = read_pairs(live_path)
    union_rows = read_pairs(union_path)
    live = set(live_rows)
    union = set(union_rows)
    h297, atlas_count = derive_h297(signatures_path, live)
    h297_set = set(h297)
    k597_source = sorted(pair for pair in live_rows if pair[1] >= 597)
    k597_band = sorted(k597_source, key=lambda pair: (-pair[1], pair[0]))
    h_union = sorted(h297_set & union)
    h_k597 = sorted(h297_set & set(k597_source))

    require(atlas_count == 24223, "signature atlas count changed")
    require(len(live_rows) == 22647, "live residual count changed")
    require(len(h297) == len(h297_set) == 42, "H297 count/uniqueness changed")
    require(max(r for _, r in h297) == 550, "H297 maximum endpoint changed")
    require(sha256_rows(h297) ==
            "ebe58d57e62a60be8325c7618146cb3e18cf9a20b51ee7270615290c109f06e6",
            "H297 normalized SHA-256 changed")
    require(fnv_rows(h297) == "211843ee21a19170", "H297 FNV changed")
    require(len(union_rows) == 1324, "old typed union count changed")
    require(len(k597_source) == 354, "K597 count changed")
    require(sha256_rows(k597_source) ==
            "e653857c9bfe1ef50e7724cfad05232b3695f88534ca289844a46d914ca52df5",
            "K597 source-order SHA-256 changed")
    require(sha256_rows(k597_band) ==
            "0574428f9b5b50be4d473ebc5403f81d22e5fe8cfdbe85747bf56152352c77e8",
            "K597 band-order SHA-256 changed")
    require(not h_union, "H297 intersects old typed union")
    require(not h_k597, "H297 intersects K597")

    print("INDEX297_INDEPENDENT_ROWSET_OVERLAP_AUDIT_V1")
    print(
        f"INPUT SIGNATURES_SHA256 {sha256_file(signatures_path)} "
        f"LIVE_SHA256 {sha256_file(live_path)} UNION_SHA256 {sha256_file(union_path)}"
    )
    print(
        f"H297 COUNT {len(h297)} MAX_ENDPOINT {max(r for _, r in h297)} "
        f"SHA256 {sha256_rows(h297)} FNV {fnv_rows(h297)}"
    )
    print(
        f"OLD_TYPED_UNION COUNT {len(union_rows)} SHA256 {sha256_rows(union_rows)} "
        f"FNV {fnv_rows(union_rows)} OVERLAP {len(h_union)} "
        f"OVERLAP_SHA256 {sha256_rows(h_union)} OVERLAP_FNV {fnv_rows(h_union)}"
    )
    print(
        f"K597 COUNT {len(k597_source)} SOURCE_SHA256 {sha256_rows(k597_source)} "
        f"SOURCE_FNV {fnv_rows(k597_source)} BAND_SHA256 {sha256_rows(k597_band)} "
        f"BAND_FNV {fnv_rows(k597_band)} OVERLAP {len(h_k597)} "
        f"OVERLAP_SHA256 {sha256_rows(h_k597)} OVERLAP_FNV {fnv_rows(h_k597)}"
    )
    print("VERDICT PASS H297_DISJOINT_OLD_TYPED_UNION_AND_K597")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
