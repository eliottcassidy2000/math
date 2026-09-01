#!/usr/bin/env python3
"""Byte, identity, and scope verifier for the THM-4313 packet."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

OFFSET = 0xCBF29CE484222325
PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def fnv_words(words: list[int]) -> int:
    state = OFFSET
    for word in words:
        for shift in range(0, 64, 8):
            state = ((state ^ ((word >> shift) & 0xFF)) * PRIME) & MASK64
    return state


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same(packet: Path, stem: str, suffix: str) -> None:
    left = packet / f"{stem}_O2{suffix}"
    right = packet / f"{stem}_O3{suffix}"
    assert left.read_bytes() == right.read_bytes(), f"byte mismatch: {stem}"


def read_plain_rows(path: Path) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    assert rows == sorted(set(rows))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    args = parser.parse_args()
    packet = args.packet

    manifest_path = packet / "SHA256SUMS"
    manifest: dict[str, str] = {}
    for line in manifest_path.read_text(encoding="ascii").splitlines():
        digest, relative = line.split("  ", 1)
        assert len(digest) == 64 and all(ch in "0123456789abcdef" for ch in digest)
        assert relative not in manifest
        manifest[relative] = digest
    actual_paths = {
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file() and path != manifest_path
    }
    assert set(manifest) == actual_paths
    for relative, digest in manifest.items():
        assert sha256(packet / relative) == digest, f"manifest mismatch: {relative}"

    for stem, suffix in [
        ("endpoint592_raw", ".out"),
        ("endpoint592_pair", ".csv"),
        ("endpoint592_failures", ".csv"),
        ("endpoint592_response_multiplicity", ".out"),
        ("endpoint592_obligation_multiplicity", ".csv"),
        ("verify_global_packing8", ".out"),
        ("endpoint592_cover_direct_audit", ".out"),
        ("full_universe_integer_dual_pricing", ".out"),
        ("full_universe_integer_dual_violations", ".csv"),
        ("singleton_capacity_cover43_final", ".out"),
        ("singleton467_cover43_final", ".csv"),
        ("protected467_cover43_final", ".csv"),
        ("safe_old_stats467_cover43_final", ".csv"),
        ("exchange43_final_raw", ".out"),
        ("exchange43_final_pair", ".csv"),
        ("exchange43_final_failures", ".csv"),
    ]:
        same(packet, stem, suffix)

    with (packet / "endpoint592_failures_O3.csv").open(
            newline="", encoding="ascii") as handle:
        failures = list(csv.DictReader(handle))
    assert len(failures) == 2468
    assert Counter(int(row["q"]) for row in failures) == {
        96: 13, 100: 3, 105: 48, 192: 1,
        210: 288, 256: 2101, 416: 14,
    }
    failure_words: list[int] = []
    for row in failures:
        failure_words += [int(row["q"]), int(row["r"]),
                          int(row["body_hex"], 16)]
    assert fnv_words(failure_words) == 0x2209B8D6760280CC

    with (packet / "endpoint592_obligation_multiplicity_O3.csv").open(
            newline="", encoding="ascii") as handle:
        multiplicities = list(csv.DictReader(handle))
    assert len(multiplicities) == 2468
    totals = [int(row["total_responders"]) for row in multiplicities]
    assert min(totals) == 423 and max(totals) == 34654
    assert sum(int(row["rank8_responders"]) == 0
               for row in multiplicities) == 6
    multiplicity_words: list[int] = []
    for row in multiplicities:
        multiplicity_words += [
            int(row["q"]), int(row["r"]), int(row["body_hex"], 16),
            int(row["rank8_responders"]), int(row["rank9_responders"]),
        ]
    assert fnv_words(multiplicity_words) == 0x1A9370CC4B793710

    packing = (packet / "verify_global_packing8_O3.out").read_text(
        encoding="ascii")
    assert "CERTIFICATE 8 FNV 52d324c3e51cc66b" in packing
    assert "PAIR_AUDIT Q256_PAIRS 21 CROSS_PAIRS 7 TOTAL 28 CO_RESPONSES 0" in packing
    assert "LOWER_BOUND 8" in packing and "VERDICT PASS" in packing

    cover_path = packet / "cover43.csv"
    with cover_path.open(newline="", encoding="ascii") as handle:
        cover_rows = list(csv.DictReader(handle))
    cover = [int(row["mask_hex"], 16) for row in cover_rows]
    assert len(cover) == len(set(cover)) == 43
    assert fnv_words(cover) == 0xCA3CB80F471F2E7E
    assert sha256(cover_path) == (
        "de7938a803d1d0f1bc99df6e1f5cb3ee6e90cf02c91fe3594e3d18f2061e3a8e")
    cover_verify = json.loads(
        (packet / "cover43_verify.json").read_text(encoding="ascii"))
    assert cover_verify["selected"] == 43
    assert cover_verify["covered_obligations"] == 2468
    assert cover_verify["missing_obligations"] == 0
    with (packet / "retained_pool_masks.csv").open(
            newline="", encoding="ascii") as handle:
        pool_rows = list(csv.DictReader(handle))
    pool_masks = [int(row["mask_hex"], 16) for row in pool_rows]
    assert len(pool_masks) == len(set(pool_masks)) == 37497
    assert all(mask.bit_count() == int(row["rank"])
               and int(row["rank"]) in (8, 9)
               for mask, row in zip(pool_masks, pool_rows, strict=True))
    assert fnv_words(pool_masks) == 0x557039EEB8DB4ED4
    assert set(cover).issubset(pool_masks)
    direct = (packet / "endpoint592_cover_direct_audit_O3.out").read_text(
        encoding="ascii")
    assert "INCIDENCES 7842" in direct and "MISSING 0" in direct

    dual = json.loads(
        (packet / "pool_dual_integer_certificate.json").read_text(
            encoding="ascii"))
    assert dual["scope"] == "all 37,497 exported retained-pool columns only"
    assert dual["pool_mask_count"] == len(pool_masks)
    assert dual["pool_mask_fnv64"] == "557039eeb8db4ed4"
    assert dual["dual_numerator"] == 35261737295
    assert dual["maximum_column_numerator"] == dual["scale"] == 1_000_000_000
    assert dual["column_violations"] == 0
    pricing = (packet / "full_universe_integer_dual_pricing_O3.out").read_text(
        encoding="ascii")
    assert "RANK8_RESPONDERS 188462 OMITTED 186480 VIOLATING 508" in pricing
    assert "RANK9_RESPONDERS 2205776 OMITTED 2170261 VIOLATING 20478" in pricing
    assert "MAX_ALL_NUMERATOR 3849498499 MASK 20188724" in pricing
    assert "VIOLATION_FNV 9f5c2aabbea12034" in pricing
    with (packet / "full_universe_integer_dual_violations_O3.csv").open(
            newline="", encoding="ascii") as handle:
        assert len(list(csv.DictReader(handle))) == 20986

    deletion_path = packet / "delete43_low_activity.txt"
    deletions = [int(line, 16) for line in
                 deletion_path.read_text(encoding="ascii").splitlines()]
    assert len(deletions) == len(set(deletions)) == 43
    assert fnv_words(deletions) == 0x0DD14EF0FE3EEC62
    assert sha256(deletion_path) == (
        "0b882a04d40cfe987c7784ee2d704ce84fb893db7a7157089f47c6f984088a00")

    with (packet / "safe_old_stats467_cover43_final_O3.csv").open(
            newline="", encoding="ascii") as handle:
        safe_rows = list(csv.DictReader(handle))
    assert len(safe_rows) == 3135
    assert [int(row["mask_hex"], 16) for row in safe_rows[:43]] == deletions
    with (packet / "protected467_cover43_final_O3.csv").open(
            newline="", encoding="ascii") as handle:
        protected = list(csv.DictReader(handle))
    assert len(protected) == 412
    assert sum(row["origin"] == "carrier" for row in protected) == 369
    assert {int(row["mask_hex"], 16) for row in protected
            if row["origin"] == "addition"} == set(cover)
    with (packet / "singleton467_cover43_final_O3.csv").open(
            newline="", encoding="ascii") as handle:
        assert len(list(csv.DictReader(handle))) == 1510

    with (packet / "exchange43_final_pair_O3.csv").open(
            newline="", encoding="ascii") as handle:
        pair_rows = list(csv.DictReader(handle))
    assert len(pair_rows) == 467
    assert len({(int(row["q"]), int(row["r"])) for row in pair_rows}) == 467
    assert all(int(row["failures"]) == 0 for row in pair_rows)
    assert (packet / "exchange43_final_failures_O3.csv").read_text(
        encoding="ascii") == "q,r,body_hex\n"
    raw = (packet / "exchange43_final_raw_O3.out").read_text(encoding="ascii")
    for needle in [
        "EXCHANGED_CARRIER 3925 FNV a0d08a38c10bdab7 RANK8 3818 RANK9 107",
        "ROWS 467 ROW_FNV 2d6aa988098aa5eb BODY_TESTS 6681439050",
        "FAILURES 0 FAILURE_FNV cbf29ce484222325",
        "VERDICT PASS",
    ]:
        assert needle in raw

    union = read_plain_rows(packet / "typed/typed_union2087.csv")
    residual = read_plain_rows(packet / "typed/final_residual20560.csv")
    top = read_plain_rows(packet / "typed/residual_top591.csv")
    assert len(union) == 2087 and len(residual) == 20560 and len(top) == 13
    assert set(union).isdisjoint(residual)
    assert all(r == 591 for _, r in top)
    assert fnv_words([value for row in union for value in row]) == 0x23E4136827B770A5
    assert fnv_words([value for row in residual for value in row]) == 0x8D797592A729E0B3
    assert fnv_words([value for row in top for value in row]) == 0xFC332C0697C671C7

    print("THM4313_ENDPOINT592_PACKET_VERIFY PASS")
    print("failures=2468 packing_lb=8 cover=43 deletions=43")
    print("target_rows=467 body_tests=6681439050 exchange_failures=0")
    print("typed_union=2087 residual=20560 next_endpoint=591 next_rows=13")


if __name__ == "__main__":
    main()
