#!/usr/bin/env python3
"""Hardened verifier for the independent endpoint-588 scratch packet."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

FNV_BASIS = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnv_add(state: int, word: int) -> int:
    for byte in int(word).to_bytes(8, "little", signed=False):
        state ^= byte
        state = (state * FNV_PRIME) & ((1 << 64) - 1)
    return state


def verify_manifest(packet: Path) -> int:
    manifest = packet / "SHA256SUMS"
    need(manifest.is_file(), "missing SHA256SUMS")
    entries: dict[str, str] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        need(len(digest) == 64 and name not in entries, "malformed manifest")
        entries[name] = digest
    actual = {
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS" and path.suffix != ".exe"
    }
    need(set(entries) == actual, "manifest closure mismatch")
    for name, digest in entries.items():
        need(sha(packet / name) == digest, f"manifest hash mismatch: {name}")
    return len(entries)


def same(packet: Path, left: str, right: str) -> None:
    need((packet / left).read_bytes() == (packet / right).read_bytes(),
         f"mode disagreement: {left} / {right}")


def read_pair_ledger(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    need(len(lines) == 67 and lines[0].startswith("q,r,active,"),
         "pair ledger shape changed")
    rows = []
    for line in lines[1:]:
        fields = line.split(",")
        need(len(fields) == 12 and int(fields[1]) == 588,
             "malformed pair row")
        rows.append((int(fields[0]), int(fields[1]), int(fields[10])))
    need(rows == sorted(rows) and sum(row[2] for row in rows) == 144867 and
         sum(row[2] != 0 for row in rows) == 10,
         "pair failure census changed")


def read_failures(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    need(len(lines) == 144868 and lines[0] == "q,r,body_hex",
         "failure ledger shape changed")
    expected = {25: 9, 50: 99836, 96: 1130, 100: 14, 105: 60,
                206: 7, 210: 43799, 256: 7, 420: 1, 462: 4}
    counts: dict[int, int] = {}
    prior = (-1, -1)
    state = FNV_BASIS
    for line in lines[1:]:
        q_text, r_text, body_text = line.split(",")
        q, r, body = int(q_text), int(r_text), int(body_text, 16)
        need(r == 588 and body.bit_count() == 9 and (q, body) > prior,
             "failure ordering/rank changed")
        prior = (q, body); counts[q] = counts.get(q, 0) + 1
        for word in (q, r, body): state = fnv_add(state, word)
    need(counts == expected and state == 0x6F51FA88F3B09CDC,
         "failure counts/FNV changed")


def read_literal(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    need(len(lines) == 144868 and lines[0].startswith("q,r,ordinal,body_hex,"),
         "literal detail shape changed")
    minimum: tuple[int, int, int] | None = None
    equality = 0
    for line in lines[1:]:
        q, r, _ordinal, body, low_mass, full_mass, low_ticks, full_ticks = line.split(",")
        q_i, r_i, body_i = int(q), int(r), int(body, 16)
        low_m, full_m = int(low_mass), int(full_mass)
        low_t, full_t = int(low_ticks), int(full_ticks)
        need(r_i == 588 and low_m <= full_m and 0 < low_t <= full_t,
             "literal lower-bound inequality/positivity failed")
        equality += low_m == full_m
        item = (low_t, q_i, body_i)
        if minimum is None or item < minimum: minimum = item
    need(minimum == (706240375488648, 210, 0x06B82090) and equality == 109452,
         "literal minimum/equality census changed")


def read_carrier(path: Path) -> None:
    masks = [int(line, 16) for line in path.read_text(encoding="utf-8").splitlines()]
    need(len(masks) == len(set(masks)) == 3925, "carrier size/distinctness changed")
    need(sum(mask.bit_count() == 8 for mask in masks) == 3809 and
         sum(mask.bit_count() == 9 for mask in masks) == 116,
         "carrier rank split changed")
    state = FNV_BASIS
    for mask in masks: state = fnv_add(state, mask)
    need(state == 0xEEAE5518D84CCAC5, "carrier FNV changed")


def read_typed(packet: Path) -> None:
    expected = {
        "typed_union2207.csv": (2207, 0x18D067B5614CF47F,
            "f03c84f15d9a149b0a083b50e922118e814d5644f5fa21e7011ae1c414ff3675"),
        "final_residual20440.csv": (20440, 0x794BD808E92E27CD,
            "be2d63e98beefb062e9ae4436d490ee2f630352989bf309cf85b5a1ffc44278c"),
        "residual_top587.csv": (10, 0xF48CA5F1904D6F52,
            "e2e841bccef0773513cc71d6b60ed03aa227cc701e19bc8a4673b4b7971d2a63"),
    }
    sets = {}
    for name, (count, expected_fnv, expected_sha) in expected.items():
        left = packet / "typed_normal" / name
        right = packet / "typed_opt" / name
        need(left.read_bytes() == right.read_bytes() and sha(left) == expected_sha,
             f"typed mode/SHA mismatch: {name}")
        rows = [tuple(map(int, line.split(",")))
                for line in left.read_text(encoding="utf-8").splitlines()]
        need(rows == sorted(set(rows)) and len(rows) == count,
             f"typed count/order mismatch: {name}")
        state = FNV_BASIS
        for row in rows:
            for word in row: state = fnv_add(state, word)
        need(state == expected_fnv, f"typed FNV mismatch: {name}")
        sets[name] = set(rows)
    need(sets["typed_union2207.csv"].isdisjoint(sets["final_residual20440.csv"]),
         "typed successor overlaps")
    need(sets["residual_top587.csv"] <= sets["final_residual20440.csv"] and
         all(r == 587 for _, r in sets["residual_top587.csv"]),
         "next frontier mismatch")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--packet", type=Path, required=True)
    packet = parser.parse_args().packet.resolve()
    manifest_files = verify_manifest(packet)
    same(packet, "endpoint588_cleanroom_pair_O2.csv", "endpoint588_cleanroom_pair_O3.csv")
    same(packet, "endpoint588_cleanroom_failures_O2.csv", "endpoint588_cleanroom_failures_O3.csv")
    same(packet, "endpoint588_cleanroom_literal_O2.csv", "endpoint588_cleanroom_literal_O3.csv")
    read_carrier(packet / "inherited_carrier3925.txt")
    read_pair_ledger(packet / "endpoint588_cleanroom_pair_O2.csv")
    read_failures(packet / "endpoint588_cleanroom_failures_O2.csv")
    read_literal(packet / "endpoint588_cleanroom_literal_O2.csv")
    read_typed(packet)
    carrier_out = (packet / "endpoint588_cleanroom_carrier_O2.out").read_text()
    literal_out = (packet / "endpoint588_cleanroom_literal_O2.out").read_text()
    need("PAIR_LEDGER_FNV 204a72794170e186" in carrier_out and
         "FAILURE_FNV 6f51fa88f3b09cdc" in carrier_out,
         "carrier transcript changed")
    need("LOW_POSITIVE 144867 FULL_POSITIVE 144867" in literal_out and
         "SUMMARY_FNV c3593788a0207b62" in literal_out,
         "literal transcript changed")
    print("INDEPENDENT_ENDPOINT588_PACKET_VERIFY PASS")
    print(f"manifest_files={manifest_files}")
    print("rows=66 body_tests=944271900 failures=144867 failed_rows=10")
    print("failure_fnv=6f51fa88f3b09cdc pair_fnv=204a72794170e186")
    print("literal_positive=144867 minimum_ticks=706240375488648 minimum_q=210 minimum_body=06b82090")
    print("typed_union=2207 residual=20440 next_endpoint=587 next_rows=10")
    print("scope=FINITE_EXACT_NO_PHYSICAL_ENTRY_NO_LRC14")


if __name__ == "__main__":
    main()
