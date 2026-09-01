#!/usr/bin/env python3
"""Independent typed row-set join for old THM-4296, K597, and H297."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

Pair = tuple[int, int]
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def pair_bytes(rows: set[Pair] | list[Pair]) -> bytes:
    return "".join(f"{q},{r}\n" for q, r in sorted(rows)).encode("ascii")


def pair_fnv(rows: set[Pair] | list[Pair]) -> str:
    state = FNV_OFFSET
    for q, r in sorted(rows):
        for value in (q, r):
            for byte in value.to_bytes(8, "little"):
                state = ((state ^ byte) * FNV_PRIME) & MASK64
    return f"{state:016x}"


def identity(rows: set[Pair] | list[Pair]) -> tuple[int, str, str]:
    raw = pair_bytes(rows)
    return len(rows), pair_fnv(rows), hashlib.sha256(raw).hexdigest()


def read_pair_ledger(path: Path, expected_count: int) -> set[Pair]:
    raw = path.read_bytes()
    rows: list[Pair] = []
    for number, line in enumerate(raw.splitlines(), 1):
        fields = line.split(b",")
        require(len(fields) == 2, f"{path}:{number}: malformed pair")
        rows.append((int(fields[0]), int(fields[1])))
    require(len(rows) == expected_count, f"{path}: row count changed")
    require(rows == sorted(set(rows)), f"{path}: not sorted and unique")
    require(raw == pair_bytes(rows), f"{path}: noncanonical pair-ledger bytes")
    return set(rows)


def derive_h297(signatures: Path, universe: set[Pair]) -> set[Pair]:
    atlas_pairs: set[Pair] = set()
    h297: set[Pair] = set()
    with signatures.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle)
        require(
            reader.fieldnames
            == ["q", "r", "inactive_count", "w0", "w1", "w2", "w3", "w4", "w5", "w6"],
            "signature-atlas header changed",
        )
        for row in reader:
            pair = (int(row["q"]), int(row["r"]))
            require(pair not in atlas_pairs, f"duplicate atlas pair {pair}")
            atlas_pairs.add(pair)
            words = [int(row[f"w{i}"], 16) for i in range(7)]
            inactive_count = int(row["inactive_count"])
            require(sum(word.bit_count() for word in words) == inactive_count,
                    f"inactive-count mismatch at {pair}")
            if pair not in universe or inactive_count != 1:
                continue
            singleton_index = next(
                64 * word_index + (word & -word).bit_length() - 1
                for word_index, word in enumerate(words) if word
            )
            if singleton_index == 297:
                h297.add(pair)
    require(len(atlas_pairs) == 24_223, "signature-atlas row count changed")
    require(universe <= atlas_pairs, "universe pair absent from signature atlas")
    return h297


def describe(label: str, rows: set[Pair]) -> str:
    count, fnv, sha = identity(rows)
    return f"{label} rows={count} fnv={fnv} sha256={sha}"


def row_text(rows: set[Pair]) -> str:
    return " ".join(f"({q},{r})" for q, r in sorted(rows))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--universe", required=True, type=Path)
    parser.add_argument("--old-union", required=True, type=Path)
    parser.add_argument("--old-residual", required=True, type=Path)
    parser.add_argument("--signatures", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    universe = read_pair_ledger(args.universe, 22_647)
    old = read_pair_ledger(args.old_union, 1_324)
    maintained_old_residual = read_pair_ledger(args.old_residual, 21_323)
    require(old <= universe, "old typed union escapes universe")
    require(universe - old == maintained_old_residual,
            "maintained THM-4296 residual is not universe minus old union")
    require(identity(universe)[1:] == (
        "df5374d4aca67677",
        "14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317",
    ), "universe identity changed")
    require(identity(old)[1:] == (
        "f55ee025df29bb65",
        "57bdcd932cc2c985e81e1b2d472a469cf4c2e11c9cb21dd4b1037c4ba098562a",
    ), "old-union identity changed")
    require(identity(maintained_old_residual)[1:] == (
        "09a0dfc4515d556b",
        "a6fb3ea00f5017b2ad66cbccf75f89756aa94360652ff6a7dbf22c637a1c7656",
    ), "old-residual identity changed")

    # K597 is typed as the complete r >= 597 prefix of the inherited universe,
    # not merely the previously unproved portion of that prefix.
    k597 = {pair for pair in universe if pair[1] >= 597}
    require(identity(k597) == (
        354,
        "33b0069ca994b786",
        "e653857c9bfe1ef50e7724cfad05232b3695f88534ca289844a46d914ca52df5",
    ), "K597 complete-prefix identity changed")
    require(all(pair[1] >= 597 for pair in k597), "K597 threshold failure")
    require(all(pair in k597 for pair in universe if pair[1] >= 597),
            "K597 completeness failure")

    h297 = derive_h297(args.signatures, universe)
    require(identity(h297) == (
        42,
        "211843ee21a19170",
        "ebe58d57e62a60be8325c7618146cb3e18cf9a20b51ee7270615290c109f06e6",
    ), "H297 identity changed")

    old_k = old & k597
    old_h = old & h297
    k_h = k597 & h297
    triple = old & k597 & h297
    k_new = k597 - old
    h_new = h297 - (old | k597)
    final_union = old | k597 | h297
    final_residual = universe - final_union
    require(final_union <= universe, "final union escapes universe")
    require(final_union.isdisjoint(final_residual), "union/residual overlap")
    require(final_union | final_residual == universe, "join does not partition universe")
    require(identity(old_k) == (
        96,
        "7b59a3a3d74b98a8",
        "4baa9aa7b345383db9fed34626ded5935a2dee98a054784ec3d2d92148e9b852",
    ), "OLD/K597 overlap changed")
    require(not old_h and not k_h and not triple,
            "H297 is no longer disjoint from OLD and K597")
    require(identity(final_union) == (
        1_624,
        "11414a33ab91fef6",
        "ef9102553cd030f67ab1bdb7d6965c3efaf4b0d8aa85daa1092354c9703caf26",
    ), "final-union identity changed")
    require(identity(final_residual) == (
        21_023,
        "e93e8089e9dc58c0",
        "ce215cb53a742e5e0d0d4f16e344687da9583ce5dc1417c7ea70399fb0bf70ba",
    ), "final-residual identity changed")

    profiles: dict[str, set[Pair]] = {}
    for pair in final_union:
        label = "+".join(name for name, node in (
            ("OLD", old), ("K597", k597), ("H297", h297)
        ) if pair in node)
        profiles.setdefault(label, set()).add(pair)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ledgers = {
        "k597_complete_prefix354.csv": k597,
        "h297_ideal42_independent.csv": h297,
        "k597_strictly_added.csv": k_new,
        "h297_strictly_added.csv": h_new,
        "typed_union_old_k597_h297.csv": final_union,
        "residual_old_k597_h297.csv": final_residual,
    }
    for name, rows in ledgers.items():
        (args.out_dir / name).write_bytes(pair_bytes(rows))

    print("TYPED_GRAPH_OLD_K597_H297_INDEPENDENT_V1")
    print(describe("UNIVERSE", universe))
    print(describe("OLD", old))
    print(describe("OLD_RESIDUAL", maintained_old_residual))
    print(describe("K597_COMPLETE_PREFIX", k597))
    print(f"K597_ENDPOINT_RANGE {min(r for _, r in k597)}..{max(r for _, r in k597)}")
    print(describe("H297", h297))
    print(f"H297_ENDPOINT_RANGE {min(r for _, r in h297)}..{max(r for _, r in h297)}")
    print("PAIRWISE_AND_TRIPLE_OVERLAPS")
    for label, rows in (
        ("OLD_AND_K597", old_k),
        ("OLD_AND_H297", old_h),
        ("K597_AND_H297", k_h),
        ("OLD_AND_K597_AND_H297", triple),
    ):
        print(describe(label, rows))
    print("EXACT_MEMBERSHIP_PROFILES")
    for label in sorted(profiles):
        print(describe(f"PROFILE_{label}", profiles[label]))
    print(describe("K597_STRICTLY_ADDED", k_new))
    print(describe("H297_STRICTLY_ADDED", h_new))
    print(describe("FINAL_TYPED_UNION", final_union))
    print(describe("FINAL_RESIDUAL", final_residual))
    max_endpoint = max(r for _, r in final_residual)
    top = {pair for pair in final_residual if pair[1] == max_endpoint}
    print(f"FINAL_MAX_ENDPOINT {max_endpoint} TOP_ROWS {row_text(top)}")
    for endpoint in sorted({r for _, r in final_residual}, reverse=True)[:10]:
        layer = {pair for pair in final_residual if pair[1] == endpoint}
        print(describe(f"FINAL_LAYER_{endpoint}", layer))
    print("TYPE_GUARD typed_row_set_join_only_no_common_deck_no_LRC14")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
