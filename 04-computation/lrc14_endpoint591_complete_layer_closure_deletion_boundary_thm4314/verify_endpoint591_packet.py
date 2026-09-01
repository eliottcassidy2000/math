#!/usr/bin/env python3
"""Hardened manifest, ledger, deletion, and typed verifier for THM-4314."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import Counter, defaultdict
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


def fnv_rows(rows: list[tuple[int, int]]) -> int:
    return fnv_words([value for row in rows for value in row])


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: object) -> None:
    """Fail identically under normal Python and ``python -O``."""
    if not condition:
        raise RuntimeError(str(message))


def read_rows(path: Path, expected: int | None = None) -> list[tuple[int, int]]:
    rows = [tuple(map(int, line.split(",")))
            for line in path.read_text(encoding="ascii").splitlines()]
    require(rows == sorted(set(rows)), f"noncanonical row ledger: {path}")
    if expected is not None:
        require(len(rows) == expected,
                f"row count changed: {path}: {len(rows)} != {expected}")
    return rows


def same(left: Path, right: Path) -> None:
    require(left.read_bytes() == right.read_bytes(),
            f"byte mismatch: {left} {right}")


def read_dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--packet", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    packet = args.packet.resolve()

    manifest_path = packet / "SHA256SUMS"
    manifest: dict[str, str] = {}
    for line in manifest_path.read_text(encoding="ascii").splitlines():
        digest, relative = line.split("  ", 1)
        require(len(digest) == 64 and
                all(ch in "0123456789abcdef" for ch in digest),
                f"malformed manifest digest: {line}")
        require(relative not in manifest, f"duplicate manifest path: {relative}")
        manifest[relative] = digest
    actual_paths = {
        path.relative_to(packet).as_posix()
        for path in packet.rglob("*")
        if path.is_file() and path != manifest_path
    }
    require(set(manifest) == actual_paths, "manifest file set changed")
    for relative, digest in manifest.items():
        require(sha256(packet / relative) == digest,
                f"manifest mismatch: {relative}")

    # Compiler and interpreter variants must agree byte-for-byte.
    for stem, suffix in [
        ("endpoint591_raw", ".out"),
        ("endpoint591_pair", ".csv"),
        ("endpoint591_failures", ".csv"),
        ("endpoint591_all_witness", ".out"),
        ("endpoint591_all_witness_low", ".csv"),
        ("endpoint591_twohit", ".out"),
        ("endpoint591_twohit", ".csv"),
        ("endpoint591_twohit_pairs", ".csv"),
    ]:
        same(packet / f"{stem}_O3{suffix}",
             packet / f"{stem}_O2{suffix}")
    same(packet / "typed_endpoint591_consumer.out",
         packet / "typed_endpoint591_consumer_opt.out")
    for name in ["typed_union2100.csv", "final_residual20547.csv",
                 "residual_top590.csv"]:
        same(packet / "typed" / name,
             packet / "typed_opt" / name)

    packet4313 = repo / "05-knowledge/results/lrc14_endpoint592_fortythree_response_exchange_thm4313"
    top591_path = packet4313 / "typed/residual_top591.csv"
    require(sha256(top591_path) ==
            "b7c0cd2324bd4e7cb0819c30806faded369f5eba83f82e51d618be41e24bd211",
            "canonical top591 SHA changed")
    require(sha256(packet4313 / "cover43.csv") ==
            "de7938a803d1d0f1bc99df6e1f5cb3ee6e90cf02c91fe3594e3d18f2061e3a8e",
            "canonical cover43 SHA changed")
    require(sha256(packet4313 / "delete43_low_activity.txt") ==
            "0b882a04d40cfe987c7784ee2d704ce84fb893db7a7157089f47c6f984088a00",
            "canonical delete43 SHA changed")
    top591 = read_rows(top591_path, 13)
    require(fnv_rows(top591) == 0xFC332C0697C671C7,
            "canonical top591 FNV changed")
    require(all(r == 591 for _, r in top591),
            "canonical top591 contains another endpoint")

    pair_path = packet / "endpoint591_pair_O3.csv"
    with pair_path.open(newline="", encoding="ascii") as handle:
        pair_rows = list(csv.DictReader(handle))
    require(len(pair_rows) == 13, "pair audit row count changed")
    audited = [(int(row["q"]), int(row["r"])) for row in pair_rows]
    require(audited == top591, "pair audit target differs from top591")
    require(all(int(row["failures"]) == 0 for row in pair_rows),
            "pair audit contains failures")
    require(sum(int(row["exposed"]) for row in pair_rows) == 28791,
            "pair audit exposed-body total changed")
    require(min(int(row["minimum_hits"]) for row in pair_rows) == 2,
            "pair audit minimum hit multiplicity changed")
    require({int(row["q"]) for row in pair_rows
             if int(row["minimum_hits"]) == 2} == {96, 105, 210},
            "minimum-hit row set changed")
    pair_words: list[int] = []
    for row in pair_rows:
        pair_words.extend([
            int(row["q"]), int(row["r"]), int(row["active"]),
            int(row["active_fnv"], 16), int(row["active_joint"]),
            int(row["active_nonjoint"]), int(row["exposed"]),
            int(row["exposed_fnv"], 16), int(row["minimum_hits"]),
            int(row["maximum_hits"]), int(row["failures"]),
            int(row["failure_fnv"], 16),
        ])
    require(fnv_words(pair_words) == 0x8191899A3E142A2C,
            "pair audit FNV changed")
    failure_path = packet / "endpoint591_failures_O3.csv"
    require(failure_path.read_text(encoding="ascii") == "q,r,body_hex\n",
            "endpoint591 failure ledger is nonempty or malformed")

    raw = (packet / "endpoint591_raw_O3.out").read_text(
        encoding="ascii")
    for needle in [
        "CARRIER 3925 FNV a0d08a38c10bdab7 RANK8 3818 RANK9 107 JOINT_RETAINED 421",
        "ROWS 13 ENDPOINT 591 ROW_FNV fc332c0697c671c7",
        "BODY_TESTS 185992950",
        "SUMMARY EXPOSED 28791 HIT_INCIDENCES 1031512 FAILURES 0 FAILED_ROWS 0",
        "PAIR_LEDGER_FNV 8191899a3e142a2c",
        "NO_PHYSICAL_ENTRY_NO_LRC14",
        "VERDICT PASS",
    ]:
        require(needle in raw, f"raw transcript missing: {needle}")

    # The independent all-witness census classifies every deletion set of
    # cardinality at most two.  A body with k original witnesses can fail
    # after at most two deletions exactly when k <= 2 and every witness is
    # deleted, so the frozen low ledger is a complete certificate.
    all_text = (packet / "endpoint591_all_witness_O3.out").read_text(
        encoding="ascii")
    for needle in [
        "CARRIER 3925 FNV a0d08a38c10bdab7 ROWS 13 BODY_TESTS 185992950",
        "SUMMARY ZERO 0 ONE 2 TWO 26 LOW_FNV 8dad3ce6f3e587c2",
        "DISTINCT_TWO_PAIRS 23 MAX_PAIR_LOAD 3",
        "VERDICT PASS",
    ]:
        require(needle in all_text, f"all-witness transcript missing: {needle}")

    protected_text = (packet / "endpoint591_twohit_O3.out").read_text(
        encoding="ascii")
    for needle in [
        "CARRIER 3925 FNV a0d08a38c10bdab7 ROWS 13 EXPOSED 28791 BELOW_TWO 0",
        "TOTAL_TWOHIT 12 CENSUS_FNV 121f6e5d02874855",
        "DISTINCT_WITNESS_PAIRS 12 MAX_PAIR_LOAD 1",
        "VERDICT PASS",
    ]:
        require(needle in protected_text,
                f"protected-joint transcript missing: {needle}")

    joint_path = (repo / "05-knowledge/results/"
                  "lrc14_mixed_rank_depth_recursive_signatures_thm4296/"
                  "inputs/joint421_masks.txt")
    joint = {
        int(line, 16) for line in joint_path.read_text(encoding="ascii").splitlines()
        if line
    }
    require(len(joint) == 421, "joint deck size changed")

    low = read_dicts(packet / "endpoint591_all_witness_low_O3.csv")
    require(len(low) == 28, "low-witness ledger row count changed")
    one = [row for row in low if int(row["hits"]) == 1]
    two = [row for row in low if int(row["hits"]) == 2]
    require(len(one) == 2 and len(two) == 26,
            "one/two-witness split changed")
    require(all(int(row["r"]) == 591 and
                int(row["body_hex"], 16).bit_count() == 9
                for row in low),
            "low-witness row/body universe changed")

    critical_obligations: dict[int, list[str]] = defaultdict(list)
    for row in one:
        witness = int(row["first_mask_hex"], 16)
        body = int(row["body_hex"], 16)
        require(witness.bit_count() in (8, 9) and witness & body == 0,
                "invalid singleton witness")
        require(not row["second_mask_hex"], "singleton has a second witness")
        critical_obligations[witness].append(
            f'{row["q"]}:{row["r"]}:{row["body_hex"]}')
    critical = set(critical_obligations)
    require(critical == {0x024085D0, 0x20A09640},
            "critical singleton masks changed")
    require(critical <= joint, "critical singleton is not protected joint")

    pair_obligations: dict[tuple[int, int], list[str]] = defaultdict(list)
    pair_types: Counter[str] = Counter()
    for row in two:
        body = int(row["body_hex"], 16)
        pair = tuple(sorted((int(row["first_mask_hex"], 16),
                             int(row["second_mask_hex"], 16))))
        require(pair[0] != pair[1], "repeated two-witness mask")
        require(all(mask.bit_count() in (8, 9) and mask & body == 0
                    for mask in pair), "invalid two-witness pair")
        pair_obligations[pair].append(
            f'{row["q"]}:{row["r"]}:{row["body_hex"]}')
        joint_count = sum(mask in joint for mask in pair)
        pair_types[{0: "NN", 1: "JN", 2: "JJ"}[joint_count]] += 1
    require(len(pair_obligations) == 23, "distinct two-witness count changed")
    require(pair_types == Counter({"NN": 12, "JN": 10, "JJ": 4}),
            "two-witness pair-type counts changed")
    extra_pairs = {
        pair for pair in pair_obligations if not critical.intersection(pair)
    }
    require(len(extra_pairs) == 22, "additional unsafe-pair count changed")

    critical_rows = read_dicts(packet / "endpoint591_critical_singletons.csv")
    require({int(row["critical_mask_hex"], 16) for row in critical_rows}
            == critical and len(critical_rows) == 2,
            "derived critical-singleton ledger changed")
    pair_rows_derived = read_dicts(
        packet / "endpoint591_distinct_twohit_witness_pairs.csv")
    derived_pairs = {
        (int(row["first_mask_hex"], 16), int(row["second_mask_hex"], 16))
        for row in pair_rows_derived
    }
    require(derived_pairs == set(pair_obligations) and
            len(pair_rows_derived) == 23,
            "derived two-witness pair ledger changed")
    require(sum(int(row["additional_unsafe_pair"])
                for row in pair_rows_derived) == 22,
            "derived additional-pair flags changed")

    protected = read_dicts(packet / "endpoint591_twohit_O3.csv")
    protected_keys = {
        (int(row["q"]), int(row["r"]), row["body_hex"],
         tuple(sorted((int(row["first_mask_hex"], 16),
                       int(row["second_mask_hex"], 16)))))
        for row in protected
    }
    nn_keys = {
        (int(row["q"]), int(row["r"]), row["body_hex"],
         tuple(sorted((int(row["first_mask_hex"], 16),
                       int(row["second_mask_hex"], 16)))))
        for row in two
        if int(row["first_mask_hex"], 16) not in joint and
        int(row["second_mask_hex"], 16) not in joint
    }
    require(protected_keys == nn_keys and len(nn_keys) == 12,
            "protected-joint/nonjoint seam changed")

    total_pairs = math.comb(3925, 2)
    critical_pairs = total_pairs - math.comb(3923, 2)
    unsafe_pairs = critical_pairs + len(extra_pairs)
    require((total_pairs, critical_pairs, unsafe_pairs) ==
            (7_700_850, 7_847, 7_869),
            "deletion-pair arithmetic changed")
    boundary = (packet / "endpoint591_deletion_boundary.out").read_text(
        encoding="ascii")
    for needle in [
        "EMPTY_DELETION FAILURES 0",
        "SINGLETON_TOTAL 3925 UNSAFE 2 SAFE 3923 CRITICAL 024085d0 20a09640",
        "DOUBLE_TOTAL 7700850 UNSAFE 7869 SAFE 7692981 "
        "CRITICAL_MASK_PAIRS 7847 EXTRA_TWOHIT_PAIRS 22",
        "EXACT_QUANTIFIER EVERY_DELETION_SET_D_SUBSET_OF_FIXED_CARRIER_"
        "WITH_CARDINALITY_AT_MOST_TWO",
        "NO_CLASSIFICATION_FOR_LARGER_DELETIONS",
        "VERDICT PASS",
    ]:
        require(needle in boundary, f"deletion boundary missing: {needle}")

    identity = (packet / "endpoint591_O2_O3_identity.out").read_text(
        encoding="ascii")
    identity_specs = [
        ("ALL_TRANSCRIPT", "endpoint591_all_witness_O3.out"),
        ("ALL_LOW_LEDGER", "endpoint591_all_witness_low_O3.csv"),
        ("PROTECTED_TRANSCRIPT", "endpoint591_twohit_O3.out"),
        ("PROTECTED_BODY_LEDGER", "endpoint591_twohit_O3.csv"),
        ("PROTECTED_PAIR_LEDGER", "endpoint591_twohit_pairs_O3.csv"),
        ("CAPACITY_TRANSCRIPT", "endpoint591_raw_O3.out"),
        ("CAPACITY_PAIR_LEDGER", "endpoint591_pair_O3.csv"),
        ("CAPACITY_FAILURE_LEDGER", "endpoint591_failures_O3.csv"),
    ]
    for label, name in identity_specs:
        require(f"BYTE_IDENTICAL {label} SHA256 {sha256(packet / name)}" in
                identity, f"identity transcript mismatch: {label}")
    for label, name, rows in [
        ("DERIVED_CRITICAL_SINGLETONS",
         "endpoint591_critical_singletons.csv", " ROWS 2"),
        ("DERIVED_DISTINCT_TWOHIT_PAIRS",
         "endpoint591_distinct_twohit_witness_pairs.csv", " ROWS 23"),
        ("DERIVED_DELETION_BOUNDARY",
         "endpoint591_deletion_boundary.out", ""),
    ]:
        require(f"{label} SHA256 {sha256(packet / name)}{rows}" in identity,
                f"derived identity mismatch: {label}")
    require(identity.endswith("VERDICT PASS\n"),
            "identity transcript lacks verdict")

    old = repo / "05-knowledge/results/lrc14_mixed_rank_depth_recursive_signatures_thm4296"
    universe = read_rows(old / "inputs/current_residual22647.csv", 22647)
    prior_union = read_rows(packet4313 / "typed/typed_union2087.csv", 2087)
    prior_residual = read_rows(packet4313 / "typed/final_residual20560.csv", 20560)
    union = read_rows(packet / "typed/typed_union2100.csv", 2100)
    residual = read_rows(packet / "typed/final_residual20547.csv", 20547)
    top590 = read_rows(packet / "typed/residual_top590.csv", 13)
    require(fnv_rows(universe) == 0xDF5374D4ACA67677,
            "universe FNV changed")
    require(fnv_rows(prior_union) == 0x23E4136827B770A5,
            "prior union FNV changed")
    require(fnv_rows(prior_residual) == 0x8D797592A729E0B3,
            "prior residual FNV changed")
    require(set(prior_union).isdisjoint(prior_residual),
            "prior partition overlaps")
    require(set(prior_union) | set(prior_residual) == set(universe),
            "prior partition differs from universe")
    require(union == sorted(set(prior_union) | set(top591)),
            "successor union is not prior union plus top591")
    require(residual == sorted(set(prior_residual) - set(top591)),
            "successor residual is not prior residual minus top591")
    require(set(union).isdisjoint(residual), "successor partition overlaps")
    require(set(union) | set(residual) == set(universe),
            "successor partition differs from universe")
    require(fnv_rows(union) == 0x3B2D991DA091A7DF,
            "successor union FNV changed")
    require(fnv_rows(residual) == 0x59CA49A11D140EC5,
            "successor residual FNV changed")
    require(max(r for _, r in residual) == 590,
            "successor endpoint changed")
    require(top590 == [row for row in residual if row[1] == 590],
            "successor top590 ledger is incomplete")
    require(fnv_rows(top590) == 0x44AA8A793D162CF9,
            "successor top590 FNV changed")

    expected_hashes = {
        "endpoint591_raw_O3.out": "41d58cb62789d70b0fd53086db0201a272d330ff5b2a1a25b14a5442d54ed8c8",
        "endpoint591_pair_O3.csv": "926182f74120987aa09c5d0d7d7607afb7a9defea8d8b69322f8c7299d7f5c6f",
        "endpoint591_failures_O3.csv": "1402c9ebace944866d1f97ce9add004fef6e912c5c807a4e21de374c9372d188",
        "typed_endpoint591_consumer.out": "95ffd4a481c38b60f6fab778728a67ab116bd364744efd89a6331ecf06b94fc5",
    }
    for name, expected in expected_hashes.items():
        require(sha256(packet / name) == expected,
                f"frozen artifact SHA changed: {name}")
    require(sha256(packet / "typed/typed_union2100.csv") ==
            "1027325378caa6d4853112fcdb006796c32180b71f82a5a7db8addbec821f01c",
            "successor union SHA changed")
    require(sha256(packet / "typed/final_residual20547.csv") ==
            "b53abd545ddd28e95088231f615229a3fbfb0812510f84dc4fdfeda38e76f0c3",
            "successor residual SHA changed")
    require(sha256(packet / "typed/residual_top590.csv") ==
            "eefb84e9fec0bcd809830a3c5ed18b0bb2aea1a2f2cb8b4994fafeaf1f5d01ce",
            "successor top590 SHA changed")

    typed = (packet / "typed_endpoint591_consumer.out").read_text(
        encoding="ascii")
    for needle in [
        "NEW_UNION 2100 FNV 3b2d991da091a7df",
        "NEW_RESIDUAL 20547 FNV 59ca49a11d140ec5",
        "NEW_TOP 590 ROWS 13 FNV 44aa8a793d162cf9",
        "NOT_PHYSICAL_ENTRY_NO_LRC14",
        "VERDICT PASS",
    ]:
        require(needle in typed, f"typed transcript missing: {needle}")

    print("ENDPOINT591_PACKET_VERIFY PASS")
    print("target_rows=13 body_tests=185992950 failures=0")
    print("singleton_unsafe=2 pair_unsafe=7869 pair_safe=7692981")
    print("typed_union=2100 residual=20547 next_endpoint=590 next_rows=13")


if __name__ == "__main__":
    main()
