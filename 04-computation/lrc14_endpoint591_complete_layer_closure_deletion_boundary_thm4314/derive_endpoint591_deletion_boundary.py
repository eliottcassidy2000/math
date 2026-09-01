#!/usr/bin/env python3
"""Derive and verify the exact endpoint-591 deletion boundary artifacts.

This verifier is deliberately certificate-only: it consumes the independent
O2/O3 exhaustive census artifacts, checks byte identity and fixed semantic
invariants, and derives the critical-singleton and two-witness-pair ledgers.
It does not rerun the 185,992,950-case carrier scans.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_SHA256 = {
    "all transcript": "47c23d0257e0939380c5593ac504a3c68f1ce36b73eb6cfef028f8b2f2a360aa",
    "all low ledger": "a83bf1c001da47e878404f5a8828a51aea1a0aa0ba0fca021e9d955ae9b54301",
    "protected transcript": "b19e37e576a38312f98820ba3a1e11e24f72cdcb475feed224ca4c97a635c5ef",
    "protected body ledger": "484cb07b030d93f9b667b507e10f8cfe29314f837c0c1c1c955df3ddf74f529d",
    "protected pair ledger": "32e26f3bba3bfa1697d4113a2c4b1c3a87fad9c4ca4ff4cf506df56597553569",
    "capacity transcript": "41d58cb62789d70b0fd53086db0201a272d330ff5b2a1a25b14a5442d54ed8c8",
    "capacity pair ledger": "926182f74120987aa09c5d0d7d7607afb7a9defea8d8b69322f8c7299d7f5c6f",
    "capacity failure ledger": "1402c9ebace944866d1f97ce9add004fef6e912c5c807a4e21de374c9372d188",
}

ROWS = (72, 96, 100, 105, 153, 192, 210, 256, 260, 294, 366, 384, 520)
CRITICAL_EXPECTED = {0x024085D0, 0x20A09640}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_identical(label: str, first: Path, second: Path) -> str:
    first_bytes = first.read_bytes()
    second_bytes = second.read_bytes()
    assert first_bytes == second_bytes, f"{label}: O2/O3 byte mismatch"
    digest = hashlib.sha256(first_bytes).hexdigest()
    assert digest == EXPECTED_SHA256[label], f"{label}: frozen SHA changed"
    return digest


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--joint",
        type=Path,
        required=True,
    )
    args = parser.parse_args()
    inputs = args.input_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    identity_specs = (
        ("all transcript", inputs / "endpoint591_all_witness_O2.out",
         inputs / "endpoint591_all_witness_O3.out"),
        ("all low ledger", inputs / "endpoint591_all_witness_low_O2.csv",
         inputs / "endpoint591_all_witness_low_O3.csv"),
        ("protected transcript", inputs / "endpoint591_twohit_O2.out",
         inputs / "endpoint591_twohit_O3.out"),
        ("protected body ledger", inputs / "endpoint591_twohit_O2.csv",
         inputs / "endpoint591_twohit_O3.csv"),
        ("protected pair ledger", inputs / "endpoint591_twohit_pairs_O2.csv",
         inputs / "endpoint591_twohit_pairs_O3.csv"),
        ("capacity transcript", inputs / "endpoint591_raw_O2.out",
         inputs / "endpoint591_raw_O3.out"),
        ("capacity pair ledger", inputs / "endpoint591_pair_O2.csv",
         inputs / "endpoint591_pair_O3.csv"),
        ("capacity failure ledger", inputs / "endpoint591_failures_O2.csv",
         inputs / "endpoint591_failures_O3.csv"),
    )
    identities = [(label, check_identical(label, o2, o3))
                  for label, o2, o3 in identity_specs]

    all_text = (inputs / "endpoint591_all_witness_O3.out").read_text(encoding="utf-8")
    protected_text = (inputs / "endpoint591_twohit_O3.out").read_text(encoding="utf-8")
    capacity_text = (inputs / "endpoint591_raw_O3.out").read_text(encoding="utf-8")
    assert "CARRIER 3925 FNV a0d08a38c10bdab7 ROWS 13 BODY_TESTS 185992950" in all_text
    assert "SUMMARY ZERO 0 ONE 2 TWO 26 LOW_FNV 8dad3ce6f3e587c2 " \
           "DISTINCT_TWO_PAIRS 23 MAX_PAIR_LOAD 3 " \
           "LEAST_MAX_PAIR 002016c9,02ac0409" in all_text
    assert "CARRIER 3925 FNV a0d08a38c10bdab7 ROWS 13 EXPOSED 28791 " \
           "BELOW_TWO 0" in protected_text
    assert "TOTAL_TWOHIT 12 CENSUS_FNV 121f6e5d02874855 " \
           "DISTINCT_WITNESS_PAIRS 12 MAX_PAIR_LOAD 1" in protected_text
    assert "ROWS 13 ENDPOINT 591 ROW_FNV fc332c0697c671c7" in capacity_text
    assert "SUMMARY EXPOSED 28791 HIT_INCIDENCES 1031512 FAILURES 0 " \
           "FAILED_ROWS 0 FAILURE_FNV cbf29ce484222325 " \
           "PAIR_LEDGER_FNV 8191899a3e142a2c" in capacity_text
    assert (inputs / "endpoint591_failures_O3.csv").read_text(
        encoding="utf-8"
    ) == "q,r,body_hex\n"

    joint = {
        int(line.strip(), 16)
        for line in args.joint.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    assert len(joint) == 421

    low = read_csv(inputs / "endpoint591_all_witness_low_O3.csv")
    assert len(low) == 28
    assert all(int(row["r"]) == 591 and int(row["q"]) in ROWS for row in low)
    assert all(int(row["body_hex"], 16).bit_count() == 9 for row in low)
    assert all(int(row["hits"]) in (1, 2) for row in low)

    one = [row for row in low if int(row["hits"]) == 1]
    two = [row for row in low if int(row["hits"]) == 2]
    assert len(one) == 2 and len(two) == 26
    critical_obligations: dict[int, list[str]] = defaultdict(list)
    for row in one:
        witness = int(row["first_mask_hex"], 16)
        assert witness.bit_count() in (8, 9)
        assert not row["second_mask_hex"]
        assert witness & int(row["body_hex"], 16) == 0
        critical_obligations[witness].append(
            f'{row["q"]}:{row["r"]}:{row["body_hex"]}'
        )
    critical = set(critical_obligations)
    assert critical == CRITICAL_EXPECTED and critical <= joint

    pair_obligations: dict[tuple[int, int], list[str]] = defaultdict(list)
    pair_type_by_obligation: Counter[str] = Counter()
    for row in two:
        body = int(row["body_hex"], 16)
        pair = tuple(sorted((int(row["first_mask_hex"], 16),
                             int(row["second_mask_hex"], 16))))
        assert pair[0] != pair[1]
        assert all(mask.bit_count() in (8, 9) and mask & body == 0 for mask in pair)
        pair_obligations[pair].append(f'{row["q"]}:{row["r"]}:{row["body_hex"]}')
        joint_count = sum(mask in joint for mask in pair)
        pair_type_by_obligation[{0: "NN", 1: "JN", 2: "JJ"}[joint_count]] += 1
    assert len(pair_obligations) == 23
    assert pair_type_by_obligation == Counter({"NN": 12, "JN": 10, "JJ": 4})
    extra_pairs = {pair for pair in pair_obligations if not critical.intersection(pair)}
    assert len(extra_pairs) == 22

    protected = read_csv(inputs / "endpoint591_twohit_O3.csv")
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
        if int(row["first_mask_hex"], 16) not in joint
        and int(row["second_mask_hex"], 16) not in joint
    }
    assert protected_keys == nn_keys and len(nn_keys) == 12

    singleton_rows = [
        {
            "critical_mask_hex": f"{mask:08x}",
            "rank": mask.bit_count(),
            "protected_joint": 1,
            "obligation_count": len(critical_obligations[mask]),
            "obligations_q_r_body": ";".join(critical_obligations[mask]),
        }
        for mask in sorted(critical)
    ]
    singleton_path = output / "endpoint591_critical_singletons.csv"
    write_csv(singleton_path, list(singleton_rows[0]), singleton_rows)

    pair_rows: list[dict[str, object]] = []
    for pair, obligations in sorted(pair_obligations.items()):
        joint_count = sum(mask in joint for mask in pair)
        pair_rows.append({
            "first_mask_hex": f"{pair[0]:08x}",
            "second_mask_hex": f"{pair[1]:08x}",
            "pair_type": {0: "NN", 1: "JN", 2: "JJ"}[joint_count],
            "contains_singleton_critical_mask": int(bool(critical.intersection(pair))),
            "additional_unsafe_pair": int(not critical.intersection(pair)),
            "obligation_count": len(obligations),
            "obligations_q_r_body": ";".join(obligations),
        })
    pair_path = output / "endpoint591_distinct_twohit_witness_pairs.csv"
    write_csv(pair_path, list(pair_rows[0]), pair_rows)

    carrier_size = 3925
    total_pairs = math.comb(carrier_size, 2)
    critical_pairs = total_pairs - math.comb(carrier_size - len(critical), 2)
    unsafe_pairs = critical_pairs + len(extra_pairs)
    assert (total_pairs, critical_pairs, unsafe_pairs) == (7_700_850, 7_847, 7_869)

    boundary_path = output / "endpoint591_deletion_boundary.out"
    boundary_text = (
        "LRC14_ENDPOINT591_DELETION_BOUNDARY_V2\n"
        "UNIVERSE FIXED_ORDERED_13_ROWS_X_ALL_RANK9_BODIES "
        "BODY_TESTS 185992950 FIXED_CARRIER 3925 FNV a0d08a38c10bdab7\n"
        "WITNESS ACTIVE_CARRIER_MASK_OF_RANK8_OR9_DISJOINT_FROM_BODY\n"
        "CLOSURE_AFTER_DELETION FOR_EVERY_ROW_BODY_AT_LEAST_ONE_"
        "UNDELETED_WITNESS\n"
        "EMPTY_DELETION FAILURES 0\n"
        "SINGLETON_TOTAL 3925 UNSAFE 2 SAFE 3923 CRITICAL "
        + " ".join(f"{mask:08x}" for mask in sorted(critical)) + "\n"
        f"DOUBLE_TOTAL {total_pairs} UNSAFE {unsafe_pairs} "
        f"SAFE {total_pairs - unsafe_pairs} CRITICAL_MASK_PAIRS {critical_pairs} "
        f"EXTRA_TWOHIT_PAIRS {len(extra_pairs)}\n"
        f"TWOHIT_OBLIGATIONS {len(two)} DISTINCT_TWOHIT_PAIRS "
        f"{len(pair_obligations)} PAIR_TYPES_BY_OBLIGATION JJ:4 JN:10 NN:12\n"
        "MAX_PAIR_LOAD 3 PAIR 002016c9,02ac0409 OBLIGATIONS "
        "96:591:3512e000 96:591:3d10e000 192:591:3512e000\n"
        "EXACT_QUANTIFIER EVERY_DELETION_SET_D_SUBSET_OF_FIXED_CARRIER_"
        "WITH_CARDINALITY_AT_MOST_TWO\n"
        "NONCLAIMS NO_CLASSIFICATION_FOR_LARGER_DELETIONS NO_PHYSICAL_ENTRY "
        "NO_LRC14\n"
        "VERDICT PASS\n"
    )
    boundary_path.write_bytes(boundary_text.encode("ascii"))

    identity_path = output / "endpoint591_O2_O3_identity.out"
    identity_lines = ["LRC14_ENDPOINT591_O2_O3_IDENTITY_V1"]
    identity_lines.extend(
        f"BYTE_IDENTICAL {label.replace(' ', '_').upper()} SHA256 {digest}"
        for label, digest in identities
    )
    identity_lines.extend([
        f"DERIVED_CRITICAL_SINGLETONS SHA256 {sha256(singleton_path)} ROWS 2",
        f"DERIVED_DISTINCT_TWOHIT_PAIRS SHA256 {sha256(pair_path)} ROWS 23",
        f"DERIVED_DELETION_BOUNDARY SHA256 {sha256(boundary_path)}",
        "VERDICT PASS",
    ])
    identity_path.write_bytes(("\n".join(identity_lines) + "\n").encode("ascii"))


if __name__ == "__main__":
    main()
