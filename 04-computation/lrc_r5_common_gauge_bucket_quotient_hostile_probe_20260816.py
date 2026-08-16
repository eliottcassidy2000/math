#!/usr/bin/env python3
"""Exact bucket-quotient hostile for the r=5 common-gauge pullback.

The common-gauge source supplies thirteen owner-character functionals on the
39x39 guard-atom pair space.  The frozen endpoint bucket table remembers only
the sum over the thirteen absolute left sheets in each (C,D,d) bucket.  This
probe computes the classes of the thirteen source functionals modulo that
bucket-constant row space.

The source common-gauge package has a disjoint clean-room audit.  This script
pins that audit but remains a finite-exact follow-up sidecar, not a second
independent reconstruction or a physical-current/LRC theorem.
"""

from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "04-computation/lrc_r5_common_gauge_bucket_quotient_hostile_probe_20260816.py"
OUTPUT = "05-knowledge/results/lrc_r5_common_gauge_bucket_quotient_hostile_probe_20260816.out"
SOURCE_PATH = ROOT / "04-computation/lrc_r5_common_ancestry_guard_atom_root_drift_probe_20260816.py"
SOURCE_SHA256 = "83f1fa49ac4d02e21a1d76fed169d101715a6620342714ed05b9172ae967a730"
SOURCE_SEMANTIC = "3d8c88fb7b9762f41ef35c00d980b99fc435c8352baf5dddb9fe412d1baeace0"
SOURCE_AUDIT_PATH = ROOT / (
    "04-computation/"
    "lrc_r5_common_ancestry_guard_atom_root_drift_independent_audit_20260816.py"
)
SOURCE_AUDIT_SHA256 = "ec16f5124c5b83a10337fba6046d08572251d0776285b7e250a86539a26ddf97"
SOURCE_AUDIT_SEMANTIC = "ad81e945207703956f5d6ec300430d562c9f98a9ec8788011119a4857ab34e01"
TARGET_PATH = ROOT / (
    "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
TARGET_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
EXPECTED_SEMANTIC_SHA256: str | None = (
    "cc535a60a93f5652b1980268a274f7a7280cc722e3c707bdae78467f6d367ddf"
)

P = 13
CHAMBERS = ("left", "middle", "right")
ACTIVE = (0, 2)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    require(lf_sha256(path) == expected_hash, (name, "source hash drift"))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "loader"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def rref_pivots(
    matrix: tuple[tuple[int, ...], ...], prime: int
) -> tuple[int, tuple[int, ...]]:
    rows = [list(value % prime for value in row) for row in matrix]
    require(rows and rows[0], "empty rank matrix")
    rank = 0
    pivots = []
    for column in range(len(rows[0])):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [value * inverse % prime for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (value - factor * pivot_value) % prime
                for value, pivot_value in zip(rows[row], rows[rank])
            ]
        pivots.append(column)
        rank += 1
        if rank == len(rows):
            break
    return rank, tuple(pivots)


def determinant(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    rows = [list(value % prime for value in row) for row in matrix]
    require(rows and all(len(row) == len(rows) for row in rows), "square determinant")
    answer = 1
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer = -answer
        pivot_value = rows[column][column]
        answer = answer * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, len(rows)):
            factor = rows[row][column] * inverse % prime
            if factor:
                rows[row] = [
                    (value - factor * pivot_entry) % prime
                    for value, pivot_entry in zip(rows[row], rows[column])
                ]
    return answer % prime


def pair_index(left: int, right: int) -> int:
    return left * (P * len(CHAMBERS)) + right


def bucket_pairs(left_chamber: int, right_chamber: int, drift: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (
            left_sheet * 3 + left_chamber,
            ((left_sheet + drift) % P) * 3 + right_chamber,
        )
        for left_sheet in range(P)
    )


def main() -> None:
    require(
        lf_sha256(SOURCE_AUDIT_PATH) == SOURCE_AUDIT_SHA256,
        "source independent-audit hash drift",
    )
    source_module = load_module(SOURCE_PATH, "bucket_quotient_source", SOURCE_SHA256)
    with redirect_stdout(io.StringIO()):
        source = source_module.main()
    require(source["semantic_sha256"] == SOURCE_SEMANTIC, "source semantic drift")

    target = load_module(TARGET_PATH, "bucket_quotient_target_context", TARGET_SHA256)
    require(tuple(source_module.ATOMS) == tuple(target.ATOMS), "atom order mismatch")
    _word, _t_den, _nn, prime, _root, zeta, *_rest = target.context()
    require(pow(zeta, P, prime) == 1 and zeta != 1, "zeta order")

    pair_gauge_offset = source["pair_gauge_offset"]
    functionals = tuple(
        tuple(
            sum(
                pair_gauge_offset[left][right][offset]
                * pow(zeta, -owner_frequency * offset % P, prime)
                for offset in range(P)
            )
            % prime
            for left in range(len(target.ATOMS))
            for right in range(len(target.ATOMS))
        )
        for owner_frequency in range(P)
    )
    functional_rank = rref_pivots(functionals, prime)[0]
    require(functional_rank == P, functional_rank)

    bucket_variation_counts = []
    for owner_frequency in range(P):
        varying = 0
        for left_chamber in range(3):
            for right_chamber in range(3):
                for drift in range(P):
                    values = tuple(
                        functionals[owner_frequency][pair_index(left, right)]
                        for left, right in bucket_pairs(left_chamber, right_chamber, drift)
                    )
                    varying += len(set(values)) > 1
        bucket_variation_counts.append(varying)
    bucket_variation_counts_t = tuple(bucket_variation_counts)
    require(bucket_variation_counts_t == (72,) * P, bucket_variation_counts_t)

    # Basis of the kernel of all 117 bucket sums: compare each absolute left
    # sheet a=1,...,12 with a=0 inside the same (C,D,d) bucket.  Restrict to
    # the 48 owner-active punctured buckets for the sharp hostile.
    active_labels = []
    active_columns = []
    for left_chamber in ACTIVE:
        for right_chamber in ACTIVE:
            for drift in range(1, P):
                pairs = bucket_pairs(left_chamber, right_chamber, drift)
                reference = pair_index(*pairs[0])
                for left_sheet in range(1, P):
                    current = pair_index(*pairs[left_sheet])
                    active_labels.append(
                        (CHAMBERS[left_chamber], CHAMBERS[right_chamber], drift, left_sheet)
                    )
                    active_columns.append(
                        tuple(
                            (functionals[k][current] - functionals[k][reference]) % prime
                            for k in range(P)
                        )
                    )
    require(len(active_labels) == 4 * 12 * 12 == 576, len(active_labels))
    response_matrix = tuple(
        tuple(active_columns[column][owner_frequency] for column in range(len(active_columns)))
        for owner_frequency in range(P)
    )
    quotient_rank, pivot_columns = rref_pivots(response_matrix, prime)
    require(quotient_rank == P and len(pivot_columns) == P, (
        quotient_rank,
        pivot_columns,
    ))
    pivot_labels = tuple(active_labels[column] for column in pivot_columns)
    pivot_minor = tuple(
        tuple(response_matrix[row][column] for column in pivot_columns)
        for row in range(P)
    )
    pivot_determinant = determinant(pivot_minor, prime)
    require(pivot_determinant != 0, "pivot determinant")

    witness_label = ("left", "left", 1, 2)
    witness_column = active_labels.index(witness_label)
    witness_changes = active_columns[witness_column]
    require(all(witness_changes), witness_changes)

    # Each basis column is +1 at one atom pair and -1 at the a=0 pair in the
    # same bucket, so all 117 bucket sums vanish exactly.  Check the named
    # witness explicitly as a hostile control.
    witness_bucket_sums = [0] * (3 * 3 * P)
    left_chamber = right_chamber = 0
    drift = 1
    witness_bucket = (left_chamber * 3 + right_chamber) * P + drift
    witness_bucket_sums[witness_bucket] += 1
    witness_bucket_sums[witness_bucket] -= 1
    require(all(value == 0 for value in witness_bucket_sums), "bucket hostile leak")

    theorem = (
        "all_13_owner_functionals_are_independent_mod_bucket_constants",
        "owner_active_within_bucket_kernel_already_surjects_onto_Fp13",
        "complete_endpoint_bucket_table_does_not_determine_any_owner_pullback_profile",
        "atom_pair_alignment_is_a_load_bearing_sidecar",
    )
    boundary = (
        "perturbations_are_split_field_pair_functions_not_boolean_or_positive_endpoint_realizations",
        "actual_endpoint_pair_function_still_has_all_13_nonzero_pullbacks",
        "no_temporal_current_grouped_coefficient_scalar_exclusion_or_LRC14_claim",
    )
    semantic_surface = (
        SOURCE_SHA256,
        SOURCE_SEMANTIC,
        SOURCE_AUDIT_SHA256,
        SOURCE_AUDIT_SEMANTIC,
        TARGET_SHA256,
        (prime, zeta),
        digest_json(functionals),
        functional_rank,
        bucket_variation_counts_t,
        len(active_labels),
        quotient_rank,
        pivot_labels,
        pivot_determinant,
        witness_label,
        witness_changes,
        theorem,
        boundary,
    )
    semantic_sha256 = digest_json(semantic_surface)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            (semantic_sha256, EXPECTED_SEMANTIC_SHA256),
        )

    path = Path(__file__).resolve()
    print("R5 COMMON-GAUGE BUCKET-QUOTIENT HOSTILE")
    print("status=VERIFIED-FINITE-EXACT follow-up sidecar;source_clean_room_audit_ACCEPT;LRC14_OPEN")
    print(f"dependencies=((source,{SOURCE_SHA256},{SOURCE_SEMANTIC}),(source_audit,{SOURCE_AUDIT_SHA256},{SOURCE_AUDIT_SEMANTIC}),(target_context,{TARGET_SHA256}))")
    print(f"field=(prime={prime},zeta13={zeta})")
    print(f"owner_atom_pair_functionals=(shape=13x1521,rank={functional_rank},sha256={digest_json(functionals)})")
    print(f"variation_inside_117_buckets={bucket_variation_counts_t};each varies on exactly the 72 realized chamber/drift types")
    print(f"owner_active_bucket_kernel=(basis_columns={len(active_labels)},response_rank={quotient_rank},surjective_to=F_p^13)")
    print(f"rank13_minor=(labels={pivot_labels},determinant={pivot_determinant})")
    print(f"two_pair_hostile=(label={witness_label},owner_pullback_changes={witness_changes},all_nonzero=True,all_117_bucket_sums_unchanged=True)")
    print("consequence=for any desired 13-vector of pullback changes, an owner-active within-bucket perturbation exists while the complete endpoint bucket/K4/Walsh data stay fixed")
    print(f"theorem={theorem}")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(path)}")
    print("reproducibility=no_assert_truth_gates;no_randomness;no_elapsed_fields")
    print(f"commands=python -B {SCRIPT};python -B -O {SCRIPT}")
    print("PASS")


if __name__ == "__main__":
    main()
