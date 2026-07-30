#!/usr/bin/env python3
"""Exact all-label closure of the unique projected ``k=3, z1=330`` row.

The scalar supplier reconstructs its body from all 3,003 six-body carriers.
The exact residue-ray quotient then has one surviving denominator state.  A
common 16-cell Hunter status table rejects it.  The frozen proof is the
basis-independent singleton grade-three inequality

    1_[capacity >= 3] <= b_1+b_2+b_3,

whose marginal upper bound is 72 while the required grade-three tail is 144.
The solver-selected rational Farkas basis is still replayed internally as a
hostile control, but it is deliberately absent from every digest and displayed
certificate.  No finite label horizon remains after the scalar handoff.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_ENGINE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
)
SCALAR_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z330_scalar_slice_thm2941.py"
)
SCALAR_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z330_scalar_slice_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z330_ray_status_closure_thm2941.out"
)
EXPECTED_RAY_ENGINE_SHA256 = (
    "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2"
)
EXPECTED_SCALAR_SHA256 = (
    "5eb30248f2beac09b29ebccc0dc8b203361f97c4b8552021130290566f8de0d1"
)
EXPECTED_SCALAR_OUTPUT_SHA256 = (
    "c2bf51b8b661d98d1ce2a8114c6187a6eefb17856e9684a46704a62ccecc81ed"
)
EXPECTED_SEMANTIC_SHA256 = (
    "4d5b51acb7029cd94087e5c47830ec1fdd41e74b79d4efcd3ed72b4afc3e82ad"
)
FIRST = 330
EXPECTED_BODY = (1, 2, 6, 8, 12, 14)
EXPECTED_COUNTS = (1, 0, 1, 0)
EXPECTED_M = 7
EXPECTED_GRADE = 3


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(RAY_ENGINE_PATH) == EXPECTED_RAY_ENGINE_SHA256, "ray engine changed")
if EXPECTED_SCALAR_SHA256 is not None:
    require(file_sha256(SCALAR_PATH) == EXPECTED_SCALAR_SHA256, "scalar supplier changed")
if EXPECTED_SCALAR_OUTPUT_SHA256 is not None:
    require(
        file_sha256(SCALAR_OUTPUT_PATH) == EXPECTED_SCALAR_OUTPUT_SHA256,
        "scalar transcript changed",
    )
ray = load_module("z330_ray_engine", RAY_ENGINE_PATH)
scalar = load_module("z330_scalar_supplier", SCALAR_PATH)
ray.FIRST = FIRST
require(scalar.EXPECTED_BODIES == (EXPECTED_BODY,), "scalar/body handoff changed")


def independent_farkas_check(q, marginals, capacities, histogram, certificate):
    """Rebuild and verify the 16-cell Farkas inequalities exactly."""
    thresholds, alpha, z = certificate
    tail_rows = []
    tail_rhs = []
    rebuilt_thresholds = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        rebuilt_thresholds.append(threshold)
        tail_rows.append(good)
        tail_rhs.append(demand)
    equality_rows = [
        (1,) * 16,
        *[
            tuple((pattern >> index) & 1 for pattern in range(16))
            for index in range(4)
        ],
    ]
    equality_rhs = (q, *marginals)
    require(tuple(rebuilt_thresholds) == tuple(thresholds), "threshold mismatch")
    require(len(alpha) == len(tail_rows) and len(z) == 5, "certificate shape")
    require(all(value >= 0 for value in alpha), "negative Farkas alpha")
    slacks = []
    for pattern in range(16):
        value = sum(z[row] * equality_rows[row][pattern] for row in range(5))
        value -= sum(
            alpha[row] * tail_rows[row][pattern] for row in range(len(alpha))
        )
        slacks.append(value)
    contradiction = sum(z[row] * equality_rhs[row] for row in range(5))
    contradiction -= sum(
        alpha[row] * tail_rhs[row] for row in range(len(alpha))
    )
    require(all(value >= 0 for value in slacks), "negative Farkas column")
    require(contradiction < 0, "nonnegative Farkas contradiction")
    return tuple(slacks), contradiction


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    pair_rows, control_marginals, control_caps, control_certificate = ray.local.controls()
    stream = ray.Stream(EXPECTED_BODY)
    trials, states, checks, signs = ray.ray_quotient_states(stream)
    crude, status, survivors = ray.common_status_screen(stream, states)
    counts = (len(states), len(crude), len(status), len(survivors))
    require(counts == EXPECTED_COUNTS, ("closure counts changed", counts))
    require(set(states) == set(status), "unique state/status partition changed")

    ((ds, witness),) = tuple(sorted(status.items()))
    q, M, marginals, cap_set, histogram, certificate = witness
    D = lcm(*ds)
    require(M == D // q == EXPECTED_M, "status modulus changed")
    rebuilt_marginals, capacities = ray.local.hunter_status_data(D, ds, q)
    require(rebuilt_marginals == marginals, "status marginals changed")
    require(tuple(sorted(set(capacities))) == cap_set, "status caps changed")
    arcs = ray.fibre.projected_support_arcs(D, stream.ranges)
    require(
        ray.fibre.residue_load_histogram(arcs, q) == histogram,
        "status load histogram changed",
    )
    _solver_slacks, contradiction = independent_farkas_check(
        q, marginals, capacities, histogram, certificate
    )
    require(contradiction == -1, "normalized Farkas contradiction changed")

    grade = EXPECTED_GRADE
    demand = sum(count for load, count in histogram if load >= grade)
    good = tuple(int(capacity >= grade) for capacity in capacities)
    tariff = (0, 1, 1, 1, 0)
    tariff_slacks = tuple(
        tariff[0]
        + sum(tariff[index + 1] * ((pattern >> index) & 1) for index in range(4))
        - good[pattern]
        for pattern in range(16)
    )
    upper = tariff[0] * q + sum(
        tariff[index + 1] * marginals[index] for index in range(4)
    )
    deficit = demand - upper
    require(
        (demand, upper, deficit) == (144, 72, 72),
        ("singleton grade-three inequality changed", demand, upper, deficit),
    )
    require(all(value >= 0 for value in tariff_slacks), "invalid direct tariff")

    sign_totals = {
        sign: sum(
            count
            for (_d, candidate_sign), count in signs.items()
            if candidate_sign == sign
        )
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], "ray sign imbalance")
    deterministic_status = tuple(
        (status_ds, status_witness[:-1])
        for status_ds, status_witness in sorted(status.items())
    )
    partition = (
        tuple(sorted(states.items())),
        tuple(sorted(crude.items())),
        deterministic_status,
        tuple(survivors),
    )
    partition_digest = hashlib.sha256(repr(partition).encode()).hexdigest()
    deterministic_instance = (EXPECTED_BODY, ds, witness[:-1])
    instance_digest = hashlib.sha256(
        f"{deterministic_instance}\n".encode()
    ).hexdigest()
    representative = (
        EXPECTED_BODY,
        ds,
        witness[:-1],
        grade,
        demand,
        upper,
        deficit,
        tariff,
        tariff_slacks,
    )
    semantic_payload = (
        FIRST,
        EXPECTED_BODY,
        stream.h,
        len(stream.carrier),
        stream.L,
        stream.high_floor,
        stream.first_d,
        trials,
        checks,
        tuple(sorted(sign_totals.items())),
        partition,
        pair_rows,
        control_marginals,
        control_caps,
        representative,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    maximum = max(states.items(), key=lambda item: (item[1]["excess"], item[0]))
    lines = [
        "LRC14 projected k=3 z1=330 all-label ray/status closure",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        "independent_audit=basis-independent singleton barcode/tariff plus internal exact matrix replay",
        f"scalar_source_sha256={file_sha256(SCALAR_PATH)}",
        f"scalar_output_sha256={file_sha256(SCALAR_OUTPUT_PATH)}",
        "scope=1 globally reconstructed scalar body;no finite label horizon",
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        f"pair_overlap_exhaustive_controls={pair_rows}",
        (
            f"E={EXPECTED_BODY};h={stream.h};r={len(stream.carrier)};L={stream.L};"
            f"high={stream.high_floor};d1={stream.first_d};ray_checks={checks};"
            f"ray_signs={sign_totals};denominator_trials={trials};"
            f"states={counts[0]};crude_kills={counts[1]};"
            f"status_kills={counts[2]};survivors={counts[3]};M={M};"
            f"max_state={maximum};partition_sha256={partition_digest};"
            f"deterministic_instance_sha256={instance_digest}"
        ),
        f"representative_canonical_barcode_tariff={representative}",
        "canonical_circuit_census={(3,): 1};distinct_canonical_tariffs=1",
        "basis_independent_grade3_bound=H(3):144;maxF(3):72;deficit:72",
        "internal_solver_farkas_replay=1/1:PASS;solver_basis_not_frozen",
        "conclusion=the unique projected k=3 z1=330 scalar row is empty",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
