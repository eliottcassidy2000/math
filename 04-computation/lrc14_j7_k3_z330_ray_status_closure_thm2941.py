#!/usr/bin/env python3
"""Exact all-label closure of the unique projected ``k=3, z1=330`` row.

The scalar supplier reconstructs its body from all 3,003 six-body carriers.
The exact residue-ray quotient then has one surviving denominator state.  A
common 16-cell Hunter status table rejects it by an exact rational Farkas
certificate, replayed below by a second exact matrix checker independent of
the status solver.  No finite label horizon remains after the scalar handoff.
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
    "6ff8676255d51d818d7c24102a8fc755e673544f0ac6b99be4bfc262c892df1e"
)
EXPECTED_SCALAR_SHA256 = (
    "6d2a1a48a9b3f48e11f067431c9afb5da2a18ad9da2b69dc4e8580d6378d5082"
)
EXPECTED_SCALAR_OUTPUT_SHA256 = (
    "c2bf51b8b661d98d1ce2a8114c6187a6eefb17856e9684a46704a62ccecc81ed"
)
EXPECTED_SEMANTIC_SHA256 = (
    "2e2fbfa8e85c04d658c7ed04f44d7021a8f08abf7038be0ed56211980b0fe927"
)
FIRST = 330
EXPECTED_BODY = (1, 2, 6, 8, 12, 14)
EXPECTED_COUNTS = (1, 0, 1, 0)
EXPECTED_M = 7


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
    slacks, contradiction = independent_farkas_check(
        q, marginals, capacities, histogram, certificate
    )
    require(contradiction == -1, "normalized Farkas contradiction changed")

    sign_totals = {
        sign: sum(
            count
            for (_d, candidate_sign), count in signs.items()
            if candidate_sign == sign
        )
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], "ray sign imbalance")
    partition = (
        tuple(sorted(states.items())),
        tuple(sorted(crude.items())),
        tuple(sorted(status.items())),
        tuple(survivors),
    )
    partition_digest = hashlib.sha256(repr(partition).encode()).hexdigest()
    certificate_digest = hashlib.sha256(
        f"{EXPECTED_BODY}|{ds}|{witness}\n".encode()
    ).hexdigest()
    representative = (EXPECTED_BODY, ds, witness, slacks, contradiction)
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
        control_certificate,
        representative,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    maximum = max(states.items(), key=lambda item: (item[1]["excess"], item[0]))
    lines = [
        "LRC14 projected k=3 z1=330 all-label ray/status closure",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        "independent_audit=local exact 16-cell matrix rebuild",
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
            f"certificate_sha256={certificate_digest}"
        ),
        f"representative_exact_farkas={representative}",
        "independent_exact_farkas_checks=1/1:PASS",
        "conclusion=the unique projected k=3 z1=330 scalar row is empty",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
