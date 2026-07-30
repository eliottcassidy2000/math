#!/usr/bin/env python3
"""Exact all-label ray/status frontier of the 45 projected ``k=3,z1=324`` rows.

The lossless scalar body atlas supplies exactly 45 bodies at first drift 324.
For each body, the periodic residue-ray law computes the attained all-label
maximum in every denominator multiset.  Crude all-divisor fibre overload and
the common 16-cell Hunter table reject all but one of the resulting 2,346
states.  Every Hunter rejection is independently replayed below with exact
rational matrix arithmetic.

The sole residual is

    E=(2,8,10,11,12,14),  (d1,d2,d3,d4)=(3920,4620,10780,10780).

This file deliberately stops at that exact frontier; the adjacent projected
closure retains the unit direction erased by the denominator quotient and
closes the residual by a two-cell antipodal phase certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_ENGINE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
)
BODY_ATLAS_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
)
BODY_ATLAS_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z324_ray_status_frontier_thm2941.out"
)
EXPECTED_RAY_ENGINE_SHA256 = (
    "6ff8676255d51d818d7c24102a8fc755e673544f0ac6b99be4bfc262c892df1e"
)
EXPECTED_BODY_ATLAS_SHA256 = (
    "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded"
)
EXPECTED_BODY_ATLAS_OUTPUT_SHA256 = (
    "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
)
EXPECTED_SEMANTIC_SHA256 = (
    "10e91b0cb91f4c9f8e36cac8d46b6bed7bd62619e4b703d13f6cfeeb6f52c09f"
)

FIRST = 324
EXPECTED_BODY_COUNT = 45
EXPECTED_TOTALS = (2_346, 702, 1_643, 1)
EXPECTED_M_HISTOGRAM = (
    (2, 24),
    (3, 262),
    (4, 71),
    (5, 598),
    (6, 220),
    (7, 457),
    (12, 2),
    (14, 4),
    (21, 4),
    (55, 1),
)
RESIDUAL_BODY = (2, 8, 10, 11, 12, 14)
RESIDUAL_DS = (3920, 4620, 10780, 10780)


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
require(file_sha256(BODY_ATLAS_PATH) == EXPECTED_BODY_ATLAS_SHA256, "body atlas changed")
require(
    file_sha256(BODY_ATLAS_OUTPUT_PATH) == EXPECTED_BODY_ATLAS_OUTPUT_SHA256,
    "body atlas transcript changed",
)
ray = load_module("z324_ray_status_engine", RAY_ENGINE_PATH)
ray.FIRST = FIRST


def bodies_from_atlas():
    bodies = []
    for line in BODY_ATLAS_OUTPUT_PATH.read_text().splitlines():
        if not line.startswith("row=E=") or f";z1={FIRST};" not in line:
            continue
        body_text = line[6:].split(";", 1)[0]
        bodies.append(tuple(map(int, body_text.split(","))))
    require(
        len(bodies) == len(set(bodies)) == EXPECTED_BODY_COUNT,
        ("z324 body ledger changed", len(bodies), len(set(bodies))),
    )
    require(tuple(sorted(bodies)) == tuple(bodies), "z324 body order changed")
    return tuple(bodies)


def independent_farkas_check(q, marginals, capacities, histogram, certificate):
    """Rebuild the common-table inequalities and verify the Farkas row exactly."""
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
    require(all(value >= 0 for value in alpha), "negative Farkas multiplier")
    for pattern in range(16):
        slack = sum(z[row] * equality_rows[row][pattern] for row in range(5))
        slack -= sum(
            alpha[row] * tail_rows[row][pattern] for row in range(len(alpha))
        )
        require(slack >= 0, ("negative Farkas column", pattern, slack))
    contradiction = sum(z[row] * equality_rhs[row] for row in range(5))
    contradiction -= sum(
        alpha[row] * tail_rhs[row] for row in range(len(alpha))
    )
    require(contradiction < 0, ("nonnegative Farkas contradiction", contradiction))
    return contradiction


def evaluate_body(body):
    stream = ray.Stream(body)
    trials, states, checks, signs = ray.ray_quotient_states(stream)
    crude, status, survivors = ray.common_status_screen(stream, states)
    require(set(states) == set(crude) | set(status) | set(survivors), "bad partition")
    require(not (set(crude) & set(status)), "kill stages overlap")

    for ds, witness in crude.items():
        gap, q, M, target, capacity = witness
        D = lcm(*ds)
        require(M == D // q, "crude cofactor changed")
        require(target == dict(stream.target_data(D))[q], "crude target changed")
        require(
            capacity == sum(ray.local.fibre_cap(D, d, q) for d in ds),
            "crude capacity changed",
        )
        require(gap == target - capacity and gap > 0, "invalid crude witness")

    certificate_digest = hashlib.sha256()
    contradictions = []
    for ds, witness in sorted(status.items()):
        q, M, marginals, cap_set, histogram, certificate = witness
        D = lcm(*ds)
        require(M == D // q, "status cofactor changed")
        rebuilt_marginals, capacities = ray.local.hunter_status_data(D, ds, q)
        require(rebuilt_marginals == marginals, "status marginals changed")
        require(tuple(sorted(set(capacities))) == cap_set, "status caps changed")
        arcs = ray.fibre.projected_support_arcs(D, stream.ranges)
        require(
            ray.fibre.residue_load_histogram(arcs, q) == histogram,
            "status load histogram changed",
        )
        contradictions.append(
            independent_farkas_check(q, marginals, capacities, histogram, certificate)
        )
        certificate_digest.update(f"{body}|{ds}|{witness}\n".encode())

    sign_totals = {
        sign: sum(
            count
            for (_d, candidate_sign), count in signs.items()
            if candidate_sign == sign
        )
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], "ray sign imbalance")
    return (
        body,
        stream.L,
        checks,
        tuple(sorted(sign_totals.items())),
        trials,
        tuple(sorted(states)),
        tuple(sorted(crude)),
        tuple(sorted(status)),
        tuple(survivors),
        tuple(sorted(Counter(witness[1] for witness in status.values()).items())),
        certificate_digest.hexdigest(),
        min(contradictions, default=None),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=max(1, mp.cpu_count() // 2))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")
    bodies = bodies_from_atlas()
    pair_rows, control_marginals, control_caps, control_certificate = ray.local.controls()
    if args.processes == 1:
        records = tuple(evaluate_body(body) for body in bodies)
    else:
        with mp.get_context("spawn").Pool(args.processes) as pool:
            records = tuple(sorted(pool.imap_unordered(evaluate_body, bodies, chunksize=1)))

    totals = tuple(
        sum(len(record[index]) for record in records) for index in (5, 6, 7, 8)
    )
    require(totals == EXPECTED_TOTALS, ("global counts changed", totals))
    global_m = Counter()
    for record in records:
        global_m.update(dict(record[9]))
    require(tuple(sorted(global_m.items())) == EXPECTED_M_HISTOGRAM, "M histogram changed")
    residuals = tuple(
        (record[0], ds) for record in records for ds in record[8]
    )
    require(residuals == ((RESIDUAL_BODY, RESIDUAL_DS),), ("residual changed", residuals))
    semantic_hash = hashlib.sha256(repr(records).encode()).hexdigest()
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=324 exact all-label ray/status frontier",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        f"body_atlas_source_sha256={file_sha256(BODY_ATLAS_PATH)}",
        f"body_atlas_output_sha256={file_sha256(BODY_ATLAS_OUTPUT_PATH)}",
        "scope=45 lossless scalar-atlas bodies;no finite label horizon",
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        f"pair_overlap_exhaustive_controls={pair_rows}",
        (
            f"totals=states:{totals[0]};crude_kills:{totals[1]};"
            f"status_kills:{totals[2]};survivors:{totals[3]}"
        ),
        f"status_M_histogram={dict(EXPECTED_M_HISTOGRAM)}",
    ]
    for record in records:
        body, L, checks, signs, trials, states, crude, status, survivors, mhist, cdigest, minimum = record
        partition_digest = hashlib.sha256(
            repr((states, crude, status, survivors)).encode()
        ).hexdigest()
        lines.append(
            f"E={body};L={L};ray_checks={checks};ray_signs={dict(signs)};"
            f"denominator_trials={trials};states={len(states)};"
            f"crude_kills={len(crude)};status_kills={len(status)};"
            f"survivors={survivors};status_M={dict(mhist)};"
            f"partition_sha256={partition_digest};certificate_sha256={cdigest};"
            f"min_exact_farkas_contradiction={minimum}"
        )
    lines.extend(
        (
            f"sole_residual=E:{RESIDUAL_BODY};denominators:{RESIDUAL_DS}",
            "independent_exact_farkas_checks=1643/1643:PASS",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
