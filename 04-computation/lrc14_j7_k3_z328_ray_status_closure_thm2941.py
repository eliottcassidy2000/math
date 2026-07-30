#!/usr/bin/env python3
"""Exact all-label closure of the nine projected ``k=3, z1=328`` rows.

The scalar supplier reconstructs the nine bodies globally from all 3,003
six-body carriers.  The residue-ray/reversal law then computes the exact
infinite-label maximum for every denominator class.  Every remaining state is
rejected by either a crude fibre overload or the common 16-cell Hunter status
table.  Each status rejection is replayed below by a second exact rational
matrix checker independent of the status solver.  No finite label horizon
remains after the scalar handoff.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from collections import Counter
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_ENGINE_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
)
SCALAR_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z328_scalar_slice_thm2941.py"
)
SCALAR_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z328_scalar_slice_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z328_ray_status_closure_thm2941.out"
)
EXPECTED_RAY_ENGINE_SHA256 = (
    "6ff8676255d51d818d7c24102a8fc755e673544f0ac6b99be4bfc262c892df1e"
)
EXPECTED_SCALAR_SHA256 = (
    "917a05cdee397747edd94a0e926ff746b960a9ef16120fd7d8f139b42bc0f5c6"
)
EXPECTED_SCALAR_OUTPUT_SHA256 = (
    "f5c11f364a626141af181d84f39d48030ef91a8ddf7d74b9602bb15cd7eb626e"
)
EXPECTED_SEMANTIC_SHA256 = (
    "3b2ecca5cb0130f69a658ad70b7e1701e01bc22c677a76781377d82d5d20daea"
)
FIRST = 328
EXPECTED_COUNTS = {
    (1, 2, 6, 8, 12, 14): (1, 0, 1),
    (1, 2, 8, 10, 12, 14): (3, 0, 3),
    (1, 4, 6, 8, 12, 14): (1, 0, 1),
    (1, 4, 8, 10, 12, 14): (35, 16, 19),
    (1, 6, 8, 10, 12, 14): (2, 0, 2),
    (2, 4, 6, 8, 12, 14): (1, 0, 1),
    (2, 4, 8, 10, 12, 14): (8, 2, 6),
    (2, 6, 8, 10, 12, 14): (8, 0, 8),
    (2, 8, 10, 11, 12, 14): (26, 18, 8),
}
EXPECTED_TOTALS = (85, 36, 49, 0)
EXPECTED_M_HISTOGRAM = ((4, 2), (5, 22), (6, 4), (7, 21))


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
require(file_sha256(SCALAR_PATH) == EXPECTED_SCALAR_SHA256, "scalar supplier changed")
if EXPECTED_SCALAR_OUTPUT_SHA256 is not None:
    require(
        file_sha256(SCALAR_OUTPUT_PATH) == EXPECTED_SCALAR_OUTPUT_SHA256,
        "scalar transcript changed",
    )
ray = load_module("z328_ray_engine", RAY_ENGINE_PATH)
scalar = load_module("z328_scalar_supplier", SCALAR_PATH)
ray.FIRST = FIRST
BODIES = scalar.EXPECTED_BODIES
require(tuple(EXPECTED_COUNTS) == BODIES, "body/count order mismatch")


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


def evaluate_body(body):
    stream = ray.Stream(body)
    trials, states, checks, signs = ray.ray_quotient_states(stream)
    crude, status, survivors = ray.common_status_screen(stream, states)
    require(
        (len(states), len(crude), len(status)) == EXPECTED_COUNTS[body],
        ("body counts changed", body, len(states), len(crude), len(status)),
    )
    require(set(states) == set(crude) | set(status) | set(survivors), "bad partition")
    require(not (set(crude) & set(status)), "kill stages overlap")

    for ds, witness in crude.items():
        gap, q, M, target, capacity = witness
        D = lcm(*ds)
        require(M == D // q, "crude modulus changed")
        require(target == dict(stream.target_data(D))[q], "crude target changed")
        require(
            capacity == sum(ray.local.fibre_cap(D, d, q) for d in ds),
            "crude capacity changed",
        )
        require(gap == target - capacity and gap > 0, "invalid crude witness")

    contradictions = []
    certificate_digest = hashlib.sha256()
    representative = None
    for ds, witness in sorted(status.items()):
        q, M, marginals, cap_set, histogram, certificate = witness
        D = lcm(*ds)
        require(M == D // q, "status modulus changed")
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
        contradictions.append(contradiction)
        if representative is None:
            representative = (body, ds, witness, slacks, contradiction)
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
    partition = (
        tuple(sorted(states.items())),
        tuple(sorted(crude.items())),
        tuple(sorted(status.items())),
        tuple(survivors),
    )
    return (
        body,
        stream.h,
        len(stream.carrier),
        stream.L,
        stream.high_floor,
        stream.first_d,
        checks,
        tuple(sorted(sign_totals.items())),
        len(ray.support.divisors(stream.L)) - 1,
        trials,
        partition,
        tuple(sorted(Counter(w[1] for w in status.values()).items())),
        max(states.items(), key=lambda item: (item[1]["excess"], item[0]), default=None),
        hashlib.sha256(repr(partition).encode()).hexdigest(),
        certificate_digest.hexdigest(),
        min(contradictions, default=None),
        representative,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    pair_rows, control_marginals, control_caps, control_certificate = ray.local.controls()
    records = tuple(evaluate_body(body) for body in BODIES)
    totals = tuple(
        sum(len(record[10][index]) for record in records) for index in range(4)
    )
    require(totals == EXPECTED_TOTALS, ("global counts changed", totals))
    global_m = Counter()
    for record in records:
        global_m.update(dict(record[11]))
    require(tuple(sorted(global_m.items())) == EXPECTED_M_HISTOGRAM, "M histogram changed")
    representative = next(record[16] for record in records if record[16] is not None)
    semantic_payload = (
        FIRST,
        BODIES,
        records,
        totals,
        tuple(sorted(global_m.items())),
        pair_rows,
        control_marginals,
        control_caps,
        control_certificate,
        representative,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=328 all-label ray/status closure",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        "independent_audit=local exact 16-cell matrix rebuild",
        f"scalar_source_sha256={file_sha256(SCALAR_PATH)}",
        f"scalar_output_sha256={file_sha256(SCALAR_OUTPUT_PATH)}",
        "scope=9 globally reconstructed scalar bodies;no finite label horizon",
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        f"pair_overlap_exhaustive_controls={pair_rows}",
        (
            f"totals=states:{totals[0]};crude_kills:{totals[1]};"
            f"status_kills:{totals[2]};survivors:{totals[3]}"
        ),
        f"status_M_histogram={dict(sorted(global_m.items()))}",
    ]
    for record in records:
        (
            body,
            h,
            components,
            L,
            high_floor,
            first_d,
            checks,
            signs,
            denominator_classes,
            trials,
            partition,
            status_histogram,
            maximum,
            partition_digest,
            certificate_digest,
            min_contradiction,
            _representative,
        ) = record
        states, crude, status, survivors = partition
        lines.append(
            f"E={body};h={h};r={components};L={L};high={high_floor};d1={first_d};"
            f"ray_checks={checks};ray_signs={dict(signs)};"
            f"denominator_classes={denominator_classes};denominator_trials={trials};"
            f"states={len(states)};crude_kills={len(crude)};"
            f"status_kills={len(status)};survivors={len(survivors)};"
            f"status_M={dict(status_histogram)};max_state={maximum};"
            f"partition_sha256={partition_digest};certificate_sha256={certificate_digest};"
            f"min_exact_farkas_contradiction={min_contradiction}"
        )
    lines.extend(
        (
            f"representative_exact_farkas={representative}",
            "independent_exact_farkas_checks=49/49:PASS",
            "conclusion=all 9 projected k=3 z1=328 scalar rows are empty",
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
