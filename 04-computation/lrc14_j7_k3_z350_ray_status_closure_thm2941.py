#!/usr/bin/env python3
"""All-label closure of the 53 projected k=3 scalar rows at z1=350.

The exact residue-ray law and its reversal sidecar turn each denominator
class into finitely many attained nonnegative hyperbolic rays.  A two-state
max-plus merge therefore computes the exact all-label scalar maximum while
retaining the projected high-label obligation.  Every surviving denominator
state is then rejected either by a crude exact fibre overload or by one
common 16-cell Hunter status table.  Every status rejection carries an exact
rational Farkas certificate, replayed here by a second checker.

No finite label horizon or omitted-tail estimate occurs after the imported
53-body scalar candidate supplier.
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
SCALAR_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z350_scalar_slice_thm2941.py"
)
SCALAR_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z350_scalar_slice_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z350_ray_status_closure_thm2941.out"
)
EXPECTED_RAY_ENGINE_SHA256 = (
    "6ff8676255d51d818d7c24102a8fc755e673544f0ac6b99be4bfc262c892df1e"
)
EXPECTED_SCALAR_SHA256 = (
    "5aaa5a70e71c4408bae13e2bd3fc1dffb7f6b271fd59cc8feb2ef4aea37205cf"
)
EXPECTED_SCALAR_OUTPUT_SHA256 = (
    "ad00a09b68620d495b991cd0662863407808d794cfcf4ed22917bcb1b790affa"
)
EXPECTED_SEMANTIC_SHA256 = (
    "5a47df48d852157d55cb7d7aa11d69d9ebd04493f0d905020f8bd0a335c1c3d8"
)
FIRST = 350
EXPECTED_COUNTS = {
    (1, 2, 3, 8, 12, 14): (1, 0, 1),
    (1, 2, 4, 8, 10, 14): (1, 0, 1),
    (1, 2, 4, 8, 12, 14): (5, 0, 5),
    (1, 2, 4, 10, 12, 14): (2, 0, 2),
    (1, 2, 5, 8, 12, 14): (1, 1, 0),
    (1, 2, 6, 8, 10, 14): (2, 1, 1),
    (1, 2, 6, 8, 12, 14): (12, 5, 7),
    (1, 2, 6, 9, 12, 14): (5, 0, 5),
    (1, 2, 6, 10, 12, 14): (23, 1, 22),
    (1, 2, 8, 10, 12, 14): (88, 21, 67),
    (1, 3, 4, 8, 12, 14): (1, 0, 1),
    (1, 4, 5, 8, 12, 14): (17, 0, 17),
    (1, 4, 6, 8, 10, 14): (4, 3, 1),
    (1, 4, 6, 8, 12, 14): (18, 8, 10),
    (1, 4, 6, 9, 12, 14): (1, 0, 1),
    (1, 4, 6, 10, 12, 14): (17, 2, 15),
    (1, 4, 8, 10, 12, 14): (259, 60, 199),
    (1, 4, 10, 11, 12, 14): (0, 0, 0),
    (1, 5, 6, 8, 12, 14): (6, 6, 0),
    (1, 6, 8, 9, 12, 14): (5, 2, 3),
    (1, 6, 8, 10, 12, 14): (115, 51, 64),
    (1, 6, 8, 10, 13, 14): (0, 0, 0),
    (1, 8, 10, 11, 12, 14): (119, 0, 119),
    (2, 3, 4, 8, 12, 14): (3, 0, 3),
    (2, 3, 5, 8, 12, 14): (1, 1, 0),
    (2, 3, 6, 8, 12, 14): (2, 0, 2),
    (2, 3, 8, 10, 12, 14): (47, 19, 28),
    (2, 3, 8, 11, 12, 14): (6, 0, 6),
    (2, 4, 5, 8, 10, 14): (1, 0, 1),
    (2, 4, 5, 8, 12, 14): (4, 0, 4),
    (2, 4, 6, 8, 10, 14): (24, 1, 23),
    (2, 4, 6, 8, 12, 14): (21, 3, 18),
    (2, 4, 6, 9, 12, 14): (1, 0, 1),
    (2, 4, 6, 10, 12, 14): (19, 2, 17),
    (2, 4, 8, 10, 12, 14): (259, 35, 224),
    (2, 4, 10, 11, 12, 14): (0, 0, 0),
    (2, 5, 6, 8, 10, 14): (3, 1, 2),
    (2, 5, 6, 8, 12, 14): (52, 32, 20),
    (2, 5, 8, 10, 12, 14): (13, 6, 7),
    (2, 5, 8, 11, 12, 14): (119, 0, 119),
    (2, 6, 8, 9, 12, 14): (29, 13, 16),
    (2, 6, 8, 10, 12, 14): (256, 51, 205),
    (2, 6, 8, 10, 13, 14): (119, 98, 21),
    (2, 6, 8, 11, 12, 14): (42, 7, 35),
    (2, 6, 9, 10, 12, 14): (11, 0, 11),
    (2, 6, 10, 11, 12, 14): (141, 107, 34),
    (2, 8, 9, 10, 12, 14): (15, 7, 8),
    (2, 8, 10, 11, 12, 14): (1027, 597, 430),
    (2, 8, 10, 12, 13, 14): (119, 78, 41),
    (3, 4, 8, 10, 12, 14): (12, 4, 8),
    (4, 6, 8, 9, 12, 14): (1, 0, 1),
    (4, 6, 8, 10, 12, 14): (32, 12, 20),
    (4, 8, 10, 11, 12, 14): (119, 60, 59),
}
EXPECTED_TOTALS = (3_200, 1_295, 1_905, 0)
EXPECTED_M_HISTOGRAM = ((2, 68), (3, 46), (4, 125), (5, 1_274), (6, 91), (7, 301))


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
ray = load_module("z350_ray_engine", RAY_ENGINE_PATH)
scalar = load_module("z350_scalar_supplier", SCALAR_PATH)
ray.FIRST = FIRST
BODIES = scalar.EXPECTED_BODIES
require(tuple(EXPECTED_COUNTS) == BODIES, "body/count order mismatch")


def independent_farkas_check(q, marginals, capacities, histogram, certificate):
    """Rebuild and check the Farkas inequalities using exact rationals."""
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
        target_table = dict(stream.target_data(D))
        require(M == D // q and target == target_table[q], "crude target changed")
        require(
            capacity == sum(ray.local.fibre_cap(D, d, q) for d in ds),
            "crude capacity changed",
        )
        require(gap == target - capacity and gap > 0, "invalid crude witness")

    certificate_digest = hashlib.sha256()
    contradictions = []
    representative = None
    for ds, witness in sorted(status.items()):
        q, M, marginals, cap_set, histogram, certificate = witness
        D = lcm(*ds)
        require(M == D // q, "status modulus changed")
        rebuilt_marginals, capacities = ray.local.hunter_status_data(D, ds, q)
        require(rebuilt_marginals == marginals, "status marginals changed")
        require(tuple(sorted(set(capacities))) == cap_set, "status caps changed")
        arcs = ray.fibre.projected_support_arcs(D, stream.ranges)
        rebuilt_histogram = ray.fibre.residue_load_histogram(arcs, q)
        require(rebuilt_histogram == histogram, "status histogram changed")
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
    status_histogram = tuple(
        sorted(Counter(witness[1] for witness in status.values()).items())
    )
    maximum = max(
        states.items(), key=lambda item: (item[1]["excess"], item[0]), default=None
    )
    partition = (
        tuple(sorted(states.items())),
        tuple(sorted(crude.items())),
        tuple(sorted(status.items())),
        tuple(survivors),
    )
    partition_digest = hashlib.sha256(repr(partition).encode()).hexdigest()
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
        status_histogram,
        maximum,
        partition_digest,
        certificate_digest.hexdigest(),
        min(contradictions, default=None),
        representative,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=max(1, mp.cpu_count() // 2))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    pair_rows, control_marginals, control_caps, control_certificate = ray.local.controls()
    with mp.Pool(args.processes) as pool:
        records = tuple(sorted(pool.imap_unordered(evaluate_body, BODIES, chunksize=1)))

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
        tuple(records),
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
        "LRC14 projected k=3 z1=350 all-label ray/status closure",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        f"scalar_source_sha256={file_sha256(SCALAR_PATH)}",
        f"scalar_output_sha256={file_sha256(SCALAR_OUTPUT_PATH)}",
        (
            "scope=53 exact scalar candidate bodies at z1=350;"
            "three distinct later nonaligned labels;no finite horizon"
        ),
        (
            "ray_law=(z+L)delta(z+L)=zdelta(z);"
            "A(L-b)=-A(b);all denominator maxima attained"
        ),
        f"pair_overlap_exhaustive_controls={pair_rows}",
        (
            "common_table_controls=loads(9,9):FEASIBLE;"
            "loads(14,8):EXACT_FARKAS_INFEASIBLE"
        ),
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
            f"partition_sha256={partition_digest};"
            f"certificate_sha256={certificate_digest};"
            f"min_exact_farkas_contradiction={min_contradiction}"
        )
    lines.extend(
        (
            f"representative_exact_farkas={representative}",
            "independent_exact_farkas_checks=1905/1905:PASS",
            "conclusion=all 53 projected k=3 z1=350 scalar rows are empty",
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
