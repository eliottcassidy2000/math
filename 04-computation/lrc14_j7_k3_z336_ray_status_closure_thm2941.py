#!/usr/bin/env python3
"""All-label closure of the eight projected ``k=3, z1=336`` rows.

The imported scalar supplier reconstructs the eight bodies globally from all
3,003 six-body carriers.  The exact residue-ray law and reversal sidecar then
replace every infinite denominator class by finitely many attained rays.  A
two-state max-plus merge computes the exact all-label scalar maximum.  Every
remaining denominator state is rejected either by a crude fibre overload or
by one common 16-cell Hunter status table.  All status rejections carry exact
rational Farkas certificates, replayed below by an independent checker.

Thus no finite label horizon or omitted-tail estimate is used after the scalar
candidate handoff.
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
    ROOT / "04-computation" / "lrc14_j7_k3_z336_scalar_slice_thm2941.py"
)
SCALAR_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z336_scalar_slice_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z336_ray_status_closure_thm2941.out"
)
EXPECTED_RAY_ENGINE_SHA256 = (
    "6ff8676255d51d818d7c24102a8fc755e673544f0ac6b99be4bfc262c892df1e"
)
EXPECTED_SCALAR_SHA256 = (
    "8ca30b3f6e8336e87b558c5e24007153aa0f0331e3b00dfe66c93cd4c93dd789"
)
EXPECTED_SCALAR_OUTPUT_SHA256 = (
    "e3de2a3d6c65c0e8eb62d400ff47deafe48c8aa6a2db45788eab45e8f50892b9"
)
EXPECTED_SEMANTIC_SHA256 = (
    "6ec3a326ec7fa66053150954b1f2a9d3910fbe044032de10c6ca57ff0e53d228"
)
FIRST = 336
EXPECTED_COUNTS = {
    (1, 2, 6, 8, 12, 14): (1, 1, 0),
    (1, 2, 6, 10, 12, 14): (1, 1, 0),
    (1, 4, 6, 8, 12, 14): (1, 1, 0),
    (1, 4, 8, 10, 12, 14): (13, 8, 5),
    (2, 4, 6, 8, 12, 14): (2, 1, 1),
    (2, 4, 8, 10, 12, 14): (8, 3, 5),
    (2, 6, 8, 10, 12, 14): (12, 5, 7),
    (2, 8, 10, 11, 12, 14): (71, 51, 20),
}
EXPECTED_TOTALS = (109, 71, 38, 0)
EXPECTED_M_HISTOGRAM = ((4, 8), (5, 18), (6, 3), (7, 9))


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
ray = load_module("z336_ray_engine", RAY_ENGINE_PATH)
scalar = load_module("z336_scalar_supplier", SCALAR_PATH)
ray.FIRST = FIRST
BODIES = scalar.EXPECTED_BODIES
require(tuple(EXPECTED_COUNTS) == BODIES, "body/count order mismatch")


def independent_farkas_check(q, marginals, capacities, histogram, certificate):
    """Rebuild and verify the status-table Farkas inequalities exactly."""
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
        "LRC14 projected k=3 z1=336 all-label ray/status closure",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        f"scalar_source_sha256={file_sha256(SCALAR_PATH)}",
        f"scalar_output_sha256={file_sha256(SCALAR_OUTPUT_PATH)}",
        (
            "scope=8 globally reconstructed scalar candidate bodies at z1=336;"
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
            "independent_exact_farkas_checks=38/38:PASS",
            "conclusion=all 8 projected k=3 z1=336 scalar rows are empty",
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
