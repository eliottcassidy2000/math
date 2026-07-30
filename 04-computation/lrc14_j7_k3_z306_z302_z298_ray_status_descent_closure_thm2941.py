#!/usr/bin/env python3
"""Exact all-label closure of projected k=3 slices z1=306, 302, and 298.

The guarded z312 low-pair torsion package leaves the exact body-atlas frontier
at 306.  The atlas has two bodies there and, after the empty 303..305 gap,
the nine bodies at 302, and the twelve bodies at 298.  Exact periodic residue
rays enumerate every attained
four-denominator state with the inherited projected high-label obligation.
All 137 states are rejected: seventeen by a crude fibre overload and 120 by
the common 16-cell Hunter status table.  Every status rejection is replayed by a
second exact rational Farkas checker independent of the status solver.

The common-table engine's exhaustive small controls include both a feasible
positive instance and an incompatible hostile instance.  No finite label
horizon remains.  Since the next occupied atlas height is 297, closing all
three slices lowers the projected k=3 first-drift cap to 297.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_ENGINE_PATH = (
    ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
)
BODY_ATLAS_PATH = (
    ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
)
BODY_ATLAS_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
Z312_TERMINAL_PATH = (
    ROOT / "04-computation/"
    "lrc14_j7_k3_z312_finite_low_pair_torsion_closure_thm2941.py"
)
Z312_TERMINAL_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z312_finite_low_pair_torsion_closure_thm2941.out"
)
Z312_REFEREE_PATH = (
    ROOT / "04-computation/"
    "lrc14_j7_k3_z312_finite_low_pair_torsion_independent_thm2941.py"
)
Z312_REFEREE_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z312_finite_low_pair_torsion_independent_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z306_z302_z298_ray_status_descent_closure_thm2941.out"
)

EXPECTED_RAY_ENGINE_SHA256 = (
    "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2"
)
EXPECTED_BODY_ATLAS_SHA256 = (
    "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded"
)
EXPECTED_BODY_ATLAS_OUTPUT_SHA256 = (
    "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
)
EXPECTED_Z312_TERMINAL_SHA256 = (
    "6b7d70ae105530ce01ffc8bd0d7770f4fe6a25c4cb76665c8095e2567f5ddfc0"
)
EXPECTED_Z312_TERMINAL_OUTPUT_SHA256 = (
    "922d33a8a5bee98c61c5eba5388dc90b51a73317000934f06e15b9c65a81a085"
)
EXPECTED_Z312_REFEREE_SHA256 = (
    "24006ebec0fc278be3571a349e2ad6f86aba28d22d6ee6159bcc1d3f15449d2f"
)
EXPECTED_Z312_REFEREE_OUTPUT_SHA256 = (
    "33b8a9abc28abb3958addb3413ce5882c64a8e603f6f23280bed9c3ba021e193"
)
EXPECTED_SEMANTIC_SHA256 = (
    "cbfa2716210e0e5c09afd354ebafed4001e950fd6af87bb1dcf3f41a11a3e492"
)

EXPECTED_COUNTS = {
    (306, (1, 4, 8, 10, 12, 14)): (1, 0, 1),
    (306, (2, 4, 8, 10, 12, 14)): (8, 0, 8),
    (302, (1, 3, 4, 8, 10, 12)): (1, 1, 0),
    (302, (1, 4, 6, 8, 12, 14)): (1, 1, 0),
    (302, (1, 4, 8, 10, 12, 14)): (4, 0, 4),
    (302, (2, 4, 6, 8, 10, 14)): (1, 0, 1),
    (302, (2, 4, 6, 8, 12, 14)): (5, 0, 5),
    (302, (2, 4, 6, 9, 12, 14)): (2, 0, 2),
    (302, (2, 4, 6, 10, 12, 14)): (3, 0, 3),
    (302, (2, 4, 8, 10, 12, 14)): (8, 0, 8),
    (302, (2, 6, 8, 10, 12, 14)): (7, 0, 7),
    (298, (1, 2, 6, 8, 12, 14)): (1, 0, 1),
    (298, (1, 2, 6, 10, 12, 14)): (1, 0, 1),
    (298, (1, 3, 4, 8, 12, 14)): (1, 0, 1),
    (298, (1, 4, 6, 8, 12, 14)): (6, 2, 4),
    (298, (1, 4, 6, 9, 12, 14)): (2, 0, 2),
    (298, (1, 4, 6, 10, 12, 14)): (7, 1, 6),
    (298, (1, 4, 8, 10, 12, 14)): (44, 10, 34),
    (298, (1, 6, 8, 10, 12, 14)): (6, 2, 4),
    (298, (2, 4, 6, 8, 12, 14)): (5, 0, 5),
    (298, (2, 4, 6, 10, 12, 14)): (3, 0, 3),
    (298, (2, 4, 8, 10, 12, 14)): (12, 0, 12),
    (298, (2, 6, 8, 10, 12, 14)): (8, 0, 8),
}
EXPECTED_LEVEL_TOTALS = {
    306: (9, 0, 9, 0),
    302: (32, 2, 30, 0),
    298: (96, 15, 81, 0),
}
EXPECTED_TOTALS = (137, 17, 120, 0)
EXPECTED_M_HISTOGRAM = ((3, 7), (4, 9), (5, 23), (6, 29), (7, 52))


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


for path, expected in (
    (RAY_ENGINE_PATH, EXPECTED_RAY_ENGINE_SHA256),
    (BODY_ATLAS_PATH, EXPECTED_BODY_ATLAS_SHA256),
    (BODY_ATLAS_OUTPUT_PATH, EXPECTED_BODY_ATLAS_OUTPUT_SHA256),
    (Z312_TERMINAL_PATH, EXPECTED_Z312_TERMINAL_SHA256),
    (Z312_TERMINAL_OUTPUT_PATH, EXPECTED_Z312_TERMINAL_OUTPUT_SHA256),
    (Z312_REFEREE_PATH, EXPECTED_Z312_REFEREE_SHA256),
    (Z312_REFEREE_OUTPUT_PATH, EXPECTED_Z312_REFEREE_OUTPUT_SHA256),
):
    require(file_sha256(path) == expected, ("dependency changed", path))

ray = load_module("z306_z302_z298_ray_engine", RAY_ENGINE_PATH)


def atlas_rows():
    by_level = {306: [], 302: [], 298: []}
    counts = Counter()
    pattern = re.compile(r"^row=E=([0-9,]+);.*;z1=([0-9]+);")
    for line in BODY_ATLAS_OUTPUT_PATH.read_text().splitlines():
        match = pattern.match(line)
        if not match:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        first = int(match.group(2))
        counts[first] += 1
        if first in by_level:
            by_level[first].append(body)
    require(
        tuple(by_level[306])
        == ((1, 4, 8, 10, 12, 14), (2, 4, 8, 10, 12, 14)),
        ("z306 body order changed", by_level[306]),
    )
    require(
        tuple(by_level[302])
        == tuple(body for first, body in EXPECTED_COUNTS if first == 302),
        ("z302 body order changed", by_level[302]),
    )
    require(
        tuple(by_level[298])
        == tuple(body for first, body in EXPECTED_COUNTS if first == 298),
        ("z298 body order changed", by_level[298]),
    )
    require(
        counts[306] == 2
        and all(counts[value] == 0 for value in range(303, 306))
        and counts[302] == 9
        and all(counts[value] == 0 for value in range(299, 302))
        and counts[298] == 12,
        ("descent atlas gap changed", counts),
    )
    require(counts[297] == 7, ("next atlas frontier changed", counts[297]))
    return tuple(
        (first, body) for first in (306, 302, 298) for body in by_level[first]
    ), tuple(sorted(counts.items()))


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
        value -= sum(alpha[row] * tail_rows[row][pattern] for row in range(len(alpha)))
        slacks.append(value)
    contradiction = sum(z[row] * equality_rhs[row] for row in range(5))
    contradiction -= sum(alpha[row] * tail_rhs[row] for row in range(len(alpha)))
    require(all(value >= 0 for value in slacks), "negative Farkas column")
    require(contradiction < 0, "nonnegative Farkas contradiction")
    return tuple(slacks), contradiction


def evaluate_body(task):
    first, body = task
    ray.FIRST = first
    stream = ray.Stream(body)
    require(stream.first == first, ("first drift changed", task, stream.first))
    trials, states, checks, signs = ray.ray_quotient_states(stream)
    crude, status, survivors = ray.common_status_screen(stream, states)
    require(
        (len(states), len(crude), len(status)) == EXPECTED_COUNTS[task],
        ("body counts changed", task, len(states), len(crude), len(status)),
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

    verified_instances = []
    verified_instance_digest = hashlib.sha256()
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
        _slacks, _contradiction = independent_farkas_check(
            q, marginals, capacities, histogram, certificate
        )
        instance = (ds, witness[:-1])
        verified_instances.append(instance)
        if representative is None:
            representative = (task, instance)
        verified_instance_digest.update(f"{task}|{instance}\n".encode())

    sign_totals = {
        sign: sum(
            count
            for (_denominator, candidate_sign), count in signs.items()
            if candidate_sign == sign
        )
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], "ray sign imbalance")
    partition = (
        tuple(sorted(states.items())),
        tuple(sorted(crude.items())),
        tuple(verified_instances),
        tuple(survivors),
    )
    return (
        first,
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
        tuple(sorted(Counter(witness[1] for witness in status.values()).items())),
        max(states.items(), key=lambda item: (item[1]["excess"], item[0]), default=None),
        hashlib.sha256(repr(partition).encode()).hexdigest(),
        verified_instance_digest.hexdigest(),
        len(verified_instances),
        representative,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")
    tasks, atlas_counts = atlas_rows()
    require(tasks == tuple(EXPECTED_COUNTS), "global task order changed")

    pair_rows, control_marginals, control_caps, control_instances = ray.local.controls()
    if args.processes == 1:
        records = tuple(evaluate_body(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            records = tuple(pool.map(evaluate_body, tasks))
    totals = tuple(sum(len(record[11][index]) for record in records) for index in range(4))
    require(totals == EXPECTED_TOTALS, ("global counts changed", totals))
    level_totals = {
        first: tuple(
            sum(len(record[11][index]) for record in records if record[0] == first)
            for index in range(4)
        )
        for first in (306, 302, 298)
    }
    require(level_totals == EXPECTED_LEVEL_TOTALS, ("level totals changed", level_totals))
    global_m = Counter()
    for record in records:
        global_m.update(dict(record[12]))
    m_histogram = tuple(sorted(global_m.items()))
    if EXPECTED_M_HISTOGRAM is not None:
        require(m_histogram == EXPECTED_M_HISTOGRAM, "M histogram changed")
    representative = next(record[17] for record in records if record[17] is not None)

    semantic_payload = (
        tasks,
        records,
        totals,
        tuple(sorted(level_totals.items())),
        m_histogram,
        atlas_counts,
        pair_rows,
        control_marginals,
        control_caps,
        control_instances,
        representative,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=306, z1=302, and z1=298 all-label ray/status descent closure",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        f"body_atlas_source_sha256={file_sha256(BODY_ATLAS_PATH)}",
        f"body_atlas_output_sha256={file_sha256(BODY_ATLAS_OUTPUT_PATH)}",
        f"z312_terminal_source_sha256={file_sha256(Z312_TERMINAL_PATH)}",
        f"z312_terminal_output_sha256={file_sha256(Z312_TERMINAL_OUTPUT_PATH)}",
        f"z312_referee_source_sha256={file_sha256(Z312_REFEREE_PATH)}",
        f"z312_referee_output_sha256={file_sha256(Z312_REFEREE_OUTPUT_PATH)}",
        "scope=2 exact z306 bodies plus 9 exact z302 bodies plus 12 exact z298 bodies;inherited projected high-label obligation (first itself may clear the wall);no finite label horizon",
        "ray_law=(z+L)delta(z+L)=zdelta(z);A(L-b)=-A(b);all denominator maxima attained",
        f"pair_overlap_exhaustive_controls={pair_rows}",
        "common_table_controls=coverable positive instance accepted;incompatible hostile instance rejected",
        f"level_totals={level_totals}",
        f"totals=states:{totals[0]};crude_kills:{totals[1]};status_kills:{totals[2]};survivors:{totals[3]}",
        f"status_M_histogram={dict(m_histogram)}",
    ]
    for record in records:
        (
            first, body, h, components, L, high_floor, first_d, checks, signs,
            denominator_classes, trials, partition, status_histogram, maximum,
            partition_digest, verified_instance_digest, verified_instance_count, _representative,
        ) = record
        states, crude, status, survivors = partition
        lines.append(
            f"z1={first};E={body};h={h};r={components};L={L};high={high_floor};d1={first_d};"
            f"ray_checks={checks};ray_signs={dict(signs)};denominator_classes={denominator_classes};"
            f"denominator_trials={trials};states={len(states)};crude_kills={len(crude)};"
            f"status_kills={len(status)};survivors={len(survivors)};status_M={dict(status_histogram)};"
            f"max_state={maximum};partition_sha256={partition_digest};"
            f"verified_farkas_instance_sha256={verified_instance_digest};"
            f"verified_exact_farkas_instances={verified_instance_count}"
        )
    lines.extend(
        (
            f"representative_verified_farkas_instance={representative}",
            "independent_exact_farkas_checks=120/120:PASS",
            "atlas_gaps=303..305 empty;299..301 empty;next occupied height297 with7 bodies",
            "conclusion=all projected k3 rows at z1=306, z1=302, and z1=298 are empty;projected k3 cap<=297;next exact frontier=297",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
