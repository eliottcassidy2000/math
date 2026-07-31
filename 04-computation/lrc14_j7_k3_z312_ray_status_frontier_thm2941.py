#!/usr/bin/env python3
"""Exact all-label ray/status frontier at the occupied ``k=3,z1=312`` slice.

The lossless scalar body atlas has exactly 80 bodies at first drift 312.  For
each body, the periodic residue-ray law computes the attained maximum in every
four-denominator multiset, with the projected high-label obligation retained.
Crude all-divisor fibre overload rejects 7,481 of the 18,249 attained states;
the common 16-cell Hunter table rejects another 10,698.  Every Hunter rejection
is replayed here with exact rational matrix arithmetic independent of the
status solver.  Replay hashes bind the deterministic infeasible instances,
not the noncanonical dual basis selected by the solver.

Seventy states remain, on four bodies.  For each one the script records its
attained maximizing literal packet and checks the older physical-denominator
bridge predicates directly.  All 70 pass both the expected-spike and
support-status relaxations, so the bridge gives no further reduction.  In this
residual family no denominator 3 or 4 occurs and denominator 2 occurs only
once.  Consequently the small-status allowance is exactly zero: an isolated
denominator-2 activity event has marginal 2/7 < 55/91, while an empty small
pattern contributes nothing.

This is deliberately an OPEN frontier package.  Denominator multisets erase
unit directions and literal phase locations, and the bridge also forgets those
coordinates.  The 70 maximizing packets are geometry probes, not an
exhaustive list of all scalar-admissible literal packets.  No closure and no
improvement of the inherited scalar height cap is claimed here.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter
from fractions import Fraction as F
from math import gcd, lcm
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
    "lrc14_j7_k3_z312_ray_status_frontier_thm2941.out"
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
EXPECTED_SEMANTIC_SHA256 = (
    "01738342225279bff3f3421aea6f79c598478565e4dd528c070999e29c165e30"
)

FIRST = 312
EXPECTED_BODY_COUNT = 80
EXPECTED_TOTALS = (18_249, 7_481, 10_698, 70)
EXPECTED_M_HISTOGRAM = (
    (2, 146),
    (3, 1_807),
    (4, 3_582),
    (5, 1_937),
    (6, 808),
    (7, 2_115),
    (9, 3),
    (10, 1),
    (11, 4),
    (12, 105),
    (13, 153),
    (14, 8),
    (15, 1),
    (21, 2),
    (24, 14),
    (28, 1),
    (40, 2),
    (42, 4),
    (70, 5),
)
THREE_SAFE_FLOOR = F(55, 91)

EXPECTED_RESIDUALS = {
    (1, 8, 10, 11, 12, 14): (
        (2, 5390, 10780, 18480),
        (392, 4620, 5390, 11760),
        (784, 4312, 4620, 5390),
        (784, 4620, 5390, 11760),
        (1617, 1848, 3920, 5390),
        (1617, 1848, 5390, 11760),
        (1617, 1848, 5390, 43120),
        (1617, 1848, 5390, 129360),
        (1848, 2352, 5390, 16170),
        (1848, 3920, 4312, 5390),
        (1848, 3920, 5390, 16170),
        (1848, 4312, 5390, 43120),
        (1848, 5390, 11760, 16170),
        (1960, 4620, 5390, 11760),
        (2352, 4312, 4620, 5390),
        (3920, 4312, 4620, 5390),
        (3920, 4620, 5390, 11760),
        (4312, 4620, 5390, 8624),
        (4312, 4620, 5390, 11760),
        (4312, 4620, 5390, 25872),
        (4312, 4620, 5390, 43120),
        (4312, 4620, 5390, 129360),
        (4620, 5390, 8624, 11760),
        (4620, 5390, 11760, 21560),
        (4620, 5390, 11760, 43120),
    ),
    (1, 8, 10, 12, 13, 14): (
        (490, 1040, 1911, 15288),
        (490, 1040, 7644, 15288),
        (490, 1911, 3120, 15288),
        (490, 1911, 3640, 10192),
        (490, 1911, 3640, 30576),
        (490, 1911, 3640, 50960),
        (490, 1911, 3640, 152880),
        (490, 1911, 7280, 15288),
        (490, 1911, 10192, 10920),
        (490, 1911, 10920, 30576),
        (490, 1911, 10920, 50960),
        (490, 1911, 10920, 152880),
        (490, 1911, 15288, 21840),
        (490, 3120, 7644, 15288),
        (490, 3640, 3822, 10192),
        (490, 3640, 7644, 10192),
        (490, 3640, 7644, 30576),
        (490, 3640, 7644, 50960),
        (490, 3640, 7644, 152880),
        (490, 3640, 9555, 10192),
        (490, 3640, 10192, 15288),
        (490, 3640, 10192, 19110),
        (490, 3640, 10192, 38220),
        (490, 3640, 10192, 76440),
        (490, 3640, 15288, 50960),
        (490, 3640, 19110, 30576),
        (490, 3640, 19110, 50960),
        (490, 3640, 19110, 152880),
        (490, 3640, 30576, 38220),
        (490, 3640, 38220, 50960),
        (490, 3640, 38220, 152880),
        (490, 7280, 7644, 15288),
        (490, 7644, 10192, 10920),
        (490, 7644, 10920, 30576),
        (490, 7644, 10920, 50960),
        (490, 7644, 10920, 152880),
        (490, 7644, 15288, 21840),
        (490, 10192, 10920, 15288),
        (490, 10920, 15288, 30576),
        (490, 10920, 15288, 50960),
        (490, 10920, 15288, 152880),
    ),
    (1, 8, 11, 12, 13, 14): (
        (1078, 3822, 8008, 112112),
        (1078, 3822, 8008, 336336),
    ),
    (2, 8, 10, 11, 12, 14): (
        (3920, 4312, 4620, 5390),
        (4312, 4620, 5390, 11760),
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


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
ray = load_module("z312_ray_status_engine", RAY_ENGINE_PATH)
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
        ("z312 body ledger changed", len(bodies), len(set(bodies))),
    )
    require(tuple(sorted(bodies)) == tuple(bodies), "z312 body order changed")
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
    return contradiction < 0


def bridge_record(stream, ds):
    """Apply the concrete expected-spike/support-status bridge to one state."""
    D = lcm(*ds)
    require(stream.L % D == 0 and D % 7 == 0, ("bad resolving denominator", ds))
    q = D // 7
    arcs = ray.fibre.projected_support_arcs(D, stream.ranges)
    support_count = sum(right - left for left, right in arcs)
    histogram = ray.fibre.residue_load_histogram(arcs, q)
    c = sum(q % d != 0 for d in ds)
    N_c = sum(count for load, count in histogram if load > c)
    expected_margin = 13 * (4 - c) * q - 55 * N_c
    expected_pass = N_c == 0 or expected_margin > 0

    small_pattern = (ds.count(2), ds.count(3), ds.count(4))
    require(small_pattern[1:] == (0, 0), ("unexpected d=3/4 residual", ds))
    require(small_pattern[0] <= 1, ("multiple d=2 residual", ds))
    if small_pattern[0]:
        # The only nonempty small event has marginal 2/7, already at most the
        # 55/91 safe floor, so no positive support load can be granted to it.
        require(F(2, 7) <= THREE_SAFE_FLOOR, "small-status direction changed")
    status_allowance = 0
    large_capacity = sum(
        (D // d) * ((d + 6) // 7) for d in ds if d > 4
    )
    support_slack = large_capacity + status_allowance - support_count
    support_pass = support_slack >= 0
    return (
        D,
        q,
        c,
        N_c,
        expected_margin,
        expected_pass,
        small_pattern,
        status_allowance,
        support_count,
        large_capacity,
        support_slack,
        support_pass,
    )


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

    verified_instance_digest = hashlib.sha256()
    verified_farkas_checks = 0
    all_negative = True
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
        certificate_negative = independent_farkas_check(
            q, marginals, capacities, histogram, certificate
        )
        verified_farkas_checks += 1
        all_negative = all_negative and certificate_negative
        # The exact certificate was replayed above.  HiGHS may choose a
        # different valid dual basis, so freeze the deterministic infeasible
        # instance and never the solver-selected certificate representative.
        verified_instance_digest.update(
            f"{body}|{ds}|{witness[:-1]}\n".encode()
        )

    status_instances = tuple(
        (ds, witness[:-1]) for ds, witness in sorted(status.items())
    )

    sign_totals = {
        sign: sum(
            count
            for (_d, candidate_sign), count in signs.items()
            if candidate_sign == sign
        )
        for sign in (-1, 0, 1)
    }
    require(sign_totals[-1] == sign_totals[1], "ray sign imbalance")

    packet_rows = []
    bridge_rows = []
    for ds in survivors:
        state = states[ds]
        labels = state["labels"]
        label_denominators = tuple(
            sorted(stream.L // gcd(stream.L, label) for label in labels)
        )
        require(label_denominators == ds, ("packet denominator mismatch", ds, labels))
        require(labels[0] == FIRST and len(set(labels)) == 4, ("packet identity changed", labels))
        high_labels = tuple(label for label in labels if label >= stream.high_floor)
        require(len(high_labels) == 1, ("high-slot geometry changed", ds, labels))
        packet_rows.append((ds, labels, state["total"], state["excess"], high_labels))
        bridge_rows.append((ds, bridge_record(stream, ds)))

    return (
        body,
        stream.L,
        stream.high_floor,
        stream.first_d,
        checks,
        tuple(sorted(sign_totals.items())),
        trials,
        tuple(sorted(states)),
        tuple(sorted(crude)),
        status_instances,
        tuple(packet_rows),
        tuple(bridge_rows),
        tuple(sorted(Counter(witness[1] for witness in status.values()).items())),
        verified_instance_digest.hexdigest(),
        verified_farkas_checks,
        all_negative,
    )


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=max(1, mp.cpu_count() // 2))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")
    bodies = bodies_from_atlas()
    pair_rows, control_marginals, control_caps, control_instances = ray.local.controls()
    if args.processes == 1:
        records = tuple(evaluate_body(body) for body in bodies)
    else:
        with mp.get_context("spawn").Pool(args.processes) as pool:
            records = tuple(sorted(pool.imap_unordered(evaluate_body, bodies, chunksize=1)))

    totals = tuple(
        sum(len(record[index]) for record in records) for index in (7, 8, 9, 10)
    )
    require(totals == EXPECTED_TOTALS, ("global counts changed", totals))
    require(
        sum(record[14] for record in records) == totals[2],
        "exact Farkas check count changed",
    )
    require(all(record[15] for record in records), "nonnegative Farkas replay")
    global_m = Counter()
    for record in records:
        global_m.update(dict(record[12]))
    require(tuple(sorted(global_m.items())) == EXPECTED_M_HISTOGRAM, "M histogram changed")

    observed_residuals = {
        record[0]: tuple(packet[0] for packet in record[10])
        for record in records
        if record[10]
    }
    require(observed_residuals == EXPECTED_RESIDUALS, ("residual ledger changed", observed_residuals))
    all_packets = tuple(
        (record[0], packet) for record in records for packet in record[10]
    )
    all_bridges = tuple(
        (record[0], ds, bridge) for record in records for ds, bridge in record[11]
    )
    require(len(all_packets) == len(all_bridges) == 70, "residual cardinality changed")
    require(all(bridge[5] and bridge[11] for _body, _ds, bridge in all_bridges), "bridge killed residual")
    require(sum(bridge[6][0] for _body, _ds, bridge in all_bridges) == 1, "d=2 ledger changed")

    semantic_payload = (
        FIRST,
        records,
        tuple(sorted(global_m.items())),
        pair_rows,
        control_marginals,
        control_caps,
        control_instances,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    c_histogram = Counter(bridge[2] for _body, _ds, bridge in all_bridges)
    D_histogram = Counter(bridge[0] for _body, _ds, bridge in all_bridges)
    positive_expected_margins = tuple(
        bridge[4] for _body, _ds, bridge in all_bridges if bridge[3] > 0
    )
    support_slacks = tuple(bridge[10] for _body, _ds, bridge in all_bridges)
    lines = [
        "LRC14 projected k=3 z1=312 exact all-label ray/status OPEN frontier",
        "status=FINITE-EXACT through ray/status/bridge;OPEN at 70 literal-phase residuals",
        f"ray_engine_sha256={file_sha256(RAY_ENGINE_PATH)}",
        f"body_atlas_source_sha256={file_sha256(BODY_ATLAS_PATH)}",
        f"body_atlas_output_sha256={file_sha256(BODY_ATLAS_OUTPUT_PATH)}",
        "scope=80 lossless scalar-atlas bodies;no finite label horizon",
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
        (
            body, L, high_floor, first_d, checks, signs, trials, states, crude,
            status, packets, _bridges, mhist, cdigest, verified_checks,
            all_negative,
        ) = record
        partition_digest = hashlib.sha256(
            repr((states, crude, status, tuple(packet[0] for packet in packets))).encode()
        ).hexdigest()
        lines.append(
            f"E={body};L={L};high_floor={high_floor};d1={first_d};"
            f"ray_checks={checks};ray_signs={dict(signs)};denominator_trials={trials};"
            f"states={len(states)};crude_kills={len(crude)};"
            f"status_kills={len(status)};survivors={len(packets)};"
            f"status_M={dict(mhist)};partition_sha256={partition_digest};"
            f"verified_farkas_instance_sha256={cdigest};"
            f"verified_farkas_checks={verified_checks};"
            f"all_farkas_contradictions_negative={int(all_negative)}"
        )

    lines.extend(
        (
            f"survivor_body_counts={{{', '.join(f'{body}: {len(rows)}' for body, rows in EXPECTED_RESIDUALS.items())}}}",
            "packet_geometry=attained maximizing packet per residual;first label 312;four distinct labels;exactly one label at/above body high wall",
            "packet_high_slot_histogram={1: 70}",
            "bridge_scope=necessary expected-spike AND one-threshold support-status relaxation;not sufficient for literal realization",
            "expected_spike_predicate=N_c==0 OR 55*N_c<13*(4-c)*q",
            "support_status_predicate=large_capacity+small_allowance>=projected_support;small_allowance=0 on all residuals",
            f"bridge_c_histogram={dict(sorted(c_histogram.items()))}",
            f"bridge_resolving_D_histogram={dict(sorted(D_histogram.items()))}",
            f"bridge_positive_N_min_expected_margin={min(positive_expected_margins, default=None)}",
            f"bridge_min_support_slack={min(support_slacks)}",
            "bridge_admissible=expected_spike:70/70;support_status:70/70;joint:70/70;reduction:0",
        )
    )
    for body, packet in all_packets:
        ds, labels, total, excess, high_labels = packet
        bridge = next(
            row for candidate_body, candidate_ds, row in all_bridges
            if candidate_body == body and candidate_ds == ds
        )
        lines.append(
            f"residual=E:{body};ds:{ds};max_labels:{labels};high_labels:{high_labels};"
            f"total:{ftext(total)};excess:{ftext(excess)};D:{bridge[0]};q:{bridge[1]};"
            f"c:{bridge[2]};N_c:{bridge[3]};expected_margin:{bridge[4]};"
            f"small_pattern:{bridge[6]};support:{bridge[8]};"
            f"large_capacity:{bridge[9]};support_slack:{bridge[10]}"
        )
    lines.extend(
        (
            "independent_exact_farkas_checks=10698/10698:PASS",
            f"semantic_sha256={semantic_hash}",
            "consequence=OPEN residual 70;no z1 cap improvement claimed",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
