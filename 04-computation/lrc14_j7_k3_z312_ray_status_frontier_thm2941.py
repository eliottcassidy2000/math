#!/usr/bin/env python3
"""Exact all-label ray/status frontier for projected ``k=3,z1=312``.

The lossless scalar body atlas contains exactly 80 bodies at first drift 312.
For each body this verifier rebuilds the attained, horizon-free denominator
quotient, applies the inherited crude all-divisor screen, and independently
checks every common-status Farkas certificate.  It leaves 70 denominator
states on four bodies for the literal projected-geometry addendum.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAY_PATH = ROOT / "04-computation" / "lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
ATLAS_PATH = ROOT / "04-computation" / "lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
ATLAS_OUTPUT_PATH = ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
DEFAULT_OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_k3_z312_ray_status_frontier_thm2941.out"
EXPECTED_RAY_SHA256 = "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2"
EXPECTED_ATLAS_SHA256 = "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded"
EXPECTED_ATLAS_OUTPUT_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
EXPECTED_SEMANTIC_SHA256 = "89e226970a503b8ec331ec30236fa39a7e80d1ec2f9b238dc61520daedc32902"
FIRST = 312
EXPECTED_BODY_COUNT = 80
EXPECTED_TOTALS = (18249, 7481, 10698, 70)
EXPECTED_RESIDUAL_COUNTS = {
    (1, 8, 10, 11, 12, 14): 25,
    (1, 8, 10, 12, 13, 14): 41,
    (1, 8, 11, 12, 13, 14): 2,
    (2, 8, 10, 11, 12, 14): 2,
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(RAY_PATH) == EXPECTED_RAY_SHA256, "ray engine changed")
require(file_sha256(ATLAS_PATH) == EXPECTED_ATLAS_SHA256, "body atlas changed")
require(file_sha256(ATLAS_OUTPUT_PATH) == EXPECTED_ATLAS_OUTPUT_SHA256, "body atlas transcript changed")
ray = load("z312_frontier_ray", RAY_PATH)
ray.FIRST = FIRST


def bodies_from_atlas():
    bodies = []
    for line in ATLAS_OUTPUT_PATH.read_text().splitlines():
        if not line.startswith("row=E=") or ";z1=312;" not in line:
            continue
        body_text = line[6:].split(";", 1)[0]
        bodies.append(tuple(map(int, body_text.split(","))))
    require(
        len(bodies) == len(set(bodies)) == EXPECTED_BODY_COUNT,
        ("z312 body slice changed", len(bodies), len(set(bodies))),
    )
    return tuple(bodies)


def independent_farkas_check(q, marginals, capacities, histogram, certificate):
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
        slack -= sum(alpha[row] * tail_rows[row][pattern] for row in range(len(alpha)))
        require(slack >= 0, ("negative Farkas column", pattern, slack))
    contradiction = sum(z[row] * equality_rhs[row] for row in range(5))
    contradiction -= sum(alpha[row] * tail_rhs[row] for row in range(len(alpha)))
    require(contradiction < 0, ("nonnegative Farkas contradiction", contradiction))
    return contradiction


def digest(value):
    return hashlib.sha256(repr(value).encode()).hexdigest()


def evaluate_body(body):
    stream = ray.Stream(body)
    trials, states, checks, signs = ray.ray_quotient_states(stream)
    crude, status, survivors = ray.common_status_screen(stream, states)
    contradictions = []
    for ds, witness in sorted(status.items()):
        q, _M, marginals, _cap_set, histogram, certificate = witness
        D = __import__("math").lcm(*ds)
        rebuilt_marginals, capacities = ray.local.hunter_status_data(D, ds, q)
        require(rebuilt_marginals == marginals, "status marginals changed")
        contradictions.append(
            independent_farkas_check(q, marginals, capacities, histogram, certificate)
        )
    sign_totals = tuple(
        (sign, sum(count for (_d, candidate), count in signs.items() if candidate == sign))
        for sign in (-1, 0, 1)
    )
    require(sign_totals[0][1] == sign_totals[2][1], "ray sign imbalance")
    residual_states = tuple((ds, states[ds]) for ds in survivors)
    record = (
        body,
        stream.L,
        stream.high_floor,
        checks,
        sign_totals,
        trials,
        len(states),
        len(crude),
        len(status),
        residual_states,
        digest(tuple(sorted(states.items()))),
        digest(tuple(sorted(crude.items()))),
        digest(tuple(sorted(status.items()))),
        min(contradictions, default=None),
    )
    print(
        f"DONE E={body} states={len(states)} crude={len(crude)} "
        f"status={len(status)} survivors={len(survivors)}",
        flush=True,
    )
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=6)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    bodies = bodies_from_atlas()
    controls = ray.local.controls()
    with mp.Pool(args.processes) as pool:
        records = tuple(sorted(pool.imap_unordered(evaluate_body, bodies, chunksize=1)))
    totals = tuple(sum(record[index] for record in records) for index in (6, 7, 8))
    residual_total = sum(len(record[9]) for record in records)
    require((*totals, residual_total) == EXPECTED_TOTALS, ("frontier totals changed", totals, residual_total))
    residual_counts = {
        record[0]: len(record[9]) for record in records if record[9]
    }
    require(residual_counts == EXPECTED_RESIDUAL_COUNTS, ("residual body ledger changed", residual_counts))
    farkas_count = sum(record[8] for record in records)
    require(farkas_count == EXPECTED_TOTALS[2], "Farkas audit count changed")
    semantic_payload = (bodies, controls, records, EXPECTED_TOTALS)
    semantic_hash = digest(semantic_payload)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=3 z1=312 exact all-label ray/status frontier",
        f"ray_engine_sha256={file_sha256(RAY_PATH)}",
        f"body_atlas_source_sha256={file_sha256(ATLAS_PATH)}",
        f"body_atlas_output_sha256={file_sha256(ATLAS_OUTPUT_PATH)}",
        "scope=80 lossless scalar-atlas bodies;no finite label horizon",
        "ray_law=(z+L)delta(z+L)=zdelta(z);A(L-b)=-A(b);all denominator maxima attained",
        f"pair_overlap_exhaustive_controls={controls[0]}",
        "totals=states:18249;crude_kills:7481;status_kills:10698;survivors:70",
        "residual_bodies=25+41+2+2 states on four bodies",
    ]
    for record in records:
        body, L, high, checks, signs, trials, states_n, crude_n, status_n, residual_states, state_hash, crude_hash, status_hash, min_contra = record
        lines.append(
            f"E={body};L={L};high={high};ray_checks={checks};ray_signs={dict(signs)};"
            f"denominator_trials={trials};states={states_n};crude_kills={crude_n};"
            f"status_kills={status_n};survivors={len(residual_states)};"
            f"state_sha256={state_hash};crude_sha256={crude_hash};"
            f"status_sha256={status_hash};min_exact_farkas_contradiction={min_contra}"
        )
        if residual_states:
            lines.append(
                f"  residual_denominators={tuple(ds for ds, _state in residual_states)}"
            )
            lines.append(f"  residual_states={residual_states}")
    lines.extend(
        (
            f"independent_exact_farkas_checks={farkas_count}/{farkas_count}:PASS",
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
