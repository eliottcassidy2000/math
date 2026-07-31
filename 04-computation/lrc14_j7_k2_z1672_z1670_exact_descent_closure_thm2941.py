#!/usr/bin/env python3
"""Exact projected k=2 descent from the THM-2972 cap through z1=1670.

An all-3003-body scalar replay on 1668..1679 proves that the only occupied
heights are 1668, 1670, and 1672.  The unique z1=1672 row and both z1=1670
rows then close under the unrestricted ray/status/projected engine.  Status
dual bases are checked over Q but removed from canonical replay digests.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import subprocess
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, lcm
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
THM2972_COMMIT = "7d96fd42099a"
THM2972_SOURCE_REL = "04-computation/lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.py"
THM2972_OUTPUT_REL = "05-knowledge/results/lrc14_j7_k2_all_four_high_punctured_packet_torsion_closure_thm2972.out"
ATLAS_ENGINE_SOURCE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1810_1835_thm2941.py"
STATUS_ENGINE_SOURCE = ROOT / "04-computation/lrc14_j7_k2_z1724_exact_descent_closure_thm2941.py"
STATUS_ENGINE_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_z1724_exact_descent_closure_thm2941.out"
OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k2_z1672_z1670_exact_descent_closure_thm2941.out"
    if HERE.parent.name == "04-computation"
    else HERE.with_suffix(".out")
)

EXPECTED_THM2972_SOURCE_SHA256 = "b92f7c6cf295f29a3ec2c64730f6bd4486073b9197f322202df98df497369d97"
EXPECTED_THM2972_OUTPUT_SHA256 = "f6b34a08c30964f2c6ebf71b5e4f580fb6c75a4cf62762cb530ef892f9d4af25"
EXPECTED_ATLAS_ENGINE_SHA256 = "a09b13e994ad6ab35e3324bd336e773b435b07859e6a3c924b84ec77f3e2aced"
EXPECTED_STATUS_ENGINE_SHA256 = "11d4c0836cc23304a57cd72fc5d8eb4a14479cb9ca00395c64d45f184c0ba313"
EXPECTED_STATUS_OUTPUT_SHA256 = "ef38b8230cdba56391877d175c80b8a8e8d804d299497270bc7fdb59b0598245"
EXPECTED_ATLAS_PROFILE_SHA256 = "10b07e2eb300d3ddbd1c4b99340a6a52ee6ef3a54d27cbd30c1a2138196f602f"
EXPECTED_ATLAS_SURVIVOR_SHA256 = "bb45977b91e5b064be778e65bd5cb45f020607716ff72902674721a675a32f64"
EXPECTED_CLOSURE_PROFILE_SHA256 = "c010bfb55df4401529dcdd48370c003834a438732e7dc2f874db22ad9631aad4"
EXPECTED_SEMANTIC_SHA256 = "6f3d3bd793c6ef513390cf15601ccbf3ee370cdfb385d543aed6df51b7c6493a"

START = 1668
END = 1679
EXPECTED_CANDIDATE_ROWS = 36_036
EXPECTED_HEIGHTS = ((1668, 9), (1670, 2), (1672, 1))
EXPECTED_ATLAS_KEYS = tuple(
    sorted(
        (
            (1668, (1, 2, 10, 11, 12, 14)),
            (1668, (1, 4, 8, 10, 12, 14)),
            (1668, (1, 6, 8, 10, 12, 14)),
            (1668, (1, 6, 9, 10, 12, 14)),
            (1668, (1, 6, 10, 11, 12, 14)),
            (1668, (1, 8, 10, 11, 12, 14)),
            (1668, (2, 4, 8, 10, 12, 14)),
            (1668, (2, 6, 8, 10, 12, 14)),
            (1668, (2, 8, 10, 11, 12, 14)),
            (1670, (1, 4, 8, 10, 12, 14)),
            (1670, (1, 6, 8, 10, 12, 14)),
            (1672, (1, 4, 8, 10, 12, 14)),
        )
    )
)
EXPECTED_HIGH_KEYS = tuple(
    sorted(
        (
            (1668, (1, 2, 10, 11, 12, 14)),
            (1668, (1, 6, 10, 11, 12, 14)),
            (1668, (1, 8, 10, 11, 12, 14)),
            (1668, (2, 8, 10, 11, 12, 14)),
        )
    )
)
CLOSURE_CASES = (
    (1670, (1, 4, 8, 10, 12, 14)),
    (1670, (1, 6, 8, 10, 12, 14)),
    (1672, (1, 4, 8, 10, 12, 14)),
)
EXPECTED_COUNTS = {
    CLOSURE_CASES[0]: (170, 24, 146, 0, 0),
    CLOSURE_CASES[1]: (1, 0, 0, 1, 1),
    CLOSURE_CASES[2]: (88, 8, 80, 0, 0),
}
EMPTY_DIGEST = sha256(b"()").hexdigest()
EXPECTED_EXACT_AUDIT = {
    CLOSURE_CASES[0]: (
        (
            "be0032ac5897299cb2995847c6d644a5efa28f57756a7ff5c361f7bb0ef97bff",
            "067e2a739064ede6614151128cd9bb8a37c6cd0f67fecdcc15ee9ce031cd2f8f",
            "43f3c899f26879481586355b18a9d924db38c2ad8ccf22bac4e561aaaeb6ffa7",
            EMPTY_DIGEST,
        ),
        None, 0, None, EMPTY_DIGEST,
    ),
    CLOSURE_CASES[1]: (
        (
            "8e2845a17f9c354f392028998bfad183a155c4567f0382fd012dd434729e048c",
            EMPTY_DIGEST,
            EMPTY_DIGEST,
            "4015991e031b38f4086e451350af5e8059e22c05f140fff2cffa181cc8eb58be",
        ),
        F(8929, 15197), 1, F(1),
        "2891878ad27c5aa660c82c54f2513b5425ab201b1954c3d458f824a507345277",
    ),
    CLOSURE_CASES[2]: (
        (
            "77d72cf86fce6add9fb33c01e0ae7e35f0424619c49b5271199755744035be49",
            "0860083380b840b591ec034d1166cc464144657c8f7dbc28ff92f924f7f29efe",
            "bac229cc1f2e7e5e7365efd5547a2fee9fc0edeb7eeb95a38ae16db270a2c0f2",
            EMPTY_DIGEST,
        ),
        None, 0, None, EMPTY_DIGEST,
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def dependency_bytes(relative):
    path = ROOT / relative
    if path.exists():
        return path.read_bytes()
    return subprocess.check_output(
        ("git", "show", f"{THM2972_COMMIT}:{relative}"), cwd=ROOT
    )


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(sha256(dependency_bytes(THM2972_SOURCE_REL)).hexdigest() == EXPECTED_THM2972_SOURCE_SHA256, "THM-2972 source changed")
require(sha256(dependency_bytes(THM2972_OUTPUT_REL)).hexdigest() == EXPECTED_THM2972_OUTPUT_SHA256, "THM-2972 output changed")
require(file_sha256(ATLAS_ENGINE_SOURCE) == EXPECTED_ATLAS_ENGINE_SHA256, "atlas engine changed")
require(file_sha256(STATUS_ENGINE_SOURCE) == EXPECTED_STATUS_ENGINE_SHA256, "status engine changed")
require(file_sha256(STATUS_ENGINE_OUTPUT) == EXPECTED_STATUS_OUTPUT_SHA256, "status output changed")

ASPEC = spec_from_file_location("k2_z1670_atlas_engine", ATLAS_ENGINE_SOURCE)
require(ASPEC is not None and ASPEC.loader is not None, "cannot load atlas engine")
ATLAS = module_from_spec(ASPEC)
ASPEC.loader.exec_module(ATLAS)
SSPEC = spec_from_file_location("k2_z1670_status_engine", STATUS_ENGINE_SOURCE)
require(SSPEC is not None and SSPEC.loader is not None, "cannot load status engine")
STATUS_PARENT = module_from_spec(SSPEC)
SSPEC.loader.exec_module(STATUS_PARENT)
ENGINE = STATUS_PARENT.ENGINE


def configure_atlas():
    ATLAS.B.START = START
    ATLAS.B.END = END
    ATLAS.B.LABELS = np.arange(START, ATLAS.B.HORIZON + 1, dtype=np.int64)


def atlas_profile(body):
    return ATLAS.B.profile(body)


def exact_status_audit(case, L, status_kills, stage_digest):
    canonical = tuple(
        (ds, upper, labels, witness[:-1])
        for ds, upper, labels, witness in status_kills
    )
    require(sha256(repr(canonical).encode()).hexdigest() == stage_digest, (case, "noncanonical status digest"))
    require(
        not status_kills or sha256(repr(tuple(status_kills)).encode()).hexdigest() != stage_digest,
        (case, "raw LP certificate entered status digest"),
    )
    actual_L, ranges = ENGINE.U.support.safe_cell_ranges(case[1])
    require(actual_L == L, (case, "safe-cell ruler changed"))
    arcs_cache = {}
    for ds, _upper, _labels, witness in status_kills:
        q, M, marginals, capacity_values, histogram, certificate = witness
        D = lcm(*ds)
        require(D == q * M, (case, ds, "status factorization"))
        if D not in arcs_cache:
            arcs_cache[D] = ENGINE.U.fibre.projected_support_arcs(D, ranges)
        actual_marginals, capacities = ENGINE.hunter_status_data5(D, ds, q)
        actual_histogram = ENGINE.U.fibre.residue_load_histogram(arcs_cache[D], q)
        require(actual_marginals == marginals, (case, ds, "status marginals"))
        require(tuple(sorted(set(capacities))) == capacity_values, (case, ds, "status capacities"))
        require(actual_histogram == histogram, (case, ds, "status histogram"))
        thresholds, alpha, z = certificate
        tail_rows = []
        tail_rhs = []
        for threshold in thresholds:
            demand = sum(count for load, count in histogram if load >= threshold)
            good = tuple(int(capacity >= threshold) for capacity in capacities)
            require(not all(good), (case, ds, threshold, "inactive dual row"))
            tail_rows.append(good)
            tail_rhs.append(demand)
        equality_rows = [
            (1,) * 32,
            *[tuple((pattern >> index) & 1 for pattern in range(32)) for index in range(5)],
        ]
        equality_rhs = (q, *marginals)
        require(len(alpha) == len(tail_rows) and len(z) == 6, (case, ds, "dual dimensions"))
        slacks = tuple(
            sum(z[row] * equality_rows[row][pattern] for row in range(6))
            - sum(alpha[row] * tail_rows[row][pattern] for row in range(len(tail_rows)))
            for pattern in range(32)
        )
        contradiction = (
            sum(z[row] * equality_rhs[row] for row in range(6))
            - sum(alpha[row] * tail_rhs[row] for row in range(len(tail_rows)))
        )
        require(all(value >= 0 for value in alpha), (case, ds, "negative multiplier"))
        require(all(value >= 0 for value in slacks), (case, ds, "negative exact slack"))
        require(contradiction < 0, (case, ds, "nonnegative contradiction"))
    return len(status_kills)


def closure_profile(case):
    first, body = case
    carrier = ENGINE.U.suffix.A.carrier_for(body)
    require(ENGINE.P.A.carrier_for(body) == carrier, (case, "carrier disagreement"))
    h = F(sum(right - left for left, right in carrier), ENGINE.U.suffix.A.RULER)
    lower = h * ENGINE.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    result = ENGINE.ray_and_status(first, body, carrier, h, lower, L)
    certificate_count = exact_status_audit(case, L, result[8], result[10][2])
    if result[9]:
        projected = ENGINE.projected_packets(
            first, body, carrier, h, lower, L, result[0], result[4], result[9]
        )
    else:
        projected = (0, 0, 0, None, 0, None, None, EMPTY_DIGEST, ())
    counts = (len(result[6]), len(result[7]), len(result[8]), len(result[9]), projected[1])
    require(counts == EXPECTED_COUNTS[case], (case, "counts changed", counts))
    require(certificate_count == counts[2], (case, "certificate count"))
    require(projected[1] == projected[2], (case, "projected survivor"))
    audit = (result[10], projected[3], projected[4], projected[6], projected[7])
    require(audit == EXPECTED_EXACT_AUDIT[case], (case, "audit changed", audit))
    return (
        first, body, h, len(carrier), L, lower, result[4], result[5], result[1],
        result[2], result[3], counts, result[10], *projected[:-1], certificate_count,
    )


def render(atlas_profiles, closure_profiles):
    require(len(atlas_profiles) == comb(14, 6) == 3003, "atlas body universe changed")
    candidate_rows = sum(row[7] for row in atlas_profiles)
    survivors = tuple(
        sorted(
            (survivor for row in atlas_profiles for survivor in row[10]),
            key=lambda row: (row[1], row[0]),
        )
    )
    keys = tuple((row[1], row[0]) for row in survivors)
    heights = tuple(sorted(Counter(row[1] for row in survivors).items()))
    high_keys = tuple(
        (row[1], row[0])
        for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate ledger changed")
    require(keys == EXPECTED_ATLAS_KEYS, ("atlas keys changed", keys))
    require(heights == EXPECTED_HEIGHTS, ("atlas heights changed", heights))
    require(high_keys == EXPECTED_HIGH_KEYS, ("HIGH-TAIL keys changed", high_keys))
    require(all(case in keys and case not in high_keys for case in CLOSURE_CASES), "closure routing changed")
    require(not any(first in {1669, 1671} or 1673 <= first <= 1679 for first, _body in keys), "descent gap occupied")
    atlas_profile_hash = sha256(repr(tuple(atlas_profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    if EXPECTED_ATLAS_PROFILE_SHA256 is not None:
        require(atlas_profile_hash == EXPECTED_ATLAS_PROFILE_SHA256, "atlas profile digest changed")
    if EXPECTED_ATLAS_SURVIVOR_SHA256 is not None:
        require(survivor_hash == EXPECTED_ATLAS_SURVIVOR_SHA256, "atlas survivor digest changed")

    require(len(closure_profiles) == len(CLOSURE_CASES) == 3, "closure universe changed")
    totals = tuple(sum(row[11][index] for row in closure_profiles) for index in range(5))
    kills = sum(row[15] for row in closure_profiles)
    certificate_count = sum(row[21] for row in closure_profiles)
    require(totals == (259, 32, 226, 1, 1), "closure ledger changed")
    require(kills == totals[-1] == 1, "terminal packet survived")
    require(certificate_count == totals[2] == 226, "certificate ledger changed")
    margins = tuple(row[16] for row in closure_profiles if row[16] is not None)
    require(margins == (F(8929, 15197),), "projected margin changed")
    closure_hash = sha256(repr(tuple(closure_profiles)).encode()).hexdigest()
    if EXPECTED_CLOSURE_PROFILE_SHA256 is not None:
        require(closure_hash == EXPECTED_CLOSURE_PROFILE_SHA256, "closure profile digest changed")
    semantic_payload = (
        EXPECTED_THM2972_SOURCE_SHA256,
        EXPECTED_THM2972_OUTPUT_SHA256,
        EXPECTED_ATLAS_ENGINE_SHA256,
        EXPECTED_STATUS_ENGINE_SHA256,
        START,
        END,
        candidate_rows,
        heights,
        high_keys,
        CLOSURE_CASES,
        totals,
        kills,
        certificate_count,
        1668,
        atlas_profile_hash,
        survivor_hash,
        closure_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact descent through z1=1672 and z1=1670",
        f"thm2972_source_sha256={EXPECTED_THM2972_SOURCE_SHA256}",
        f"thm2972_output_sha256={EXPECTED_THM2972_OUTPUT_SHA256}",
        f"atlas_engine_sha256={file_sha256(ATLAS_ENGINE_SOURCE)}",
        f"status_engine_sha256={file_sha256(STATUS_ENGINE_SOURCE)}",
        f"status_output_sha256={file_sha256(STATUS_ENGINE_OUTPUT)}",
        "inherited_cap=1679;dependency=THM-2972",
        "universe=six_body_roots:3003;body_labels:1..14;ordered distinct later nonaligned labels;scalar necessary condition",
        f"atlas_band=1668..1679;candidate_rows={candidate_rows};survivors={len(survivors)};height_counts={heights};high_tail_rows={len(high_keys)}",
        "empty_first_bands=1669,1671,1673..1679;closed_heights=1670,1672;next_occupied_height=1668",
        "status_replay=226 infeasibility certificates verified exactly over Q;solver-selected dual bases excluded from stage and semantic digests",
        f"closure_counts=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};status_survivors:{totals[3]};literal_packets:{totals[4]};projected_kills:{kills};survivors:0",
        "global_minimum_margin=8929/15197;maximum_prefix_cells=1",
    ]
    for first, body in keys:
        lines.append(
            f"ATLAS;z1={first};E={','.join(map(str, body))};"
            f"type={'HIGH-TAIL' if (first, body) in high_keys else 'EXACT'}"
        )
    for row in closure_profiles:
        counts = row[11]
        lines.append(
            f"CASE;z1={row[0]};E={','.join(map(str, row[1]))};h={ftext(row[2])};"
            f"carrier_components={row[3]};L={row[4]};lower={ftext(row[5])};"
            f"delta1={ftext(row[6])};first_d={row[7]};ray_sha256={row[8]};"
            f"divisors={row[9]};trials={row[10]};scalar={counts[0]};crude={counts[1]};"
            f"status={counts[2]};status_survivors={counts[3]};stage_sha256={'/'.join(row[12])};"
            f"cells={row[13]};packets={row[14]};projected_kills={row[15]};"
            f"minimum_margin={ftext(row[16])};maximum_prefix_cells={row[17]};"
            f"minimum_row={row[18]};direct_mass={ftext(row[19])};state_sha256={row[20]};"
            f"exact_status_certificates={row[21]};conclusion=EMPTY"
        )
    lines.extend((
        "consequence=projected k=2 first drift label z1<=1668",
        f"atlas_profile_sha256={atlas_profile_hash}",
        f"atlas_survivor_sha256={survivor_hash}",
        f"closure_profile_sha256={closure_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    os.environ["PYTHONHASHSEED"] = "0"
    configure_atlas()
    bodies = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        atlas_profiles = [atlas_profile(body) for body in bodies]
    else:
        with mp.get_context("spawn").Pool(
            args.workers, initializer=configure_atlas
        ) as pool:
            atlas_profiles = list(pool.imap(atlas_profile, bodies, chunksize=4))
    atlas_profiles.sort(key=lambda row: row[0])
    closure_profiles = tuple(closure_profile(case) for case in CLOSURE_CASES)
    payload = render(tuple(atlas_profiles), closure_profiles)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
