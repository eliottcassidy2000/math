#!/usr/bin/env python3
"""Uniform hybrid closure of all five projected k=2 rows at z1=1780.

Four rows close over every choice of four distinct later nonaligned labels by
the exact ray/status/projected pipeline.  The remaining atlas row is marked
HIGH-TAIL: its exact forced-high optimum (unrestricted top three singleton
excesses plus the best wall-eligible singleton) is already below the scalar
threshold.  Both routes are horizon-free.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENT_SOURCE = (
    ROOT / "04-computation" / "lrc14_j7_k2_z1784_all_label_closure_thm2941.py"
)
PARENT_OUTPUT = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1784_all_label_closure_thm2941.out"
)
OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" / "lrc14_j7_k2_z1780_hybrid_closure_thm2941.out"
)
EXPECTED_PARENT_SOURCE_SHA256 = (
    "65d7fbc94b76f7438fe9fad635687aa69cdfca59fff2c6a2e16f00e203834914"
)
EXPECTED_PARENT_OUTPUT_SHA256 = (
    "17da2b42eb6a88fa9fcec2c4db8e0e06fbb1dfa925df1a4ce30aa6b6df0c65ad"
)
EXPECTED_PROFILE_SHA256 = (
    "ec75c60b591057fb7c864fe0feeb0b871435300d0f273fb8fa0a74c62279d883"
)
EXPECTED_SEMANTIC_SHA256 = (
    "27f5db376a16c44b2e840c51fbedc31bff225ab690fb17ed629d032936fca953"
)

FIRST = 1780
QUANTIFIER = "distinct later nonaligned labels"
HIGH_OBLIGATION = "one later label at or beyond the projected wall"
PROJECTED_WALL_RATIO = F(13, 150)
SUFFIX_SLOTS = 4
EXACT_CASES = tuple(
    (FIRST, body)
    for body in (
        (1, 4, 8, 10, 12, 14),
        (2, 4, 8, 10, 12, 14),
        (2, 6, 8, 10, 12, 14),
        (4, 6, 8, 10, 12, 14),
    )
)
HIGH_CASE = (FIRST, (1, 8, 10, 12, 13, 14))
ALL_CASES = tuple(sorted((*EXACT_CASES, HIGH_CASE)))
EXPECTED_HIGH_GAP = F(-36718081, 292931597140)
EXPECTED_COUNTS = {
    EXACT_CASES[0]: (95, 28, 67, 0, 0),
    EXACT_CASES[1]: (10, 1, 0, 9, 9),
    EXACT_CASES[2]: (9, 2, 0, 7, 10),
    EXACT_CASES[3]: (1, 0, 0, 1, 1),
}
EXPECTED_AUDIT = {
    EXACT_CASES[0]: (
        (
            "f6c19fe22f3d6fa7e5b2f0158917a29429faeefa5cf19affd7c27920a96e414c",
            "1365feaa13a57f7db6b7db84216628acc251e2a551ba237e1e48a66ea46fc6f0",
            "7b27019bda3bf7703f863b3e9403944e6f5ae08265a4a9aab169413d5a492bb6",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        ),
        None,
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    EXACT_CASES[1]: (
        (
            "bf233bafdbac5ed79ebe9e9953af4f89cb0e8763eb6198e6c55ebb0e172138cd",
            "9d0f755a5b11cdfe9273a5856bd9105a67f7fc2fd39ab2eb1d799acd5080e13d",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "d466abc5d28346d15af55b3e77b8b618876317862bc69bc9fd19609d4667098a",
        ),
        F(66, 91),
        1,
        "5180ac3cd8979a342a76143dbfb205fef8c47ce5dcd37920b0f936051428347c",
    ),
    EXACT_CASES[2]: (
        (
            "dc602dcdee55e43431e8f1c522946f06fc5858d69c81ab7cfd53f53b49eb90aa",
            "162b54a872a66a6373f0552a186d06586a358c1eeda088c946a88488e8ac0864",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "e000c6e7794722a5a80f39b3fc69b68a2cfb51bb6acb82fef17fa17f0b7555fd",
        ),
        F(1896001, 2948309),
        2,
        "9d80397dc3faa56312da927071d32e6965d39fb9b6a9e2462ffd67442335a9aa",
    ),
    EXACT_CASES[3]: (
        (
            "0a82481108e1b7c0d6752946fe12b7d3f50d1eb00ec5bdffe9f6d673d66a1ca2",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "3c5f15a43ba398627a2e6d45f29063b77954a6b100ae446ab3d3606d8ea11643",
        ),
        F(8215, 16471),
        1,
        "26c61a0305ecd8aa79d8fe111723c383e30885371d467a243b8f4cdd0a8e4223",
    ),
}
CHECK_EXPECTED_AUDIT = True


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(PARENT_SOURCE) == EXPECTED_PARENT_SOURCE_SHA256, "parent source changed")
require(file_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256, "parent output changed")
SPEC = spec_from_file_location("k2_z1780_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1784 parent")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)
E = P.E
H = P.P.P.H
H.EXPECTED_GAPS[HIGH_CASE] = EXPECTED_HIGH_GAP
ATLAS_KEYS = P.ATLAS_KEYS
ATLAS_HIGH_KEYS = P.P.P.ATLAS_HIGH_KEYS
require(tuple(key for key in ATLAS_KEYS if key[0] == FIRST) == ALL_CASES, "z1780 rows changed")
require(
    tuple(key for key in ATLAS_HIGH_KEYS if key[0] == FIRST) == (HIGH_CASE,),
    "z1780 forced-high atlas dependency changed",
)
require(
    H.HIGH_WALL_RATIO == E.U.suffix.PROJECTED_RATIOS[2] == PROJECTED_WALL_RATIO,
    "projected-wall ratio changed",
)
require(H.SUFFIX_SLOTS == SUFFIX_SLOTS, "forced-high suffix arity changed")
require(max(first for first, _body in ATLAS_KEYS if first < FIRST) == 1776, "next height changed")


def exact_profile(case):
    first, body = case
    carrier = E.U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), E.U.suffix.A.RULER)
    lower = h * E.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    result = E.ray_and_status(first, body, carrier, h, lower, L)
    amplitudes, ray_digest, divisor_count, trials, first_delta, first_d, scalar, crude, status, states, stages = result
    projected = E.projected_packets(first, body, carrier, h, lower, L, amplitudes, first_delta, states) if states else (0, 0, 0, None, 0, None, None, sha256(b"()").hexdigest(), ())
    counts = (len(scalar), len(crude), len(status), len(states), projected[1])
    require(counts == EXPECTED_COUNTS[case], (case, "counts changed", counts))
    require(projected[1] == projected[2], (case, "packet survived"))
    audit = (stages, projected[3], projected[4], projected[7])
    if CHECK_EXPECTED_AUDIT:
        require(audit == EXPECTED_AUDIT[case], (case, "audit changed", audit))
    return ("ALL-LABEL", first, body, h, len(carrier), L, lower, first_delta, first_d, ray_digest, divisor_count, trials, counts, stages, *projected[:-1])


def profile(task):
    route, first, body = task
    if route == "ALL-LABEL":
        return exact_profile((first, body))
    require(route == "HIGH-SCALAR", (task, "unknown route"))
    return (route, *H.profile((first, body)))


TASKS = tuple([("ALL-LABEL", *case) for case in EXACT_CASES] + [("HIGH-SCALAR", *HIGH_CASE)])


def render(profiles):
    exact_rows = tuple(row for row in profiles if row[0] == "ALL-LABEL")
    high_rows = tuple(row for row in profiles if row[0] == "HIGH-SCALAR")
    require(len(exact_rows) == 4 and len(high_rows) == 1, "route split changed")
    totals = tuple(sum(row[12][index] for row in exact_rows) for index in range(5))
    kills = sum(row[16] for row in exact_rows)
    require(totals == (115, 31, 67, 17, 20) and kills == 20, "ledger changed")
    require(high_rows[0][-1] == EXPECTED_HIGH_GAP < 0, "high gap changed")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (EXPECTED_PARENT_SOURCE_SHA256, EXPECTED_PARENT_OUTPUT_SHA256, EXACT_CASES, HIGH_CASE, QUANTIFIER, HIGH_OBLIGATION, PROJECTED_WALL_RATIO, SUFFIX_SLOTS, totals, kills, EXPECTED_HIGH_GAP, profile_hash)
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 hybrid exact closure at z1=1780",
        f"parent_source_sha256={file_sha256(PARENT_SOURCE)}",
        f"parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        f"scope=four rows over all {QUANTIFIER};one row under {HIGH_OBLIGATION};no finite label horizon",
        "routes=four all-label ray/status/projected closures;one exact forced-high scalar maximum",
        "forced_high_dependency=atlas HIGH-TAIL;projected wall=floor(13L/150)+1;four suffix slots",
        "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);two-aligned cap=25/91",
        f"all_label_counts=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};status_survivors:{totals[3]};literal_packets:{totals[4]};projected_kills:{kills};survivors:0",
    ]
    for row in profiles:
        if row[0] == "ALL-LABEL":
            (_route, first, body, h, components, L, lower, first_delta, first_d, ray_digest, divisor_count, trials, counts, stages, cells, packets, packet_kills, margin, max_cells, minimum_row, direct_mass, state_digest) = row
            lines.append(f"CASE;route=ALL-LABEL;z1={first};E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};lower={ftext(lower)};delta1={ftext(first_delta)};first_d={first_d};ray_sha256={ray_digest};denominators={divisor_count};trials={trials};counts={counts};stage_sha256={stages};body_cells={cells};packets={packets};kills={packet_kills};min_margin={ftext(margin)};max_prefix_cells={max_cells};direct_control_mass={ftext(direct_mass)};state_sha256={state_digest};minimum_row={minimum_row};conclusion=EMPTY")
        else:
            (_route, first, body, h, components, L, high_floor, first_delta, lower, positive, negative, zero, ray_digest, rank4, omitted, best_high, branch, constrained, upper, gap) = row
            lines.append(f"CASE;route=HIGH-SCALAR;z1={first};E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};high_floor={high_floor};delta1={ftext(first_delta)};lower={ftext(lower)};ray_signs=+{positive}/-{negative}/0:{zero};ray_sha256={ray_digest};unrestricted_top4={rank4};first_omitted={omitted};best_high={best_high};branch={branch};constrained={constrained};upper={ftext(upper)};gap={ftext(gap)};conclusion=SCALAR-EMPTY")
    lines.extend(("global_z1780_rows=5;empty=5;survivors=0", "atlas=next occupied height:1776", "consequence=projected k=2 first drift label z1<=1776", f"profile_sha256={profile_hash}", f"semantic_sha256={semantic_hash}", "all_exact_controls=PASS"))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(len(TASKS), mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1 and args.hash_seed >= 0, "invalid execution controls")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        profiles = [profile(task) for task in TASKS]
    else:
        with mp.get_context("spawn").Pool(min(args.workers, len(TASKS))) as pool:
            profiles = list(pool.imap(profile, TASKS))
    order = {task: index for index, task in enumerate(TASKS)}
    profiles.sort(key=lambda row: order[(row[0], row[1], row[2])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
