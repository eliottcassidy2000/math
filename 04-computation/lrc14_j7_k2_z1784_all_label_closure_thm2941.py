#!/usr/bin/env python3
"""Uniform all-label closure of all six projected k=2 rows at z1=1784.

The exact 1750..1799 atlas has six rows at 1784 and no row on 1785..1789.
One row was routed HIGH-TAIL by the scalar atlas, but the stronger all-label
ray/status/projected pipeline closes it without using that obligation.  Thus
this verifier treats every row uniformly over every choice of four distinct
later nonaligned labels, with no finite label horizon.
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
    ROOT / "04-computation" / "lrc14_j7_k2_z1788_ray_status_closure_thm2941.py"
)
PARENT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1788_ray_status_closure_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_z1784_all_label_closure_thm2941.out"
)
EXPECTED_PARENT_SOURCE_SHA256 = (
    "ecd2f8e8179545f971ed6bd4a33b816c0f34f0693aeb0dda37c3a5c445a70435"
)
EXPECTED_PARENT_OUTPUT_SHA256 = (
    "c1090a83e242f2ded17e8d8aabd74983f0925466f55a85c33bfccfbe385f57b6"
)
EXPECTED_PROFILE_SHA256 = (
    "8edd0da2514cbdb502303523bfcaf547662b447780d0e8c8efb58386177332ab"
)
EXPECTED_SEMANTIC_SHA256 = (
    "e49cc640ba9654b5a0331b387368414242029632c09f674babd08f0068a12137"
)

FIRST = 1784
QUANTIFIER = "distinct later nonaligned labels"
CASES = tuple(
    (FIRST, body)
    for body in (
        (1, 4, 8, 10, 12, 14),
        (1, 6, 8, 10, 12, 14),
        (2, 4, 8, 10, 12, 14),
        (2, 6, 8, 10, 12, 14),
        (2, 8, 9, 10, 12, 14),
        (4, 6, 8, 10, 12, 14),
    )
)
EXPECTED_COUNTS = {
    CASES[0]: (96, 28, 68, 0, 0),
    CASES[1]: (5, 4, 0, 1, 3),
    CASES[2]: (10, 1, 0, 9, 9),
    CASES[3]: (20, 3, 0, 17, 26),
    CASES[4]: (7, 0, 0, 7, 7),
    CASES[5]: (2, 1, 0, 1, 1),
}
EXPECTED_AUDIT = {
    CASES[0]: (
        (
            "7219e50cd5bd12d3fe822e554f5176ccd6890e0ad966fa717f939fd08d500825",
            "93487db134709cd1233f8fd1d89ac0f859552b52ab2db997ccdb314b876d29d7",
            "9a2f883b742c2833d995d7168d5fe9401a18e5f0e45a27a32ecb4a8674d317ad",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
        ),
        None,
        0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    CASES[1]: (
        (
            "2f928984e637f3fecd26fabd31b9a2844dcd9b59f74d7295b6a71457f5aa5330",
            "9ffe30f1f92b78705718308fb2e681eca8e868b8ba2aa049c146a974312b19f0",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "22565beb6b5168d1809013afbbe31abb8392866cc9aa90a10fb0945144f0b597",
        ),
        F(1026, 16471),
        2,
        "a24458c835a52f562af8121e2929d149a32f7d37d7640b2a461500dd0733b4d3",
    ),
    CASES[2]: (
        (
            "75d22d7f852e1d4d8a3b8fcc30de02d8549524103b0a6aa0560ddc1c53eb0c05",
            "ce9e968587584ebf882fc705a4947526bdeb6d74bedc7aa927975e04dd0b4146",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "ca475eba02a678f80dec0a980a06c07815380c475711eeb358fbaa28338f318c",
        ),
        F(66, 91),
        1,
        "4e59788066e4b8df7cff4df040b9ff4b6a70da52e8963b23f57c1ee0d2f6b4bc",
    ),
    CASES[3]: (
        (
            "d2c0d6e30a0549e8fb063ad6d40c6c41b8014ba96e70b47570091301eecaa0fa",
            "a8d9384d26874fecc4a527ec15b553629c95edd50ce03e4383a1e4be643a4ec7",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "a29400326b044d2e2a41c0644e92cefc5e7d547c22916b7c43d676aac80dc890",
        ),
        F(723886, 3632447),
        2,
        "e3157706a0a603126ff4ea078aca3a7f85560d7eaee5dee6695ee19fb722a68f",
    ),
    CASES[4]: (
        (
            "4cfd43465d0df7330dbcacde22ad98d2dd3e5d7fb110799f5a9c7890f25d80b3",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "ff3c26600496b16cd00ac5c6feae2386fab7c462de6bc2bd1e983cae0529b917",
        ),
        F(66, 91),
        4,
        "31ea9cf985f7a7035c315de9dd85998f693f832ef1fbb8b4924559b465077d75",
    ),
    CASES[5]: (
        (
            "d60decbd1f10aacfad435a71880034f95b51e07fd1fe9c8c372353df92dc3303",
            "78a41a5c8e2a59263764c087031bb3e69766b1660dfae20f7757f2a1e1690df9",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "cbd90e87bb307dc83b35231258f3e1f8cc73aae7bacb05de762d5958b9f26b1f",
        ),
        F(3980, 20293),
        1,
        "6eef0e115b007ec5d2dc94d8ea223f04fc06aaa8bafce7806e1dd8a2e61c016e",
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
SPEC = spec_from_file_location("k2_z1784_parent", PARENT_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1788 parent")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)
E = P.P.E
ATLAS_KEYS = P.P.ATLAS_KEYS
require(tuple(key for key in ATLAS_KEYS if key[0] == FIRST) == CASES, "z1784 atlas rows changed")
require(
    not any(first in (1785, 1786, 1787, 1789) for first, _body in ATLAS_KEYS),
    "unproved upper atlas gap changed",
)
require(max(first for first, _body in ATLAS_KEYS if first < FIRST) == 1780, "next height changed")


def profile(case):
    first, body = case
    carrier = E.U.suffix.A.carrier_for(body)
    require(E.P.A.carrier_for(body) == carrier, (case, "carrier engines disagree"))
    h = F(sum(right - left for left, right in carrier), E.U.suffix.A.RULER)
    lower = h * E.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    (
        amplitudes,
        ray_digest,
        divisor_count,
        trials,
        first_delta,
        first_d,
        scalar,
        crude_kills,
        status_kills,
        states,
        stage_digests,
    ) = E.ray_and_status(first, body, carrier, h, lower, L)
    if states:
        projected = E.projected_packets(first, body, carrier, h, lower, L, amplitudes, first_delta, states)
    else:
        projected = (0, 0, 0, None, 0, None, None, sha256(b"()").hexdigest(), ())
    counts = (len(scalar), len(crude_kills), len(status_kills), len(states), projected[1])
    require(counts == EXPECTED_COUNTS[case], (case, "counts changed", counts))
    require(projected[1] == projected[2], (case, "a projected packet survived"))
    audit = (stage_digests, projected[3], projected[4], projected[7])
    if CHECK_EXPECTED_AUDIT:
        require(audit == EXPECTED_AUDIT[case], (case, "audit changed", audit))
    return (
        first, body, h, len(carrier), L, lower, first_delta, first_d,
        ray_digest, divisor_count, trials, counts, stage_digests, *projected[:-1],
    )


def render(profiles):
    require(tuple((row[0], row[1]) for row in profiles) == CASES, "case order changed")
    totals = tuple(sum(row[11][index] for row in profiles) for index in range(5))
    projected_kills = sum(row[15] for row in profiles)
    require(totals == (140, 37, 68, 35, 46), "global ledger changed")
    require(projected_kills == totals[-1] == 46, "global packet survived")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_PARENT_SOURCE_SHA256, EXPECTED_PARENT_OUTPUT_SHA256, CASES,
        QUANTIFIER, totals, projected_kills, profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 all-label exact closure at z1=1784",
        f"parent_source_sha256={file_sha256(PARENT_SOURCE)}",
        f"parent_output_sha256={file_sha256(PARENT_OUTPUT)}",
        f"scope=six z1784 atlas rows;all {QUANTIFIER};no finite label horizon",
        "pipeline=residue rays;denominator multisets;all-divisor crude;common-K5 status;scalar-slack literal packets;projected residual",
        "projected_identity=P_(E,Z)=phi_L(C_E minus union_z D_z);two-aligned cap=25/91",
        f"global_counts=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};status_survivors:{totals[3]};literal_packets:{totals[4]};projected_kills:{projected_kills};survivors:0",
    ]
    for row in profiles:
        (
            first, body, h, components, L, lower, first_delta, first_d,
            ray_digest, divisor_count, trials, counts, stage_digests, body_cells,
            packets, kills, minimum_margin, maximum_cells, minimum_row,
            direct_mass, state_digest,
        ) = row
        lines.append(
            f"CASE;z1={first};E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};"
            f"lower={ftext(lower)};delta1={ftext(first_delta)};first_d={first_d};"
            f"ray_sha256={ray_digest};denominators={divisor_count};trials={trials};"
            f"counts={counts};stage_sha256={stage_digests};body_cells={body_cells};"
            f"packets={packets};kills={kills};min_margin={ftext(minimum_margin)};"
            f"max_prefix_cells={maximum_cells};direct_control_mass={ftext(direct_mass)};"
            f"state_sha256={state_digest};minimum_row={minimum_row};conclusion=EMPTY"
        )
    lines.extend((
        "global_z1784_rows=6;empty=6;survivors=0",
        "atlas=empty:1785..1790 after inherited closures;next occupied height:1780",
        "consequence=projected k=2 first drift label z1<=1780",
        f"profile_sha256={profile_hash}", f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(len(CASES), mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        profiles = [profile(case) for case in CASES]
    else:
        with mp.get_context("spawn").Pool(min(args.workers, len(CASES))) as pool:
            profiles = list(pool.imap(profile, CASES))
    order = {case: index for index, case in enumerate(CASES)}
    profiles.sort(key=lambda row: order[(row[0], row[1])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
