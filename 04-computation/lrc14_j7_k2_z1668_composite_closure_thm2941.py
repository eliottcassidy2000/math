#!/usr/bin/env python3
"""Close the projected k=2 shells at z1=1668 and z1=1665 exactly."""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
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
BASE_NAME = "lrc14_j7_k2_z1672_z1670_exact_descent_closure_thm2941"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
if not BASE_SOURCE.exists():
    BASE_SOURCE = HERE.parent / f"{BASE_NAME}.py"
if not BASE_OUTPUT.exists():
    BASE_OUTPUT = HERE.parent / f"{BASE_NAME}.out"
TORSION_SOURCE = ROOT / "04-computation/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.py"
TORSION_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_exact_one_high_located_torsion_closure_thm2970.out"
OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k2_z1668_composite_closure_thm2941.out"
    if HERE.parent.name == "04-computation"
    else HERE.with_suffix(".out")
)

EXPECTED_BASE_SOURCE_SHA256 = "b87b7e0e73207e026323269cb2493cf0793456b8292eb6effe94cdd823bb11e4"
EXPECTED_BASE_OUTPUT_SHA256 = "bfcf783f84037f9b1417b10c4a1a6ec71a2aeec54c8f20383a1752a3d3149134"
EXPECTED_TORSION_SOURCE_SHA256 = "eeba4fa6968365269b776832cbcbd3753affcac10915726eb28fb7a817c9eec4"
EXPECTED_TORSION_OUTPUT_SHA256 = "8f1dd4702b64caef2b94cfdf8e474bba137900ef48200bf5d9aabb2e440e454f"
EXPECTED_ATLAS_PROFILE_SHA256 = "1fcd59d389411cc61911bf0ce4e927b1b34d16b1047acf011d0089e5a297e44a"
EXPECTED_ATLAS_SURVIVOR_SHA256 = "5c84b42ba96b8a04ce868c4fbad686c73118d063dac478a79f9cc1ec54c7aa00"
EXPECTED_CLOSURE_PROFILE_SHA256 = "9a68432607b516375a650df2d6fc05b197cd6dae26b077d2cd55c27c474c580b"
EXPECTED_SEMANTIC_SHA256 = "17f0eaf3fbe0a3e6aab565b5d862d5307160f13b47a56b89f55882826d95d1f4"

START = 1660
END = 1668
EXPECTED_CANDIDATE_ROWS = 27_027
EXPECTED_HEIGHTS = ((1660, 1), (1665, 1), (1668, 9))
NEXT_CASE = (1660, (1, 4, 8, 10, 12, 14))
DESCENT_CASE = (1665, (1, 4, 8, 10, 12, 14))
HIGH_CASES = tuple(
    (1668, body)
    for body in (
        (1, 2, 10, 11, 12, 14),
        (1, 6, 10, 11, 12, 14),
        (1, 8, 10, 11, 12, 14),
        (2, 8, 10, 11, 12, 14),
    )
)
ORDINARY_CASES = tuple(
    case
    for case in (
        DESCENT_CASE,
        (1668, (1, 4, 8, 10, 12, 14)),
        (1668, (1, 6, 8, 10, 12, 14)),
        (1668, (1, 6, 9, 10, 12, 14)),
        (1668, (2, 4, 8, 10, 12, 14)),
        (1668, (2, 6, 8, 10, 12, 14)),
    )
)
ALL_CLOSURE_CASES = tuple(sorted((*HIGH_CASES, *ORDINARY_CASES)))
EXPECTED_ATLAS_KEYS = tuple(sorted((NEXT_CASE, *ALL_CLOSURE_CASES)))
EXPECTED_COUNTS = {
    ORDINARY_CASES[0]: (1, 0, 0, 1, 1),
    ORDINARY_CASES[1]: (252, 46, 206, 0, 0),
    ORDINARY_CASES[2]: (11, 2, 9, 0, 0),
    ORDINARY_CASES[3]: (1, 1, 0, 0, 0),
    ORDINARY_CASES[4]: (32, 0, 14, 18, 18),
    ORDINARY_CASES[5]: (1, 0, 0, 1, 1),
}
EMPTY_DIGEST = sha256(b"()").hexdigest()
EXPECTED_ORDINARY_AUDIT = {
    ORDINARY_CASES[0]: (
        (
            "bf7c7af434b2dff91276da70d1ce104d5ee8aab769b430b0672514da910b04d5",
            EMPTY_DIGEST,
            EMPTY_DIGEST,
            "327308388a3a66f064a2bf276d08c1b551f105225395672b250fc70c9b2dea80",
        ), F(66, 91), 1, F(1),
        "c57692cbbb649e1342b6ef1c9cd034e087193accb4c6539e6c939ed022be288c",
    ),
    ORDINARY_CASES[1]: (
        (
            "4a4df5edb79e2f4846fc1e34edd21de9ed5ab9245e490c48f8d3f02edb77ed0b",
            "cdc2ba58f0a3262e634ad8ec54a255001953b0cd49358fbbf6ebd8332984c3d4",
            "2e8c226bb476ff8b19a686f793214f5012f5e0726c40dadbd438b87428f3b532",
            EMPTY_DIGEST,
        ), None, 0, None, EMPTY_DIGEST,
    ),
    ORDINARY_CASES[2]: (
        (
            "28cdcb12a94f51b6f9f30ae291b3df149c436759c94d2d6ed72221abe87720d8",
            "956d8580bbf269b5a1dc9527398d2d422fd4f2e9d6f6550cf34c449719fd3928",
            "c56b261cb4a4893653a53851a853f4423d9c1b41f4a13eb36fea7d484b70c1d7",
            EMPTY_DIGEST,
        ), None, 0, None, EMPTY_DIGEST,
    ),
    ORDINARY_CASES[3]: (
        (
            "764ff30e5f14bdd26a1ebe5ee7d333638212c3d6b4094bbf128b84eecf8b4c5a",
            "c564f17ee425880e1332d8d35afde9ffdd729283e154e73b0b81e45d0703ea72",
            EMPTY_DIGEST,
            EMPTY_DIGEST,
        ), None, 0, None, EMPTY_DIGEST,
    ),
    ORDINARY_CASES[4]: (
        (
            "d4424abb552b85732244f590248511ddc5993ccb3eaa201777f9ef4e8ae8946d",
            EMPTY_DIGEST,
            "627e0fe6ee6c0cc3867194441a9e2d6ea272fd6479a1e67f4359895fb339c72b",
            "dc064f2526a44d71b229ea5f479dbe18a54944af141bea8d912fc8791e1fdba4",
        ), F(1296, 5551), 1, F(1),
        "bb4a53b80ebec6e1bcb005b285ba7f33c0f2df169ab084435b1fafa5b7fdd405",
    ),
    ORDINARY_CASES[5]: (
        (
            "f00cc91711fcb09fc8f6aa6453943d556befa3950d2713e7180fbc365ddddc4b",
            EMPTY_DIGEST,
            EMPTY_DIGEST,
            "5aa5a322c166423a5da036808cd19c393814a09fb2c7a9f534a93c458bca88e1",
        ), F(66, 91), 1, F(1),
        "38e26ec7f85617357208b879b19b06a3e9513157e15967b84c8250fe33ffde6e",
    ),
}
EXPECTED_HIGH_GAPS = {
    HIGH_CASES[0]: (
        F(87451119983, 135501948281700),
        F(81911847151, 677509741408500),
        (1800, 1836, 1750), 4620, 6062,
    ),
    HIGH_CASES[1]: (
        F(43898289181, 60518299341500),
        F(5759211506579, 27777899397748500),
        (1750, 1722, 1836), 4620, 6062,
    ),
    HIGH_CASES[2]: (
        F(11241200339, 13397139929628),
        F(9319641068953, 45684247160031480),
        (1836, 1708, 1736), 2156, 11580,
    ),
    HIGH_CASES[3]: (
        F(7716527477, 9048177294520),
        F(3136055955749, 15228082386677160),
        (1708, 1736, 1836), 2156, 11580,
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return "NONE" if value is None else f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1672/z1670 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1672/z1670 output changed")
require(file_sha256(TORSION_SOURCE) == EXPECTED_TORSION_SOURCE_SHA256, "torsion source changed")
require(file_sha256(TORSION_OUTPUT) == EXPECTED_TORSION_OUTPUT_SHA256, "torsion output changed")

BSPEC = spec_from_file_location("k2_z1668_base", BASE_SOURCE)
require(BSPEC is not None and BSPEC.loader is not None, "cannot load z1672/z1670 base")
BASE = module_from_spec(BSPEC)
BSPEC.loader.exec_module(BASE)
TSPEC = spec_from_file_location("k2_z1668_torsion", TORSION_SOURCE)
require(TSPEC is not None and TSPEC.loader is not None, "cannot load torsion engine")
TORSION = module_from_spec(TSPEC)
TSPEC.loader.exec_module(TORSION)
ENGINE = BASE.ENGINE


def configure_atlas():
    BASE.ATLAS.B.START = START
    BASE.ATLAS.B.END = END
    BASE.ATLAS.B.LABELS = np.arange(START, BASE.ATLAS.B.HORIZON + 1, dtype=np.int64)


def atlas_profile(body):
    return BASE.ATLAS.B.profile(body)


def ordinary_profile(case):
    first, body = case
    carrier = ENGINE.U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), ENGINE.U.suffix.A.RULER)
    lower = h * ENGINE.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    result = ENGINE.ray_and_status(first, body, carrier, h, lower, L)
    certificate_count = BASE.exact_status_audit(case, L, result[8], result[10][2])
    projected = (
        ENGINE.projected_packets(first, body, carrier, h, lower, L, result[0], result[4], result[9])
        if result[9]
        else (0, 0, 0, None, 0, None, None, EMPTY_DIGEST, ())
    )
    counts = (len(result[6]), len(result[7]), len(result[8]), len(result[9]), projected[1])
    require(counts == EXPECTED_COUNTS[case], (case, "counts changed", counts))
    require(certificate_count == counts[2], (case, "certificate count"))
    require(projected[1] == projected[2], (case, "projected survivor"))
    audit = (result[10], projected[3], projected[4], projected[6], projected[7])
    require(audit == EXPECTED_ORDINARY_AUDIT[case], (case, "ordinary audit changed", audit))
    return (
        "ALL-LABEL", first, body, h, len(carrier), L, lower, result[4], result[5],
        result[1], result[2], result[3], counts, result[10], *projected[:-1], certificate_count,
    )


def high_profile(case):
    first, body = case
    modulus = 14 * lcm(*body)
    carrier, mass_numerator, amplitudes = TORSION.amplitude_data(body, modulus)
    h = F(mass_numerator, TORSION.RULER)
    first_delta = TORSION.delta(amplitudes, first)
    lower = h / 91
    high_floor = (13 * modulus) // 150 + 1
    low_rows = sorted(
        ((TORSION.delta(amplitudes, label), label) for label in range(first + 1, high_floor)),
        key=lambda item: (-item[0], item[1]),
    )
    required = lower - first_delta
    maxima = TORSION.high_ray_maxima(amplitudes, modulus, high_floor)
    best_d, best_high = max(maxima.items(), key=lambda item: item[1][0])
    two_gap = required - (
        2 * best_high[0]
        + max(best_high[0], low_rows[0][0])
        + max(best_high[0], low_rows[1][0])
    )
    exact_gap = required - sum((value for value, _label in low_rows[:3]), F()) - best_high[0]
    audit = (two_gap, exact_gap, tuple(label for _value, label in low_rows[:3]), best_d, best_high[1])
    require(audit == EXPECTED_HIGH_GAPS[case], (case, "high audit changed", audit))
    require(two_gap > 0 and exact_gap > 0, (case, "forced-high row survived"))
    amplitude_hash = sha256(repr(tuple(int(value) for value in amplitudes)).encode()).hexdigest()
    return (
        "HIGH-EXACT-RAY", first, body, h, len(carrier), modulus, high_floor,
        first_delta, lower, two_gap, exact_gap, audit[2], best_d, best_high[1], amplitude_hash,
    )


def closure_profile(case):
    return ordinary_profile(case) if case in ORDINARY_CASES else high_profile(case)


def render(atlas_profiles, closure_profiles):
    require(len(atlas_profiles) == comb(14, 6) == 3003, "atlas universe changed")
    candidate_rows = sum(row[7] for row in atlas_profiles)
    survivors = tuple(
        sorted((survivor for row in atlas_profiles for survivor in row[10]), key=lambda row: (row[1], row[0]))
    )
    keys = tuple((row[1], row[0]) for row in survivors)
    heights = tuple(sorted(Counter(row[1] for row in survivors).items()))
    high_keys = tuple(
        (row[1], row[0]) for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate ledger changed")
    require(keys == EXPECTED_ATLAS_KEYS, ("atlas keys changed", keys))
    require(heights == EXPECTED_HEIGHTS, ("atlas heights changed", heights))
    require(high_keys == HIGH_CASES, ("HIGH partition changed", high_keys))
    require(
        not any(1661 <= first <= 1664 or first in (1666, 1667) for first, _body in keys),
        "descent gap occupied",
    )
    atlas_hash = sha256(repr(tuple(atlas_profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    if EXPECTED_ATLAS_PROFILE_SHA256 is not None:
        require(atlas_hash == EXPECTED_ATLAS_PROFILE_SHA256, "atlas digest changed")
    if EXPECTED_ATLAS_SURVIVOR_SHA256 is not None:
        require(survivor_hash == EXPECTED_ATLAS_SURVIVOR_SHA256, "survivor digest changed")

    ordinary = tuple(row for row in closure_profiles if row[0] == "ALL-LABEL")
    high = tuple(row for row in closure_profiles if row[0] == "HIGH-EXACT-RAY")
    require((len(ordinary), len(high)) == (6, 4), "closure partition changed")
    totals = tuple(sum(row[12][index] for row in ordinary) for index in range(5))
    kills = sum(row[16] for row in ordinary)
    certificates = sum(row[22] for row in ordinary)
    require(totals == (298, 49, 229, 20, 20), "ordinary ledger changed")
    require(kills == totals[-1] == 20, "ordinary packet survived")
    require(certificates == totals[2] == 229, "certificate ledger changed")
    require(min(row[17] for row in ordinary if row[17] is not None) == F(1296, 5551), "margin changed")
    closure_hash = sha256(repr(tuple(closure_profiles)).encode()).hexdigest()
    if EXPECTED_CLOSURE_PROFILE_SHA256 is not None:
        require(closure_hash == EXPECTED_CLOSURE_PROFILE_SHA256, "closure digest changed")
    semantic_hash = sha256(repr((
        EXPECTED_BASE_SOURCE_SHA256, EXPECTED_BASE_OUTPUT_SHA256,
        EXPECTED_TORSION_SOURCE_SHA256, START, END, candidate_rows, heights,
        high_keys, totals, kills, certificates, 1660, atlas_hash, survivor_hash, closure_hash,
    )).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=2 composite closure of the z1=1668 and z1=1665 shells",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        f"torsion_source_sha256={file_sha256(TORSION_SOURCE)}",
        f"torsion_output_sha256={file_sha256(TORSION_OUTPUT)}",
        "inherited_cap=1668;scope=all nine z1668 rows and the sole z1665 row;all distinct later nonaligned labels;no finite label horizon",
        f"atlas_band=1660..1668;candidate_rows={candidate_rows};height_counts={heights};empty_firsts=1661..1664,1666,1667;next_occupied_height=1660",
        "partition=4 forced-high exact-ray scalar-empty rows;6 unrestricted all-label rows",
        f"ordinary_counts=scalar:{totals[0]};crude:{totals[1]};status:{totals[2]};status_survivors:{totals[3]};literal_packets:{totals[4]};projected_kills:{kills};survivors:0",
        "global_minimum_margin=1296/5551;maximum_prefix_cells=1",
        "status_replay=229 infeasibility certificates verified exactly over Q;solver-selected dual bases excluded from digests",
    ]
    for row in closure_profiles:
        if row[0] == "HIGH-EXACT-RAY":
            lines.append(
                f"CASE;route=HIGH-EXACT-RAY;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};high_floor={row[6]};"
                f"delta1={ftext(row[7])};lower={ftext(row[8])};two_high_gap={ftext(row[9])};"
                f"exact_one_high_gap={ftext(row[10])};top_lows={row[11]};best_high_d={row[12]};"
                f"best_high_label={row[13]};amplitude_sha256={row[14]};conclusion=SCALAR-EMPTY"
            )
        else:
            counts = row[12]
            lines.append(
                f"CASE;route=ALL-LABEL;z1={row[1]};E={','.join(map(str,row[2]))};"
                f"h={ftext(row[3])};r={row[4]};L={row[5]};lower={ftext(row[6])};"
                f"delta1={ftext(row[7])};first_d={row[8]};ray_sha256={row[9]};"
                f"divisors={row[10]};trials={row[11]};counts={counts};stage_sha256={row[13]};"
                f"cells={row[14]};packets={row[15]};kills={row[16]};minimum_margin={ftext(row[17])};"
                f"maximum_prefix_cells={row[18]};direct_mass={ftext(row[20])};state_sha256={row[21]};"
                f"exact_status_certificates={row[22]};conclusion=EMPTY"
            )
    lines.extend((
        "consequence=projected k=2 first drift label z1<=1660",
        f"atlas_profile_sha256={atlas_hash}",
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
        with mp.get_context("spawn").Pool(args.workers, initializer=configure_atlas) as pool:
            atlas_profiles = list(pool.imap(atlas_profile, bodies, chunksize=4))
    atlas_profiles.sort(key=lambda row: row[0])
    if args.workers == 1:
        closure_profiles = tuple(closure_profile(case) for case in ALL_CLOSURE_CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            closure_profiles = tuple(pool.imap(closure_profile, ALL_CLOSURE_CASES))
    payload = render(tuple(atlas_profiles), closure_profiles)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
