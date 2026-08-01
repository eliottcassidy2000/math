#!/usr/bin/env python3
"""Exact closure of all five projected k=2 scalar rows at z1=1656.

The canonical THM-2941 lower atlas leaves exactly five rows at its boundary.
Four have ruler 11760.  For those rows this referee reuses the canonical
all-label chain

    residue rays -> all-divisor capacity -> common K5 status
    -> finite literal packets -> lossless projected residual.

Every returned Farkas vector is rebuilt and checked over Q, while the replay
digest retains only the deterministic infeasible instance, never the
solver-selected dual basis.  The fifth row has ruler 152880 and lies below
the forced-high wall.  Its unrestricted top four pass the scalar inequality,
but replacing the fourth term by the exact best wall-eligible ray point makes
the sharp constrained maximum negative.  Thus all five rows are empty and
the projected k=2 first-drift bound improves from 1656 to 1655.

This is a finite-exact boundary input for reserved THM-2980.  It does not
prove the full 1600..1679 band theorem or LRC(14).
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "04-computation").is_dir())
BASE_NAME = "lrc14_j7_k2_z1660_status_descent_closure_thm2941"
BASE_SOURCE = ROOT / f"04-computation/{BASE_NAME}.py"
BASE_OUTPUT = ROOT / f"05-knowledge/results/{BASE_NAME}.out"
OUTPUT_PATH = (
    ROOT / "05-knowledge/results/lrc14_j7_k2_z1656_boundary_closure_thm2980.out"
    if HERE.parent.name == "04-computation"
    else HERE.with_suffix(".out")
)

EXPECTED_BASE_SOURCE_SHA256 = "b7d3344888034c1abc4cc5da72b73c8a81992e287ba7fc843b35cfbc6d5b305a"
EXPECTED_BASE_OUTPUT_SHA256 = "1a1a85083ebb382d478139ce5206a20131ed412c5e47c61b26b58975cb2fed2e"
EXPECTED_PROFILE_SHA256 = "0928270d68caccbd55a4769b86be50a2527daf4b403f8fd66e7292eae713d35c"
EXPECTED_SEMANTIC_SHA256 = "10c0562fa613c8e953bde494d132940dfc39dfb74c5b81b467dfdcb888365321"

FIRST = 1656
QUANTIFIER = "distinct later nonaligned labels"
ALIGNED_TWO_UNION_CAP = F(25, 91)
HIGH_WALL_RATIO = F(13, 150)
SCALAR_ETA = F(1, 91)
ORDINARY_CASES = (
    (FIRST, (1, 2, 8, 10, 12, 14)),
    (FIRST, (1, 4, 8, 10, 12, 14)),
    (FIRST, (1, 6, 8, 10, 12, 14)),
    (FIRST, (2, 4, 8, 10, 12, 14)),
)
EXCEPTIONAL_CASE = (FIRST, (1, 8, 10, 12, 13, 14))
BOUNDARY_CASES = tuple(sorted((*ORDINARY_CASES, EXCEPTIONAL_CASE)))
NEXT_PARENT_CASE = (1660, (1, 4, 8, 10, 12, 14))

EXPECTED_ORDINARY = {
    ORDINARY_CASES[0]: (
        (28, 0, 12, 16),
        (
            "c3fe8c97838f623ea6ef50db1ad7ebe9d7a3f1298cb96fde7ebfbaa751a0ce4d",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "b88dad0fafd9fe5ae07459a84b9face2b812c1f84eb8cd1028e8108635710dbe",
            "ae199c050f256ad5a44643cd74db510305b699a2c7335bffe2b4639c4cf0a482",
        ),
        (
            4196,
            19,
            19,
            F(4390, 320229),
            4,
            (F(4390, 320229), (1656, 1736, 1836, 1992, 2142), F(1015, 3519), 3),
            F(1),
            "84c2ed9e18d334b34141815b6bf237e8f9779e3c492402fc7365f94df642cf14",
        ),
    ),
    ORDINARY_CASES[1]: (
        (474, 82, 391, 1),
        (
            "e9b13f81d111507f657f883b12ea503631fa8aabd832a83f89526d83d696e349",
            "dde62f8b68da08f939f96cfd6c21e559825ab381b6172dfe7a6087bdbc4d26bf",
            "2e87b0765565f78e3b2be1f6d9b4255cbb15f73d6bddc053857cce4e89ade871",
            "f9ad5639ad74dcfb46bebdd44286cf910844e0fef7f33465bef26dcf523b864a",
        ),
        (
            4196,
            1,
            1,
            F(66, 91),
            1,
            (F(66, 91), (1656, 1708, 1743, 1836, 2060), F(1), 1),
            F(1),
            "3f3b2d8cd64378c70429a4f8f5e3dccc89f8433cf1a2824250c36e8b2e7f638d",
        ),
    ),
    ORDINARY_CASES[2]: (
        (45, 9, 28, 8),
        (
            "6b663f793a2e7d5ef6e05b1c031330d9606c4ddb276bbda88cfa3a3461209c92",
            "818a619261a5e892f7f65cbbaea0110abda3caa38fc86db9eb00e209e50f3c0b",
            "117e844889f93436562eed805e766f4bcca37a6f89201a4b717a80a295be7421",
            "b7abc8a2e30718802eee1a2a538565945c63af6e2e4fc3b3469c413fdcebbea9",
        ),
        (
            4196,
            8,
            8,
            F(109909, 667849),
            4,
            (F(109909, 667849), (1656, 1722, 1790, 1810, 1836), F(3224, 7339), 1),
            F(1),
            "dbba10575ec67e33c94d5e9b6e28d3add12fe2876ed5fc4be91cbddf7744b430",
        ),
    ),
    ORDINARY_CASES[3]: (
        (11, 0, 3, 8),
        (
            "4780a1b6ed9cb5246a71e911af998a59bcaafd9ea947d80e82b995a8dbe91686",
            "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
            "ced34003250544629a5babb727c68fcd620ef451a4caa81b85852466e5163796",
            "c8dc21a6c3fbc39a60cd73fd911984addc8baf67f2123960f69922a572066a59",
        ),
        (
            4496,
            8,
            8,
            F(1296, 5551),
            1,
            (F(1296, 5551), (1656, 1708, 1750, 1810, 1836), F(31, 61), 1),
            F(1),
            "cc2ba4dde0d751568e47bb1106128b809a29f9eec8c8a7f7ae835df81da20831",
        ),
    ),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE_SOURCE) == EXPECTED_BASE_SOURCE_SHA256, "z1660 source changed")
require(file_sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256, "z1660 output changed")
SPEC = spec_from_file_location("k2_z1656_boundary_base", BASE_SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load z1660 base")
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
ENGINE = BASE.ENGINE
EXACT_STATUS_AUDIT = BASE.BASE.BASE.exact_status_audit

require(
    tuple(sorted(BASE.EXPECTED_KEYS)) == tuple(sorted((*BOUNDARY_CASES, NEXT_PARENT_CASE))),
    "canonical z1656 boundary universe changed",
)
BASE_TEXT = BASE_OUTPUT.read_text(encoding="utf-8")
require(
    "atlas_band=1656..1660;candidate_rows=15015;"
    "height_counts=((1656, 5), (1660, 1));empty_firsts=1657..1659;"
    "next_occupied_height=1656" in BASE_TEXT
    and "consequence=projected k=2 first drift label z1<=1656" in BASE_TEXT,
    "canonical z1656 handoff changed",
)


def ordinary_profile(case):
    first, body = case
    carrier = ENGINE.U.suffix.A.carrier_for(body)
    require(ENGINE.P.A.carrier_for(body) == carrier, (case, "carrier engines disagree"))
    h = F(sum(right - left for left, right in carrier), ENGINE.U.suffix.A.RULER)
    lower = h * ENGINE.U.suffix.ETAS[2]
    L = 14 * lcm(*body)
    high_wall = ENGINE.U.suffix.PROJECTED_RATIOS[2] * L
    high_floor = max(15, high_wall.numerator // high_wall.denominator + 1)
    require(L == 11760 and first >= high_floor, (case, "ordinary boundary type changed"))
    result = ENGINE.ray_and_status(first, body, carrier, h, lower, L)
    counts = tuple(len(result[index]) for index in (6, 7, 8, 9))
    expected_counts, expected_stages, expected_projected = EXPECTED_ORDINARY[case]
    require(counts == expected_counts, (case, "status counts changed", counts))
    require(result[10] == expected_stages, (case, "canonical stage digests changed"))
    certificate_count = EXACT_STATUS_AUDIT(case, L, result[8], result[10][2])
    require(certificate_count == counts[2], (case, "exact certificate ledger changed"))
    require(counts[3] > 0, (case, "projection route became vacuous"))
    projected = ENGINE.projected_packets(
        first, body, carrier, h, lower, L, result[0], result[4], result[9]
    )
    require(projected[:8] == expected_projected, (case, "projected closure changed", projected[:8]))
    require(
        projected[1] == projected[2] > 0
        and projected[3] > 0
        and projected[6] == 1,
        (case, "projected consequence control failed", projected[:8]),
    )
    return (
        first,
        body,
        h,
        len(carrier),
        L,
        lower,
        result[4],
        result[5],
        result[1],
        result[2],
        result[3],
        counts,
        result[10],
        certificate_count,
        *projected[:8],
    )


def exceptional_high_profile():
    first, body = EXCEPTIONAL_CASE
    U = ENGINE.U
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * SCALAR_ETA
    L = 14 * lcm(*body)
    high_wall = HIGH_WALL_RATIO * L
    high_floor = max(15, high_wall.numerator // high_wall.denominator + 1)
    require(L == 152880 and first < high_floor, "exceptional high-wall type changed")

    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + L) * U.delta(carrier, h, residue + L) == amplitude,
            (EXCEPTIONAL_CASE, "ray recurrence", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(amplitudes[L - residue] == -amplitudes[residue] for residue in range(1, L))
        and L * U.delta(carrier, h, L) == 0
        and signs[1] == signs[-1]
        and sum(signs.values()) == L - 1,
        (EXCEPTIONAL_CASE, "ray symmetry changed", signs),
    )

    arbitrary = []
    high = []
    omitted = []
    for residue in range(1, L):
        amplitude = amplitudes[residue]
        if amplitude <= 0:
            continue
        first_label = residue
        if first_label <= first:
            first_label += ((first + 1 - first_label + L - 1) // L) * L
        for offset in range(4):
            label = first_label + offset * L
            arbitrary.append((amplitude / label, label, residue, offset))
        fifth_label = first_label + 4 * L
        omitted.append((amplitude / fifth_label, fifth_label, residue))
        high_label = residue
        if high_label < high_floor:
            high_label += ((high_floor - high_label + L - 1) // L) * L
        high.append((amplitude / high_label, high_label, residue))

    rank4 = tuple(sorted(arbitrary, key=lambda row: (-row[0], row[1:]))[:4])
    omitted_max = min(omitted, key=lambda row: (-row[0], row[1:]))
    best_high = min(high, key=lambda row: (-row[0], row[1:]))
    require(
        len(rank4) == 4
        and omitted_max[0] <= rank4[-1][0]
        and all(label < high_floor for _value, label, _residue, _offset in rank4),
        (EXCEPTIONAL_CASE, "top-ray truncation or branch changed"),
    )
    constrained = (*rank4[:3], (best_high[0], best_high[1], best_high[2], 0))
    require(
        len({row[1] for row in constrained}) == 4
        and any(row[1] >= high_floor for row in constrained),
        (EXCEPTIONAL_CASE, "forced-high maximizer malformed"),
    )
    first_delta = U.delta(carrier, h, first)
    unrestricted_upper = first_delta + sum((row[0] for row in rank4), F())
    unrestricted_gap = unrestricted_upper - lower
    upper = first_delta + sum((row[0] for row in constrained), F())
    gap = upper - lower
    require(
        unrestricted_gap > 0 and gap < 0,
        (EXCEPTIONAL_CASE, "sharp forced-high boundary changed", unrestricted_gap, gap),
    )
    for value, label, residue, _offset in constrained:
        require(
            U.delta(carrier, h, label) == value
            and label > first
            and label % L == residue % L,
            (EXCEPTIONAL_CASE, "chosen singleton control", label),
        )
    profile = (
        first,
        body,
        h,
        len(carrier),
        L,
        high_floor,
        first_delta,
        L // gcd(L, first),
        lower,
        signs[1],
        signs[-1],
        signs[0],
        sha256(repr(tuple(amplitudes)).encode()).hexdigest(),
        rank4,
        omitted_max,
        best_high,
        constrained,
        unrestricted_upper,
        unrestricted_gap,
        upper,
        gap,
    )
    expected = (
        FIRST,
        EXCEPTIONAL_CASE[1],
        F(811, 2548),
        38,
        152880,
        13250,
        F(6373, 9230130),
        6370,
        F(811, 231868),
        76434,
        76434,
        11,
        "83725ca3db416fdd7c1adf56b57a472128377ac328e2da4f8dee69b1f9a4f5ef",
        (
            (F(149, 160524), 2340, 2340, 0),
            (F(1877, 2233959), 2004, 2004, 0),
            (F(1081, 1364454), 1836, 1836, 0),
            (F(215, 276458), 1736, 1736, 0),
        ),
        (F(439, 110895330), 646620, 35100),
        (F(149, 909636), 13260, 13260),
        (
            (F(149, 160524), 2340, 2340, 0),
            (F(1877, 2233959), 2004, 2004, 0),
            (F(1081, 1364454), 1836, 1836, 0),
            (F(149, 909636), 13260, 13260, 0),
        ),
        None,
        None,
        F(44743106, 13102169535),
        F(-56372503, 681312815820),
    )
    require(profile[:-4] == expected[:-4], (EXCEPTIONAL_CASE, "exact profile changed"))
    require(profile[-2:] == expected[-2:], (EXCEPTIONAL_CASE, "forced-high gap changed"))
    return profile


def render(ordinary_profiles, high_profile):
    require(tuple((row[0], row[1]) for row in ordinary_profiles) == ORDINARY_CASES, "ordinary universe changed")
    scalar_total = sum(row[11][0] for row in ordinary_profiles)
    crude_total = sum(row[11][1] for row in ordinary_profiles)
    status_total = sum(row[11][2] for row in ordinary_profiles)
    residual_total = sum(row[11][3] for row in ordinary_profiles)
    certificate_total = sum(row[13] for row in ordinary_profiles)
    packet_total = sum(row[15] for row in ordinary_profiles)
    kill_total = sum(row[16] for row in ordinary_profiles)
    global_margin = min(row[17] for row in ordinary_profiles)
    require(
        (scalar_total, crude_total, status_total, residual_total, certificate_total, packet_total, kill_total)
        == (558, 91, 434, 33, 434, 36, 36),
        "global ordinary ledger changed",
    )
    require(global_margin == F(4390, 320229), "global projected margin changed")
    profile_payload = (ordinary_profiles, high_profile)
    profile_hash = sha256(repr(profile_payload).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "boundary profile digest changed")
    semantic_payload = (
        EXPECTED_BASE_SOURCE_SHA256,
        EXPECTED_BASE_OUTPUT_SHA256,
        BOUNDARY_CASES,
        QUANTIFIER,
        ALIGNED_TWO_UNION_CAP,
        HIGH_WALL_RATIO,
        scalar_total,
        crude_total,
        status_total,
        residual_total,
        certificate_total,
        packet_total,
        kill_total,
        global_margin,
        high_profile[-3],
        high_profile[-1],
        1655,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "boundary semantic digest changed")

    lines = [
        "LRC14 projected k=2 exact z1=1656 boundary closure",
        f"referee_source_sha256={file_sha256(HERE)}",
        f"base_source_sha256={file_sha256(BASE_SOURCE)}",
        f"base_output_sha256={file_sha256(BASE_OUTPUT)}",
        (
            "inherited_cap=1656;universe=exactly five scalar-atlas rows at z1=1656;"
            f"all {QUANTIFIER};no finite label horizon"
        ),
        "boundary_partition=ordinary_all_label:4;forced_high_scalar:1;empty:5;survivors:0",
        (
            f"ordinary_ledger=scalar:{scalar_total};crude:{crude_total};status:{status_total};"
            f"status_survivors:{residual_total};exact_status_certificates:{certificate_total};"
            f"literal_packets:{packet_total};projected_kills:{kill_total}"
        ),
        "status_replay=all rational Farkas certificates independently checked over Q;solver-selected dual bases excluded from digests",
    ]
    for row in ordinary_profiles:
        lines.append(
            f"CASE;route=ALL-LABEL-PROJECTED;z1={row[0]};E={','.join(map(str, row[1]))};"
            f"h={ftext(row[2])};r={row[3]};L={row[4]};lower={ftext(row[5])};"
            f"delta1={ftext(row[6])};first_d={row[7]};ray_sha256={row[8]};"
            f"denominators={row[9]};trials={row[10]};counts={row[11]};"
            f"stage_sha256={'/'.join(row[12])};exact_status_certificates={row[13]};"
            f"cells={row[14]};packets={row[15]};projected_kills={row[16]};"
            f"minimum_margin={ftext(row[17])};maximum_prefix_cells={row[18]};"
            f"minimum_row={row[19]};direct_mass={ftext(row[20])};state_sha256={row[21]};"
            "conclusion=EMPTY"
        )
    lines.append(
        f"CASE;route=FORCED-HIGH-SCALAR;z1={high_profile[0]};"
        f"E={','.join(map(str, high_profile[1]))};h={ftext(high_profile[2])};"
        f"r={high_profile[3]};L={high_profile[4]};high_floor={high_profile[5]};"
        f"delta1={ftext(high_profile[6])};first_d={high_profile[7]};"
        f"lower={ftext(high_profile[8])};ray_signs=+{high_profile[9]}/-{high_profile[10]}/0:{high_profile[11]};"
        f"ray_sha256={high_profile[12]};rank4={high_profile[13]};"
        f"first_omitted={high_profile[14]};best_high={high_profile[15]};"
        f"constrained={high_profile[16]};unrestricted_upper={ftext(high_profile[17])};"
        f"unrestricted_gap={ftext(high_profile[18])};upper={ftext(high_profile[19])};"
        f"gap={ftext(high_profile[20])};conclusion=SCALAR-EMPTY"
    )
    lines.extend(
        (
            "hostile_controls=all 33 ordinary status survivors require projected closure;exceptional unrestricted top four has positive gap but exact forced-high maximum has negative gap",
            "consequence=projected k=2 first drift label z1<=1655",
            "scope_boundary=finite-exact THM-2980 input only;full 1600..1679 band and LRC(14) remain open",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(4, mp.cpu_count() or 1))
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    os.environ["PYTHONHASHSEED"] = "0"
    high_profile = exceptional_high_profile()
    if args.workers == 1:
        ordinary_profiles = tuple(ordinary_profile(case) for case in ORDINARY_CASES)
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            ordinary_profiles = tuple(pool.map(ordinary_profile, ORDINARY_CASES))
    payload = render(ordinary_profiles, high_profile)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
