#!/usr/bin/env python3
"""Exact forced-high ray closure for the six tail rows on 1800..1835.

For a fixed six-body carrier with ruler L, periodicity gives

    delta(r+mL) = A(r)/(r+mL),  1 <= r < L.

Thus every positive residue ray is strictly decreasing.  With four later
distinct nonaligned drift slots, only the first four labels of each positive
ray can enter the unrestricted top four; the best label meeting the forced
projected wall is the first wall-eligible label on some positive ray.  If the
unrestricted top four miss the wall, the exact constrained optimum is
therefore the unrestricted top three plus the best wall-eligible point.

The two scalar atlases leave exactly six such rows.  This verifier checks
the ray law and antipodes at every residue and proves the exact constrained
scalar maximum is strictly below h/91 in all six cases.  No label horizon
is used.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UNIFORM = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
BAND_1810 = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1810_1835_thm2941.out"
)
BAND_1800 = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1800_1809_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_high_wall_descent_1800_1810_closure_thm2941.out"
)
EXPECTED_UNIFORM_SHA256 = (
    "34ab29162ed33d90093e6d2bf781def36c420a1cd6596158b5d6579a3a8f3f46"
)
EXPECTED_BAND_1810_SHA256 = (
    "b29e1ccd9c3406c14bcfe2a53d5b6cee990c3d7a5e876bd5badcb27e2b506f0e"
)
EXPECTED_BAND_1800_SHA256 = (
    "ef9c9a2042c8d13d2f94acf8463fff1207bc6d03e8f1a3cb8589e8fb7d0db918"
)
EXPECTED_PROFILE_SHA256 = (
    "56f1cccd78211b8e1fa3640bb40c8f41ebefcc73ad21ce3b561bbee595921114"
)
EXPECTED_SEMANTIC_SHA256 = (
    "f456f0956986e6797d648fd3ba18579268a5102e121346650a4a08f7635d57e4"
)

QUANTIFIER = "distinct later nonaligned labels"

CASES = (
    (1810, (2, 6, 8, 9, 10, 14)),
    (1810, (2, 8, 9, 10, 12, 14)),
    (1807, (1, 2, 10, 12, 13, 14)),
    (1807, (1, 4, 10, 12, 13, 14)),
    (1807, (1, 8, 10, 12, 13, 14)),
    (1800, (1, 2, 10, 12, 13, 14)),
)
EXPECTED_GAPS = {
    CASES[0]: F(-58501985687, 2792778698922030),
    CASES[1]: F(-18923911, 6064666012860),
    CASES[2]: F(-29632076147, 106231479134508),
    CASES[3]: F(-71176511873, 265578697836270),
    CASES[4]: F(-17302995, 182999964056),
    CASES[5]: F(-7673254783, 38212762278600),
}
SUFFIX_SLOTS = 4
HIGH_WALL_RATIO = F(13, 150)
SCALAR_ETA = F(1, 91)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(UNIFORM) == EXPECTED_UNIFORM_SHA256, "uniform engine changed")
require(file_sha256(BAND_1810) == EXPECTED_BAND_1810_SHA256, "1810 atlas changed")
require(file_sha256(BAND_1800) == EXPECTED_BAND_1800_SHA256, "1800 atlas changed")
SPEC = spec_from_file_location("k2_high_wall_uniform", UNIFORM)
require(SPEC is not None and SPEC.loader is not None, "cannot load uniform engine")
U = module_from_spec(SPEC)
SPEC.loader.exec_module(U)


def profile(case):
    first, body = case
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * SCALAR_ETA
    L = 14 * lcm(*body)
    require(U.suffix.A.RULER % L == 0, (case, "body ruler left master ruler"))
    high_wall = HIGH_WALL_RATIO * L
    high_floor = max(15, high_wall.numerator // high_wall.denominator + 1)
    require(first < high_floor, (case, "case does not force a later high label"))

    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require(
            (residue + L) * U.delta(carrier, h, residue + L) == amplitude,
            (case, "ray recurrence", residue),
        )
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(
        all(amplitudes[L - residue] == -amplitudes[residue] for residue in range(1, L)),
        (case, "ray antipode"),
    )
    require(L * U.delta(carrier, h, L) == 0, (case, "aligned ray is nonzero"))
    require(
        signs[1] == signs[-1] and sum(signs.values()) == L - 1,
        (case, "ray sign census", signs),
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
        for offset in range(SUFFIX_SLOTS):
            label = first_label + offset * L
            arbitrary.append((amplitude / label, label, residue, offset))
        fifth_label = first_label + SUFFIX_SLOTS * L
        omitted.append((amplitude / fifth_label, fifth_label, residue))
        high_label = residue
        if high_label < high_floor:
            high_label += ((high_floor - high_label + L - 1) // L) * L
        high.append((amplitude / high_label, high_label, residue))

    rank4 = tuple(sorted(arbitrary, key=lambda row: (-row[0], row[1:]))[:4])
    require(len(rank4) == 4 and rank4[-1][0] > 0, (case, "positive top four missing"))
    omitted_max = min(omitted, key=lambda row: (-row[0], row[1:]))
    require(
        omitted_max[0] <= rank4[-1][0],
        (case, "first-four-per-ray truncation failed", omitted_max, rank4[-1]),
    )
    best_high = min(high, key=lambda row: (-row[0], row[1:]))
    if any(label >= high_floor for _value, label, _residue, _offset in rank4):
        constrained = rank4
        branch = "UNRESTRICTED_TOP4_ALREADY_HIGH"
    else:
        constrained = (*rank4[:3], (best_high[0], best_high[1], best_high[2], 0))
        branch = "GLOBAL_TOP3_PLUS_BEST_HIGH"
    require(
        len({row[1] for row in constrained}) == SUFFIX_SLOTS
        and any(row[1] >= high_floor for row in constrained),
        (case, "constrained suffix malformed", constrained),
    )
    first_delta = U.delta(carrier, h, first)
    upper = first_delta + sum((row[0] for row in constrained), F())
    gap = upper - lower
    require(gap < 0, (case, "forced-high scalar row survived", gap))
    expected_gap = EXPECTED_GAPS[case]
    if expected_gap is not None:
        require(gap == expected_gap, (case, "exact gap changed", gap))
    for value, label, residue, _offset in constrained:
        require(
            U.delta(carrier, h, label) == value
            and label > first
            and label % L == residue % L,
            (case, "chosen singleton control", label),
        )
    return (
        first,
        body,
        h,
        len(carrier),
        L,
        high_floor,
        first_delta,
        lower,
        signs[1],
        signs[-1],
        signs[0],
        sha256(repr(tuple(amplitudes)).encode()).hexdigest(),
        rank4,
        omitted_max,
        best_high,
        branch,
        tuple(constrained),
        upper,
        gap,
    )


def render(profiles):
    require(tuple((row[0], row[1]) for row in profiles) == CASES, "case universe changed")
    require(all(row[-1] < 0 for row in profiles), "a high-wall row survived")
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        CASES,
        QUANTIFIER,
        SUFFIX_SLOTS,
        HIGH_WALL_RATIO,
        SCALAR_ETA,
        tuple(row[-1] for row in profiles),
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 forced-high exact-ray descent closure 1800..1810",
        f"uniform_engine_sha256={file_sha256(UNIFORM)}",
        f"scalar_band_1810_sha256={file_sha256(BAND_1810)}",
        f"scalar_band_1800_sha256={file_sha256(BAND_1800)}",
        f"scope=six inherited HIGH-TAIL rows;all {QUANTIFIER};no label horizon",
        (
            "ray_law=delta(r+mL)=A(r)/(r+mL);search=first four points of every positive "
            "ray plus the first wall-eligible point"
        ),
        "forced_high_wall=floor(13L/150)+1;scalar_necessary=sum(delta)>=h/91",
    ]
    for row in profiles:
        (
            first,
            body,
            h,
            components,
            L,
            high_floor,
            first_delta,
            lower,
            positive,
            negative,
            zero,
            ray_digest,
            rank4,
            omitted_max,
            best_high,
            branch,
            constrained,
            upper,
            gap,
        ) = row
        lines.append(
            f"CASE;z1={first};E={','.join(map(str, body))};h={ftext(h)};"
            f"r={components};L={L};high_floor={high_floor};"
            f"delta1={ftext(first_delta)};lower={ftext(lower)};"
            f"ray_signs=+{positive}/-{negative}/0:{zero};ray_sha256={ray_digest};"
            f"unrestricted_top4={rank4};first_omitted={omitted_max};"
            f"best_high={best_high};branch={branch};constrained={constrained};"
            f"upper={ftext(upper)};gap={ftext(gap)};conclusion=SCALAR-EMPTY"
        )
    lines.extend(
        (
            "global_high_wall_rows=6;scalar_empty=6;survivors=0",
            f"profile_sha256={profile_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
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
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, CASES))
    order = {case: index for index, case in enumerate(CASES)}
    profiles.sort(key=lambda row: order[(row[0], row[1])])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
