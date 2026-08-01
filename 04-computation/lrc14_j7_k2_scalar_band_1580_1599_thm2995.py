#!/usr/bin/env python3
"""Exact projected-k2 scalar atlas on 1580 <= z1 <= 1599."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1800_1809_thm2941.py"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1580_1599_thm2995.out"
EXPECTED_BASE_SHA256 = "000b89e4e6745bcfece0f78b7fbd386300dd4162ac22d3da3eaa69aef16cd495"
EXPECTED_PROFILE_SHA256 = "0f02d797f9a3bc39871d4a4ccf72c42525ca3ca4109773469b8fe95b780780d3"
EXPECTED_SURVIVOR_SHA256 = "aab25b18f0bb9e2cb33eaa378176e720e1375c19891be9c97c67371f72ae2636"
EXPECTED_ROWS = 60_060
EXPECTED_SURVIVORS = 26
START = 1580
END = 1599


def ftext(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    base_hash = sha256(BASE.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")).hexdigest()
    if base_hash != EXPECTED_BASE_SHA256:
        raise RuntimeError(("direct scalar base changed", base_hash))

    spec = spec_from_file_location("k2_band_base_for_1580_1599", BASE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load base band")
    band = module_from_spec(spec)
    spec.loader.exec_module(band)
    band.D.B.START = START
    band.D.B.END = END
    band.D.B.LABELS = np.arange(START, band.D.B.HORIZON + 1, dtype=np.int64)

    bodies = tuple(combinations(range(1, 15), 6))
    profiles = [band.D.B.profile(body) for body in bodies]
    survivors = tuple(
        sorted(
            (survivor for profile in profiles for survivor in profile[10]),
            key=lambda row: (row[1], row[0]),
        )
    )
    heights = {}
    for row in survivors:
        heights[row[1]] = heights.get(row[1], 0) + 1
    profile_digest = sha256(repr(tuple(profiles)).encode()).hexdigest()
    survivor_digest = sha256(repr(survivors).encode()).hexdigest()
    candidate_rows = sum(profile[7] for profile in profiles)
    if candidate_rows != EXPECTED_ROWS or len(survivors) != EXPECTED_SURVIVORS:
        raise RuntimeError((candidate_rows, len(survivors)))
    if profile_digest != EXPECTED_PROFILE_SHA256:
        raise RuntimeError(("profile digest", profile_digest))
    if survivor_digest != EXPECTED_SURVIVOR_SHA256:
        raise RuntimeError(("survivor digest", survivor_digest))

    lines = [
        f"LRC14 projected k=2 scalar scout {START}..{END}",
        f"direct_base_lf_sha256={EXPECTED_BASE_SHA256}",
        f"candidate_rows={candidate_rows};survivors={len(survivors)}",
        f"height_counts={tuple(sorted(heights.items()))}",
    ]
    for row in survivors:
        body, first, first_delta, chosen, upper, lower, gap, h, components, modulus, high_floor = row
        lines.append(
            f"SURVIVOR;E={','.join(map(str, body))};h={ftext(h)};r={components};L={modulus};"
            f"largest_floor={high_floor};z1={first};delta1={ftext(first_delta)};suffix="
            + ",".join(f"{kind}:{label}:{ftext(value)}" for value, label, kind in chosen)
            + f";lower={ftext(lower)};upper={ftext(upper)};gap={ftext(gap)}"
        )
    lines.extend(
        (
            f"profile_sha256={profile_digest}",
            f"survivor_sha256={survivor_digest}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
