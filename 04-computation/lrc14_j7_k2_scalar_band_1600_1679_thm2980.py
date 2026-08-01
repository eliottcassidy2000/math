#!/usr/bin/env python3
"""Canonical exact projected-k2 scalar atlas on 1600 <= z1 <= 1679."""

from __future__ import annotations

import argparse
import multiprocessing as mp
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_k2_scalar_band_1800_1809_thm2941.py"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1600_1679_thm2980.out"
EXPECTED_BASE_SHA256 = "000b89e4e6745bcfece0f78b7fbd386300dd4162ac22d3da3eaa69aef16cd495"
EXPECTED_PROFILE_SHA256 = "8ef2eead4eb68eeff33e38152c09648032e5652886594c71214c54bd55e39d83"
EXPECTED_SURVIVOR_SHA256 = "b94748a736bc0069efb3b091b112527ad42a146b7eb549bd15f101658d2249d6"
EXPECTED_ROWS = 240_240
EXPECTED_SURVIVORS = 68
START = 1600
END = 1679


def ftext(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    args = parser.parse_args()
    if sha256(BASE.read_bytes().replace(b"\r\n", b"\n")).hexdigest() != EXPECTED_BASE_SHA256:
        raise RuntimeError("direct scalar base changed")

    spec = spec_from_file_location("k2_band_base_for_low_scout", BASE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load base band")
    band = module_from_spec(spec)
    spec.loader.exec_module(band)
    band.D.B.START = START
    band.D.B.END = END
    band.D.B.LABELS = np.arange(START, band.D.B.HORIZON + 1, dtype=np.int64)

    def profile(body):
        return band.D.B.profile(body)

    bodies = tuple(combinations(range(1, 15), 6))
    # A serial path keeps this scratch driver import-simple on Windows.  The
    # vectorized body engine is the expensive layer and remains exact.
    profiles = [profile(body) for body in bodies]
    survivors = tuple(
        sorted(
            (survivor for row in profiles for survivor in row[10]),
            key=lambda row: (row[1], row[0]),
        )
    )
    heights = {}
    for row in survivors:
        heights[row[1]] = heights.get(row[1], 0) + 1
    profile_digest = sha256(repr(tuple(profiles)).encode()).hexdigest()
    survivor_digest = sha256(repr(survivors).encode()).hexdigest()
    candidate_rows = sum(row[7] for row in profiles)
    if candidate_rows != EXPECTED_ROWS or len(survivors) != EXPECTED_SURVIVORS:
        raise RuntimeError((candidate_rows, len(survivors)))
    if profile_digest != EXPECTED_PROFILE_SHA256 or survivor_digest != EXPECTED_SURVIVOR_SHA256:
        raise RuntimeError((profile_digest, survivor_digest))
    lines = [
        f"LRC14 projected k=2 scalar scout {START}..{END}",
        f"direct_base_lf_sha256={EXPECTED_BASE_SHA256}",
        f"candidate_rows={sum(row[7] for row in profiles)};survivors={len(survivors)}",
        f"height_counts={tuple(sorted(heights.items()))}",
    ]
    for row in survivors:
        body, first, first_delta, chosen, upper, lower, gap, h, components, L, high_floor = row
        lines.append(
            f"SURVIVOR;E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};"
            f"largest_floor={high_floor};z1={first};delta1={ftext(first_delta)};suffix="
            + ",".join(
                f"{kind}:{label}:{ftext(value)}" for value, label, kind in chosen
            )
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
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
