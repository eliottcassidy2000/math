#!/usr/bin/env python3
"""Exact THM-2941 scalar band for projected k=2, 1790 <= z1 <= 1799.

It reuses the pinned 1800..1809 wrapper, then changes only the first-label
interval before evaluating all 3003 bodies.
The inherited engine retains exact integration through label 7000 and its
rigorous omitted-tail upper bound.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_k2_scalar_band_1800_1809_thm2941.py"
START = 1790
END = 1799
EXPECTED_CANDIDATE_ROWS = 30030
EXPECTED_KEYS = (
    ((1, 4, 8, 10, 12, 14), 1790),
    ((1, 6, 8, 10, 12, 14), 1790),
    ((1, 8, 10, 12, 13, 14), 1790),
    ((2, 4, 8, 10, 12, 14), 1790),
    ((2, 6, 8, 9, 10, 14), 1790),
    ((2, 6, 8, 10, 12, 14), 1790),
    ((2, 8, 9, 10, 12, 14), 1790),
    ((4, 6, 8, 10, 12, 14), 1790),
)
EXPECTED_HIGH_TAIL_KEYS = (
    ((1, 8, 10, 12, 13, 14), 1790),
    ((2, 6, 8, 9, 10, 14), 1790),
    ((2, 8, 9, 10, 12, 14), 1790),
)
EXPECTED_PROFILE_SHA256 = "08ca435286a7e364e251a7467d7b9ed7d2dd577d8d6f4f5b6edbafb1a4caf1dd"
EXPECTED_SURVIVOR_SHA256 = "11ea44dbcba74f9fbdc5c35d1618cae876d977d62bd17f97c9a2558a27452a4a"
EXPECTED_SEMANTIC_SHA256 = "262260bc37789b023ddf4f3f0d46887707ea93265f9ae7f074c01da41ec2e833"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ftext(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


SPEC = spec_from_file_location("k2_band_1800_base_for_1790_scout", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load base band")
BAND = module_from_spec(SPEC)
SPEC.loader.exec_module(BAND)
BAND.D.B.START = START
BAND.D.B.END = END
BAND.D.B.LABELS = np.arange(START, BAND.D.B.HORIZON + 1, dtype=np.int64)


def profile(body):
    return BAND.D.B.profile(body)


def render(profiles) -> str:
    require(len(profiles) == 3003, "body universe changed")
    candidate_rows = sum(row[7] for row in profiles)
    survivors = tuple(
        sorted(
            (survivor for row in profiles for survivor in row[10]),
            key=lambda row: (row[0], row[1]),
        )
    )
    keys = tuple((row[0], row[1]) for row in survivors)
    high_tail_keys = tuple(
        (row[0], row[1])
        for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    semantic = (
        START,
        END,
        BAND.D.B.HORIZON,
        BAND.D.B.K,
        BAND.D.B.SUFFIX_SLOTS,
        BAND.D.B.PROJECTED_RATIO,
        candidate_rows,
        keys,
        high_tail_keys,
        survivor_hash,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic).encode()).hexdigest()
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate-row census changed")
    require(keys == EXPECTED_KEYS, "survivor key census changed")
    require(high_tail_keys == EXPECTED_HIGH_TAIL_KEYS, "high-tail key census changed")
    require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    require(survivor_hash == EXPECTED_SURVIVOR_SHA256, "survivor digest changed")
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 scalar band 1790..1799",
        f"candidate_rows={candidate_rows};survivors={len(survivors)};high_tail={len(high_tail_keys)}",
        f"keys={keys}",
        f"high_tail_keys={high_tail_keys}",
    ]
    for row in survivors:
        body, first, first_delta, chosen, upper, lower, gap, h, components, L, high_floor = row
        lines.append(
            f"SURVIVOR;E={','.join(map(str, body))};h={ftext(h)};r={components};L={L};"
            f"largest_floor={high_floor};z1={first};delta1={ftext(first_delta)};suffix="
            + ",".join(f"{kind}:{label}:{ftext(value)}" for value, label, kind in chosen)
            + f";lower={ftext(lower)};upper={ftext(upper)};gap={ftext(gap)}"
        )
    lines.extend(
        (
            f"profile_sha256={profile_hash}",
            f"survivor_sha256={survivor_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    bodies = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        profiles = [profile(body) for body in bodies]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            profiles = list(pool.imap(profile, bodies, chunksize=4))
    profiles.sort(key=lambda row: row[0])
    output = render(profiles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
