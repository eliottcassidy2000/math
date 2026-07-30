#!/usr/bin/env python3
"""Exact THM-2941 scalar band for projected k=2, 1770 <= z1 <= 1779."""

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
START, END = 1770, 1779
EXPECTED_CANDIDATE_ROWS = 30030
EXPECTED_KEYS = (((1, 4, 8, 10, 12, 14), 1776),)
EXPECTED_PROFILE_SHA256 = "7de16804fe4a588d9b0c722dc53be9c154f02b517cc46abc8e80d343eb62b1d0"
EXPECTED_SURVIVOR_SHA256 = "bf328ce173825df510e8b7c63957b8182a5d094edcf82b3c3d1fd035f9e3e16d"
EXPECTED_SEMANTIC_SHA256 = "31c260435a37ee8ac0f472690f93cd593a86e22cb386b765ee4b717e6d6537ab"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ftext(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


SPEC = spec_from_file_location("k2_band_1800_base_for_1770", BASE)
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
    survivors = tuple(sorted(
        (survivor for row in profiles for survivor in row[10]),
        key=lambda row: (row[0], row[1]),
    ))
    keys = tuple((row[0], row[1]) for row in survivors)
    high_tail_keys = tuple(
        (row[0], row[1]) for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    semantic = (
        START, END, BAND.D.B.HORIZON, BAND.D.B.K, BAND.D.B.SUFFIX_SLOTS,
        BAND.D.B.PROJECTED_RATIO, candidate_rows, keys, high_tail_keys,
        survivor_hash, profile_hash,
    )
    semantic_hash = sha256(repr(semantic).encode()).hexdigest()
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate-row census changed")
    require(keys == EXPECTED_KEYS, "survivor key census changed")
    require(high_tail_keys == (), "unexpected high-tail survivor")
    require(profile_hash == EXPECTED_PROFILE_SHA256, "profile digest changed")
    require(survivor_hash == EXPECTED_SURVIVOR_SHA256, "survivor digest changed")
    require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 scalar band 1770..1779",
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
    lines.extend((
        f"profile_sha256={profile_hash}",
        f"survivor_sha256={survivor_hash}",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ))
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
