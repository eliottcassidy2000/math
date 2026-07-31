#!/usr/bin/env python3
"""Global projected k=2 scalar atlas on 1810 <= z1 <= 1835.

This is the reusable descent step immediately below the proved z1=1836
slice.  It runs the pinned all-body scalar-envelope engine on all
C(14,6)=3003 six-body roots.  Labels through 7000 are integrated exactly;
each later slot beyond that horizon is bounded by 6r/(49z), so the scan is
uniform over all later labels rather than a finite-label experiment.

Exactly ten rows survive the scalar necessary condition: eight at z1=1810,
one at z1=1812, and one at z1=1824.  In particular z1=1811 and the open
bands 1813..1823 and 1825..1835 are empty.  Eight surviving rows have exact
suffix maximizers; the two remaining z1=1810 rows retain an explicit
HIGH-TAIL upper-bound placeholder and are routed to exact rays separately.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_scalar_band_1836_1836_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1810_1835_thm2941.out"
)
EXPECTED_BASE_SHA256 = (
    "4afbd91c68d6dbfd792f7cccede15bc0d0f5ae70ec6af30670edcc13466c4c3b"
)
EXPECTED_PROFILE_SHA256 = (
    "6492e0dce3d8289674e72aa8cb42fdd2f69e96efe4f956967ea04cec480f487f"
)
EXPECTED_SURVIVOR_SHA256 = (
    "f23f76bb3a656916ba9ef3dd063aef151dd8cea2059c7c7789327b68b0be8240"
)
EXPECTED_SEMANTIC_SHA256 = (
    "790e1a9fe294be6ca6d47204374dcc0f585f3d9532aec2d4b5d689e3df0e21da"
)

START = 1810
END = 1835
EXPECTED_CANDIDATE_ROWS = 78_078
EXPECTED_KEYS = (
    ((1, 2, 8, 10, 12, 14), 1810),
    ((1, 4, 8, 10, 12, 14), 1810),
    ((1, 4, 8, 10, 12, 14), 1812),
    ((1, 4, 8, 10, 12, 14), 1824),
    ((1, 6, 8, 10, 12, 14), 1810),
    ((2, 4, 8, 10, 12, 14), 1810),
    ((2, 6, 8, 9, 10, 14), 1810),
    ((2, 6, 8, 10, 12, 14), 1810),
    ((2, 8, 9, 10, 12, 14), 1810),
    ((4, 6, 8, 10, 12, 14), 1810),
)
EXPECTED_HIGH_TAIL_KEYS = (
    ((2, 6, 8, 9, 10, 14), 1810),
    ((2, 8, 9, 10, 12, 14), 1810),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(BASE) == EXPECTED_BASE_SHA256, "scalar-atlas base changed")
SPEC = spec_from_file_location("k2_descent_scalar_base", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load scalar-atlas base")
B = module_from_spec(SPEC)
SPEC.loader.exec_module(B)

# The base engine is parameterized by these module globals.  Reset all three
# together before any worker is spawned; a spawn worker imports this module
# and performs the same deterministic reset.
B.START = START
B.END = END
B.LABELS = np.arange(START, B.HORIZON + 1, dtype=np.int64)


def profile(body):
    return B.profile(body)


def render(profiles):
    require(len(profiles) == comb(14, 6) == 3003, "body universe changed")
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
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate-row count changed")
    require(keys == EXPECTED_KEYS, ("survivor keys changed", keys))
    require(high_tail_keys == EXPECTED_HIGH_TAIL_KEYS, "HIGH-TAIL routing changed")
    require(
        all(
            all(
                (label is not None and kind == "EXACT")
                or (label is None and kind in ("TAIL-0", "HIGH-TAIL"))
                for _value, label, kind in row[3]
            )
            for row in survivors
        ),
        "a suffix certificate is malformed",
    )
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    survivor_hash = sha256(repr(survivors).encode()).hexdigest()
    require(profile_hash == EXPECTED_PROFILE_SHA256, "global profile digest changed")
    if EXPECTED_SURVIVOR_SHA256 is not None:
        require(survivor_hash == EXPECTED_SURVIVOR_SHA256, "survivor digest changed")
    semantic_payload = (
        EXPECTED_BASE_SHA256,
        START,
        END,
        B.HORIZON,
        B.K,
        B.SUFFIX_SLOTS,
        B.PROJECTED_RATIO,
        candidate_rows,
        keys,
        high_tail_keys,
        survivor_hash,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 global projected k=2 scalar band 1810..1835",
        f"scalar_atlas_base_sha256={file_sha256(BASE)}",
        (
            "universe=six_body_roots:3003;body_labels:1..14;"
            "ordered distinct later nonaligned labels;scalar necessary condition"
        ),
        (
            f"first_band={START}..{END};first_candidate_rows={candidate_rows};"
            f"exact_horizon={B.HORIZON};omitted_slot_bound=6r/[49(H+1)]"
        ),
        (
            "projected_wall=max(15,floor(13L/150)+1);"
            "forced-high omitted slot starts max(H+1,projected_wall)"
        ),
        (
            f"global_survivors={len(survivors)};surviving_firsts="
            f"{tuple(row[1] for row in survivors)};high_tail_rows={high_tail_keys}"
        ),
        "empty_firsts=1811;empty_first_bands=1813..1823,1825..1835",
    ]
    for row in survivors:
        (
            body,
            first,
            first_delta,
            chosen,
            upper,
            lower,
            gap,
            h,
            components,
            L,
            high_floor,
        ) = row
        lines.append(
            f"SURVIVOR;E={','.join(map(str, body))};h={ftext(h)};"
            f"r={components};L={L};largest_floor={high_floor};"
            f"z1={first};delta1={ftext(first_delta)};suffix="
            + ",".join(
                f"{kind}:{label}:{ftext(value)}" for value, label, kind in chosen
            )
            + f";lower={ftext(lower)};upper={ftext(upper)};gap={ftext(gap)}"
        )
    lines.extend(
        (
            f"global_profile_sha256={profile_hash}",
            f"survivor_sha256={survivor_hash}",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
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
