#!/usr/bin/env python3
"""Global projected k=2 scalar atlas on 1750 <= z1 <= 1799.

This is the next exact all-body descent band below the proved empty interval
1800..1835.  The guarded parent engine is rerun on all C(14,6)=3003 bodies.
Singleton excesses through label 7000 are integrated exactly; each omitted
later nonaligned slot is bounded by 6r/(49z).  The result is therefore uniform
over all distinct later nonaligned labels, not an exhaustive finite-label
search.

Exactly 35 scalar rows survive.  The highest occupied height is z1=1790 with
eight bodies; z1=1791..1799 is empty.  Five of the eight frontier rows have
exact suffix maximizers and three carry the explicit HIGH-TAIL placeholder
for a separately certified forced-high ray maximum.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PARENT = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_scalar_band_1800_1809_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1750_1799_thm2941.out"
)
EXPECTED_PARENT_SHA256 = (
    "000b89e4e6745bcfece0f78b7fbd386300dd4162ac22d3da3eaa69aef16cd495"
)
EXPECTED_PROFILE_SHA256 = (
    "3bcb4c6e57756834f44ac6a4d0090a66e980ba1b1e24b1d47ecbd6d8764a014d"
)
EXPECTED_SURVIVOR_SHA256 = (
    "9b7742f3f306f3cec85436de20754dedb7b35e72cce9d231221440914a6a700c"
)
EXPECTED_SEMANTIC_SHA256 = (
    "0b1af1c4808d398529b779388e0b17590d7e4a45a779714763c7bedc462c166a"
)

QUANTIFIER = "ordered distinct later nonaligned labels"

START = 1750
END = 1799
EXPECTED_CANDIDATE_ROWS = 150_143
EXPECTED_HEIGHTS = (
    (1750, 12),
    (1758, 1),
    (1768, 1),
    (1776, 1),
    (1780, 5),
    (1784, 6),
    (1788, 1),
    (1790, 8),
)
EXPECTED_HIGH_KEYS = (
    (1750, (1, 2, 10, 11, 12, 14)),
    (1750, (1, 2, 10, 12, 13, 14)),
    (1750, (1, 4, 10, 12, 13, 14)),
    (1750, (1, 6, 10, 12, 13, 14)),
    (1768, (1, 8, 10, 12, 13, 14)),
    (1780, (1, 8, 10, 12, 13, 14)),
    (1784, (2, 8, 9, 10, 12, 14)),
    (1790, (1, 8, 10, 12, 13, 14)),
    (1790, (2, 6, 8, 9, 10, 14)),
    (1790, (2, 8, 9, 10, 12, 14)),
)
EMPTY_FIRST_BANDS = (
    "1751..1757,1759..1767,1769..1775,1777..1779,"
    "1781..1783,1785..1787,1789,1791..1799"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(PARENT) == EXPECTED_PARENT_SHA256, "1800..1809 parent changed")
SPEC = spec_from_file_location("k2_descent_1750_parent", PARENT)
require(SPEC is not None and SPEC.loader is not None, "cannot load parent atlas")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)

P.D.B.START = START
P.D.B.END = END
P.D.B.LABELS = np.arange(START, P.D.B.HORIZON + 1, dtype=np.int64)


def profile(body):
    return P.D.B.profile(body)


def render(profiles):
    require(len(profiles) == comb(14, 6) == 3003, "body universe changed")
    candidate_rows = sum(row[7] for row in profiles)
    survivors = tuple(
        sorted(
            (survivor for row in profiles for survivor in row[10]),
            key=lambda row: (row[1], row[0]),
        )
    )
    height_counts = tuple(sorted(Counter(row[1] for row in survivors).items()))
    high_keys = tuple(
        (row[1], row[0])
        for row in survivors
        if any(kind == "HIGH-TAIL" for _value, _label, kind in row[3])
    )
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate count changed")
    require(len(survivors) == 35, "survivor count changed")
    require(height_counts == EXPECTED_HEIGHTS, "surviving-height profile changed")
    require(high_keys == EXPECTED_HIGH_KEYS, "HIGH-TAIL routing changed")
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
    require(survivor_hash == EXPECTED_SURVIVOR_SHA256, "survivor digest changed")
    semantic_payload = (
        EXPECTED_PARENT_SHA256,
        QUANTIFIER,
        START,
        END,
        P.D.B.HORIZON,
        P.D.B.K,
        P.D.B.SUFFIX_SLOTS,
        P.D.B.PROJECTED_RATIO,
        candidate_rows,
        height_counts,
        high_keys,
        survivor_hash,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 global projected k=2 scalar band 1750..1799",
        f"parent_atlas_source_sha256={file_sha256(PARENT)}",
        (
            "universe=six_body_roots:3003;body_labels:1..14;"
            f"{QUANTIFIER};scalar necessary condition"
        ),
        (
            f"first_band={START}..{END};first_candidate_rows={candidate_rows};"
            f"exact_horizon={P.D.B.HORIZON};omitted_slot_bound=6r/[49(H+1)]"
        ),
        (
            "projected_wall=max(15,floor(13L/150)+1);"
            "forced-high omitted slot starts max(H+1,projected_wall)"
        ),
        (
            f"global_survivors={len(survivors)};height_counts={height_counts};"
            f"high_tail_rows={high_keys}"
        ),
        f"empty_first_bands={EMPTY_FIRST_BANDS}",
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
