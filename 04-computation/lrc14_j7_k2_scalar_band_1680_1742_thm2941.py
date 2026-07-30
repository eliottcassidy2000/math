#!/usr/bin/env python3
"""Global projected k=2 scalar atlas on 1680 <= z1 <= 1742.

This is the next exact all-body descent band below the inherited closure
through z1=1743.  The guarded parent engine is rerun on all C(14,6)=3003
bodies.  Singleton excesses through label 7000 are integrated exactly;
each omitted later nonaligned slot is bounded by 6r/(49z).  The result is
therefore uniform over all distinct later nonaligned labels, not a finite
label search.

Exactly 58 scalar rows survive.  The highest occupied height is z1=1736,
with fifteen bodies, while every height z1=1737..1742 is scalar-empty.
The survivors split into 39 exact-suffix rows and 19 explicit HIGH-TAIL
rows for separate exact-ray treatment.
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
    / "lrc14_j7_k2_scalar_band_1750_1799_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
)
EXPECTED_PARENT_SHA256 = (
    "9ad7c58575b25c79b41bcf4226710c934720bf3e1394abf573a647eec5af87ff"
)
EXPECTED_PROFILE_SHA256 = (
    "3726e31959c7c3b1ba9d96aafd906fdfed2d613adec49eb963afa074e39c2de7"
)
EXPECTED_SURVIVOR_SHA256 = (
    "eddc5446fa044d40d8471da3ed682cb9797b4b183452e80a1ca717cab997f068"
)
EXPECTED_SEMANTIC_SHA256 = (
    "6d2e875e16397448840bf65073c2d203d13effb134830f60bdf7a1c454dba770"
)

QUANTIFIER = "ordered distinct later nonaligned labels"
START = 1680
END = 1742
EXPECTED_CANDIDATE_ROWS = 189_105
EXPECTED_HEIGHTS = (
    (1683, 1),
    (1694, 10),
    (1702, 3),
    (1708, 14),
    (1722, 11),
    (1724, 2),
    (1732, 2),
    (1736, 15),
)
EXPECTED_HIGH_KEYS = (
    (1683, (1, 2, 10, 11, 12, 14)),
    (1694, (1, 2, 10, 12, 13, 14)),
    (1694, (1, 4, 10, 12, 13, 14)),
    (1694, (2, 8, 9, 10, 12, 14)),
    (1708, (1, 2, 10, 11, 12, 14)),
    (1708, (1, 4, 10, 11, 12, 14)),
    (1708, (2, 8, 9, 10, 12, 14)),
    (1708, (2, 8, 10, 11, 12, 14)),
    (1708, (4, 8, 10, 11, 12, 14)),
    (1722, (1, 2, 10, 11, 12, 14)),
    (1722, (1, 6, 10, 11, 12, 14)),
    (1722, (2, 3, 10, 11, 12, 14)),
    (1722, (2, 8, 9, 10, 12, 14)),
    (1736, (1, 2, 10, 11, 12, 14)),
    (1736, (1, 2, 10, 12, 13, 14)),
    (1736, (1, 4, 10, 12, 13, 14)),
    (1736, (1, 8, 10, 12, 13, 14)),
    (1736, (2, 3, 10, 11, 12, 14)),
    (1736, (2, 8, 9, 10, 12, 14)),
)
EMPTY_FIRST_BANDS = (
    "1680..1682,1684..1693,1695..1701,1703..1707,1709..1721,1723,"
    "1725..1731,1733..1735,1737..1742"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(file_sha256(PARENT) == EXPECTED_PARENT_SHA256, "parent atlas changed")
SPEC = spec_from_file_location("k2_descent_1680_parent", PARENT)
require(SPEC is not None and SPEC.loader is not None, "cannot load parent atlas")
P = module_from_spec(SPEC)
SPEC.loader.exec_module(P)

P.P.D.B.START = START
P.P.D.B.END = END
P.P.D.B.LABELS = np.arange(START, P.P.D.B.HORIZON + 1, dtype=np.int64)


def profile(body):
    return P.P.D.B.profile(body)


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
    require(len(survivors) == 58, "survivor count changed")
    require(height_counts == EXPECTED_HEIGHTS, "surviving heights changed")
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
        P.P.D.B.HORIZON,
        P.P.D.B.K,
        P.P.D.B.SUFFIX_SLOTS,
        P.P.D.B.PROJECTED_RATIO,
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
        "LRC14 global projected k=2 scalar band 1680..1742",
        f"parent_atlas_source_sha256={file_sha256(PARENT)}",
        (
            "universe=six_body_roots:3003;body_labels:1..14;"
            f"{QUANTIFIER};scalar necessary condition"
        ),
        (
            f"first_band={START}..{END};first_candidate_rows={candidate_rows};"
            f"exact_horizon={P.P.D.B.HORIZON};omitted_slot_bound=6r/[49(H+1)]"
        ),
        (
            "projected_wall=max(15,floor(13L/150)+1);"
            "forced-high omitted slot starts max(H+1,projected_wall)"
        ),
        (
            f"global_survivors={len(survivors)};height_counts={height_counts};"
            f"ordinary_rows={len(survivors)-len(high_keys)};"
            f"high_tail_rows={len(high_keys)}:{high_keys}"
        ),
        f"empty_first_bands={EMPTY_FIRST_BANDS}",
        "consequence=all scalar-eligible projected k=2 rows with 1737<=z1<=1742 are empty",
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
