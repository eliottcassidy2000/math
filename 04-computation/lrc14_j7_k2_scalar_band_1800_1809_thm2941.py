#!/usr/bin/env python3
"""Global projected k=2 scalar atlas on 1800 <= z1 <= 1809.

The pinned all-body atlas is rerun on all C(14,6)=3003 six-body roots,
with exact singleton integration through label 7000 and the rigorous
6r/(49z) bound on every omitted later slot.  Thus its six survivors are an
all-label upper-envelope result, not a finite search.

The only survivors are the hostile L=11760 body at z1=1800 and 1805, and
four explicitly routed HIGH-TAIL rows (one at z1=1800 and three at z1=1807).
All other first labels in this band are scalar-empty.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DESCENT_BASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k2_scalar_band_1810_1835_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1800_1809_thm2941.out"
)
EXPECTED_DESCENT_BASE_SHA256 = (
    "a09b13e994ad6ab35e3324bd336e773b435b07859e6a3c924b84ec77f3e2aced"
)
EXPECTED_PROFILE_SHA256 = (
    "30c39d884aee43f34396c6f9d2ee8b1bfb034bddc3e235311555ccb55b96d67e"
)
EXPECTED_SURVIVOR_SHA256 = (
    "b1b7263821e64795ebaac43b9c1f310b3e2ec6337d947b7c37e3b76d172dcc20"
)
EXPECTED_SEMANTIC_SHA256 = (
    "babab337409ef839d6a20caa495a645969b8fc5bc5c4e05fb3de18486f7f9fbc"
)

START = 1800
END = 1809
EXPECTED_CANDIDATE_ROWS = 30_030
EXPECTED_KEYS = (
    ((1, 2, 10, 12, 13, 14), 1800),
    ((1, 2, 10, 12, 13, 14), 1807),
    ((1, 4, 8, 10, 12, 14), 1800),
    ((1, 4, 8, 10, 12, 14), 1805),
    ((1, 4, 10, 12, 13, 14), 1807),
    ((1, 8, 10, 12, 13, 14), 1807),
)
EXPECTED_HIGH_TAIL_KEYS = (
    ((1, 2, 10, 12, 13, 14), 1800),
    ((1, 2, 10, 12, 13, 14), 1807),
    ((1, 4, 10, 12, 13, 14), 1807),
    ((1, 8, 10, 12, 13, 14), 1807),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(
    file_sha256(DESCENT_BASE) == EXPECTED_DESCENT_BASE_SHA256,
    "1810..1835 descent-atlas dependency changed",
)
SPEC = spec_from_file_location("k2_descent_1810_base", DESCENT_BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load descent atlas")
D = module_from_spec(SPEC)
SPEC.loader.exec_module(D)

D.B.START = START
D.B.END = END
D.B.LABELS = np.arange(START, D.B.HORIZON + 1, dtype=np.int64)


def profile(body):
    return D.B.profile(body)


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
    require(candidate_rows == EXPECTED_CANDIDATE_ROWS, "candidate count changed")
    require(keys == EXPECTED_KEYS, ("survivor keys changed", keys))
    require(high_tail_keys == EXPECTED_HIGH_TAIL_KEYS, "HIGH-TAIL routes changed")
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
        EXPECTED_DESCENT_BASE_SHA256,
        START,
        END,
        D.B.HORIZON,
        D.B.K,
        D.B.SUFFIX_SLOTS,
        D.B.PROJECTED_RATIO,
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
        "LRC14 global projected k=2 scalar band 1800..1809",
        f"descent_atlas_source_sha256={file_sha256(DESCENT_BASE)}",
        (
            "universe=six_body_roots:3003;body_labels:1..14;"
            "ordered distinct later nonaligned labels;scalar necessary condition"
        ),
        (
            f"first_band={START}..{END};first_candidate_rows={candidate_rows};"
            f"exact_horizon={D.B.HORIZON};omitted_slot_bound=6r/[49(H+1)]"
        ),
        (
            "projected_wall=max(15,floor(13L/150)+1);"
            "forced-high omitted slot starts max(H+1,projected_wall)"
        ),
        (
            f"global_survivors={len(survivors)};surviving_firsts="
            f"{tuple(row[1] for row in survivors)};high_tail_rows={high_tail_keys}"
        ),
        "empty_firsts=1801,1802,1803,1804,1806,1808,1809",
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
