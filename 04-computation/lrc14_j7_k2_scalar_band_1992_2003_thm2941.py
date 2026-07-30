#!/usr/bin/env python3
"""Global one-pass projected k=2 scalar atlas on 1992 <= z1 <= 2003.

The universe is all C(14,6)=3003 body roots.  For every admissible first
nonaligned drift in the band, compute the rigorous top-four suffix envelope
over all distinct later labels.  Labels through 7000 are integrated exactly;
each omitted slot uses the proved delta(z)<=6r/(49z) bound.  When the
projected-safe-arc wall forces a large final drift, its exact/tail candidate
is selected in the same single-high-slot envelope as the parent THM-2941
scan.  Thus no finite horizon is treated as exhaustive.

Only one scalar necessary row survives globally: the hostile body
E=(1,4,8,10,12,14) at z1=1992.  Its chosen suffix is exact
(no tail placeholders).  Transport closure is a separate step.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
UNIFORM = (
    ROOT
    / "04-computation"
    / "lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k2_scalar_band_1992_2003_thm2941.out"
)
EXPECTED_UNIFORM_SHA256 = (
    "dfa4788297b8c31fc9b5dce1afadf29d20b267cb4159fa95dadb9346b1980b36"
)
EXPECTED_PROFILE_SHA256 = (
    "ba117ecaaf23dfa859f39d0f414c1391272e90e35c0754a8f3a8f28f60929f6e"
)
EXPECTED_SEMANTIC_SHA256 = (
    "3c434cf63235754627f9c288d78cbb61d6dd3b54425ceeb4e6e09ef65e82587e"
)
HOSTILE_BODY = (1, 4, 8, 10, 12, 14)
START = 1992
END = 2003
HORIZON = 7000
K = 2
SUFFIX_SLOTS = 4
PROJECTED_RATIO = F(13, 150)
EXPECTED_SURVIVORS = (
    (
        HOSTILE_BODY,
        1992,
        F(121, 170814),
        (
            (F(1867, 2119740), 2060, "EXACT"),
            (F(263, 310415), 2172, "EXACT"),
            (F(821, 1049580), 2142, "EXACT"),
            (F(2911, 3724980), 2534, "EXACT"),
        ),
        F(6496500403, 1624087555020),
        F(1049, 267540),
        F(417952777, 5278284553815),
        F(1049, 2940),
        34,
        11760,
        1020,
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


require(
    file_sha256(UNIFORM) == EXPECTED_UNIFORM_SHA256,
    "uniform ray/status dependency changed",
)
spec = spec_from_file_location("k2_global_band_uniform_support", UNIFORM)
require(spec is not None and spec.loader is not None, "cannot load uniform")
U = module_from_spec(spec)
spec.loader.exec_module(U)


LABELS = np.arange(START, HORIZON + 1, dtype=np.int64)
RULER = U.suffix.A.RULER
ONE_FOURTEENTH = RULER // 14
ONE_SEVENTH = RULER // 7
THIRTEEN_FOURTEENTHS = 13 * RULER // 14


def primitive_vector(numerators):
    cycles = numerators // RULER
    remainders = numerators % RULER
    return (
        cycles * ONE_SEVENTH
        + np.minimum(remainders, ONE_FOURTEENTH)
        + np.maximum(0, remainders - THIRTEEN_FOURTEENTHS)
    )


def exact_top_indices(amplitudes, eligible, limit):
    """Float-discover, then exact-expand and exact-sort the requested head."""
    indices = np.flatnonzero(eligible)
    require(len(indices) >= limit, "exact head underfilled")
    scores = amplitudes[indices].astype(np.float64) / LABELS[indices]
    seed_local = np.argpartition(scores, -limit)[-limit:]
    seed = indices[seed_local]
    threshold = min(
        seed,
        key=lambda index: F(int(amplitudes[index]), int(LABELS[index])),
    )
    threshold_a = int(amplitudes[threshold])
    threshold_z = int(LABELS[threshold])
    comparisons = (
        amplitudes[indices] * threshold_z
        - threshold_a * LABELS[indices]
    )
    expanded = indices[comparisons >= 0]
    ordered = sorted(
        (int(index) for index in expanded),
        key=lambda index: (
            -F(int(amplitudes[index]), int(LABELS[index])),
            int(LABELS[index]),
        ),
    )
    answer = tuple(ordered[:limit])
    floor_index = answer[-1]
    floor_a = int(amplitudes[floor_index])
    floor_z = int(LABELS[floor_index])
    answer_set = set(answer)
    outside = np.asarray(
        [index for index in indices if int(index) not in answer_set],
        dtype=np.int64,
    )
    require(
        np.all(
            amplitudes[outside] * floor_z
            <= floor_a * LABELS[outside]
        ),
        "float discovery omitted an exact head value",
    )
    return answer


def top_insert(rows, item, limit=SUFFIX_SLOTS):
    rows.append(item)
    rows.sort(
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        )
    )
    del rows[limit:]


def suffix_upper(arbitrary_top, high_top, tail, high_tail, constrained):
    arbitrary = list(arbitrary_top)
    arbitrary.extend(
        (tail, None, f"TAIL-{index}") for index in range(SUFFIX_SLOTS)
    )
    arbitrary.sort(
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        )
    )
    if not constrained:
        return tuple(arbitrary[:SUFFIX_SLOTS])

    high = list(high_top)
    high.append((high_tail, None, "HIGH-TAIL"))
    high.sort(
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        )
    )
    selected_high = high[0]
    rest = [
        row
        for row in arbitrary
        if selected_high[1] is None or row[1] != selected_high[1]
    ]
    answer = (selected_high, *rest[: SUFFIX_SLOTS - 1])
    require(len(answer) == SUFFIX_SLOTS, "constrained suffix underfilled")
    return answer


def profile(body):
    carrier = U.suffix.A.carrier_for(body)
    mass_numerator = sum(right - left for left, right in carrier)
    components = len(carrier)
    h = F(mass_numerator, RULER)
    lower = h * U.suffix.ETAS[K]
    L = 14 * lcm(*body)
    wall = PROJECTED_RATIO * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)

    coverage_numerators = np.zeros(len(LABELS), dtype=np.int64)
    for left, right in carrier:
        coverage_numerators += (
            primitive_vector(LABELS * right)
            - primitive_vector(LABELS * left)
        )
    amplitudes = 7 * coverage_numerators - mass_numerator * LABELS
    max_product = int(np.max(np.abs(amplitudes))) * HORIZON
    require(max_product < 2**62, ("exact comparison overflow risk", body))
    control_labels = (
        START,
        END,
        HORIZON,
        START + sum(body) % (END - START + 1),
    )
    require(
        all(
            F(
                int(amplitudes[label - START]),
                7 * RULER * label,
            )
            == U.delta(carrier, h, label)
            for label in control_labels
        ),
        ("vector/scalar singleton control failed", body, control_labels),
    )

    admissible = LABELS % L != 0
    base_mask = admissible & (LABELS > END)
    base_indices = exact_top_indices(
        amplitudes, base_mask, SUFFIX_SLOTS
    )
    arbitrary_top = [
        (
            F(int(amplitudes[index]), 7 * RULER * int(LABELS[index])),
            int(LABELS[index]),
            "EXACT",
        )
        for index in base_indices
    ]
    arbitrary_top.sort(key=lambda row: (-row[0], row[1], row[2]))

    high_top = []
    if high_floor <= HORIZON:
        high_mask = admissible & (LABELS >= high_floor)
        high_index = exact_top_indices(amplitudes, high_mask, 1)[0]
        high_top.append(
            (
                F(
                    int(amplitudes[high_index]),
                    7 * RULER * int(LABELS[high_index]),
                ),
                int(LABELS[high_index]),
                "EXACT",
            )
        )

    tail = F(6 * components, 49 * (HORIZON + 1))
    high_tail = F(
        6 * components,
        49 * max(HORIZON + 1, high_floor),
    )
    survivors = []
    candidate_count = 0
    best_gap = None
    best_first = None
    for first in range(END, START - 1, -1):
        index = first - START
        if first % L != 0:
            candidate_count += 1
            first_delta = F(
                int(amplitudes[index]),
                7 * RULER * first,
            )
            chosen = suffix_upper(
                arbitrary_top,
                high_top,
                tail,
                high_tail,
                first < high_floor,
            )
            upper = first_delta + sum(
                (value for value, _label, _kind in chosen), F(0)
            )
            gap = upper - lower
            if best_gap is None or gap > best_gap:
                best_gap = gap
                best_first = first
            if gap >= 0:
                survivors.append(
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
                    )
                )
        if first % L != 0:
            top_insert(
                arbitrary_top,
                (
                    F(
                        int(amplitudes[index]),
                        7 * RULER * first,
                    ),
                    first,
                    "EXACT",
                ),
            )

    require(best_gap is not None and best_first is not None, "empty band")
    return (
        body,
        h,
        components,
        L,
        high_floor,
        tail,
        high_tail,
        candidate_count,
        best_first,
        best_gap,
        tuple(survivors),
    )


def render(profiles):
    require(len(profiles) == comb(14, 6) == 3003, "body universe changed")
    candidate_rows = sum(row[7] for row in profiles)
    survivors = tuple(
        sorted(
            (survivor for row in profiles for survivor in row[10]),
            key=lambda row: (row[0], row[1]),
        )
    )
    require(survivors == EXPECTED_SURVIVORS, ("global band changed", survivors))
    require(
        all(
            all(label is not None and kind == "EXACT" for _value, label, kind in row[3])
            for row in survivors
        ),
        "a surviving suffix uses an analytic tail placeholder",
    )
    profile_hash = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(
            profile_hash == EXPECTED_PROFILE_SHA256,
            "global profile digest changed",
        )
    semantic_payload = (
        START,
        END,
        HORIZON,
        K,
        SUFFIX_SLOTS,
        PROJECTED_RATIO,
        candidate_rows,
        survivors,
        profile_hash,
    )
    semantic_hash = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            "global scalar-band semantic digest changed",
        )

    lines = [
        "LRC14 global projected k=2 scalar band 1992..2003",
        f"uniform_ray_status_source_sha256={file_sha256(UNIFORM)}",
        (
            "universe=six_body_roots:3003;body_labels:1..14;"
            "ordered distinct later nonaligned labels;scalar necessary condition"
        ),
        (
            f"first_band={START}..{END};first_candidate_rows={candidate_rows};"
            f"exact_horizon={HORIZON};omitted_slot_bound=6r/[49(H+1)]"
        ),
        (
            "projected_wall=max(15,floor(13L/150)+1);"
            "forced-high omitted slot starts max(H+1,projected_wall)"
        ),
        (
            f"global_survivors={len(survivors)};"
            f"surviving_firsts={tuple(row[1] for row in survivors)};"
            f"surviving_bodies={tuple(row[0] for row in survivors)}"
        ),
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
                f"{label}:{ftext(value)}" for value, label, _kind in chosen
            )
            + f";lower={ftext(lower)};upper={ftext(upper)};gap={ftext(gap)}"
        )
    lines.extend(
        (
            "interval_1993_2003_global_survivors=()",
            "interval_1992_2003_global_survivors=(1992,)",
            f"global_profile_sha256={profile_hash}",
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
