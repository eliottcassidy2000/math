#!/usr/bin/env python3
"""One-pass exact-envelope atlas for projected ``k=3`` scalar rows.

This reconstructs every projected scalar row with first drift in
``200 <= z1 <= 378`` over the 3,003 six-body subsets of ``{1,...,14}``.
For a fixed body the singleton excesses through the inherited exact horizon
7,000 are computed once.  A descending suffix scan then answers all 179
first-drift queries from the same top-three data, including the projected
high-label obligation.  Labels beyond the horizon use the proved THM-1094
bound ``delta(z) <= 6r/(49z)``; the horizon is never treated as exhaustive.

The transcript records the complete per-height census and the explicit
25-body ``z1=364`` slice consumed by the all-label ray/status referee.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUFFIX_PATH = (
    ROOT
    / "04-computation"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.py"
)
SUFFIX_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.out"
)
OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_projected_scalar_atlas_thm2941.out"
)
EXPECTED_SUFFIX_SHA256 = (
    "a003d287f618eb301edf6974d0b67dc128c4f380a169e7809ed5b5754e8b8303"
)
EXPECTED_SUFFIX_OUTPUT_SHA256 = (
    "61e16aab8a368881c574047e576645e6b41837dc9f804f7a78d37230d843612b"
)
EXPECTED_TOTAL_ROWS = 6_060
EXPECTED_Z364_ROWS = 25
EXPECTED_UPPER_SPIKES = ((336, 8), (350, 53), (364, 25), (378, 9))
EXPECTED_Z364_BODIES = (
    (1, 2, 6, 8, 12, 14),
    (1, 2, 6, 10, 12, 14),
    (1, 2, 8, 10, 12, 14),
    (1, 4, 6, 8, 12, 14),
    (1, 4, 6, 10, 12, 14),
    (1, 4, 8, 10, 12, 14),
    (1, 6, 8, 9, 12, 14),
    (1, 6, 8, 10, 12, 14),
    (2, 3, 8, 10, 12, 14),
    (2, 3, 8, 11, 12, 14),
    (2, 4, 6, 8, 12, 14),
    (2, 4, 6, 10, 12, 14),
    (2, 4, 8, 10, 12, 14),
    (2, 6, 8, 9, 12, 14),
    (2, 6, 8, 10, 12, 14),
    (2, 6, 8, 10, 13, 14),
    (2, 6, 8, 11, 12, 14),
    (2, 6, 9, 10, 12, 14),
    (2, 6, 10, 11, 12, 14),
    (2, 8, 9, 10, 12, 14),
    (2, 8, 10, 11, 12, 14),
    (2, 8, 10, 12, 13, 14),
    (3, 4, 8, 10, 12, 14),
    (4, 6, 8, 10, 12, 14),
    (4, 8, 10, 11, 12, 14),
)
EXPECTED_SEMANTIC_SHA256 = (
    "46535b2b075e15244ec87101f44d74a8034b093244fd8559fe11b5b70425fda8"
)

FIRST_MIN = 200
FIRST_MAX = 378


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot load", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(
    file_sha256(SUFFIX_PATH) == EXPECTED_SUFFIX_SHA256,
    "THM-2941 suffix dependency changed",
)
require(
    file_sha256(SUFFIX_OUTPUT_PATH) == EXPECTED_SUFFIX_OUTPUT_SHA256,
    "THM-2941 suffix transcript changed",
)
suffix = load_module("thm2941_projected_suffix", SUFFIX_PATH)
require(
    suffix.HORIZON == 7_000
    and suffix.ETAS[3] == F(3, 91)
    and suffix.PROJECTED_RATIOS[3] == F(13, 132),
    "projected k=3 constants changed",
)


def row_text(row: tuple[object, ...]) -> str:
    body, h, components, L, high_floor, first, first_delta, chosen, lower, upper = row
    chosen_text = ",".join(
        f"{label if label is not None else source}:{value.numerator}/{value.denominator}"
        for value, label, source in chosen
    )
    gap = upper - lower
    return (
        f"E={','.join(map(str, body))};h={h.numerator}/{h.denominator};"
        f"r={components};L={L};high={high_floor};z1={first};"
        f"delta1={first_delta.numerator}/{first_delta.denominator};"
        f"suffix={chosen_text};lower={lower.numerator}/{lower.denominator};"
        f"upper={upper.numerator}/{upper.denominator};"
        f"gap={gap.numerator}/{gap.denominator}"
    )


def body_rows(body: tuple[int, ...]) -> tuple[tuple[object, ...], ...]:
    """Return every surviving atlas row for one body in one suffix pass."""

    carrier = suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), suffix.A.RULER)
    components = len(carrier)
    L = 14 * math.lcm(*body)
    require(h > 0 and components > 0, (body, "empty carrier"))

    lower = h * suffix.ETAS[3]
    wall = suffix.PROJECTED_RATIOS[3] * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    ordinary_tail = F(6 * components, 49 * (suffix.HORIZON + 1))
    high_tail = F(
        6 * components,
        49 * max(suffix.HORIZON + 1, high_floor),
    )

    # Computing the full exact suffix once is the key one-pass economy.
    delta = {
        label: suffix.A.singleton_coverage(carrier, label) - h / 7
        for label in range(FIRST_MIN, suffix.HORIZON + 1)
        if label % L
    }
    arbitrary_top: list[tuple[F, int | None, str]] = []
    high_top: list[tuple[F, int | None, str]] = []
    rows = []
    for first in range(suffix.HORIZON, FIRST_MIN - 1, -1):
        if first % L == 0:
            continue
        first_delta = delta[first]
        if first <= FIRST_MAX:
            chosen = suffix.suffix_upper(
                arbitrary_exact=arbitrary_top,
                high_exact=high_top,
                need=3,
                tail=ordinary_tail,
                high_tail=high_tail,
                constrained=first < high_floor,
            )
            upper = first_delta + sum((item[0] for item in chosen), F(0))
            if upper >= lower:
                rows.append(
                    (
                        body,
                        h,
                        components,
                        L,
                        high_floor,
                        first,
                        first_delta,
                        chosen,
                        lower,
                        upper,
                    )
                )

        item = (first_delta, first, "EXACT")
        suffix.top_insert(arbitrary_top, item, limit=3)
        if first >= high_floor:
            suffix.top_insert(high_top, item, limit=1)
    return tuple(rows)


def atlas(processes: int) -> tuple[tuple[object, ...], ...]:
    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == 3_003, "six-body universe changed")
    if processes == 1:
        nested = map(body_rows, bodies)
        rows = tuple(row for body_result in nested for row in body_result)
    else:
        with mp.Pool(processes) as pool:
            nested = pool.imap_unordered(body_rows, bodies, chunksize=4)
            rows = tuple(row for body_result in nested for row in body_result)
    return tuple(sorted(rows, key=lambda row: (row[5], row[0])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--processes",
        type=int,
        default=max(1, mp.cpu_count() // 2),
    )
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")

    rows = atlas(args.processes)
    counts = Counter(row[5] for row in rows)
    z364_bodies = tuple(row[0] for row in rows if row[5] == 364)
    require(len(rows) == EXPECTED_TOTAL_ROWS, ("atlas total changed", len(rows)))
    require(
        len(z364_bodies) == EXPECTED_Z364_ROWS,
        ("z1=364 slice changed", len(z364_bodies)),
    )
    require(
        z364_bodies == EXPECTED_Z364_BODIES,
        ("z1=364 body ledger changed", z364_bodies),
    )
    upper_counts = tuple(
        (first, counts[first])
        for first in range(336, FIRST_MAX + 1)
        if counts[first]
    )
    require(
        upper_counts == EXPECTED_UPPER_SPIKES,
        ("upper-bank spike ledger changed", upper_counts),
    )

    semantic = hashlib.sha256()
    for row in rows:
        semantic.update((row_text(row) + "\n").encode())
    semantic_sha256 = semantic.hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "projected scalar atlas semantic digest changed",
        )

    lines = [
        "LRC14 projected k=3 one-pass scalar atlas",
        f"suffix_source_sha256={EXPECTED_SUFFIX_SHA256}",
        f"suffix_output_sha256={EXPECTED_SUFFIX_OUTPUT_SHA256}",
        (
            "universe=six_body_roots=3003;first_drift=200..378;"
            "three_distinct_later_drifts;nonzero_mod_L;exact_horizon=7000"
        ),
        (
            "omitted_tail=delta(z)<=6r/(49z);"
            "projected_high_wall=max(15,floor(13L/132)+1)"
        ),
        "method=one exact descending suffix pass per body for all 179 first drifts",
        f"total_projected_scalar_rows={len(rows)}",
        f"counts_by_first={tuple(sorted(counts.items()))}",
        f"upper_bank_spikes={upper_counts};all_interstitial_counts=0",
        f"z1_364_bodies={z364_bodies}",
    ]
    lines.extend(
        f"z1_364_row={row_text(row)}" for row in rows if row[5] == 364
    )
    lines.extend(
        (
            f"semantic_sha256={semantic_sha256}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
