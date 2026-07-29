#!/usr/bin/env python3
"""Close the projected k=3 first-drift row z1=379 in THM-2941.

The inherited projected suffix census has maximum first drift 380.  The
unique row there was closed by the packet/fibre referee.  This computation
tests the next integer first drift, 379, over all 3,003 literal six-body
carriers.

For each body it reconstructs the carrier with the independent integer-ruler
engine, integrates every allowed suffix label through 7,000 exactly, and
retains the exact largest three-label excess envelope subject to the
projected high-label wall.  Omitted labels receive the proved THM-1094 bound

    delta(z) <= 6 r / (49 z).

No row survives.  A live positive control independently reconstructs the
known unique z1=380 frontier body and its exact scalar envelope.  This is a
one-row scalar closure, not a census of any z1<=378 row.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUFFIX = (
    ROOT
    / "04-computation"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.py"
)
SUFFIX_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_next_frontier_scalar_closure_thm2941.out"
)
EXPECTED_SUFFIX_SHA256 = (
    "a003d287f618eb301edf6974d0b67dc128c4f380a169e7809ed5b5754e8b8303"
)
EXPECTED_SUFFIX_OUTPUT_SHA256 = (
    "61e16aab8a368881c574047e576645e6b41837dc9f804f7a78d37230d843612b"
)
EXPECTED_SEMANTIC_SHA256 = (
    "48ab29334a93fd0087d9645513be14f884a30bd014c2f05329c1f7d0c295d4ee"
)

FIRST = 379
EXPECTED_CRUDE_ROWS = 2_579
EXPECTED_POSITIVE = (
    (1, 4, 8, 10, 12, 14),
    F(1049, 2940),
    34,
    11760,
    1159,
    380,
    F(1531, 391020),
    (
        (F(237, 171500), 1500, "EXACT"),
        (F(2981, 843780), 492, "EXACT"),
        (F(2687, 843780), 410, "EXACT"),
    ),
    F(1049, 89180),
    F(401287, 33399625),
    F(437649, 1736780500),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(SUFFIX) == EXPECTED_SUFFIX_SHA256,
    "THM-2941 suffix dependency changed",
)
require(
    file_sha256(SUFFIX_OUTPUT) == EXPECTED_SUFFIX_OUTPUT_SHA256,
    "THM-2941 suffix transcript dependency changed",
)
spec = importlib.util.spec_from_file_location("thm2941_suffix", SUFFIX)
require(spec is not None and spec.loader is not None, f"cannot import {SUFFIX}")
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)
require(
    S.HORIZON == 7_000
    and S.ETAS[3] == F(3, 91)
    and S.PROJECTED_RATIOS[3] == F(13, 132),
    "projected k=3 constants changed",
)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def evaluate(args: tuple[tuple[int, ...], int]):
    """Return ``(passed_crude, exact_record_or_none)`` for one body."""

    body, first = args
    carrier = S.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), S.A.RULER)
    components = len(carrier)
    L = 14 * math.lcm(*body)
    require(h > 0 and components > 0, (body, "empty carrier"))
    if first % L == 0:
        return False, None

    lower = h * S.ETAS[3]
    first_delta = S.A.singleton_coverage(carrier, first) - h / 7
    wall = S.PROJECTED_RATIOS[3] * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    ordinary_tail = F(6 * components, 49 * (S.HORIZON + 1))
    high_tail = F(
        6 * components,
        49 * max(S.HORIZON + 1, high_floor),
    )

    # This preliminary all-tail upper bound is only a speed optimization.
    C = F(6 * components, 49)
    next_label = first + 1
    if first < high_floor:
        crude = (
            first_delta
            + 2 * C / next_label
            + C / max(next_label, high_floor)
        )
    else:
        crude = first_delta + 3 * C / next_label
    if crude < lower:
        return False, None

    arbitrary: list[tuple[F, int | None, str]] = []
    high: list[tuple[F, int | None, str]] = []
    for label in range(first + 1, S.HORIZON + 1):
        if label % L == 0:
            continue
        delta = S.A.singleton_coverage(carrier, label) - h / 7
        S.top_insert(arbitrary, (delta, label, "EXACT"), 3)
        if label >= high_floor:
            S.top_insert(high, (delta, label, "EXACT"), 1)

    chosen = S.suffix_upper(
        arbitrary_exact=arbitrary,
        high_exact=high,
        need=3,
        tail=ordinary_tail,
        high_tail=high_tail,
        constrained=first < high_floor,
    )
    upper = first_delta + sum((value for value, _, _ in chosen), F(0))
    return True, (
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
        upper - lower,
    )


def record_text(record) -> str:
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
        gap,
    ) = record
    suffix = ",".join(
        f"{label if label is not None else source}:{ftext(value)}"
        for value, label, source in chosen
    )
    return (
        f"E={','.join(map(str, body))};h={ftext(h)};r={components};"
        f"L={L};high={high_floor};z1={first};"
        f"delta1={ftext(first_delta)};suffix={suffix};"
        f"lower={ftext(lower)};upper={ftext(upper)};gap={ftext(gap)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--processes",
        type=int,
        default=max(1, mp.cpu_count() // 2),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    bodies = tuple(combinations(range(1, 15), 6))
    with mp.Pool(args.processes) as pool:
        evaluated = tuple(
            pool.imap_unordered(
                evaluate,
                ((body, FIRST) for body in bodies),
                chunksize=8,
            )
        )
    crude_count = sum(passed for passed, _record in evaluated)
    records = tuple(
        sorted(record for _passed, record in evaluated if record is not None)
    )
    survivors = tuple(record for record in records if record[-1] >= 0)
    require(len(bodies) == 3_003, "body universe changed")
    require(crude_count == EXPECTED_CRUDE_ROWS, "crude row count changed")
    require(len(records) == EXPECTED_CRUDE_ROWS, "exact row ledger incomplete")
    require(not survivors, ("z1=379 scalar row survived", survivors[:3]))

    closest = max(records, key=lambda record: (record[-1], tuple(-x for x in record[0])))
    positive = evaluate((EXPECTED_POSITIVE[0], 380))
    require(positive == (True, EXPECTED_POSITIVE), "z1=380 positive control changed")

    semantic = hashlib.sha256()
    for record in records:
        semantic.update((record_text(record) + "\n").encode())
    semantic_sha256 = semantic.hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "z1=379 semantic ledger changed",
        )

    lines = [
        "LRC14 projected k=3 next-frontier scalar closure",
        f"suffix_source_sha256={EXPECTED_SUFFIX_SHA256}",
        f"suffix_output_sha256={EXPECTED_SUFFIX_OUTPUT_SHA256}",
        (
            "universe=six_body_roots=3003;first_drift=379;"
            "three_distinct_later_drifts;nonzero_mod_L;exact_horizon=7000"
        ),
        (
            "omitted_tail=delta(z)<=6r/(49z);"
            "projected_high_wall=max(15,floor(13L/132)+1)"
        ),
        f"crude_all_tail_rows={crude_count}",
        f"exact_projected_scalar_rows={len(survivors)}",
        f"closest_rejection={record_text(closest)}",
        f"positive_control={record_text(EXPECTED_POSITIVE)}",
        "conclusion=projected k=3 first-drift cap improves 379->378",
        f"semantic_sha256={semantic_sha256}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
