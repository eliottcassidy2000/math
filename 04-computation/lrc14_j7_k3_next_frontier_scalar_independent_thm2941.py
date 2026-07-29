#!/usr/bin/env python3
"""Independent rational-interval audit of the THM-2941 z1=379 closure.

The primary next-frontier referee uses the integer-ruler carrier and scalar
primitive.  This audit instead reconstructs each body carrier through the
THM-2923 rational interval engine and evaluates its singleton masses with the
guarded NumPy vector primitive.  It must reproduce the complete 2,579-row
semantic ledger, not only the empty survivor count.
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
PRIMARY_ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_critical_scalar_wall_balanced_boundary_thm2941.py"
)
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_k3_next_frontier_scalar_independent_thm2941.out"
)
EXPECTED_PRIMARY_ENGINE_SHA256 = (
    "37a39f7a7409848cc637c9689c519e4452eda71b83d9042bfb1fc98b3ffd8efc"
)
EXPECTED_SEMANTIC_SHA256 = (
    "48ab29334a93fd0087d9645513be14f884a30bd014c2f05329c1f7d0c295d4ee"
)

FIRST = 379
ETA = F(3, 91)
PROJECTED_RATIO = F(13, 132)
HORIZON = 7_000
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
    file_sha256(PRIMARY_ENGINE) == EXPECTED_PRIMARY_ENGINE_SHA256,
    "THM-2941 primary rational engine changed",
)
spec = importlib.util.spec_from_file_location("thm2941_primary", PRIMARY_ENGINE)
require(
    spec is not None and spec.loader is not None,
    f"cannot import {PRIMARY_ENGINE}",
)
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def choose_suffix(
    ranked: list[tuple[F, int]],
    *,
    first: int,
    high_floor: int,
    ordinary_tail: F,
    high_tail: F,
) -> tuple[tuple[F, int | None, str], ...]:
    arbitrary = [(value, label, "EXACT") for value, label in ranked[:3]]
    arbitrary.extend(
        (ordinary_tail, None, f"TAIL-{index}") for index in range(3)
    )
    arbitrary.sort(
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        )
    )
    if first >= high_floor:
        return tuple(arbitrary[:3])

    exact_high = next(
        (
            (value, label, "EXACT")
            for value, label in ranked
            if label >= high_floor
        ),
        None,
    )
    high_candidates = [(high_tail, None, "HIGH-TAIL")]
    if exact_high is not None:
        high_candidates.append(exact_high)
    selected_high = min(
        high_candidates,
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        ),
    )
    rest = [
        row
        for row in arbitrary
        if selected_high[1] is None or row[1] != selected_high[1]
    ]
    return (selected_high, *rest[:2])


def evaluate(args: tuple[tuple[int, ...], int]):
    body, first = args
    carrier, components, h = P.S.direct_carrier(body)
    require(carrier and components == len(carrier) and h > 0, "empty carrier")
    L = 14 * math.lcm(*body)
    if first % L == 0:
        return False, None
    lower = h * ETA
    first_delta = P.S.H.exact_coverages(carrier, [first])[0][0] - h / 7
    wall = PROJECTED_RATIO * L
    high_floor = max(15, wall.numerator // wall.denominator + 1)
    C = F(6 * components, 49)

    if first < high_floor:
        crude = (
            first_delta
            + 2 * C / (first + 1)
            + C / max(first + 1, high_floor)
        )
    else:
        crude = first_delta + 3 * C / (first + 1)
    if crude < lower:
        return False, None

    labels = [
        label
        for label in range(first + 1, HORIZON + 1)
        if label % L
    ]
    exact = P.S.H.exact_coverages(carrier, labels)
    ranked = sorted(
        ((value - h / 7, label) for value, label in exact),
        key=lambda item: (-item[0], item[1]),
    )
    ordinary_tail = F(6 * components, 49 * (HORIZON + 1))
    high_tail = F(
        6 * components,
        49 * max(HORIZON + 1, high_floor),
    )
    chosen = choose_suffix(
        ranked,
        first=first,
        high_floor=high_floor,
        ordinary_tail=ordinary_tail,
        high_tail=high_tail,
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
            "independent z1=379 semantic ledger changed",
        )

    lines = [
        "LRC14 projected k=3 next-frontier independent scalar audit",
        f"primary_engine_sha256={EXPECTED_PRIMARY_ENGINE_SHA256}",
        (
            "universe=six_body_roots=3003;first_drift=379;"
            "three_distinct_later_drifts;nonzero_mod_L;exact_horizon=7000"
        ),
        (
            "engine=rational_intervals+guarded_numpy_vector_primitive;"
            "omitted_tail=delta(z)<=6r/(49z)"
        ),
        f"crude_all_tail_rows={crude_count}",
        f"exact_projected_scalar_rows={len(survivors)}",
        f"closest_rejection={record_text(closest)}",
        f"positive_control={record_text(EXPECTED_POSITIVE)}",
        f"semantic_sha256={semantic_sha256}",
        "independent_agreement=PASS",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
