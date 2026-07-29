#!/usr/bin/env python3
"""Exact geometric largest-drift filters for the k-aligned suffix bank.

This composes the THM-2941 exact singleton-excess bank with two necessary
geometric largest-drift conditions.

For six body speeds, k aligned tails and d=7-k drift tails, settled LRC
through 13 applied to the 13-k body/drift speeds gives a safe arc of radius

    k / (14 (14-k) M),

where M is its maximum speed.  If this radius is greater than one body
grid cell, the arc contains a full 1/L cell.  The k aligned combs cannot
cover that cell because their normalized union has mass at most k/7.
Thus a countercover requires

    M >= floor(k L / (14 (14-k))) + 1.

For the direct six-body boundary all body speeds are at most 14 and every
external drift is at least 15.  The script imposes the maximum of 15 and
this lower bound on the largest drift, and recomputes the rigorous
exact-suffix upper envelope for `k=2,...,6`.  It separately applies the
stronger projected-safe-arc wall

    M > k L / (7 (14-k) (1-u_k)),

where `u_k` is the proved uniform `k`-comb safe floor.  Labels through
HORIZON are integrated exactly.  Every omitted label is controlled by
delta(w)<=6r/(49w), starting at the larger of HORIZON+1 and the required
largest-drift floor.  No finite horizon is treated as exhaustive.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import math
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INDEPENDENT = (
    ROOT
    / "04-computation"
    / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
)
EXPECTED_SUPPORT_SHA256 = (
    "5482e10635ecf72840bc0c083360fd7ddad65c2885d743820061bcba58cd5609"
)
EXPECTED_INDEPENDENT_SHA256 = (
    "5d25a955fe184d6c1a3d8b632b4bbf901dc996ee46ad67c5748836fcc7134404"
)
HORIZON = 7_000
SAFE_FLOORS = {
    2: F(66, 91),
    3: F(55, 91),
    4: F(558, 1183),
    5: F(478, 1365),
    6: F(61, 273),
}
ETAS = {
    k: safe_floor - F(7 - k, 7)
    for k, safe_floor in SAFE_FLOORS.items()
}
PROJECTED_RATIOS = {
    k: F(k, 7 * (14 - k)) / (1 - safe_floor)
    for k, safe_floor in SAFE_FLOORS.items()
}
WALL_RATIOS = {
    "fullcell": {
        k: F(k, 14 * (14 - k)) for k in ETAS
    },
    "projected": PROJECTED_RATIOS,
}
DEFAULT_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_aligned_projected_arc_suffix_thm2941.out"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def normalized_sha256(path: Path) -> str:
    payload = path.read_bytes()
    require(
        b"\r" not in payload.replace(b"\r\n", b""),
        f"{path.name}: lone carriage return",
    )
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


spec = importlib.util.spec_from_file_location("thm2941_fullcell_audit", INDEPENDENT)
require(spec is not None and spec.loader is not None, "cannot load audit engine")
A = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
SUPPORT_NAMES = (
    "merge_intervals",
    "danger_intervals",
    "carrier_for",
    "danger_primitive",
    "singleton_coverage",
)
support_payload = "\n".join(
    inspect.getsource(getattr(A, name)) for name in SUPPORT_NAMES
).encode()
require(
    hashlib.sha256(support_payload).hexdigest() == EXPECTED_SUPPORT_SHA256,
    "independent THM-2941 carrier/singleton support changed",
)
require(
    normalized_sha256(INDEPENDENT) == EXPECTED_INDEPENDENT_SHA256,
    "independent THM-2941 source changed",
)
require(
    A.BASE_LABEL == 15 and A.RULER == 5_045_040,
    "independent THM-2941 constants changed",
)
require(
    ETAS
    == {
        2: F(1, 91),
        3: F(3, 91),
        4: F(51, 1183),
        5: F(88, 1365),
        6: F(22, 273),
    },
    "aligned safe-surplus constants changed",
)
require(
    PROJECTED_RATIOS
    == {
        2: F(13, 150),
        3: F(13, 132),
        4: F(2366, 21875),
        5: F(2275, 18627),
        6: F(819, 5936),
    },
    "projected-safe-arc ratios changed",
)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def top_insert(
    rows: list[tuple[F, int | None, str]],
    item: tuple[F, int | None, str],
    limit: int,
) -> None:
    rows.append(item)
    rows.sort(
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        )
    )
    del rows[limit:]


def suffix_upper(
    *,
    arbitrary_exact: list[tuple[F, int | None, str]],
    high_exact: list[tuple[F, int | None, str]],
    need: int,
    tail: F,
    high_tail: F,
    constrained: bool,
) -> tuple[tuple[F, int | None, str], ...]:
    """Largest rigorous suffix envelope, with one high item if constrained."""

    arbitrary = list(arbitrary_exact[:need])
    arbitrary.extend((tail, None, f"TAIL-{index}") for index in range(need))
    arbitrary.sort(
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        )
    )
    if not constrained:
        return tuple(arbitrary[:need])

    high = list(high_exact[:1])
    high.append((high_tail, None, "HIGH-TAIL"))
    high.sort(
        key=lambda row: (
            -row[0],
            HORIZON + 1 if row[1] is None else row[1],
            row[2],
        )
    )
    selected_high = high[0]

    # Exact labels cannot be reused.  Tail entries represent distinct omitted
    # labels, so an arbitrary tail remains available after selecting HIGH-TAIL.
    rest = [
        row
        for row in arbitrary
        if selected_high[1] is None or row[1] != selected_high[1]
    ]
    answer = (selected_high, *rest[: need - 1])
    require(len(answer) == need, "suffix envelope underfilled")
    return answer


def profile(body: tuple[int, ...]) -> dict[str, object]:
    carrier = A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), A.RULER)
    components = len(carrier)
    canonical_l = 14 * math.lcm(*body)
    require(h > 0 and components > 0, f"{body}: empty carrier")

    delta = {
        label: A.singleton_coverage(carrier, label) - h / 7
        for label in range(A.BASE_LABEL, HORIZON + 1)
        if label % canonical_l
    }
    ordinary_tail = F(6 * components, 49 * (HORIZON + 1))

    caps: dict[int, int] = {}
    high_floors: dict[str, dict[int, int]] = {
        mode: {} for mode in WALL_RATIOS
    }
    for k, eta in ETAS.items():
        drift_count = 7 - k
        bound = F(6 * drift_count * components, 49) / (h * eta)
        caps[k] = bound.numerator // bound.denominator
        require(caps[k] < HORIZON, f"{body};k={k}: cap reaches horizon")
        for mode, ratios in WALL_RATIOS.items():
            wall = ratios[k] * canonical_l
            high_floors[mode][k] = max(
                15,
                wall.numerator // wall.denominator + 1,
            )

    candidate_counts = {k: 0 for k in ETAS}
    baseline_counts = {k: 0 for k in ETAS}
    constrained_counts = {
        mode: {k: 0 for k in ETAS} for mode in WALL_RATIOS
    }
    baseline_frontier: dict[int, dict[str, object] | None] = {
        k: None for k in ETAS
    }
    constrained_frontier = {
        mode: {k: None for k in ETAS} for mode in WALL_RATIOS
    }

    arbitrary_top: list[tuple[F, int | None, str]] = []
    high_top = {
        mode: {k: [] for k in ETAS} for mode in WALL_RATIOS
    }

    for first in range(HORIZON, A.BASE_LABEL - 1, -1):
        if first % canonical_l == 0:
            continue
        first_delta = delta[first]
        for k, eta in ETAS.items():
            if first > caps[k]:
                continue
            candidate_counts[k] += 1
            need = 6 - k
            lower = h * eta

            baseline_suffix = suffix_upper(
                arbitrary_exact=arbitrary_top,
                high_exact=[],
                need=need,
                tail=ordinary_tail,
                high_tail=ordinary_tail,
                constrained=False,
            )
            baseline_upper = first_delta + sum(
                (value for value, _, _ in baseline_suffix), F(0)
            )
            if baseline_upper >= lower:
                baseline_counts[k] += 1
                row = {
                    "body": body,
                    "k": k,
                    "h": h,
                    "components": components,
                    "canonical_l": canonical_l,
                    "analytic_cap": caps[k],
                    "high_floor": None,
                    "mode": "baseline",
                    "first": first,
                    "first_delta": first_delta,
                    "chosen": baseline_suffix,
                    "upper": baseline_upper,
                    "lower": lower,
                    "gap": baseline_upper - lower,
                }
                if (
                    baseline_frontier[k] is None
                    or first > baseline_frontier[k]["first"]
                ):
                    baseline_frontier[k] = row

            for mode in WALL_RATIOS:
                high_floor = high_floors[mode][k]
                if need == 0:
                    if first < high_floor:
                        continue
                    constrained_suffix = ()
                else:
                    constrained = first < high_floor
                    high_tail_start = max(HORIZON + 1, high_floor)
                    high_tail = F(6 * components, 49 * high_tail_start)
                    constrained_suffix = suffix_upper(
                        arbitrary_exact=arbitrary_top,
                        high_exact=high_top[mode][k],
                        need=need,
                        tail=ordinary_tail,
                        high_tail=high_tail,
                        constrained=constrained,
                    )
                constrained_upper = first_delta + sum(
                    (value for value, _, _ in constrained_suffix), F(0)
                )
                if constrained_upper >= lower:
                    constrained_counts[mode][k] += 1
                    row = {
                        "body": body,
                        "k": k,
                        "h": h,
                        "components": components,
                        "canonical_l": canonical_l,
                        "analytic_cap": caps[k],
                        "high_floor": high_floor,
                        "mode": mode,
                        "first": first,
                        "first_delta": first_delta,
                        "chosen": constrained_suffix,
                        "upper": constrained_upper,
                        "lower": lower,
                        "gap": constrained_upper - lower,
                    }
                    if (
                        constrained_frontier[mode][k] is None
                        or first > constrained_frontier[mode][k]["first"]
                    ):
                        constrained_frontier[mode][k] = row

        item = (first_delta, first, "EXACT")
        top_insert(arbitrary_top, item, limit=4)
        for mode in WALL_RATIOS:
            for k in ETAS:
                if first >= high_floors[mode][k]:
                    top_insert(high_top[mode][k], item, limit=4)

    return {
        "body": body,
        "caps": caps,
        "high_floors": high_floors,
        "candidate_counts": candidate_counts,
        "baseline_counts": baseline_counts,
        "constrained_counts": constrained_counts,
        "baseline_frontier": baseline_frontier,
        "constrained_frontier": constrained_frontier,
    }


def row_tuple(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["body"],
        row["k"],
        row["h"],
        row["components"],
        row["canonical_l"],
        row["analytic_cap"],
        row["mode"],
        row["high_floor"],
        row["first"],
        row["first_delta"],
        row["chosen"],
        row["upper"],
        row["lower"],
    )


def chosen_text(
    chosen: tuple[tuple[F, int | None, str], ...],
) -> str:
    return ",".join(
        (
            f"{label}:{ftext(value)}"
            if label is not None
            else f"{kind}:{ftext(value)}"
        )
        for value, label, kind in chosen
    )


def render(rows: list[dict[str, object]]) -> str:
    require(len(rows) == math.comb(14, 6) == 3_003, "root universe changed")
    lines = [
        "LRC14 k-aligned geometric largest-drift constrained suffix scan",
        f"independent_source_sha256={normalized_sha256(INDEPENDENT)}",
        f"support_sha256={hashlib.sha256(support_payload).hexdigest()}",
        "universe=(six_body_roots=3003,body_labels=1..14,"
        "external_drifts>=15,drift_mod_L_nonzero,k=2..6)",
        f"exact_horizon={HORIZON}",
        "fullcell_wall=max(15,floor(kL/[14(14-k)])+1)",
        "projected_wall=max(15,floor(alpha_k L)+1);"
        "alpha=(k2:13/150,k3:13/132,k4:2366/21875,"
        "k5:2275/18627,k6:819/5936)",
        "omitted_tail=delta(w)<=6r/(49w);ordinary omitted slots start H+1;"
        "forced-high slot starts max(H+1,largest_floor)",
    ]

    digest_rows: list[tuple[object, ...]] = []
    expected_baseline = {
        2: (9_787_617, 2_537_423, 2_340),
        3: (2_578_195, 498_150, 432),
        4: (1_460_123, 229_638, 260),
        5: (626_787, 36_301, 130),
        6: (224_618, 1_602, 44),
    }
    expected_constrained = {
        "fullcell": {
            2: (2_463_393, 2_340, 1),
            3: (419_380, 410, 2),
            4: (112_301, 208, 2),
            5: (5_942, 78, 4),
            6: (0, None, 0),
        },
        "projected": {
            2: (2_239_853, 2_142, 1),
            3: (376_020, 380, 1),
            4: (87_975, 182, 76),
            5: (4_702, 66, 1),
            6: (0, None, 0),
        },
    }
    for k in ETAS:
        candidates = sum(row["candidate_counts"][k] for row in rows)
        baseline = sum(row["baseline_counts"][k] for row in rows)
        baseline_rows = [
            row["baseline_frontier"][k]
            for row in rows
            if row["baseline_frontier"][k] is not None
        ]
        baseline_max = max(row["first"] for row in baseline_rows)
        require(
            (candidates, baseline, baseline_max) == expected_baseline[k],
            f"k={k}: baseline cross-check changed",
        )
        lines.append(
            f"k={k};drifts={7-k};candidate_first_rows={candidates};"
            f"baseline_surviving_first_rows={baseline};"
            f"baseline_max_first={baseline_max}"
        )
        for mode in WALL_RATIOS:
            constrained = sum(
                row["constrained_counts"][mode][k] for row in rows
            )
            constrained_rows = [
                row["constrained_frontier"][mode][k]
                for row in rows
                if row["constrained_frontier"][mode][k] is not None
            ]
            if constrained_rows:
                constrained_max: int | None = max(
                    row["first"] for row in constrained_rows
                )
                frontiers = sorted(
                    (
                        row
                        for row in constrained_rows
                        if row["first"] == constrained_max
                    ),
                    key=row_tuple,
                )
            else:
                constrained_max = None
                frontiers = []
            require(
                (constrained, constrained_max, len(frontiers))
                == expected_constrained[mode][k],
                f"k={k};mode={mode}: constrained aggregate changed",
            )
            digest_rows.extend(row_tuple(row) for row in frontiers)
            lines.append(
                f"  mode={mode};surviving_first_rows={constrained};"
                f"max_first={constrained_max};frontier_rows={len(frontiers)}"
            )
            for row in frontiers:
                lines.append(
                    "    FRONTIER;"
                    f"E={','.join(map(str,row['body']))};"
                    f"h={ftext(row['h'])};r={row['components']};"
                    f"L={row['canonical_l']};analytic_cap={row['analytic_cap']};"
                    f"largest_floor={row['high_floor']};"
                    f"z1={row['first']};delta1={ftext(row['first_delta'])};"
                    f"suffix={chosen_text(row['chosen'])};"
                    f"lower={ftext(row['lower'])};upper={ftext(row['upper'])};"
                    f"gap={ftext(row['gap'])}"
                )

    digest = hashlib.sha256(
        b"LRC14/k-aligned/geometric-constrained-suffix/v2\n"
        + repr(tuple(digest_rows)).encode()
    ).hexdigest()
    require(
        digest
        == "49bff626dfeeb87ce585f655b5a629e47acff10d9ebe28726245b81846a9ded6",
        "geometric constrained frontier digest changed",
    )
    lines.append(f"frontier_digest={digest}")
    lines.append("all_exact_controls=PASS")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count() or 1))
    parser.add_argument("--hash-seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(args.hash_seed >= 0, "hash seed must be nonnegative")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        rows = [profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(profile, roots, chunksize=1))
    rows.sort(key=lambda row: row["body"])
    output = render(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
