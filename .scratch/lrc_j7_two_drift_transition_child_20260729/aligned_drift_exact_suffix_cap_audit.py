#!/usr/bin/env python3
"""Exact suffix audit for the smallest drift on every k-aligned face.

Scratch artifact suitable for promotion after its dependency pin and expected
digest are frozen.  For each six-body root and each ``2 <= k <= 6``, put
``d=7-k``.  If ``k`` endpoint-grid-aligned tails and ``d`` drift tails cover
the literal carrier, then

    sum_q delta(z_q) >= h * eta_k,

where

    eta_2=1/91, eta_3=3/91, eta_4=51/1183,
    eta_5=88/1365, eta_6=22/273.

The first drift is ordered ``z_1 < ... < z_d``.  THM-1094 gives

    delta(w) <= 6 r / (49 w).

This first supplies a rootwise analytic cap.  Inside that cap the verifier
computes ``delta(z_1)`` exactly, retains the largest ``d-1`` distinct exact
suffix values through ``HORIZON``, and pads the suffix by ``d-1`` copies of
the rigorous omitted-label bound ``6r/(49(HORIZON+1))``.  Thus no finite
search horizon is treated as a theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import inspect
import math
import multiprocessing as mp
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INDEPENDENT = (
    ROOT
    / "04-computation"
    / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
)
EXPECTED_SUPPORT_SHA256 = (
    "5482e10635ecf72840bc0c083360fd7ddad65c2885d743820061bcba58cd5609"
)
HORIZON = 7_000
ETAS = {
    2: F(1, 91),
    3: F(3, 91),
    4: F(51, 1183),
    5: F(88, 1365),
    6: F(22, 273),
}


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


spec = importlib.util.spec_from_file_location("thm2941_suffix_audit", INDEPENDENT)
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
    A.BASE_LABEL == 15 and A.RULER == 5_045_040,
    "independent THM-2941 constants changed",
)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def top_insert(
    rows: list[tuple[F, int]],
    item: tuple[F, int],
    limit: int,
) -> None:
    rows.append(item)
    rows.sort(key=lambda row: (-row[0], row[1]))
    del rows[limit:]


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
    tail = F(6 * components, 49 * (HORIZON + 1))

    caps: dict[int, int] = {}
    for k, eta in ETAS.items():
        drift_count = 7 - k
        bound = F(6 * drift_count * components, 49) / (h * eta)
        caps[k] = bound.numerator // bound.denominator
        require(
            caps[k] < HORIZON,
            f"{body};k={k}: analytic cap reaches exact horizon",
        )

    candidate_counts = {k: 0 for k in ETAS}
    survivor_counts = {k: 0 for k in ETAS}
    frontier: dict[int, dict[str, object] | None] = {k: None for k in ETAS}
    exact_top: list[tuple[F, int]] = []

    for first in range(HORIZON, A.BASE_LABEL - 1, -1):
        if first % canonical_l == 0:
            continue
        first_delta = delta[first]
        for k, eta in ETAS.items():
            if first > caps[k]:
                continue
            candidate_counts[k] += 1
            need = 6 - k
            pool: list[tuple[F, int | None]] = [
                (value, label) for value, label in exact_top[:need]
            ]
            pool.extend((tail, None) for _ in range(need))
            pool.sort(
                key=lambda row: (
                    -row[0],
                    HORIZON + 1 if row[1] is None else row[1],
                )
            )
            chosen = tuple(pool[:need])
            upper = first_delta + sum((value for value, _ in chosen), F(0))
            lower = h * eta
            if upper < lower:
                continue
            survivor_counts[k] += 1
            row = {
                "body": body,
                "k": k,
                "drifts": 7 - k,
                "h": h,
                "components": components,
                "canonical_l": canonical_l,
                "analytic_cap": caps[k],
                "first": first,
                "first_delta": first_delta,
                "chosen": chosen,
                "upper": upper,
                "lower": lower,
                "gap": upper - lower,
                "tail": tail,
                "tail_used": any(label is None for _, label in chosen),
            }
            if frontier[k] is None or first > frontier[k]["first"]:
                frontier[k] = row
        top_insert(exact_top, (first_delta, first), limit=4)

    return {
        "body": body,
        "h": h,
        "components": components,
        "canonical_l": canonical_l,
        "caps": caps,
        "candidate_counts": candidate_counts,
        "survivor_counts": survivor_counts,
        "frontier": frontier,
    }


def row_tuple(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row["body"],
        row["k"],
        row["h"],
        row["components"],
        row["canonical_l"],
        row["analytic_cap"],
        row["first"],
        row["first_delta"],
        row["chosen"],
        row["upper"],
        row["lower"],
        row["tail"],
        row["tail_used"],
    )


def render(rows: list[dict[str, object]]) -> str:
    require(len(rows) == math.comb(14, 6) == 3_003, "root universe changed")
    lines = [
        "LRC14 k-aligned exact smallest-drift suffix audit",
        f"independent_source_sha256={normalized_sha256(INDEPENDENT)}",
        f"support_sha256={hashlib.sha256(support_payload).hexdigest()}",
        "universe=(six_body_roots=3003,body_labels=1..14,"
        "external_drifts>=15,drift_mod_L_nonzero,k=2..6)",
        f"exact_horizon={HORIZON}",
        "tail=delta(w)<=6r/(49w)<=6r/(49*(H+1)) for omitted w>H",
        "eta=(k2:1/91,k3:3/91,k4:51/1183,k5:88/1365,k6:22/273)",
    ]

    digest_rows: list[tuple[object, ...]] = []
    for k in ETAS:
        candidates = sum(row["candidate_counts"][k] for row in rows)
        survivors = sum(row["survivor_counts"][k] for row in rows)
        local_frontiers = [
            row["frontier"][k]
            for row in rows
            if row["frontier"][k] is not None
        ]
        require(local_frontiers, f"k={k}: survivor bank is empty")
        maximum = max(row["first"] for row in local_frontiers)
        frontier_rows = sorted(
            (
                row
                for row in local_frontiers
                if row["first"] == maximum
            ),
            key=lambda row: row_tuple(row),
        )
        digest_rows.extend(row_tuple(row) for row in frontier_rows)
        lines.append(
            f"k={k};drifts={7-k};candidate_first_rows={candidates};"
            f"surviving_first_rows={survivors};max_first={maximum};"
            f"frontier_rows={len(frontier_rows)};"
            f"tail_used_at_frontier={sum(row['tail_used'] for row in frontier_rows)}"
        )
        for row in frontier_rows:
            chosen = ",".join(
                (
                    f"{label}:{ftext(value)}"
                    if label is not None
                    else f"TAIL:{ftext(value)}"
                )
                for value, label in row["chosen"]
            )
            lines.append(
                "  FRONTIER;"
                f"E={','.join(map(str,row['body']))};"
                f"h={ftext(row['h'])};r={row['components']};"
                f"L={row['canonical_l']};analytic_cap={row['analytic_cap']};"
                f"z1={row['first']};delta1={ftext(row['first_delta'])};"
                f"suffix={chosen};lower={ftext(row['lower'])};"
                f"upper={ftext(row['upper'])};gap={ftext(row['gap'])}"
            )

    digest = hashlib.sha256(
        b"LRC14/THM2941/k-aligned-suffix/v1\n"
        + repr(tuple(digest_rows)).encode()
    ).hexdigest()
    lines.append(f"frontier_digest={digest}")
    lines.append("all_exact_controls=PASS")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, mp.cpu_count() or 1),
    )
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    roots = tuple(combinations(range(1, 15), 6))
    if args.workers == 1:
        rows = [profile(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(profile, roots, chunksize=1))
    rows.sort(key=lambda row: row["body"])
    print(render(rows), end="")


if __name__ == "__main__":
    main()
