#!/usr/bin/env python3
r"""Deterministic 50-carrier sample of the refined j=5 residual graph.

The sample consists of the first 50 lexicographic scalar-failure first-apex
carriers which have a positive global ``5m/7`` pair-cap margin but are not
already closed by ``2B<m``.  For each carrier this script:

* constructs ``theta=m-B`` and the finite head ``H={x:c(x)>=theta/2}``;
* enumerates every literal ``theta``-heavy edge inside ``H``; and
* retains an edge only if the exact global pair cap on its literal residual
  fails to prove positive measure.

This bounded sample estimates the usefulness and cost of the refined
feasibility graph.  It is not an all-carrier theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = (
    ROOT
    / "04-computation/lrc14_j5_pair_hostiles_recursive_audit_codex_20260729.py"
)
AUDIT_SHA256 = (
    "aa1880b84a3d54693dffd29d2a687c622d360aee3a975eef7757ada2828389bf"
)
SAMPLE_SIZE = 50
EXPECTED_SAMPLE_BOUNDARY = (
    ((1, 2, 3, 4, 5, 6, 7, 9), 3, 39),
    ((1, 2, 3, 4, 6, 7, 9, 10), 4, 39),
)
EXPECTED_SAMPLE_DIGEST = (
    "229e68ecbac2963c1fccf23bd04edb7d08f98a39a90e4ae0606fece6cb853148"
)
EXPECTED_COUNTS = (279, 10, 434, 31, 0, 0)
EXPECTED_EDGE_DIGEST = (
    "075feff178850a6d1ba7632ecaed9ddd56c324bba1945031f2c4b71792dccb5a"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_audit():
    require(
        hashlib.sha256(AUDIT_PATH.read_bytes()).hexdigest() == AUDIT_SHA256,
        "pair-hostile audit changed",
    )
    spec = importlib.util.spec_from_file_location("j5_refined_sample_audit", AUDIT_PATH)
    require(spec is not None and spec.loader is not None, "cannot load audit")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load_audit()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def collect_sample() -> list[dict[str, object]]:
    selected = []
    for body in A.P1.BODIES:
        for scalar_row in A.scalar_failure_body(body):
            pair_row = A.pair_threshold_profile(scalar_row)
            if (
                not pair_row["crossing"]
                and pair_row["pair_cap"] < F(5, 7) * pair_row["m"]
                and not pair_row["two_pair_direct"]
            ):
                selected.append(pair_row)
                if len(selected) == SAMPLE_SIZE:
                    return selected
    raise RuntimeError("sample universe exhausted")


def coarse_profile(pair_row: dict[str, object]) -> dict[str, object]:
    body = pair_row["body"]
    apex = pair_row["apex"]
    carrier = A.reconstruct_carrier(body, (apex,))
    mass = pair_row["m"]
    theta = mass - pair_row["pair_cap"]
    require(theta > 2 * mass / 7, f"nonfinite sample threshold: {body}, {apex}")
    level = theta / 2
    threshold = A.S2 * pair_row["r"] / (7 * (level - mass / 7))
    tail_first = A.ceiling(threshold)
    speeds = [
        speed
        for speed in range(A.FIRST_EXTERNAL, tail_first)
        if speed != apex
    ]
    rows = A.P1.THM2885.coverages_many(carrier, speeds)
    head = tuple(sorted(speed for value, speed in rows if value >= level))
    coarse_rows = []
    for x, y in combinations(head, 2):
        residual = A.P1.THM2883.subtract_local_multi(carrier, (x, y))
        union = mass - sum((right - left for left, right in residual), F(0))
        if union >= theta:
            coarse_rows.append(
                (
                    body,
                    pair_row["apex_rank"],
                    apex,
                    theta,
                    x,
                    y,
                    union,
                )
            )
    return {
        "body": body,
        "apex_rank": pair_row["apex_rank"],
        "apex": apex,
        "m": mass,
        "theta": theta,
        "head_tail_first": tail_first,
        "head": head,
        "coarse_rows": tuple(coarse_rows),
    }


def refine_edge(
    task: tuple[tuple[int, ...], int, int, F, int, int, F],
) -> dict[str, object]:
    body, apex_rank, apex, theta, x, y, union = task
    carrier = A.reconstruct_carrier(body, (apex,))
    residual = A.P1.THM2883.subtract_local_multi(carrier, (x, y))
    residual_mass = sum((right - left for left, right in residual), F(0))
    require(residual_mass > 0, f"sample coarse pair covers: {body}, {apex}, {x}, {y}")
    cap = A.global_pair_cap(body, (apex, x, y), residual)
    margin = residual_mass - cap["cap"]
    return {
        "body": body,
        "apex_rank": apex_rank,
        "apex": apex,
        "edge": (x, y),
        "union": union,
        "residual_mass": residual_mass,
        "residual_cap": cap["cap"],
        "margin": margin,
        "maximizer": cap["finite_pair"],
        "survives": margin <= 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    sample = collect_sample()
    profiles = [coarse_profile(row) for row in sample]
    tasks = [task for profile in profiles for task in profile["coarse_rows"]]
    context = mp.get_context("spawn")
    if args.workers == 1:
        refined = list(map(refine_edge, tasks))
    else:
        with context.Pool(args.workers) as pool:
            refined = pool.map(refine_edge, tasks, chunksize=1)

    by_root: dict[tuple[tuple[int, ...], int], list[dict[str, object]]] = {}
    for row in refined:
        by_root.setdefault((row["body"], row["apex"]), []).append(row)
    ledger_rows = []
    summaries = []
    for profile in profiles:
        rows = by_root.get((profile["body"], profile["apex"]), [])
        require(
            tuple(row["edge"] for row in rows)
            == tuple((task[4], task[5]) for task in profile["coarse_rows"]),
            f"sample parallel order changed: {profile['body']}, {profile['apex']}",
        )
        surviving = tuple(row["edge"] for row in rows if row["survives"])
        summaries.append(
            (
                profile["body"],
                profile["apex_rank"],
                profile["apex"],
                profile["head_tail_first"],
                len(profile["head"]),
                len(rows),
                len(surviving),
                surviving[0] if surviving else None,
            )
        )
        ledger_rows.extend(
            "E="
            + ",".join(map(str, row["body"]))
            + f";rank={row['apex_rank']};a={row['apex']};"
            + f"edge={row['edge'][0]},{row['edge'][1]};"
            + f"U={ftext(row['union'])};L={ftext(row['residual_mass'])};"
            + f"B2={ftext(row['residual_cap'])};"
            + f"margin={ftext(row['margin'])};"
            + f"maxpair={row['maximizer']};survives={int(row['survives'])}\n"
            for row in rows
        )
    sample_keys = tuple(
        (row["body"], row["apex_rank"], row["apex"])
        for row in sample
    )
    sample_digest = hashlib.sha256(
        (
            "LRC14/j5/refined-residual-graph-sample-keys/v1\n"
            + "".join(
                ",".join(map(str, body)) + f";{rank};{apex}\n"
                for body, rank, apex in sample_keys
            )
        ).encode()
    ).hexdigest()
    ledger_digest = hashlib.sha256(
        (
            "LRC14/j5/refined-residual-graph-sample-edges/v1\n"
            + "".join(ledger_rows)
        ).encode()
    ).hexdigest()
    nonempty = [summary for summary in summaries if summary[6]]
    observed_counts = (
        sum(summary[4] for summary in summaries),
        max(summary[4] for summary in summaries),
        len(tasks),
        max(summary[5] for summary in summaries),
        sum(summary[6] for summary in summaries),
        len(nonempty),
    )
    require(
        (sample_keys[0], sample_keys[-1]) == EXPECTED_SAMPLE_BOUNDARY
        and sample_digest == EXPECTED_SAMPLE_DIGEST,
        "deterministic sample keys changed",
    )
    require(
        observed_counts == EXPECTED_COUNTS
        and ledger_digest == EXPECTED_EDGE_DIGEST,
        "refined residual sample census changed",
    )
    print("LRC14 J=5 REFINED RESIDUAL GRAPH 50-CARRIER SAMPLE")
    print("status=FINITE-EXACT-SCOUT;deterministic_lex_sample;not_all_carriers")
    print(
        f"sample={len(sample)};"
        f"first={sample_keys[0]};last={sample_keys[-1]};"
        f"sample_digest={sample_digest}"
    )
    print(
        f"head_total={observed_counts[0]};"
        f"head_max={observed_counts[1]};"
        f"coarse_edges={observed_counts[2]};"
        f"coarse_max={observed_counts[3]};"
        f"refined_edges={observed_counts[4]};"
        f"nonempty_carriers={observed_counts[5]};"
        f"first_nonempty={nonempty[0] if nonempty else None};"
        f"edge_digest={ledger_digest}"
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
