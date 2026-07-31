#!/usr/bin/env python3
"""Exact child-singleton tightening on all 18,290 H4 pair flags.

This extends nonpositive_child_q1_scout.py from the 179 inherited-gap
failures to the complete THM-2901 exception-pair universe.  For every literal
pair residual it computes the attained global allowed singleton maximum,
tests the direct three-singleton closure, and records the discrepancy cutoff
for the resulting binary-link high core.

Discovery only.  A surviving flag is a finite child obligation, not a cover
witness.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "nonpositive_child_q1_scout.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_base():
    spec = importlib.util.spec_from_file_location("h4_all_child_q1_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load base scout")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_base()
H4 = B.H4
E = B.E


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_parent(row: dict[str, object]) -> dict[str, object]:
    carrier, components, mass = E.R.CORE.good_norm(
        tuple(sorted((*row["body"], row["apex"])))
    )
    require(
        components == row["components"] and mass == row["mass"],
        "parent carrier reconstruction changed",
    )
    core = tuple(row["core"])
    singleton = {
        label: value
        for value, label in E.T.coverages_many(carrier, list(core))
    }
    pair_ledger = E.PairLedger(carrier, singleton)
    result_rows = []
    for pair in combinations(core, 2):
        residual_mass = mass - pair_ledger.edge(*pair)
        residual = E.R.subtract_local_multi(carrier, pair)
        require(
            E.interval_mass(residual) == residual_mass,
            "literal pair residual mass mismatch",
        )
        forbidden = frozenset((*row["prefix"], *pair))
        q1_child, witness, singleton_cutoff, scanned = B.exact_singleton_cap(
            residual,
            residual_mass,
            forbidden,
        )
        link_delta = 5 * residual_mass / 7 - q1_child
        require(link_delta > 0, "exact child cap left a non-finite binary link")
        beta = link_delta / 2
        gamma = E.S2 * len(residual) / 7
        link_cutoff = E.ceiling(gamma / beta) - 1
        require(
            residual_mass / 7 + gamma / (link_cutoff + 1)
            <= (residual_mass - q1_child) / 2,
            "binary-link tail did not seal",
        )
        scalar_margin = residual_mass - 3 * q1_child
        result_rows.append(
            (
                row["body"],
                row["rank"],
                pair,
                residual_mass,
                q1_child,
                witness,
                singleton_cutoff,
                scanned,
                link_delta,
                link_cutoff,
                scalar_margin,
            )
        )
    return {
        "body": row["body"],
        "rank": row["rank"],
        "rows": tuple(result_rows),
    }


def nearest_rank(values: list[int]) -> tuple[tuple[int, int], ...]:
    ordered = sorted(values)
    require(ordered, "empty quantile population")
    return tuple(
        (
            percentile,
            ordered[
                min(
                    len(ordered) - 1,
                    percentile * (len(ordered) - 1) // 100,
                )
            ],
        )
        for percentile in (0, 25, 50, 75, 90, 95, 99, 100)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8))
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    parents = [H4.exact_row(fields) for fields in H4.parse_exceptions()]
    if args.workers == 1:
        parent_results = list(map(exact_parent, parents))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            parent_results = pool.map(exact_parent, parents)
    rows = [
        row
        for parent in parent_results
        for row in parent["rows"]
    ]
    require(len(rows) == 18_290, "H4 pair universe changed")
    singleton_cutoffs = [row[6] for row in rows]
    scans = [row[7] for row in rows]
    link_cutoffs = [row[9] for row in rows]
    margins = [row[10] for row in rows]
    counts = (
        len(rows),
        sum(margin > 0 for margin in margins),
        sum(margin == 0 for margin in margins),
        sum(scans),
        min(singleton_cutoffs),
        max(singleton_cutoffs),
        min(link_cutoffs),
        max(link_cutoffs),
        tuple(Counter(sum(row[10] <= 0 for row in parent["rows"]) for parent in parent_results).most_common()),
    )
    semantic = hashlib.sha256(
        b"LRC14/j6/paircap-exception/H4-all-child-q1/v1\n"
    )
    for row in rows:
        semantic.update(
            (
                f"{row[0]};{row[1]};{row[2]};{ftext(row[3])};"
                f"{ftext(row[4])};{row[5]};{row[6]};{row[7]};"
                f"{ftext(row[8])};{row[9]};{ftext(row[10])}\n"
            ).encode()
        )
    print("LRC14 j6 all-exception H4 exact child-q1 scout")
    print(f"counts={counts}")
    print(f"singleton_cutoff_quantiles={nearest_rank(singleton_cutoffs)}")
    print(f"link_cutoff_quantiles={nearest_rank(link_cutoffs)}")
    print(f"scan_quantiles={nearest_rank(scans)}")
    print(f"semantic_digest={semantic.hexdigest()}")
    print(
        "scope=all 18,290 H4 pair flags;exact child q1 and finite binary-link "
        "cutoffs;survivors are obligations, not cover witnesses"
    )


if __name__ == "__main__":
    main()
