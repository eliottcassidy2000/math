#!/usr/bin/env python3
"""Tighten the 179 inherited-nonfinite H4 pair links by exact child q1.

The parent singleton cap q1 is valid on every literal pair residual but can
be much too large.  This scout identifies the H4 pairs for which

    5*h_L/7 - q1_parent <= 0

and computes the attained global singleton maximum on the literal residual
using an exact finite head plus the THM-735 discrepancy tail.  It then
retests binary-link finiteness with q1_child.

Discovery only: a nonpositive result or a failed scalar three-cover test is
not a cover witness.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
H4_PATH = (
    ROOT
    / "04-computation/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.py"
)
H4_SHA256 = "63a80908a6380a877345f0cc4aba7a5e0ef2bb3d59b1b10d58367444ed406b75"
INITIAL_STOP = 512


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_h4():
    require(file_sha256(H4_PATH) == H4_SHA256, "H4 membership engine changed")
    spec = importlib.util.spec_from_file_location("h4_child_q1", H4_PATH)
    require(spec is not None and spec.loader is not None, "cannot load H4 engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H4 = load_h4()
E = H4.E


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def allowed_labels(stop_exclusive: int, forbidden: frozenset[int]) -> list[int]:
    return [
        label
        for label in range(E.FIRST_EXTERNAL, stop_exclusive)
        if label not in forbidden
    ]


def exact_singleton_cap(
    residual: list[tuple[F, F]],
    residual_mass: F,
    forbidden: frozenset[int],
) -> tuple[F, int, int, int]:
    components = len(residual)
    gamma = E.S2 * components / 7
    initial_labels = allowed_labels(INITIAL_STOP, forbidden)
    require(initial_labels, "empty initial singleton head")
    initial_rows = E.T.coverages_many(residual, initial_labels)
    best, witness = max(initial_rows, key=lambda item: (item[0], -item[1]))
    gap = best - residual_mass / 7
    require(gap > 0, "initial singleton head does not beat the mean floor")
    cutoff = E.ceiling(gamma / gap)

    if cutoff > INITIAL_STOP:
        rows = initial_rows + E.T.coverages_many(
            residual,
            allowed_labels(cutoff, forbidden)[len(initial_labels) :],
        )
        best, witness = max(rows, key=lambda item: (item[0], -item[1]))
    else:
        rows = initial_rows

    final_gap = best - residual_mass / 7
    final_cutoff = E.ceiling(gamma / final_gap)
    require(final_cutoff <= max(cutoff, INITIAL_STOP), "singleton head did not seal")
    require(
        E.R.coverage(residual, witness) == best,
        "singleton vector/scalar witness mismatch",
    )
    # Every allowed w>=final_cutoff has strict discrepancy coverage <best.
    return best, witness, final_cutoff, len(rows)


def main() -> None:
    rows = [H4.exact_row(fields) for fields in H4.parse_exceptions()]
    totals = Counter()
    semantic = hashlib.sha256(
        b"LRC14/j6/paircap-exception/H4-nonpositive-child-q1/v1\n"
    )
    closest_positive = None
    closest_nonpositive = None
    max_scan = None

    for row in rows:
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

        for pair in combinations(core, 2):
            residual_mass = mass - pair_ledger.edge(*pair)
            inherited_delta = 5 * residual_mass / 7 - row["q1"]
            if inherited_delta > 0:
                continue
            totals["input_nonpositive"] += 1
            residual = E.R.subtract_local_multi(carrier, pair)
            require(
                E.interval_mass(residual) == residual_mass,
                "literal pair residual mass mismatch",
            )
            forbidden = frozenset((*row["prefix"], *pair))
            q1_child, witness, cutoff, scanned = exact_singleton_cap(
                residual,
                residual_mass,
                forbidden,
            )
            child_delta = 5 * residual_mass / 7 - q1_child
            scalar_margin = residual_mass - 3 * q1_child
            totals["child_finite"] += child_delta > 0
            totals["child_equality"] += child_delta == 0
            totals["scalar_closed"] += scalar_margin > 0
            totals["scanned"] += scanned
            record = (
                abs(child_delta),
                row["body"],
                row["rank"],
                pair,
                residual_mass,
                q1_child,
                witness,
                cutoff,
            )
            if child_delta > 0:
                if closest_positive is None or record < closest_positive:
                    closest_positive = record
            elif closest_nonpositive is None or record < closest_nonpositive:
                closest_nonpositive = record
            scan_record = (
                cutoff,
                scanned,
                row["body"],
                row["rank"],
                pair,
                witness,
            )
            if max_scan is None or scan_record > max_scan:
                max_scan = scan_record
            semantic.update(
                (
                    f"{row['body']};{row['rank']};{pair};"
                    f"{ftext(residual_mass)};{ftext(row['q1'])};"
                    f"{ftext(q1_child)};{witness};{cutoff};{scanned};"
                    f"{ftext(child_delta)};{ftext(scalar_margin)}\n"
                ).encode()
            )

    require(totals["input_nonpositive"] == 179, "nonpositive pair count changed")
    print("LRC14 j6 H4 nonpositive-link exact child-q1 scout")
    print(f"counts={tuple(sorted(totals.items()))}")
    print(f"closest_positive={closest_positive}")
    print(f"closest_nonpositive={closest_nonpositive}")
    print(f"max_scan={max_scan}")
    print(f"semantic_digest={semantic.hexdigest()}")
    print(
        "scope=179 inherited-nonpositive H4 pair links;"
        "exact child singleton cap and finiteness retest only;"
        "failed tests are not cover witnesses"
    )


if __name__ == "__main__":
    main()
