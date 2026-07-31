#!/usr/bin/env python3
"""Measure parent-cap inheritance before fresh H4-pair residual profiling.

For an apex carrier C and an H4 pair L={x,y}, put

    R = C \ (D_x union D_y).

Then c_R(w) <= c_C(w) and U_R({u,v}) <= U_C({u,v}).  Consequently
the following are sound upper-bound closures of a three-cover of R:

    top3 of c_C outside the inherited exclusions and L;
    B2(C) + B1(C).

The first uses the already available apex singleton ledger.  The second
requires one exact parent pair cap per apex branch, amortized across every
H4 pair.  Only rows surviving both tests are sent to the existing fresh
child top3/B2 computation.

This is a scratch workload probe on the committed four-root battery.
"""

from __future__ import annotations

import argparse
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SUFFIX_PATH = (
    ROOT
    / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895.py"
)
MINSEED_PATH = (
    ROOT
    / "04-computation/lrc14_j6_minseed_parity_closure_thm2895_addendum.py"
)


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


Q = load(SUFFIX_PATH, "inheritance_suffix_h4")
M = load(MINSEED_PATH, "inheritance_minseed_h4")


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def top_values_outside(
    ranked: list[tuple[F, int]],
    excluded: set[int],
    count: int,
) -> tuple[F, ...]:
    values: list[F] = []
    for value, speed in ranked:
        if speed in excluded:
            continue
        values.append(value)
        if len(values) == count:
            break
    if len(values) != count:
        raise RuntimeError("parent singleton ledger has too few allowed labels")
    return tuple(values)


def profile(label: str, branches: list[dict[str, object]]) -> None:
    total_pairs = 0
    inherited_top3 = 0
    inherited_b2b1 = 0
    inherited_union = 0
    parent_cap_branches = 0
    parent_paid_pairs = 0
    fresh_profiles = 0
    fresh_top3 = 0
    fresh_b2b1 = 0
    fresh_union = 0
    fresh_top3_tail_firsts: list[int] = []
    fresh_hard: list[
        tuple[dict[str, object], tuple[int, int], dict[str, object]]
    ] = []

    for row in branches:
        ranked, _ = Q.H.globally_ranked(row["carrier"], row["excluded"])
        pair_rows: list[tuple[tuple[int, int], list[tuple[F, F]], F, F]] = []
        top3_survivors: list[
            tuple[tuple[int, int], list[tuple[F, F]], F, F]
        ] = []
        for hpair in combinations(row["H"], 2):
            residual = Q.H.R.subtract_local_multi(row["carrier"], hpair)
            residual_mass = mass(residual)
            if residual_mass <= 0:
                raise RuntimeError("H4 pair covers the apex carrier")
            parent_top3 = top_values_outside(
                ranked,
                set(hpair),
                3,
            )
            top3_margin = residual_mass - sum(parent_top3, F(0))
            item = (hpair, residual, residual_mass, top3_margin)
            pair_rows.append(item)
            if top3_margin > 0:
                inherited_top3 += 1
            else:
                top3_survivors.append(item)

        total_pairs += len(pair_rows)
        parent_cap = None
        if top3_survivors:
            parent_cap = Q.H.pair_cap(row["carrier"], row["excluded"])
            parent_cap_branches += 1
            parent_paid_pairs += parent_cap["paid"]

        branch_inherited_union = len(pair_rows) - len(top3_survivors)
        child_queue = []
        for item in top3_survivors:
            hpair, residual, residual_mass, _ = item
            if parent_cap is None:
                raise RuntimeError("missing parent pair cap")
            b2b1_margin = (
                residual_mass - parent_cap["cap"] - parent_cap["q1"]
            )
            if b2b1_margin > 0:
                inherited_b2b1 += 1
                branch_inherited_union += 1
            else:
                child_queue.append(item)
        inherited_union += branch_inherited_union

        for hpair, _, _, _ in child_queue:
            fresh = Q.pair_residual(row, hpair)
            fresh_profiles += 1
            q3 = fresh["top3"][2][0]
            threshold = (
                Q.H.S2
                * fresh["r"]
                / (7 * (q3 - fresh["m"] / 7))
            )
            fresh_top3_tail_firsts.append(
                max(Q.H.FIRST_EXTERNAL, Q.ceiling(threshold))
            )
            fresh_top3 += int(fresh["direct_margin"] > 0)
            fresh_b2b1 += int(fresh["pair_margin"] > 0)
            if fresh["adaptive_closed"]:
                fresh_union += 1
            else:
                fresh_hard.append((row, hpair, fresh))

        branch = row["branch"]
        print(
            f"{label}:E={branch['body']};a={branch['apex']};"
            f"prefix={tuple(sorted(row['excluded']))};"
            f"H={len(row['H'])};pairs={len(pair_rows)};"
            f"parent_top3={len(pair_rows)-len(top3_survivors)};"
            f"parent_union={branch_inherited_union};"
            f"fresh={len(child_queue)}"
        )

    recursive_closed = 0
    recursive_checks = 0
    for row, _, fresh in fresh_hard:
        recursive = Q.recursive_k3_close(row, fresh)
        recursive_closed += int(recursive["closed"])
        recursive_checks += recursive["checks"]

    print(
        f"{label}:TOTAL branches={len(branches)};pairs={total_pairs};"
        f"parent_top3={inherited_top3};"
        f"parent_B2B1_extra={inherited_b2b1};"
        f"parent_union={inherited_union};"
        f"parent_cap_branches={parent_cap_branches};"
        f"parent_paid_pairs={parent_paid_pairs};"
        f"fresh_profiles={fresh_profiles};"
        f"fresh_top3={fresh_top3};fresh_B2B1={fresh_b2b1};"
        f"fresh_union={fresh_union};fresh_hard={len(fresh_hard)};"
        f"recursive_closed={recursive_closed};"
        f"recursive_checks={recursive_checks}"
    )
    if fresh_top3_tail_firsts:
        ordered = sorted(fresh_top3_tail_firsts)
        quantiles = tuple(
            ordered[(len(ordered) - 1) * numerator // 4]
            for numerator in range(5)
        )
        print(
            f"{label}:FRESH_TOP3_TAILS min={ordered[0]};"
            f"quartiles={quantiles};max={ordered[-1]};"
            f"sum={sum(ordered)}"
        )
    if recursive_closed != len(fresh_hard):
        raise RuntimeError("inherited pipeline leaves an open recursive row")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--universe",
        choices=("open25", "minseed15", "both"),
        default="both",
    )
    args = parser.parse_args()
    if args.universe in ("open25", "both"):
        profile("OPEN25", Q.open_branch_profiles())
    if args.universe in ("minseed15", "both"):
        profile("MINSEED15", M.seed_branches())


if __name__ == "__main__":
    main()
