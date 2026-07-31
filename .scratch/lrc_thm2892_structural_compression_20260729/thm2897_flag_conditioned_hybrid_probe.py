#!/usr/bin/env python3
"""Exact interaction probe for THM-2897 partition caps and THM-2895 flags.

For every scalar-open branch in the proved four-root THM-2895 battery this
script compares, in increasing cost order:

1. the parent five-slot rank-selective certificate q5+2 B2<h;
2. flag-conditioned parent top-three singleton caps;
3. the actual flag union plus a reusable parent q3(L)+B2 cap;
4. a cap q3(L)+B2_L recomputed on the parent carrier after excluding the
   flag;
5. the fresh child-residual certificates used by THM-2895.

The fourth test is the exact flag-conditioned tropical convolution on the
parent carrier.  It never subtracts a child danger comb when evaluating its
caps; the actual flag union is the only residual datum required.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PRIMARY = (
    ROOT / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895.py"
)
PRIMARY_SHA = "970d77503f8d56d737e223dabb3c3562d7b19cd018ca75398e3deb054715e5f6"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load():
    require(file_sha256(PRIMARY) == PRIMARY_SHA, "THM-2895 source changed")
    spec = importlib.util.spec_from_file_location("thm2897_hybrid_primary", PRIMARY)
    require(spec is not None and spec.loader is not None, "cannot load primary")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


Q = load()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def top_outside(
    ranked: list[tuple[F, int]],
    excluded: set[int],
    count: int,
) -> tuple[tuple[F, int], ...]:
    rows = tuple(
        (value, speed)
        for value, speed in ranked
        if speed not in excluded
    )[:count]
    require(len(rows) == count, "too few conditioned singleton ranks")
    return rows


def main() -> None:
    branches = Q.open_branch_profiles()
    require(len(branches) == 25, "open branch universe changed")
    all_rows: list[tuple[dict[str, object], dict[str, object]]] = []
    conditioned_cap_calls = 0
    conditioned_paid = 0
    fresh_calls = 0
    ledger = hashlib.sha256(b"LRC14/THM2897/flag-conditioned-hybrid/v1\n")

    for branch in branches:
        carrier = branch["carrier"]
        h = mass(carrier)
        excluded = set(branch["excluded"])
        ranked, _ = Q.H.globally_ranked(carrier, excluded)
        parent = Q.H.pair_cap(carrier, excluded)
        parent_direct_margin = h - ranked[4][0] - 2 * parent["cap"]
        local: list[dict[str, object]] = []
        for flag in combinations(branch["H"], 2):
            child = Q.H.R.subtract_local_multi(carrier, flag)
            child_mass = mass(child)
            require(child_mass > 0, "flag covers the parent carrier")
            union = h - child_mass
            conditioned_top = top_outside(ranked, set(flag), 3)
            conditioned_q3 = conditioned_top[2][0]
            top3_margin = child_mass - sum(
                (value for value, _ in conditioned_top),
                F(0),
            )
            q_global_b2_global_margin = (
                child_mass - conditioned_q3 - parent["cap"]
            )

            # Pay the fully flag-conditioned parent B2 cap only if the two
            # reusable ledgers fail.  This is the strongest {1,2} tropical
            # certificate on C with the flag labels deleted.
            conditioned_margin = None
            conditioned = None
            if top3_margin <= 0 and q_global_b2_global_margin <= 0:
                conditioned = Q.H.pair_cap(
                    carrier,
                    excluded | set(flag),
                )
                conditioned_cap_calls += 1
                conditioned_paid += conditioned["paid"]
                conditioned_margin = (
                    child_mass
                    - conditioned["ranked"][2][0]
                    - conditioned["cap"]
                )

            inherited_closed = (
                top3_margin > 0
                or q_global_b2_global_margin > 0
                or (
                    conditioned_margin is not None
                    and conditioned_margin > 0
                )
            )
            fresh = None
            fresh_rank_margin = None
            if not inherited_closed:
                fresh = Q.pair_residual(branch, flag)
                fresh_calls += 1
                fresh_rank_margin = (
                    fresh["m"]
                    - fresh["top3"][2][0]
                    - fresh["cap"]["cap"]
                )
            row = {
                "flag": flag,
                "child_mass": child_mass,
                "union": union,
                "parent_direct_margin": parent_direct_margin,
                "top3_margin": top3_margin,
                "global_margin": q_global_b2_global_margin,
                "conditioned_margin": conditioned_margin,
                "inherited_closed": inherited_closed,
                "fresh_closed": (
                    None if fresh is None else fresh["adaptive_closed"]
                ),
                "fresh_rank_margin": fresh_rank_margin,
            }
            local.append(row)
            all_rows.append((branch, row))
            ledger.update(
                (
                    f"E={branch['branch']['body']};"
                    f"r={branch['branch']['rank']};"
                    f"a={branch['branch']['apex']};L={flag};"
                    f"U={ftext(union)};R={ftext(child_mass)};"
                    f"parent={ftext(parent_direct_margin)};"
                    f"T3={ftext(top3_margin)};"
                    f"GB={ftext(q_global_b2_global_margin)};"
                    f"CB={'NA' if conditioned_margin is None else ftext(conditioned_margin)};"
                    f"I={int(inherited_closed)};"
                    f"F={'NA' if fresh is None else int(fresh['adaptive_closed'])};"
                    f"FR={'NA' if fresh_rank_margin is None else ftext(fresh_rank_margin)}\n"
                ).encode()
            )

        scalar = branch["branch"]
        print(
            f"E={scalar['body']};rank={scalar['rank']};a={scalar['apex']};"
            f"pairs={len(local)};"
            f"parent_direct={int(parent_direct_margin > 0)};"
            f"parent_top3={sum(row['top3_margin'] > 0 for row in local)};"
            f"global_qB_union={sum(row['top3_margin'] > 0 or row['global_margin'] > 0 for row in local)};"
            f"conditioned_union={sum(row['inherited_closed'] for row in local)};"
            f"fresh={sum(not row['inherited_closed'] for row in local)};"
            f"branch_inherited_closed={int(all(row['inherited_closed'] for row in local))}"
        )

    parent_direct_branches = sum(
        (
            rows[0][1]["parent_direct_margin"] > 0
            if rows
            else False
        )
        for rows in (
            [
                item
                for item in all_rows
                if item[0] is branch
            ]
            for branch in branches
        )
    )
    inherited_closed = sum(row["inherited_closed"] for _, row in all_rows)
    inherited_open = len(all_rows) - inherited_closed
    branch_hybrid_closed = sum(
        all(
            row["inherited_closed"]
            for owner, row in all_rows
            if owner is branch
        )
        for branch in branches
    )
    fresh_old_failures = sum(
        row["fresh_closed"] is False for _, row in all_rows
    )
    fresh_rank_failures = sum(
        row["fresh_rank_margin"] is not None
        and row["fresh_rank_margin"] <= 0
        for _, row in all_rows
    )
    print(
        f"TOTAL branches={len(branches)};pairs={len(all_rows)};"
        f"parent_direct_branches={parent_direct_branches};"
        f"inherited_closed_pairs={inherited_closed};"
        f"inherited_open_pairs={inherited_open};"
        f"hybrid_closed_branches={branch_hybrid_closed};"
        f"conditioned_cap_calls={conditioned_cap_calls};"
        f"conditioned_paid={conditioned_paid};"
        f"fresh_calls={fresh_calls};"
        f"fresh_old_failures={fresh_old_failures};"
        f"fresh_rank_failures={fresh_rank_failures}"
    )
    print(f"ledger_sha256={ledger.hexdigest()}")
    print("scope=THM-2895 four-root/25-open-branch battery")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
