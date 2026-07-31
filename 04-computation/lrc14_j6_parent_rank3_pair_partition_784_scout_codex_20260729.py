#!/usr/bin/env python3
"""Scout the inherited-parent THM-2900 certificate on its 784-row control.

THM-2900 proves the literal-child certificate ``q3(R_L)+B2(R_L)<|R_L|``.
Its parent-carrier form instead asks

    U_C(L) + q3^C(P union L) + B2^C(P union L) < |C|.

The parent form preserves one carrier across every forced H4 pair ``L`` and
is therefore the potentially scalable route for the all-root census.  This
script compares the two sufficient certificates on the already proved
four-root universe.  It is a discovery scout until its constants are locked.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE = (
    ROOT
    / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.py"
)
ENGINE_SHA256 = (
    "9023a4042dc8def3f8781668e721549972fb1458d07f2fab89bf7ac3a6f745cc"
)
EXPECTED_COUNTS: tuple[int, ...] | None = (25, 784, 204, 580, 784, 7310)
EXPECTED_ROOT_COUNTS: tuple[tuple[tuple[int, ...], int, int], ...] | None = (
    ((2, 8, 9, 10, 11, 13, 14), 225, 48),
    ((1, 3, 9, 10, 11, 12, 14), 118, 46),
    ((2, 5, 9, 11, 12, 13, 14), 351, 83),
    ((2, 3, 4, 5, 6, 7, 8), 90, 27),
)
EXPECTED_EXTREMA: tuple[tuple[object, ...], ...] | None = (
    (
        "-135564641/3577774200",
        (2, 5, 9, 11, 12, 13, 14),
        1,
        16,
        (20, 25),
    ),
    (
        "-17489/1873980108",
        (2, 5, 9, 11, 12, 13, 14),
        3,
        17,
        (20, 38),
    ),
)
EXPECTED_LEDGER_SHA256: str | None = (
    "952b436fe47f9227de7791c7ee788bbb2b9e00b6214f8170c7d449178c51c298"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_engine():
    require(
        hashlib.sha256(ENGINE.read_bytes()).hexdigest() == ENGINE_SHA256,
        "independent engine changed",
    )
    spec = importlib.util.spec_from_file_location("parent_rank3_engine", ENGINE)
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E = load_engine()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def parent_pair_cap(
    carrier: tuple[tuple[F, F], ...],
    ranked: list[tuple[F, int]],
    tail_single: F,
) -> tuple[F, int]:
    head, _, paid = E.finite_pair_cap(carrier, ranked)
    return max(head, ranked[0][0] + tail_single), paid


def main() -> None:
    roots = [E.root_profile(body) for body in E.ROOTS]
    branches = []
    for root in roots:
        for rank in range(1, root["K"] + 1):
            branch = E.suffix_profile(root, rank)
            if not branch["closed"]:
                branches.append(E.h4_profile(branch))
    require(len(branches) == 25, "open branch universe changed")

    rows = []
    paid = 0
    for branch in branches:
        ranked, tail_single = E.sealed_ranking(
            branch["carrier"],
            branch["excluded"],
            E.PAIR_HORIZON,
            5,
        )
        for flag in combinations(branch["H"], 2):
            allowed = [
                item for item in ranked if item[1] not in set(flag)
            ]
            require(len(allowed) >= 3, "parent rank-three list exhausted")
            q3 = allowed[2][0]
            b2, paid_here = parent_pair_cap(
                branch["carrier"],
                allowed,
                tail_single,
            )
            paid += paid_here
            residual = E.subtract_many(branch["carrier"], flag)
            union_flag = branch["m"] - E.interval_mass(residual)
            parent_margin = branch["m"] - union_flag - q3 - b2

            child = E.pair_profile(branch, flag)
            child_margin = child["m"] - child["top3"][2][0] - child["B2"]
            rows.append(
                (
                    branch,
                    flag,
                    parent_margin,
                    child_margin,
                    union_flag,
                    q3,
                    b2,
                )
            )

    require(len(rows) == 784, "pair universe changed")
    parent_closed = [row for row in rows if row[2] > 0]
    child_closed = [row for row in rows if row[3] > 0]
    require(len(child_closed) == len(rows), "THM-2900 child control failed")
    sharp_parent = min(
        rows,
        key=lambda row: (
            row[2],
            row[0]["body"],
            row[0]["rank"],
            row[1],
        ),
    )
    closest_parent_failure = max(
        (row for row in rows if row[2] <= 0),
        key=lambda row: (
            row[2],
            tuple(-value for value in row[0]["body"]),
            -row[0]["rank"],
            tuple(-value for value in row[1]),
        ),
        default=None,
    )
    root_counts = tuple(
        (
            body,
            sum(row[0]["body"] == body for row in rows),
            sum(row[0]["body"] == body and row[2] > 0 for row in rows),
        )
        for body in E.ROOTS
    )
    counts = (
        len(branches),
        len(rows),
        len(parent_closed),
        len(rows) - len(parent_closed),
        len(child_closed),
        paid,
    )
    extrema = (
        (
            ftext(sharp_parent[2]),
            sharp_parent[0]["body"],
            sharp_parent[0]["rank"],
            sharp_parent[0]["apex"],
            sharp_parent[1],
        ),
        (
            ftext(closest_parent_failure[2]),
            closest_parent_failure[0]["body"],
            closest_parent_failure[0]["rank"],
            closest_parent_failure[0]["apex"],
            closest_parent_failure[1],
        ),
    )
    ledger = hashlib.sha256(b"LRC14/j6/parent-q3-plus-B2/784/v1\n")
    for branch, flag, parent_margin, child_margin, union_flag, q3, b2 in rows:
        ledger.update(
            (
                f"E={branch['body']};rank={branch['rank']};"
                f"apex={branch['apex']};L={flag};"
                f"U={ftext(union_flag)};q3={ftext(q3)};B2={ftext(b2)};"
                f"parent={ftext(parent_margin)};"
                f"child={ftext(child_margin)}\n"
            ).encode()
        )
    ledger_digest = ledger.hexdigest()
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "aggregate counts changed")
    if EXPECTED_ROOT_COUNTS is not None:
        require(root_counts == EXPECTED_ROOT_COUNTS, "root counts changed")
    if EXPECTED_EXTREMA is not None:
        require(extrema == EXPECTED_EXTREMA, "extrema changed")
    if EXPECTED_LEDGER_SHA256 is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_SHA256,
            "parent-certificate ledger changed",
        )
    print("LRC14 j6 inherited-parent q3+B2 scout")
    print(
        f"branches={counts[0]};pairs={counts[1]};"
        f"parent_closed={counts[2]};parent_open={counts[3]};"
        f"child_closed={counts[4]};paid_parent_pairs={counts[5]}"
    )
    print(f"root_counts={root_counts}")
    print(
        "minimum_parent_margin="
        f"{ftext(sharp_parent[2])};E={sharp_parent[0]['body']};"
        f"rank={sharp_parent[0]['rank']};apex={sharp_parent[0]['apex']};"
        f"H4pair={sharp_parent[1]}"
    )
    if closest_parent_failure is not None:
        print(
            "closest_parent_failure="
            f"{ftext(closest_parent_failure[2])};"
            f"E={closest_parent_failure[0]['body']};"
            f"rank={closest_parent_failure[0]['rank']};"
            f"apex={closest_parent_failure[0]['apex']};"
            f"H4pair={closest_parent_failure[1]}"
        )
    print(f"canonical_ledger_sha256={ledger_digest}")
    print(
        "mode="
        + ("DISCOVERY" if EXPECTED_LEDGER_SHA256 is None else "LOCKED")
    )
    print("scope=784 THM-2900 control rows;not all-root;not LRC14")


if __name__ == "__main__":
    main()
