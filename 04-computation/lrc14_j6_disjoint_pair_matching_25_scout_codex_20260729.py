#!/usr/bin/env python3
r"""Exact maximum-weight 2-matching scout on 25 marked j=6 suffixes.

On a scalar-hard five-slot carrier of mass ``h``, let ``q5`` be the fifth
allowed singleton order statistic and let ``B2`` be a global pair-union cap.
Put

    T = h-q5,                 L = T-B2.

If two vertex-disjoint pair edges have union weights summing to at least
``T``, each edge has weight at least ``L``.  When

    L > q1+h/7,

the strict discrepancy estimate makes the ``L``-heavy pair graph finite.
This verifier constructs that graph exactly and tests whether its maximum
weight two-edge matching reaches ``T``.

The universe is the 25 scalar-hard marked suffixes in the four-root
THM-2895 battery.  This is a scoped structural scout, not a uniform
seven-body theorem or LRC(14).
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUFFIX_PATH = (
    ROOT
    / "04-computation/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.py"
)
SUFFIX_SHA256 = (
    "6434f020c5aa4000ac81fa081881d93ac0b4190516f854fbd9d8493475baf539"
)
PAIR_PATH = (
    ROOT
    / "04-computation/lrc14_j6_h4_pair_residual_exact_kernel_codex_20260729.py"
)
PAIR_SHA256 = (
    "b82f318bf89ffd3ab4c918c87736461d068e03f25941aa25a0961d0f74b4d70a"
)

# Filled after discovery and then locked.
EXPECTED_COUNTS: tuple[int, ...] | None = (
    25,
    13,
    3,
    9,
    6,
    59,
    242,
    434,
    13,
)
EXPECTED_ROOT_COUNTS: tuple[tuple[object, ...], ...] | None = (
    ((2, 8, 9, 10, 11, 13, 14), 6, 4, 1, 2, 1),
    ((1, 3, 9, 10, 11, 12, 14), 5, 3, 2, 3, 1),
    ((2, 5, 9, 11, 12, 13, 14), 8, 2, 0, 1, 1),
    ((2, 3, 4, 5, 6, 7, 8), 6, 4, 0, 3, 3),
)
EXPECTED_EXTREMA: tuple[object, ...] | None = (
    (
        (2, 8, 9, 10, 11, 13, 14),
        6,
        29,
        "574489/206158680",
    ),
    (
        (2, 3, 4, 5, 6, 7, 8),
        5,
        36,
        "241/360360",
        (20, 44),
        (33, 78),
    ),
    (
        (2, 8, 9, 10, 11, 13, 14),
        6,
        29,
        2611,
    ),
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "ccd5ac755132dc574925345018f25833aae54535fcae2b47a204203d663df129"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path, expected_sha256: str, name: str):
    require(file_sha256(path) == expected_sha256, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load(SUFFIX_PATH, SUFFIX_SHA256, "j6_matching_suffix")
H = load(PAIR_PATH, PAIR_SHA256, "j6_matching_pair")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def profile_hard_branch(
    root: dict[str, object],
    branch: dict[str, object],
) -> dict[str, object]:
    """Build and solve the finite heavy-pair matching problem."""

    carrier = H.R.subtract_local(root["good"], branch["apex"])
    require(
        H.mass(carrier) == branch["m"]
        and len(carrier) == branch["r"],
        f"carrier mismatch: {branch['body']}, rank {branch['rank']}",
    )
    excluded = set(branch["excluded_prefix"])
    pair = H.pair_cap(carrier, excluded)
    q1 = pair["q1"]
    q5 = branch["top5"][4][0]
    require(
        pair["ranked"][4][0] == q5,
        f"q5 seal mismatch: {branch['body']}, rank {branch['rank']}",
    )
    b2 = pair["cap"]
    target = branch["m"] - q5
    level = target - b2
    gap = level - q1 - branch["m"] / 7
    base_closed = 2 * b2 < target

    heavy_edges: list[tuple[F, tuple[int, int]]] = []
    heavy_paid = 0
    direct_controls = 0
    tail_first = None
    best_matching = None
    if gap > 0:
        gamma = H.S2 * branch["r"] / 7
        tail_first = max(H.HORIZON + 1, H.ceiling(gamma / gap))
        rows = list(pair["ranked"])
        if tail_first > H.HORIZON + 1:
            rows.extend(
                H.T.coverages_many(
                    carrier,
                    [
                        speed
                        for speed in range(H.HORIZON + 1, tail_first)
                        if speed not in excluded
                    ],
                )
            )
        require(
            q1
            + branch["m"] / 7
            + H.S2 * branch["r"] / (7 * tail_first)
            <= level,
            f"heavy-pair tail did not seal: {branch['body']}, "
            f"rank {branch['rank']}",
        )
        rows.sort(key=lambda item: (-item[0], item[1]))
        first_control_done = False
        for first_index, (first_value, first) in enumerate(rows[:-1]):
            if first_value + rows[first_index + 1][0] < level:
                break
            after_first = H.R.subtract_local(carrier, first)
            for second_value, second in rows[first_index + 1 :]:
                if first_value + second_value < level:
                    break
                after_pair = H.R.subtract_local(after_first, second)
                union = branch["m"] - H.mass(after_pair)
                heavy_paid += 1
                if not first_control_done:
                    reverse = H.R.subtract_local(
                        H.R.subtract_local(carrier, second),
                        first,
                    )
                    family = tuple(
                        sorted(
                            (
                                *branch["body"],
                                branch["apex"],
                                first,
                                second,
                            )
                        )
                    )
                    direct, components, direct_mass = H.T.CORE.good_norm(
                        family
                    )
                    require(
                        after_pair == reverse == direct
                        and len(after_pair) == components
                        and H.mass(after_pair) == direct_mass,
                        f"heavy-edge reconstruction mismatch: {family}",
                    )
                    first_control_done = True
                    direct_controls += 1
                if union >= level:
                    heavy_edges.append(
                        (union, tuple(sorted((first, second))))
                    )
        require(
            first_control_done,
            f"empty finite heavy-pair candidate set: {branch['body']}, "
            f"rank {branch['rank']}",
        )
        for index, (first_weight, first_edge) in enumerate(heavy_edges):
            for second_weight, second_edge in heavy_edges[index + 1 :]:
                if set(first_edge).isdisjoint(second_edge):
                    record = (
                        first_weight + second_weight,
                        first_edge,
                        second_edge,
                    )
                    if best_matching is None or record > best_matching:
                        best_matching = record

    matching_closed = base_closed or (
        gap > 0
        and (
            best_matching is None
            or best_matching[0] < target
        )
    )
    return {
        "body": branch["body"],
        "rank": branch["rank"],
        "apex": branch["apex"],
        "prefix": branch["excluded_prefix"],
        "h": branch["m"],
        "q1": q1,
        "q5": q5,
        "B2": b2,
        "target": target,
        "level": level,
        "gap": gap,
        "eligible": gap > 0,
        "tail_first": tail_first,
        "base_closed": base_closed,
        "matching_closed": matching_closed,
        "heavy_edges": tuple(heavy_edges),
        "best_matching": best_matching,
        "pair_paid": pair["paid"],
        "heavy_paid": heavy_paid,
        "direct_controls": direct_controls,
    }


def row_line(row: dict[str, object]) -> str:
    heavy_digest = hashlib.sha256(
        (
            "LRC14/j6/disjoint-pair/heavy-edges/v1\n"
            + "".join(
                f"{edge[0]},{edge[1]}:{ftext(weight)}\n"
                for weight, edge in row["heavy_edges"]
            )
        ).encode()
    ).hexdigest()
    if row["best_matching"] is None:
        best = "none"
    else:
        best = (
            f"{ftext(row['best_matching'][0])}:"
            f"{row['best_matching'][1]}:{row['best_matching'][2]}"
        )
    return (
        f"E={row['body']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};"
        f"q1={ftext(row['q1'])};q5={ftext(row['q5'])};"
        f"B2={ftext(row['B2'])};T={ftext(row['target'])};"
        f"L={ftext(row['level'])};gap={ftext(row['gap'])};"
        f"eligible={int(row['eligible'])};"
        f"tail_first={row['tail_first']};"
        f"base={int(row['base_closed'])};"
        f"matching={int(row['matching_closed'])};"
        f"HE={len(row['heavy_edges'])};best={best};"
        f"pair_paid={row['pair_paid']};"
        f"heavy_paid={row['heavy_paid']};"
        f"direct={row['direct_controls']};"
        f"heavy_sha={heavy_digest}"
    )


def main() -> None:
    rows: list[dict[str, object]] = []
    for body in S.BATTERY_ROOTS:
        root = S.global_root(body)
        for rank in range(1, root["K"] + 1):
            branch = S.suffix_branch(root, rank)
            if not branch["closed"]:
                rows.append(profile_hard_branch(root, branch))
    require(len(rows) == 25, "scalar-hard branch universe changed")

    root_counts = tuple(
        (
            body,
            sum(row["body"] == body for row in rows),
            sum(
                row["body"] == body and row["eligible"]
                for row in rows
            ),
            sum(
                row["body"] == body and row["base_closed"]
                for row in rows
            ),
            sum(
                row["body"] == body and row["matching_closed"]
                for row in rows
            ),
            sum(
                row["body"] == body
                and not row["base_closed"]
                and row["matching_closed"]
                for row in rows
            ),
        )
        for body in S.BATTERY_ROOTS
    )
    counts = (
        len(rows),
        sum(row["eligible"] for row in rows),
        sum(row["base_closed"] for row in rows),
        sum(row["matching_closed"] for row in rows),
        sum(
            not row["base_closed"] and row["matching_closed"]
            for row in rows
        ),
        sum(len(row["heavy_edges"]) for row in rows),
        sum(row["pair_paid"] for row in rows),
        sum(row["heavy_paid"] for row in rows),
        sum(row["direct_controls"] for row in rows),
    )
    minimum_gap = min(
        (row for row in rows if row["eligible"]),
        key=lambda row: (row["gap"], row["body"], row["rank"]),
    )
    closest_tested_failure = min(
        (
            row
            for row in rows
            if row["eligible"] and not row["matching_closed"]
        ),
        key=lambda row: (
            row["best_matching"][0] - row["target"],
            row["body"],
            row["rank"],
        ),
    )
    maximum_tail = max(
        (row for row in rows if row["eligible"]),
        key=lambda row: (
            row["tail_first"],
            tuple(-x for x in row["body"]),
            -row["rank"],
        ),
    )
    extrema = (
        (
            minimum_gap["body"],
            minimum_gap["rank"],
            minimum_gap["apex"],
            ftext(minimum_gap["gap"]),
        ),
        (
            closest_tested_failure["body"],
            closest_tested_failure["rank"],
            closest_tested_failure["apex"],
            ftext(
                closest_tested_failure["best_matching"][0]
                - closest_tested_failure["target"]
            ),
            closest_tested_failure["best_matching"][1],
            closest_tested_failure["best_matching"][2],
        ),
        (
            maximum_tail["body"],
            maximum_tail["rank"],
            maximum_tail["apex"],
            maximum_tail["tail_first"],
        ),
    )
    lines = [row_line(row) for row in rows]
    ledger_digest = hashlib.sha256(
        (
            "LRC14/j6/disjoint-pair-matching-25/v1\n"
            + "\n".join(lines)
            + "\n"
        ).encode()
    ).hexdigest()

    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, f"counts changed: {counts}")
    if EXPECTED_ROOT_COUNTS is not None:
        require(
            root_counts == EXPECTED_ROOT_COUNTS,
            f"root counts changed: {root_counts}",
        )
    if EXPECTED_EXTREMA is not None:
        require(extrema == EXPECTED_EXTREMA, f"extrema changed: {extrema}")
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_DIGEST,
            "complete ledger digest changed",
        )

    print("LRC14 j6 disjoint pair matching 25-branch scout")
    print(f"counts={counts}")
    print(f"root_counts={root_counts}")
    print(f"extrema={extrema}")
    for line in lines:
        print(line)
    print(f"complete_ledger_sha256={ledger_digest}")
    if EXPECTED_COUNTS is None:
        print("mode=DISCOVERY")
    else:
        print("mode=LOCKED")
    print("graph=intrinsic undirected pair-union graph;not a tournament")
    print("scope=25 scalar-hard branches on four roots;not uniform;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
