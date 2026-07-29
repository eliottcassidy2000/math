#!/usr/bin/env python3
r"""Exact ranked-suffix scalar census on all 41,415 j=6 apex branches.

THM-2896 assigns a globally ranked hitting gate of least size ``K(E)<=25``
to each of the 3,432 seven-body roots.  Give every putative six-cover to its
unique earliest gate member ``a_r``.  THM-2893's lossless suffix refinement
then forbids all root-prefix labels ``a_1,...,a_r`` from the five remaining
slots.

For every resulting branch this verifier reconstructs the literal
first-apex residual, independently reconstructs the same direct good set,
computes the five largest allowed individual coverages on the strict suffix,
and globally seals the infinite tail by the strict THM-735(ii) discrepancy
bound.  If their sum is strictly below the residual mass, subadditivity
closes the branch.

This is the complete scalar layer after the THM-2896 gate.  Any surviving
branch still needs complement caps, heavy flags, or another nonscalar
argument.  The script makes no LRC(14) claim unless every branch closes.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GATE_PATH = (
    ROOT
    / "04-computation/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.py"
)
GATE_SHA256 = (
    "fc36f26d4c8da5b005465696b954eec700c080376eef9ee5ba74a7111def99d7"
)
GATE_OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j6_all_root_adaptive_gate_atlas_codex_20260729.out"
)
GATE_OUTPUT_SHA256 = (
    "3081b93a870faacb31d205e43f7ca87872d7a9f196f4774d8740ced6a314d80b"
)
RESIDUAL_PATH = (
    ROOT
    / "04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py"
)
RESIDUAL_SHA256 = (
    "a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e"
)

FIRST_EXTERNAL = 15
BASE_HORIZON = 1600
S2 = F(99, 70)

# Filled after discovery, then locked for ordinary and optimized replay.
EXPECTED_COUNTS: tuple[int, ...] | None = (
    3432,
    41_415,
    26_609,
    14_806,
    5,
    3427,
    912,
    13,
    165_660,
    41_415,
    2923,
    378,
)
EXPECTED_STRATUM_COUNTS: tuple[tuple[str, tuple[int, ...]], ...] | None = (
    ("low", (792, 9447, 6394, 3053, 0, 9)),
    ("one", (1848, 21_986, 14_133, 7853, 4, 11)),
    ("both", (792, 9982, 6082, 3900, 1, 10)),
)
EXPECTED_OPEN_DISTRIBUTION: tuple[tuple[int, int], ...] | None = (
    (0, 5),
    (1, 65),
    (2, 354),
    (3, 710),
    (4, 827),
    (5, 702),
    (6, 450),
    (7, 204),
    (8, 87),
    (9, 21),
    (10, 5),
    (11, 2),
)
EXPECTED_OPEN_RANK_DISTRIBUTION: tuple[tuple[int, int], ...] | None = (
    (1, 3417),
    (2, 3308),
    (3, 2853),
    (4, 2194),
    (5, 1446),
    (6, 844),
    (7, 444),
    (8, 180),
    (9, 78),
    (10, 24),
    (11, 14),
    (12, 3),
    (13, 1),
)
EXPECTED_EXTREMA: tuple[object, ...] | None = (
    ("5557/3988104120", (3, 4, 5, 7, 8, 13, 14), 4, 17),
    ("-17/119819700", (3, 5, 6, 7, 10, 12, 13), 3, 22),
    ("-27077/270270", (2, 4, 6, 8, 10, 12, 14), 1, 22),
    (
        "5028283260/1720429",
        (1, 5, 10, 11, 12, 13, 14),
        15,
        175,
        2923,
    ),
    (378, (4, 8, 10, 11, 12, 13, 14), 16, 2, "55/1323"),
)
EXPECTED_LEDGER_DIGEST: str | None = None
EXPECTED_STRATUM_DIGESTS: tuple[tuple[str, str], ...] | None = None
EXPECTED_TERMINAL_BODIES: tuple[tuple[int, ...], ...] | None = (
    (1, 2, 3, 4, 5, 6, 13),
    (1, 2, 3, 4, 6, 7, 14),
    (1, 2, 3, 4, 6, 11, 13),
    (1, 2, 3, 4, 6, 12, 13),
    (7, 8, 9, 11, 12, 13, 14),
)
EXPECTED_RANK13_OPEN: tuple[object, ...] | None = (
    (
        (1, 5, 8, 9, 11, 12, 13),
        17,
        13,
        42,
        (21, 17, 23, 19, 34, 35, 50, 28, 31, 49, 40, 46, 42),
        "8683/42042",
        34,
        (
            (37, "823/18648"),
            (20, "3713/90090"),
            (29, "4195/102312"),
            (51, "116317/2858856"),
            (30, "4769/120120"),
        ),
        "-2378297/11503321830",
    ),
)
EXPECTED_PARITY_COUNTS: tuple[int, ...] | None = None
EXPECTED_PARITY_STRATUM: tuple[tuple[object, ...], ...] | None = None
EXPECTED_PARITY_RANK: tuple[tuple[object, ...], ...] | None = None
EXPECTED_PARITY_EXTREMUM: tuple[object, ...] | None = None
EXPECTED_ROOT_PARITY_COUNTS: tuple[int, ...] | None = None
EXPECTED_ROOT_PARITY_STRATUM: tuple[tuple[object, ...], ...] | None = None
EXPECTED_ROOT_PARITY_EXTREMA: tuple[tuple[object, ...], ...] | None = None
EXPECTED_ROOT_PARITY_CUTOFF_QUANTILES: tuple[
    tuple[int, int], ...
] | None = None
THM2895_ROOTS = (
    (2, 8, 9, 10, 11, 13, 14),
    (1, 3, 9, 10, 11, 12, 14),
    (2, 5, 9, 11, 12, 13, 14),
    (2, 3, 4, 5, 6, 7, 8),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path: Path, expected_sha256: str, name: str):
    require(file_sha256(path) == expected_sha256, f"{path.name} changed")
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(
    file_sha256(GATE_OUTPUT_PATH) == GATE_OUTPUT_SHA256,
    "THM-2896 transcript changed",
)
G = load_module(GATE_PATH, GATE_SHA256, "j6_suffix_all_gate")
R = load_module(RESIDUAL_PATH, RESIDUAL_SHA256, "j6_suffix_all_residual")


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def select_controls(
    ranked: list[tuple[F, int]],
    by_speed: dict[int, F],
    tail_first: int,
) -> tuple[int, ...]:
    """Select four distinct allowed rank/boundary scalar controls."""

    candidates = (
        ranked[0][1],
        ranked[4][1],
        tail_first - 1,
        BASE_HORIZON,
        FIRST_EXTERNAL,
        FIRST_EXTERNAL + 1,
        ranked[len(ranked) // 2][1],
    )
    # A long excluded prefix can remove several named boundary controls, and
    # coincidences among the remaining named controls can leave fewer than
    # four distinct speeds.  Complete the control family deterministically
    # from the exact ranking; every ranked speed is allowed and present in
    # ``by_speed``.
    controls = tuple(
        speed
        for speed in dict.fromkeys(
            (*candidates, *(speed for _, speed in ranked))
        )
        if speed in by_speed
    )[:4]
    require(len(controls) == 4, "could not choose four scalar controls")
    return controls


def profile_branch(
    root: dict[str, object],
    rank: int,
) -> dict[str, object]:
    """Globally seal one literal suffix branch's residual top five."""

    body = root["body"]
    apex = root["top"][rank - 1][1]
    prefix = tuple(speed for _, speed in root["top"][:rank])
    excluded = set(prefix)

    literal = R.subtract_local(root["good"], apex)
    direct, components, direct_mass = G.T.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    mass = interval_mass(literal)
    require(
        literal == direct
        and len(literal) == components
        and mass == direct_mass > 0,
        f"literal/direct residual mismatch: {body}, rank {rank}, apex {apex}",
    )

    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
        if speed not in excluded
    ]
    rows = G.T.coverages_many(literal, speeds)
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q5_base = ranked[4][0]
    require(
        q5_base > mass / 7,
        f"suffix rank five misses discrepancy limit: {body}, rank {rank}",
    )
    threshold = S2 * components / (7 * (q5_base - mass / 7))
    tail_first = max(BASE_HORIZON + 1, ceiling(threshold))
    if tail_first > BASE_HORIZON + 1:
        rows.extend(
            G.T.coverages_many(
                literal,
                [
                    speed
                    for speed in range(BASE_HORIZON + 1, tail_first)
                    if speed not in excluded
                ],
            )
        )
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q5_base,
        f"suffix rank-five tail did not seal: {body}, rank {rank}",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top5 = tuple(ranked[:5])
    require(
        all(speed not in excluded for _, speed in top5),
        f"excluded root prefix survived: {body}, rank {rank}",
    )
    margin = mass - sum((value for value, _ in top5), F(0))

    by_speed = {speed: value for value, speed in rows}
    controls = select_controls(ranked, by_speed, tail_first)
    for speed in controls:
        require(
            by_speed[speed] == G.T.coverage(literal, speed),
            f"suffix vector/scalar mismatch: {body}, rank {rank}, {speed}",
        )
    return {
        "body": body,
        "stratum": root["stratum"],
        "K": root["adaptive_k"],
        "rank": rank,
        "apex": apex,
        "prefix": prefix,
        "m": mass,
        "r": components,
        "top5": top5,
        "margin": margin,
        "closed": margin > 0,
        "threshold": threshold,
        "tail_first": tail_first,
        "controls": len(controls),
        "direct_controls": 1,
    }


def profile_root(body: tuple[int, ...]) -> dict[str, object]:
    """Pickle-stable worker for one THM-2896 root and all its suffix branches."""

    root = G.profile_body(body)
    good, components, mass = G.T.CORE.good_norm(body)
    require(
        components == root["r"] and mass == root["m"],
        f"root reconstruction changed: {body}",
    )
    root["good"] = good
    root_q1, root_q1_speed = root["top"][0]
    root_parity_margin = F(2, 7) * mass - root_q1
    root_parity_cutoff = (
        ceiling(5 * S2 * components / (7 * root_parity_margin)) - 1
        if root_parity_margin > 0
        else None
    )
    branches = tuple(
        profile_branch(root, rank)
        for rank in range(1, root["adaptive_k"] + 1)
    )
    open_ranks = tuple(row["rank"] for row in branches if not row["closed"])
    closed_ranks = tuple(row["rank"] for row in branches if row["closed"])
    return {
        "body": body,
        "stratum": root["stratum"],
        "K": root["adaptive_k"],
        "gate_margin": root["adaptive_margin"],
        "m": mass,
        "r": components,
        "q1": root_q1,
        "q1_speed": root_q1_speed,
        "root_parity_margin": root_parity_margin,
        "root_parity_cutoff": root_parity_cutoff,
        "branches": branches,
        "open_ranks": open_ranks,
        "closed_ranks": closed_ranks,
        "terminal": not open_ranks,
    }


def branch_line(row: dict[str, object]) -> str:
    return (
        f"S={row['stratum']};E={','.join(map(str, row['body']))};"
        f"K={row['K']};rank={row['rank']};apex={row['apex']};"
        f"prefix={','.join(map(str, row['prefix']))};"
        f"m={ftext(row['m'])};r={row['r']};"
        f"margin={ftext(row['margin'])};closed={int(row['closed'])};"
        f"threshold={ftext(row['threshold'])};"
        f"tail_first={row['tail_first']};top5="
        + ",".join(
            f"{speed}:{ftext(value)}" for value, speed in row["top5"]
        )
        + "\n"
    )


def parity_group_summary(
    rows: list[dict[str, object]],
    key: object,
) -> tuple[object, ...]:
    """Summarize the THM-2895 p=5 singleton-gate margin on a row group."""

    require(rows, f"empty parity group: {key}")
    decorated = [
        (F(3, 7) * row["m"] - row["top5"][0][0], row)
        for row in rows
    ]
    minimum_margin, minimum_row = min(
        decorated,
        key=lambda item: (
            item[0],
            item[1]["body"],
            item[1]["rank"],
        ),
    )
    return (
        key,
        len(rows),
        sum(margin > 0 for margin, _ in decorated),
        sum(margin <= 0 for margin, _ in decorated),
        ftext(minimum_margin),
        minimum_row["body"],
        minimum_row["rank"],
        minimum_row["apex"],
    )


def root_parity_group_summary(
    roots: list[dict[str, object]],
    key: object,
) -> tuple[object, ...]:
    """Summarize the THM-2895 p=6 singleton gate on a root group."""

    require(roots, f"empty root parity group: {key}")
    minimum = min(
        roots,
        key=lambda row: (
            row["root_parity_margin"],
            row["body"],
        ),
    )
    return (
        key,
        len(roots),
        sum(row["root_parity_margin"] > 0 for row in roots),
        sum(row["root_parity_margin"] <= 0 for row in roots),
        ftext(minimum["root_parity_margin"]),
        minimum["body"],
        minimum["q1_speed"],
    )


def nearest_rank_quantiles(
    values: list[int],
    percentages: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    """Return deterministic nearest-rank order statistics."""

    require(values, "empty quantile list")
    ordered = sorted(values)
    size = len(ordered)
    return tuple(
        (
            percentage,
            ordered[
                0
                if percentage == 0
                else ceiling(F(percentage * size, 100)) - 1
            ],
        )
        for percentage in percentages
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 8),
    )
    parser.add_argument(
        "--emit-hard-ledger",
        action="store_true",
        help="emit every scalar-open row with its exact top-five order data",
    )
    args = parser.parse_args()
    require(args.workers >= 1, "workers must be positive")
    require(
        len(G.BODIES) == 3432 and len(set(G.BODIES)) == 3432,
        "seven-body universe changed",
    )

    if args.workers == 1:
        roots = list(map(profile_root, G.BODIES))
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            roots = pool.map(profile_root, G.BODIES, chunksize=1)
    require(
        tuple(row["body"] for row in roots) == G.BODIES,
        "root worker order changed",
    )
    branches = [
        branch for root in roots for branch in root["branches"]
    ]
    require(
        len(branches) == sum(root["K"] for root in roots),
        "branch universe changed",
    )
    closed = [row for row in branches if row["closed"]]
    open_rows = [row for row in branches if not row["closed"]]

    strata = ("low", "one", "both")
    stratum_counts = tuple(
        (
            name,
            (
                sum(root["stratum"] == name for root in roots),
                sum(
                    root["K"]
                    for root in roots
                    if root["stratum"] == name
                ),
                sum(row["stratum"] == name and row["closed"] for row in branches),
                sum(
                    row["stratum"] == name and not row["closed"]
                    for row in branches
                ),
                sum(
                    root["stratum"] == name and root["terminal"]
                    for root in roots
                ),
                max(
                    (
                        len(root["open_ranks"])
                        for root in roots
                        if root["stratum"] == name
                    ),
                    default=0,
                ),
            ),
        )
        for name in strata
    )
    open_distribution = tuple(
        sorted(Counter(len(root["open_ranks"]) for root in roots).items())
    )
    open_rank_distribution = tuple(
        sorted(Counter(row["rank"] for row in open_rows).items())
    )
    terminal_bodies = tuple(
        root["body"] for root in roots if root["terminal"]
    )
    rank13_open = tuple(
        (
            row["body"],
            row["K"],
            row["rank"],
            row["apex"],
            row["prefix"],
            ftext(row["m"]),
            row["r"],
            tuple(
                (speed, ftext(value)) for value, speed in row["top5"]
            ),
            ftext(row["margin"]),
        )
        for row in open_rows
        if row["rank"] == 13
    )
    parity_rows = [
        (
            F(3, 7) * row["m"] - row["top5"][0][0],
            row,
        )
        for row in open_rows
    ]
    parity_counts = (
        len(parity_rows),
        sum(margin > 0 for margin, _ in parity_rows),
        sum(margin <= 0 for margin, _ in parity_rows),
    )
    parity_stratum = tuple(
        parity_group_summary(
            [row for row in open_rows if row["stratum"] == name],
            name,
        )
        for name in ("low", "one", "both")
    )
    parity_rank = tuple(
        parity_group_summary(
            [row for row in open_rows if row["rank"] == rank],
            rank,
        )
        for rank in sorted({row["rank"] for row in open_rows})
    )
    parity_minimum_margin, parity_minimum_row = min(
        parity_rows,
        key=lambda item: (
            item[0],
            item[1]["body"],
            item[1]["rank"],
        ),
    )
    parity_extremum = (
        ftext(parity_minimum_margin),
        parity_minimum_row["body"],
        parity_minimum_row["rank"],
        parity_minimum_row["apex"],
        ftext(parity_minimum_row["m"]),
        (
            parity_minimum_row["top5"][0][1],
            ftext(parity_minimum_row["top5"][0][0]),
        ),
    )
    root_parity_counts = (
        len(roots),
        sum(root["root_parity_margin"] > 0 for root in roots),
        sum(root["root_parity_margin"] <= 0 for root in roots),
    )
    root_parity_stratum = tuple(
        root_parity_group_summary(
            [root for root in roots if root["stratum"] == name],
            name,
        )
        for name in ("low", "one", "both")
    )
    root_parity_positive = [
        root for root in roots if root["root_parity_margin"] > 0
    ]
    root_parity_failures = [
        root for root in roots if root["root_parity_margin"] <= 0
    ]
    minimum_root_parity_positive = min(
        root_parity_positive,
        key=lambda row: (
            row["root_parity_margin"],
            row["body"],
        ),
    )
    closest_root_parity_failure = max(
        root_parity_failures,
        key=lambda row: (
            row["root_parity_margin"],
            tuple(-x for x in row["body"]),
        ),
    )
    worst_root_parity_failure = min(
        root_parity_failures,
        key=lambda row: (
            row["root_parity_margin"],
            row["body"],
        ),
    )
    root_parity_extrema = tuple(
        (
            ftext(root["root_parity_margin"]),
            root["body"],
            root["q1_speed"],
            ftext(root["m"]),
            root["r"],
            root["root_parity_cutoff"],
        )
        for root in (
            minimum_root_parity_positive,
            closest_root_parity_failure,
            worst_root_parity_failure,
        )
    )
    root_parity_cutoff_quantiles = nearest_rank_quantiles(
        [
            root["root_parity_cutoff"]
            for root in root_parity_positive
        ],
        (0, 25, 50, 75, 90, 95, 99, 100),
    )
    nonmonotone_roots = sum(
        bool(root["closed_ranks"])
        and bool(root["open_ranks"])
        and min(root["closed_ranks"]) < max(root["open_ranks"])
        for root in roots
    )

    minimum_positive = min(
        closed,
        key=lambda row: (row["margin"], row["body"], row["rank"]),
    )
    closest_open = max(
        open_rows,
        key=lambda row: (
            row["margin"],
            tuple(-x for x in row["body"]),
            -row["rank"],
        ),
    )
    worst_open = min(
        open_rows,
        key=lambda row: (row["margin"], row["body"], row["rank"]),
    )
    maximum_threshold = max(
        branches,
        key=lambda row: (
            row["threshold"],
            tuple(-x for x in row["body"]),
            -row["rank"],
        ),
    )
    maximum_retained = max(
        (
            speed,
            row["body"],
            row["rank"],
            top_rank,
            value,
        )
        for row in branches
        for top_rank, (value, speed) in enumerate(row["top5"], start=1)
    )
    extrema = (
        (
            ftext(minimum_positive["margin"]),
            minimum_positive["body"],
            minimum_positive["rank"],
            minimum_positive["apex"],
        ),
        (
            ftext(closest_open["margin"]),
            closest_open["body"],
            closest_open["rank"],
            closest_open["apex"],
        ),
        (
            ftext(worst_open["margin"]),
            worst_open["body"],
            worst_open["rank"],
            worst_open["apex"],
        ),
        (
            ftext(maximum_threshold["threshold"]),
            maximum_threshold["body"],
            maximum_threshold["rank"],
            maximum_threshold["apex"],
            maximum_threshold["tail_first"],
        ),
        (
            maximum_retained[0],
            maximum_retained[1],
            maximum_retained[2],
            maximum_retained[3],
            ftext(maximum_retained[4]),
        ),
    )
    counts = (
        len(roots),
        len(branches),
        len(closed),
        len(open_rows),
        sum(root["terminal"] for root in roots),
        len(roots) - sum(root["terminal"] for root in roots),
        nonmonotone_roots,
        max((row["rank"] for row in open_rows), default=0),
        sum(row["controls"] for row in branches),
        sum(row["direct_controls"] for row in branches),
        max(row["tail_first"] for row in branches),
        maximum_retained[0],
    )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/all-root-ranked-suffix-scalar/v2\n")
    stratum_hashes = {
        name: hashlib.sha256(
            f"LRC14/j6/all-root-ranked-suffix-scalar/{name}/v2\n".encode()
        )
        for name in strata
    }
    for root in roots:
        root_line = (
            f"ROOT={','.join(map(str, root['body']))};"
            f"S={root['stratum']};K={root['K']};"
            f"gate_margin={ftext(root['gate_margin'])};"
            f"m={ftext(root['m'])};r={root['r']};"
            f"q1={root['q1_speed']}:{ftext(root['q1'])};"
            f"p6_margin={ftext(root['root_parity_margin'])};"
            f"p6_cutoff={root['root_parity_cutoff']};"
            f"closed_ranks={root['closed_ranks']};"
            f"open_ranks={root['open_ranks']}\n"
        ).encode()
        digest.update(root_line)
        stratum_hashes[root["stratum"]].update(root_line)
        for row in root["branches"]:
            encoded = branch_line(row).encode()
            digest.update(encoded)
            stratum_hashes[row["stratum"]].update(encoded)
    ledger_digest = digest.hexdigest()
    stratum_digests = tuple(
        (name, stratum_hashes[name].hexdigest()) for name in strata
    )

    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "global counts changed")
    if EXPECTED_STRATUM_COUNTS is not None:
        require(
            stratum_counts == EXPECTED_STRATUM_COUNTS,
            "stratum counts changed",
        )
    if EXPECTED_OPEN_DISTRIBUTION is not None:
        require(
            open_distribution == EXPECTED_OPEN_DISTRIBUTION,
            "open-root distribution changed",
        )
    if EXPECTED_OPEN_RANK_DISTRIBUTION is not None:
        require(
            open_rank_distribution == EXPECTED_OPEN_RANK_DISTRIBUTION,
            "open-rank distribution changed",
        )
    if EXPECTED_EXTREMA is not None:
        require(extrema == EXPECTED_EXTREMA, "extrema changed")
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_DIGEST,
            "complete ledger digest changed",
        )
    if EXPECTED_STRATUM_DIGESTS is not None:
        require(
            stratum_digests == EXPECTED_STRATUM_DIGESTS,
            "stratum digests changed",
        )
    if EXPECTED_TERMINAL_BODIES is not None:
        require(
            terminal_bodies == EXPECTED_TERMINAL_BODIES,
            "terminal root bodies changed",
        )
    if EXPECTED_RANK13_OPEN is not None:
        require(
            rank13_open == EXPECTED_RANK13_OPEN,
            "rank-thirteen survivor changed",
        )
    if EXPECTED_PARITY_COUNTS is not None:
        require(
            parity_counts == EXPECTED_PARITY_COUNTS,
            "parity-eligibility counts changed",
        )
    if EXPECTED_PARITY_STRATUM is not None:
        require(
            parity_stratum == EXPECTED_PARITY_STRATUM,
            "parity-eligibility stratum census changed",
        )
    if EXPECTED_PARITY_RANK is not None:
        require(
            parity_rank == EXPECTED_PARITY_RANK,
            "parity-eligibility rank census changed",
        )
    if EXPECTED_PARITY_EXTREMUM is not None:
        require(
            parity_extremum == EXPECTED_PARITY_EXTREMUM,
            "parity-eligibility extremum changed",
        )
    if EXPECTED_ROOT_PARITY_COUNTS is not None:
        require(
            root_parity_counts == EXPECTED_ROOT_PARITY_COUNTS,
            "root parity-eligibility counts changed",
        )
    if EXPECTED_ROOT_PARITY_STRATUM is not None:
        require(
            root_parity_stratum == EXPECTED_ROOT_PARITY_STRATUM,
            "root parity-eligibility stratum census changed",
        )
    if EXPECTED_ROOT_PARITY_EXTREMA is not None:
        require(
            root_parity_extrema == EXPECTED_ROOT_PARITY_EXTREMA,
            "root parity-eligibility extrema changed",
        )
    if EXPECTED_ROOT_PARITY_CUTOFF_QUANTILES is not None:
        require(
            root_parity_cutoff_quantiles
            == EXPECTED_ROOT_PARITY_CUTOFF_QUANTILES,
            "root parity-cutoff quantiles changed",
        )
    require(
        set(terminal_bodies).isdisjoint(THM2895_ROOTS),
        "scalar terminal roots overlap the four THM-2895 roots",
    )

    print("LRC14 j6 all-root ranked-suffix scalar census")
    if args.emit_hard_ledger:
        for row in open_rows:
            print("HARD;" + branch_line(row).rstrip())
    print(f"counts={counts}")
    print(f"stratum_counts={stratum_counts}")
    print(f"open_branches_per_root={open_distribution}")
    print(f"open_rank_distribution={open_rank_distribution}")
    print(f"terminal_bodies={terminal_bodies}")
    print(f"rank13_open={rank13_open}")
    print("terminal_bodies_disjoint_from_THM2895=PASS")
    print(f"parity_eligibility_counts={parity_counts}")
    print(f"parity_eligibility_stratum={parity_stratum}")
    print(f"parity_eligibility_rank={parity_rank}")
    print(f"parity_eligibility_extremum={parity_extremum}")
    print(f"root_parity_eligibility_counts={root_parity_counts}")
    print(f"root_parity_eligibility_stratum={root_parity_stratum}")
    print(f"root_parity_eligibility_extrema={root_parity_extrema}")
    print(
        "root_parity_cutoff_quantiles_nearest_rank="
        f"{root_parity_cutoff_quantiles}"
    )
    print(f"extrema={extrema}")
    print(f"complete_ledger_sha256={ledger_digest}")
    print(f"stratum_ledger_sha256={stratum_digests}")
    if EXPECTED_COUNTS is None:
        print("mode=DISCOVERY (expected constants are not yet locked)")
    else:
        print("mode=LOCKED")
    print("scope=all 41415 THM2896 branches;scalar suffix layer only;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
