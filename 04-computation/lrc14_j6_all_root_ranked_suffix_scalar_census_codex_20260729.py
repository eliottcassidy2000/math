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
EXPECTED_LEDGER_DIGEST: str | None = (
    "6579713e351ab89ea26c8a9e0178c6aacbd0011107d984813cb621e39b7d8ff0"
)
EXPECTED_STRATUM_DIGESTS: tuple[tuple[str, str], ...] | None = (
    (
        "low",
        "f80b53fefc96bbffba5ee6772c443474a74abb52fa8b7c13b138e2d27d9364c7",
    ),
    (
        "one",
        "616cc320e02a96e94b79f14d86efd55593eb06bca4f5237d7aa1dad8c83c3a84",
    ),
    (
        "both",
        "15a0440ec04159adac967812a424f8a6a4be6da26d06ed04a7eac666b60c7894",
    ),
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 1, 8),
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
    digest.update(b"LRC14/j6/all-root-ranked-suffix-scalar/v1\n")
    stratum_hashes = {
        name: hashlib.sha256(
            f"LRC14/j6/all-root-ranked-suffix-scalar/{name}/v1\n".encode()
        )
        for name in strata
    }
    for root in roots:
        root_line = (
            f"ROOT={','.join(map(str, root['body']))};"
            f"S={root['stratum']};K={root['K']};"
            f"gate_margin={ftext(root['gate_margin'])};"
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

    print("LRC14 j6 all-root ranked-suffix scalar census")
    print(f"counts={counts}")
    print(f"stratum_counts={stratum_counts}")
    print(f"open_branches_per_root={open_distribution}")
    print(f"open_rank_distribution={open_rank_distribution}")
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
