#!/usr/bin/env python3
r"""Exact four-root battery for the ranked first-apex suffix refinement.

THM-2893 observes that a globally ranked hitting prefix admits a lossless
partition.  If ``a_r`` is the earliest prefix label in a putative six-cover,
then after subtracting ``D_(a_r)`` all five remaining labels lie in the
strict root-order suffix: every ``a_1,...,a_r`` is forbidden.

This verifier exercises that refinement on four seven-body roots: the three
hostile roots from the earlier j=6 complement-cap battery and one low-stratum
control.  It reconstructs each global top-thirty root order, finds its least
one-hit gate ``K``, and profiles every one of the resulting 73 suffix
branches.  On each branch it:

* reconstructs the literal first-apex residual in two independent ways;
* excludes the complete root prefix through the selected apex;
* computes the exact globally sealed top five allowed residual coverages;
* closes the branch when their sum is strictly below the residual mass; and
* checks four distinct vector entries by the scalar primitive.

The strict THM-735(ii) discrepancy bound seals all infinite tails.  This is
a scoped feasibility battery, not a uniform seven-body result or LRC(14).
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM2885_PATH = (
    ROOT
    / "04-computation/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.py"
)
THM2885_SHA256 = (
    "dff97f67b1104c25589802a6a2f216b6e7bfedd58eebfa1bcce615d59c1e872f"
)
THM2885_OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_thm2885_eight_body_top15_hitting_gate_codex_20260729.out"
)
THM2885_OUTPUT_SHA256 = (
    "21a89f15fb144c406936ff62eaf039c0643e36d82ed99a9d28495181fa13e402"
)
THM2883_PATH = (
    ROOT
    / "04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py"
)
THM2883_SHA256 = (
    "a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e"
)

FIRST_EXTERNAL = 15
BASE_HORIZON = 1600
ROOT_TOP_COUNT = 30
ROOT_COVER_SIZE = 6
S2 = F(99, 70)

BATTERY_ROOTS = (
    (2, 8, 9, 10, 11, 13, 14),
    (1, 3, 9, 10, 11, 12, 14),
    (2, 5, 9, 11, 12, 13, 14),
    (2, 3, 4, 5, 6, 7, 8),
)
NONMONOTONE_BODY = (1, 3, 9, 10, 11, 12, 14)

# Filled after discovery, then locked for ordinary and optimized replay.
EXPECTED_ROOT_COUNTS: tuple[tuple[tuple[int, ...], int, int, int], ...] | None = (
    (
        ((2, 8, 9, 10, 11, 13, 14), 19, 13, 6),
        ((1, 3, 9, 10, 11, 12, 14), 20, 15, 5),
        ((2, 5, 9, 11, 12, 13, 14), 21, 13, 8),
        ((2, 3, 4, 5, 6, 7, 8), 13, 7, 6),
    )
)
EXPECTED_TOTAL_COUNTS: tuple[int, ...] | None = (
    4,
    73,
    48,
    25,
    16,
    292,
    73,
    1601,
    2769,
)
EXPECTED_EXTREMA: tuple[object, ...] | None = (
    (
        "42457/149557408",
        (1, 3, 9, 10, 11, 12, 14),
        9,
        78,
    ),
    (
        "-101/105105",
        (2, 3, 4, 5, 6, 7, 8),
        4,
        52,
    ),
    (
        "-2021405909/38818159380",
        (2, 8, 9, 10, 11, 13, 14),
        1,
        19,
    ),
    (
        "228324096/168781",
        (2, 8, 9, 10, 11, 13, 14),
    ),
    (
        "1300178880/469673",
        (2, 8, 9, 10, 11, 13, 14),
        13,
        168,
    ),
)
EXPECTED_NONMONOTONE: tuple[object, ...] | None = (
    (
        5,
        16,
        "5683263/4164400240",
        (15, 34, 46, 26, 38),
    ),
    (
        6,
        46,
        "-2298437/841260420",
        (26, 32, 58, 120, 78),
    ),
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "6510bc4bc977f4c6815fe5024e6f032c0935cf2ce6ec6f9bfdbd17ad6ed91ac7"
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
    file_sha256(THM2885_OUTPUT_PATH) == THM2885_OUTPUT_SHA256,
    "THM-2885 transcript changed",
)
T = load_module(THM2885_PATH, THM2885_SHA256, "j6_suffix_thm2885")
R = load_module(THM2883_PATH, THM2883_SHA256, "j6_suffix_thm2883")


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
    """Select four distinct scanned speeds, including rank and tail controls."""

    candidates = (
        ranked[0][1],
        ranked[4][1],
        tail_first - 1,
        BASE_HORIZON,
        FIRST_EXTERNAL,
        FIRST_EXTERNAL + 1,
        ranked[len(ranked) // 2][1],
    )
    controls = tuple(
        speed
        for speed in dict.fromkeys(candidates)
        if speed in by_speed
    )[:4]
    require(len(controls) == 4, "could not select four distinct controls")
    return controls


def global_root(body: tuple[int, ...]) -> dict[str, object]:
    """Globally seal the root top thirty and find the least one-hit gate."""

    good, components, mass = T.CORE.good_norm(body)
    rows = T.coverages_many(
        good,
        range(FIRST_EXTERNAL, BASE_HORIZON + 1),
    )
    base_ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q30 = base_ranked[ROOT_TOP_COUNT - 1][0]
    require(q30 > mass / 7, f"root q30 misses discrepancy limit: {body}")
    threshold = S2 * components / (7 * (q30 - mass / 7))
    tail_first = max(BASE_HORIZON + 1, ceiling(threshold))
    if tail_first > BASE_HORIZON + 1:
        rows.extend(
            T.coverages_many(
                good,
                range(BASE_HORIZON + 1, tail_first),
            )
        )
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q30,
        f"root tail did not seal: {body}",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top = tuple(ranked[:ROOT_TOP_COUNT])

    adaptive_k: int | None = None
    gate_margin: F | None = None
    for k in range(ROOT_TOP_COUNT - ROOT_COVER_SIZE + 1):
        margin = mass - sum(
            (value for value, _ in top[k : k + ROOT_COVER_SIZE]),
            F(0),
        )
        if margin > 0:
            adaptive_k = k
            gate_margin = margin
            break
    require(
        adaptive_k is not None and gate_margin is not None,
        f"root top30 has no gate: {body}",
    )
    by_speed = {speed: value for value, speed in rows}
    controls = select_controls(ranked, by_speed, tail_first)
    for speed in controls:
        require(
            by_speed[speed] == T.coverage(good, speed),
            f"root vector/scalar mismatch: {body}, {speed}",
        )
    return {
        "body": body,
        "good": good,
        "m": mass,
        "r": components,
        "top": top,
        "K": adaptive_k,
        "gate_margin": gate_margin,
        "threshold": threshold,
        "tail_first": tail_first,
        "controls": len(controls),
    }


def suffix_branch(root: dict[str, object], rank: int) -> dict[str, object]:
    """Profile one unique-earliest-apex branch on its strict root suffix."""

    body = root["body"]
    apex = root["top"][rank - 1][1]
    excluded_prefix = tuple(
        speed for _, speed in root["top"][:rank]
    )
    excluded = set(excluded_prefix)
    literal = R.subtract_local(root["good"], apex)
    direct, direct_r, direct_m = T.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    mass = interval_mass(literal)
    require(
        literal == direct
        and len(literal) == direct_r
        and mass == direct_m > 0,
        f"literal/direct first-apex mismatch: {body}, rank {rank}",
    )

    speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
        if speed not in excluded
    ]
    rows = T.coverages_many(literal, speeds)
    base_ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    q5 = base_ranked[4][0]
    require(
        q5 > mass / 7,
        f"suffix q5 misses discrepancy limit: {body}, rank {rank}",
    )
    threshold = S2 * direct_r / (7 * (q5 - mass / 7))
    tail_first = max(BASE_HORIZON + 1, ceiling(threshold))
    if tail_first > BASE_HORIZON + 1:
        rows.extend(
            T.coverages_many(
                literal,
                [
                    speed
                    for speed in range(BASE_HORIZON + 1, tail_first)
                    if speed not in excluded
                ],
            )
        )
    require(
        mass / 7 + S2 * direct_r / (7 * tail_first) <= q5,
        f"suffix tail did not seal: {body}, rank {rank}",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top5 = tuple(ranked[:5])
    require(
        all(speed not in excluded for _, speed in top5),
        f"forbidden prefix label survived: {body}, rank {rank}",
    )
    margin = mass - sum((value for value, _ in top5), F(0))

    by_speed = {speed: value for value, speed in rows}
    controls = select_controls(ranked, by_speed, tail_first)
    for speed in controls:
        require(
            by_speed[speed] == T.coverage(literal, speed),
            f"suffix vector/scalar mismatch: {body}, rank {rank}, {speed}",
        )
    return {
        "body": body,
        "rank": rank,
        "apex": apex,
        "excluded_prefix": excluded_prefix,
        "m": mass,
        "r": direct_r,
        "top5": top5,
        "margin": margin,
        "closed": margin > 0,
        "threshold": threshold,
        "tail_first": tail_first,
        "controls": len(controls),
        "direct_controls": 1,
    }


def branch_line(row: dict[str, object]) -> str:
    return (
        f"E={','.join(map(str, row['body']))};rank={row['rank']};"
        f"apex={row['apex']};"
        f"prefix={','.join(map(str, row['excluded_prefix']))};"
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
    require(
        len(BATTERY_ROOTS) == 4 and len(set(BATTERY_ROOTS)) == 4,
        "suffix battery universe changed",
    )
    roots = [global_root(body) for body in BATTERY_ROOTS]
    branches = [
        suffix_branch(root, rank)
        for root in roots
        for rank in range(1, root["K"] + 1)
    ]
    root_counts = tuple(
        (
            root["body"],
            root["K"],
            sum(
                row["closed"]
                for row in branches
                if row["body"] == root["body"]
            ),
            sum(
                not row["closed"]
                for row in branches
                if row["body"] == root["body"]
            ),
        )
        for root in roots
    )
    closed = [row for row in branches if row["closed"]]
    open_rows = [row for row in branches if not row["closed"]]
    minimum_positive = min(
        closed,
        key=lambda row: (row["margin"], row["body"], row["rank"]),
    )
    closest_open = max(
        open_rows,
        key=lambda row: (row["margin"], tuple(-x for x in row["body"]), -row["rank"]),
    )
    worst_open = min(
        open_rows,
        key=lambda row: (row["margin"], row["body"], row["rank"]),
    )
    maximum_root_threshold = max(
        roots,
        key=lambda row: (row["threshold"], tuple(-x for x in row["body"])),
    )
    maximum_suffix_threshold = max(
        branches,
        key=lambda row: (
            row["threshold"],
            tuple(-x for x in row["body"]),
            -row["rank"],
        ),
    )
    total_counts = (
        len(roots),
        len(branches),
        len(closed),
        len(open_rows),
        sum(root["controls"] for root in roots),
        sum(row["controls"] for row in branches),
        sum(row["direct_controls"] for row in branches),
        max(root["tail_first"] for root in roots),
        max(row["tail_first"] for row in branches),
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
            ftext(maximum_root_threshold["threshold"]),
            maximum_root_threshold["body"],
        ),
        (
            ftext(maximum_suffix_threshold["threshold"]),
            maximum_suffix_threshold["body"],
            maximum_suffix_threshold["rank"],
            maximum_suffix_threshold["apex"],
        ),
    )
    nonmonotone_rows = {
        row["rank"]: row
        for row in branches
        if row["body"] == NONMONOTONE_BODY and row["rank"] in (5, 6)
    }
    require(
        set(nonmonotone_rows) == {5, 6}
        and nonmonotone_rows[5]["closed"]
        and not nonmonotone_rows[6]["closed"],
        "rank-five/rank-six nonmonotonicity changed",
    )
    nonmonotone = tuple(
        (
            rank,
            nonmonotone_rows[rank]["apex"],
            ftext(nonmonotone_rows[rank]["margin"]),
            tuple(speed for _, speed in nonmonotone_rows[rank]["top5"]),
        )
        for rank in (5, 6)
    )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/rank-first-suffix-scalar-battery/v1\n")
    for root in roots:
        digest.update(
            (
                f"ROOT={','.join(map(str, root['body']))};K={root['K']};"
                f"m={ftext(root['m'])};r={root['r']};"
                f"gate_margin={ftext(root['gate_margin'])};"
                f"threshold={ftext(root['threshold'])};"
                f"tail_first={root['tail_first']};top30="
                + ",".join(
                    f"{speed}:{ftext(value)}"
                    for value, speed in root["top"]
                )
                + "\n"
            ).encode()
        )
    for row in branches:
        digest.update(branch_line(row).encode())
    ledger_digest = digest.hexdigest()

    if EXPECTED_ROOT_COUNTS is not None:
        require(root_counts == EXPECTED_ROOT_COUNTS, "root counts changed")
    if EXPECTED_TOTAL_COUNTS is not None:
        require(total_counts == EXPECTED_TOTAL_COUNTS, "total counts changed")
    if EXPECTED_EXTREMA is not None:
        require(extrema == EXPECTED_EXTREMA, "extrema changed")
    if EXPECTED_NONMONOTONE is not None:
        require(
            nonmonotone == EXPECTED_NONMONOTONE,
            "rank-five/rank-six witness changed",
        )
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_DIGEST,
            "complete ledger digest changed",
        )

    print("LRC14 j6 ranked first-apex suffix scalar battery")
    print(f"root_counts={root_counts}")
    print(f"total_counts={total_counts}")
    print(f"extrema={extrema}")
    print(f"rank5_rank6_nonmonotone={nonmonotone}")
    for body, gate_k, closed_count, open_count in root_counts:
        body_rows = [row for row in branches if row["body"] == body]
        print(
            f"E={body};K={gate_k};closed={closed_count};open={open_count};"
            f"closed_ranks={[row['rank'] for row in body_rows if row['closed']]};"
            f"open_ranks={[row['rank'] for row in body_rows if not row['closed']]}"
        )
    print(f"complete_ledger_sha256={ledger_digest}")
    if EXPECTED_TOTAL_COUNTS is None:
        print("mode=DISCOVERY (expected constants are not yet locked)")
    else:
        print("mode=LOCKED")
    print("scope=4/3432 roots;73 suffix branches;not uniform;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
