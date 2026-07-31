#!/usr/bin/env python3
"""Independent exact audit of the full fixed-order four-root j=6 result (THM-2895).

This program deliberately imports no repository LRC implementation.  It
reconstructs danger sets as rational interval unions, evaluates all comb
coverages from a fresh periodic antiderivative, and independently rebuilds:

* the four root rankings and all 73 unique-earliest-apex suffix branches;
* the 25 scalar-open branches and their finite H4 cores;
* all 784 H4-pair residuals, their global top-three and B2+q1 certificates;
* the five rows left by the adaptive union; and
* the recursive H2/heavy-edge/longest-component singleton closure.

The intended scope is exactly the fixed four-root battery, not all seven-body
roots and not LRC(14).
"""

from __future__ import annotations

import hashlib
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import lcm

import numpy as np


FIRST_EXTERNAL = 15
ROOT_HORIZON = 1600
PAIR_HORIZON = 2500
ROOT_TOP = 30
ROOT_COVER_SIZE = 6
S2 = F(99, 70)

ROOTS = (
    (2, 8, 9, 10, 11, 13, 14),
    (1, 3, 9, 10, 11, 12, 14),
    (2, 5, 9, 11, 12, 13, 14),
    (2, 3, 4, 5, 6, 7, 8),
)

EXPECTED_ROOT_COUNTS = (
    ((2, 8, 9, 10, 11, 13, 14), 19, 13, 6),
    ((1, 3, 9, 10, 11, 12, 14), 20, 15, 5),
    ((2, 5, 9, 11, 12, 13, 14), 21, 13, 8),
    ((2, 3, 4, 5, 6, 7, 8), 13, 7, 6),
)
EXPECTED_PAIR_TOTALS = (784, 771, 773, 779, 5, 7551)
EXPECTED_STRATA = (
    ("rank1", 361, 348, 350, 356, 5),
    ("nonrank1", 423, 423, 423, 423, 0),
)
EXPECTED_RECURSIVE = (
    (
        (2, 8, 9, 10, 11, 13, 14),
        1,
        19,
        (37, 108),
        (17, 23, 46),
        ((17, 23),),
        ("121/34776",),
        (41,),
        23,
    ),
    (
        (2, 8, 9, 10, 11, 13, 14),
        1,
        19,
        (37, 125),
        (17, 23, 46),
        ((17, 23),),
        ("13/3500",),
        (38,),
        20,
    ),
    (
        (2, 8, 9, 10, 11, 13, 14),
        1,
        19,
        (108, 125),
        (17, 23, 46),
        ((17, 23),),
        ("125/28728",),
        (32,),
        15,
    ),
    (
        (2, 5, 9, 11, 12, 13, 14),
        1,
        16,
        (23, 25),
        (20, 37),
        ((20, 37),),
        ("33/7252",),
        (31,),
        13,
    ),
    (
        (2, 5, 9, 11, 12, 13, 14),
        1,
        16,
        (25, 34),
        (20, 23),
        ((20, 23),),
        ("19/4508",),
        (33,),
        15,
    ),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(intervals: tuple[tuple[F, F], ...]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def merge(
    intervals: list[tuple[F, F]] | tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    out: list[list[F]] = []
    for left, right in sorted(intervals):
        require(F(0) <= left <= right <= F(1), "interval left [0,1]")
        if left == right:
            continue
        if out and left <= out[-1][1]:
            out[-1][1] = max(out[-1][1], right)
        else:
            out.append([left, right])
    return tuple((left, right) for left, right in out)


@lru_cache(maxsize=None)
def danger(speed: int) -> tuple[tuple[F, F], ...]:
    """The closed danger comb D_speed on the unit circle, cut at zero."""

    require(speed >= 1, "nonpositive speed")
    radius = F(1, 14 * speed)
    pieces: list[tuple[F, F]] = []
    for residue in range(speed):
        center = F(residue, speed)
        left = center - radius
        right = center + radius
        if left < 0:
            pieces.append((F(0), right))
            pieces.append((left + 1, F(1)))
        elif right > 1:
            pieces.append((left, F(1)))
            pieces.append((F(0), right - 1))
        else:
            pieces.append((left, right))
    return merge(pieces)


def subtract_one(
    carrier: tuple[tuple[F, F], ...],
    speed: int,
) -> tuple[tuple[F, F], ...]:
    """Subtract one danger comb by a two-pointer interval sweep."""

    bad = danger(speed)
    out: list[tuple[F, F]] = []
    bad_start = 0
    for left, right in carrier:
        cursor = left
        while bad_start < len(bad) and bad[bad_start][1] <= left:
            bad_start += 1
        index = bad_start
        while index < len(bad) and bad[index][0] < right:
            bad_left, bad_right = bad[index]
            if cursor < bad_left:
                out.append((cursor, min(right, bad_left)))
            cursor = max(cursor, bad_right)
            if cursor >= right:
                break
            index += 1
        if cursor < right:
            out.append((cursor, right))
    return tuple(out)


def subtract_many(
    carrier: tuple[tuple[F, F], ...],
    speeds: tuple[int, ...],
) -> tuple[tuple[F, F], ...]:
    out = carrier
    for speed in speeds:
        out = subtract_one(out, speed)
    return out


def good_set(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    return subtract_many(((F(0), F(1)),), speeds)


def tooth_antiderivative(value: F) -> F:
    integer = value.numerator // value.denominator
    residue = value - integer
    return (
        F(integer, 7)
        + min(residue, F(1, 14))
        + max(F(0), residue - F(13, 14))
    )


def scalar_coverage(
    carrier: tuple[tuple[F, F], ...],
    speed: int,
) -> F:
    return sum(
        (
            tooth_antiderivative(speed * right)
            - tooth_antiderivative(speed * left)
        )
        / speed
        for left, right in carrier
    )


def coverages(
    carrier: tuple[tuple[F, F], ...],
    speeds: list[int] | range,
) -> list[tuple[F, int]]:
    """Fresh exact vector engine for the periodic antiderivative."""

    speed_list = list(speeds)
    if not speed_list:
        return []
    require(carrier, "empty carrier passed to coverage engine")
    denominator = lcm(
        *(endpoint.denominator for piece in carrier for endpoint in piece)
    )
    max_speed = max(speed_list)
    require(
        14 * denominator * max_speed < 2**63,
        "int64 exact primitive guard failed",
    )
    endpoints = np.array(
        [
            endpoint.numerator * (denominator // endpoint.denominator)
            for piece in carrier
            for endpoint in piece
        ],
        dtype=np.int64,
    )
    out: list[tuple[F, int]] = []
    for start in range(0, len(speed_list), 5000):
        batch = speed_list[start : start + 5000]
        axes = np.array(batch, dtype=np.int64)
        scaled = endpoints[:, None] * axes[None, :]
        quotient, remainder = np.divmod(scaled, denominator)
        primitive = (
            2 * denominator * quotient
            + np.minimum(14 * remainder, denominator)
            + np.maximum(0, 14 * remainder - 13 * denominator)
        )
        numerators = (primitive[1::2] - primitive[0::2]).sum(axis=0)
        out.extend(
            (F(int(numerator), 14 * denominator * speed), speed)
            for numerator, speed in zip(numerators, batch)
        )
    return out


def discrepancy_tail(
    carrier: tuple[tuple[F, F], ...],
    first_speed: int,
) -> F:
    return (
        interval_mass(carrier) / 7
        + S2 * len(carrier) / (7 * first_speed)
    )


def sealed_ranking(
    carrier: tuple[tuple[F, F], ...],
    excluded: set[int],
    horizon: int,
    needed_rank: int,
) -> tuple[list[tuple[F, int]], F]:
    rows = coverages(
        carrier,
        [
            speed
            for speed in range(FIRST_EXTERNAL, horizon + 1)
            if speed not in excluded
        ],
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    tail = discrepancy_tail(carrier, horizon + 1)
    require(
        len(ranked) >= needed_rank and tail <= ranked[needed_rank - 1][0],
        "finite ranking did not seal",
    )
    controls = tuple(
        dict.fromkeys(
            (
                ranked[0][1],
                ranked[needed_rank - 1][1],
                FIRST_EXTERNAL
                if FIRST_EXTERNAL not in excluded
                else FIRST_EXTERNAL + 1,
                horizon,
            )
        )
    )
    for speed in controls:
        if speed not in excluded:
            by_speed = dict((speed0, value) for value, speed0 in rows)
            require(
                by_speed[speed] == scalar_coverage(carrier, speed),
                f"scalar/vector mismatch at {speed}",
            )
    return ranked, tail


def root_profile(body: tuple[int, ...]) -> dict[str, object]:
    carrier = good_set(body)
    mass = interval_mass(carrier)
    rows = coverages(carrier, range(FIRST_EXTERNAL, ROOT_HORIZON + 1))
    base = sorted(rows, key=lambda item: (-item[0], item[1]))
    q30 = base[ROOT_TOP - 1][0]
    require(q30 > mass / 7, "root rank thirty misses limiting density")
    threshold = S2 * len(carrier) / (7 * (q30 - mass / 7))
    tail_first = max(ROOT_HORIZON + 1, ceiling(threshold))
    if tail_first > ROOT_HORIZON + 1:
        rows.extend(
            coverages(
                carrier,
                range(ROOT_HORIZON + 1, tail_first),
            )
        )
    require(
        discrepancy_tail(carrier, tail_first) <= q30,
        "root discrepancy tail failed",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top = tuple(ranked[:ROOT_TOP])
    gate = None
    gate_margin = None
    for index in range(ROOT_TOP - ROOT_COVER_SIZE + 1):
        margin = mass - sum(
            (value for value, _ in top[index : index + ROOT_COVER_SIZE]),
            F(0),
        )
        if margin > 0:
            gate = index
            gate_margin = margin
            break
    require(gate is not None and gate_margin is not None, "root gate absent")
    require(
        scalar_coverage(carrier, top[0][1]) == top[0][0],
        "root scalar control failed",
    )
    return {
        "body": body,
        "carrier": carrier,
        "m": mass,
        "top": top,
        "K": gate,
        "gate_margin": gate_margin,
        "tail_first": tail_first,
    }


def suffix_profile(root: dict[str, object], rank: int) -> dict[str, object]:
    body = root["body"]
    apex = root["top"][rank - 1][1]
    prefix = tuple(speed for _, speed in root["top"][:rank])
    excluded = set(prefix)
    literal = subtract_one(root["carrier"], apex)
    direct = good_set(tuple(sorted((*body, apex))))
    require(literal == direct, "first-apex literal/direct mismatch")
    mass = interval_mass(literal)
    rows = coverages(
        literal,
        [
            speed
            for speed in range(FIRST_EXTERNAL, ROOT_HORIZON + 1)
            if speed not in excluded
        ],
    )
    base = sorted(rows, key=lambda item: (-item[0], item[1]))
    q5 = base[4][0]
    require(q5 > mass / 7, "suffix rank five misses limiting density")
    threshold = S2 * len(literal) / (7 * (q5 - mass / 7))
    tail_first = max(ROOT_HORIZON + 1, ceiling(threshold))
    if tail_first > ROOT_HORIZON + 1:
        rows.extend(
            coverages(
                literal,
                [
                    speed
                    for speed in range(ROOT_HORIZON + 1, tail_first)
                    if speed not in excluded
                ],
            )
        )
    require(
        discrepancy_tail(literal, tail_first) <= q5,
        "suffix discrepancy tail failed",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top5 = tuple(ranked[:5])
    margin = mass - sum((value for value, _ in top5), F(0))
    require(
        scalar_coverage(literal, top5[0][1]) == top5[0][0],
        "suffix scalar control failed",
    )
    return {
        "body": body,
        "rank": rank,
        "apex": apex,
        "prefix": prefix,
        "excluded": excluded,
        "carrier": literal,
        "m": mass,
        "r": len(literal),
        "top5": top5,
        "margin": margin,
        "closed": margin > 0,
        "tail_first": tail_first,
    }


def h4_profile(branch: dict[str, object]) -> dict[str, object]:
    ranked, tail = sealed_ranking(
        branch["carrier"],
        branch["excluded"],
        PAIR_HORIZON,
        1,
    )
    q1 = ranked[0][0]
    singleton_margin = F(3, 7) * branch["m"] - q1
    require(singleton_margin > 0, "H4 singleton hypothesis failed")
    level = (branch["m"] - q1) / 4
    require(level > branch["m"] / 7, "H4 threshold not finite")
    tail_first = max(
        FIRST_EXTERNAL,
        ceiling(
            S2
            * branch["r"]
            / (7 * (level - branch["m"] / 7))
        ),
    )
    require(
        tail_first <= PAIR_HORIZON + 1
        and discrepancy_tail(branch["carrier"], tail_first) <= level,
        "H4 tail seal failed",
    )
    core = tuple(
        sorted(
            speed
            for value, speed in ranked
            if value >= level
        )
    )
    require(
        all(speed not in branch["excluded"] for speed in core),
        "excluded prefix entered H4",
    )
    return {
        **branch,
        "q1": q1,
        "single_tail": tail,
        "singleton_margin": singleton_margin,
        "level": level,
        "Htail": tail_first,
        "H": core,
    }


def finite_pair_cap(
    carrier: tuple[tuple[F, F], ...],
    ranked: list[tuple[F, int]],
) -> tuple[F, tuple[int, int], int]:
    """Maximize finite-head pair union with singleton-sum pruning."""

    mass = interval_mass(carrier)
    head = F(0)
    maximizer = None
    paid = 0
    for first_index, (first_value, first) in enumerate(ranked[:-1]):
        if first_value + ranked[first_index + 1][0] <= head:
            break
        after_first = subtract_one(carrier, first)
        for second_value, second in ranked[first_index + 1 :]:
            if first_value + second_value <= head:
                break
            survivor = subtract_one(after_first, second)
            union = mass - interval_mass(survivor)
            paid += 1
            if union > head:
                head = union
                maximizer = tuple(sorted((first, second)))
    require(maximizer is not None and paid > 0, "empty finite pair search")
    return head, maximizer, paid


def pair_profile(
    branch: dict[str, object],
    pair: tuple[int, int],
) -> dict[str, object]:
    residual = subtract_many(branch["carrier"], pair)
    residual_mass = interval_mass(residual)
    require(residual_mass > 0, "H4 pair covered its carrier")
    direct = good_set(
        tuple(sorted((*branch["body"], branch["apex"], *pair)))
    )
    require(
        residual == direct and residual_mass == interval_mass(direct),
        "H4-pair literal/direct mismatch",
    )
    excluded = set(branch["excluded"]) | set(pair)
    ranked, tail_single = sealed_ranking(
        residual,
        excluded,
        PAIR_HORIZON,
        3,
    )
    top3 = tuple(ranked[:3])
    direct_margin = residual_mass - sum(
        (value for value, _ in top3),
        F(0),
    )
    head, maximizer, paid = finite_pair_cap(residual, ranked)
    pair_cap = max(head, ranked[0][0] + tail_single)
    pair_margin = residual_mass - pair_cap - ranked[0][0]
    return {
        "hpair": pair,
        "residual": residual,
        "m": residual_mass,
        "r": len(residual),
        "excluded": excluded,
        "ranked": ranked,
        "tail_single": tail_single,
        "q1": ranked[0][0],
        "top3": top3,
        "direct_margin": direct_margin,
        "B2": pair_cap,
        "B2head": head,
        "B2max": maximizer,
        "pair_margin": pair_margin,
        "paid": paid,
        "adaptive_closed": direct_margin > 0 or pair_margin > 0,
    }


def recursively_close(
    branch: dict[str, object],
    pair: dict[str, object],
) -> dict[str, object]:
    residual = pair["residual"]
    mass = pair["m"]
    q1 = pair["q1"]
    singleton_margin = F(5, 7) * mass - q1
    require(singleton_margin > 0, "recursive singleton hypothesis failed")
    theta = mass - q1
    level = theta / 2
    require(level > mass / 7, "recursive H2 threshold not finite")
    tail_first = max(
        FIRST_EXTERNAL,
        ceiling(S2 * pair["r"] / (7 * (level - mass / 7))),
    )
    require(
        tail_first <= PAIR_HORIZON + 1
        and discrepancy_tail(residual, tail_first) <= level,
        "recursive H2 tail failed",
    )
    core = tuple(
        sorted(
            speed
            for value, speed in pair["ranked"]
            if value >= level
        )
    )
    heavy: list[tuple[int, int]] = []
    tails: list[tuple[tuple[F, F], ...]] = []
    for edge in combinations(core, 2):
        after = subtract_many(residual, edge)
        union = mass - interval_mass(after)
        if union >= theta:
            heavy.append(edge)
            tails.append(after)
    longest: list[F] = []
    horizons: list[int] = []
    checks = 0
    covers: list[tuple[tuple[int, int], int]] = []
    for edge, after in zip(heavy, tails):
        require(interval_mass(after) > 0, "heavy edge covered residual")
        direct = good_set(
            tuple(
                sorted(
                    (
                        *branch["body"],
                        branch["apex"],
                        *pair["hpair"],
                        *edge,
                    )
                )
            )
        )
        require(after == direct, "recursive literal/direct mismatch")
        ell = max(right - left for left, right in after)
        horizon_fraction = F(1, 7) / ell
        horizon = (
            horizon_fraction.numerator // horizon_fraction.denominator
        )
        longest.append(ell)
        horizons.append(horizon)
        excluded = set(pair["excluded"]) | set(edge)
        for speed in range(FIRST_EXTERNAL, horizon + 1):
            if speed in excluded:
                continue
            checks += 1
            scalar = scalar_coverage(after, speed)
            survivor = subtract_one(after, speed)
            require(
                interval_mass(survivor)
                == interval_mass(after) - scalar,
                "recursive scalar/subtraction mismatch",
            )
            if not survivor:
                covers.append((edge, speed))
        require(
            F(1, 7 * (horizon + 1)) < ell,
            "longest-component horizon did not seal",
        )
    return {
        "singleton_margin": singleton_margin,
        "H": core,
        "heavy": tuple(heavy),
        "longest": tuple(longest),
        "horizons": tuple(horizons),
        "checks": checks,
        "covers": tuple(covers),
        "closed": not covers,
    }


def main() -> None:
    roots = [root_profile(body) for body in ROOTS]
    branches = [
        suffix_profile(root, rank)
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
    require(root_counts == EXPECTED_ROOT_COUNTS, "root/suffix counts changed")
    open_branches = [
        h4_profile(branch) for branch in branches if not branch["closed"]
    ]
    require(
        (len(branches), len(open_branches)) == (73, 25),
        "suffix universe changed",
    )

    pair_rows: list[
        tuple[dict[str, object], dict[str, object]]
    ] = []
    for branch in open_branches:
        for pair in combinations(branch["H"], 2):
            pair_rows.append((branch, pair_profile(branch, pair)))

    failures = [
        (branch, pair)
        for branch, pair in pair_rows
        if not pair["adaptive_closed"]
    ]
    totals = (
        len(pair_rows),
        sum(pair["direct_margin"] > 0 for _, pair in pair_rows),
        sum(pair["pair_margin"] > 0 for _, pair in pair_rows),
        len(pair_rows) - len(failures),
        len(failures),
        sum(pair["paid"] for _, pair in pair_rows),
    )
    require(totals == EXPECTED_PAIR_TOTALS, "pair totals changed")
    strata = tuple(
        (
            name,
            len(rows),
            sum(pair["direct_margin"] > 0 for _, pair in rows),
            sum(pair["pair_margin"] > 0 for _, pair in rows),
            sum(pair["adaptive_closed"] for _, pair in rows),
            sum(not pair["adaptive_closed"] for _, pair in rows),
        )
        for name, rows in (
            (
                "rank1",
                [
                    row
                    for row in pair_rows
                    if row[0]["rank"] == 1
                ],
            ),
            (
                "nonrank1",
                [
                    row
                    for row in pair_rows
                    if row[0]["rank"] != 1
                ],
            ),
        )
    )
    require(strata == EXPECTED_STRATA, "rank stratification changed")

    recursive_rows = [
        (branch, pair, recursively_close(branch, pair))
        for branch, pair in failures
    ]
    recursive_data = tuple(
        (
            branch["body"],
            branch["rank"],
            branch["apex"],
            pair["hpair"],
            recursive["H"],
            recursive["heavy"],
            tuple(ftext(value) for value in recursive["longest"]),
            recursive["horizons"],
            recursive["checks"],
        )
        for branch, pair, recursive in recursive_rows
    )
    require(recursive_data == EXPECTED_RECURSIVE, "recursive rows changed")
    require(
        all(recursive["closed"] for _, _, recursive in recursive_rows),
        "recursive singleton closure failed",
    )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/full-independent-four-root-audit/v1\n")
    for branch in open_branches:
        digest.update(
            (
                f"B={','.join(map(str, branch['body']))};"
                f"r={branch['rank']};a={branch['apex']};"
                f"H={','.join(map(str, branch['H']))};"
                f"d1={ftext(branch['singleton_margin'])}\n"
            ).encode()
        )
    for branch, pair in pair_rows:
        digest.update(
            (
                f"P={','.join(map(str, branch['body']))};"
                f"r={branch['rank']};a={branch['apex']};"
                f"e={pair['hpair'][0]},{pair['hpair'][1]};"
                f"d3={ftext(pair['direct_margin'])};"
                f"d21={ftext(pair['pair_margin'])};"
                f"b2={pair['B2max'][0]},{pair['B2max'][1]}\n"
            ).encode()
        )

    print("FULL J6 FOUR-ROOT INDEPENDENT AUDIT")
    print("engine=fresh rational interval sweep + periodic antiderivative")
    print(f"root_counts={root_counts}")
    print(
        f"suffix_counts=all:{len(branches)},scalar_closed:"
        f"{sum(branch['closed'] for branch in branches)},"
        f"scalar_open:{len(open_branches)}"
    )
    print(
        "H4_sizes="
        + repr(
            tuple(
                (
                    branch["body"],
                    branch["rank"],
                    branch["apex"],
                    len(branch["H"]),
                )
                for branch in open_branches
            )
        )
    )
    print(f"pair_totals={totals}")
    print(f"rank_strata={strata}")
    print(
        "recursive_rows="
        + repr(
            tuple(
                (
                    branch["body"],
                    branch["rank"],
                    branch["apex"],
                    pair["hpair"],
                    recursive["H"],
                    recursive["heavy"],
                    recursive["horizons"],
                    recursive["checks"],
                    recursive["covers"],
                )
                for branch, pair, recursive in recursive_rows
            )
        )
    )
    print(f"independent_ledger_sha256={digest.hexdigest()}")
    print("scope=4 roots;73 suffix branches;25 open;not uniform LRC14")
    print("all_independent_controls=PASS")


if __name__ == "__main__":
    main()
