#!/usr/bin/env python3
r"""Exact minimum heavy-seed census for four ranked-suffix LRC gates.

Let ``A`` be a finite first-apex hitting gate.  For ``a in A`` and
``P subset A\{a}``, define ``T_a(P)`` to mean that the five largest
individual coverages of the literal ``a``-residual, after excluding
``P union {a}``, sum to less than the residual mass.  Then ``T_a`` is
monotone in ``P``.

Treat a seed ``H subset A`` as the apex branches reserved for nonscalar
cap/flag work.  Starting from ``P=H``, repeatedly activate every remaining
``a`` satisfying ``T_a(P)``.  Activated apices may be scalar-peeled in any
compatible order, and the final closure is the least fixed point of this
monotone process.  If the closure is all of ``A``, only the seed branches
need nonscalar work.  Declaring a seed active is not itself a proof of that
branch: after ordering the seed, each seed branch must be certified on the
suffix excluding only the earlier seeds.

This verifier computes the exact minimum seed size on the same four roots
as the ranked-suffix scalar battery.  All smaller seed subsets are exhausted.
Each activation predicate is globally exact: for each apex it combines the
exact coverages of the other gate labels with a globally tail-sealed top five
over the complement of the entire gate.

This is a scoped target-set-selection battery, not a uniform theorem or
LRC(14).
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCALAR_PATH = (
    ROOT
    / "04-computation/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.py"
)
SCALAR_SHA256 = (
    "6434f020c5aa4000ac81fa081881d93ac0b4190516f854fbd9d8493475baf539"
)
SCALAR_OUTPUT_PATH = (
    ROOT
    / "05-knowledge/results/lrc14_j6_rank_first_suffix_scalar_battery_codex_20260729.out"
)
SCALAR_OUTPUT_SHA256 = (
    "d03f4be7ead138135447b4d720e91d215b1668c2182a099da92129289e605ee9"
)

FIRST_EXTERNAL = 15
BASE_HORIZON = 1600
S2 = F(99, 70)

# Filled after discovery, then locked for ordinary and optimized replay.
EXPECTED_ROOT_RESULTS: tuple[object, ...] | None = (
    (
        (2, 8, 9, 10, 11, 13, 14),
        19,
        6,
        0,
        4,
        (17, 23, 46, 24),
        (2, 3, 4, 5),
        2,
        "7d94fbcef8a10fde2ca3c401fe30322d728a4099dc9119a4c3780613ce4d46cd",
        19,
        0,
        None,
        (19, 17, 23, 46, 24),
        (2, 10, 1, 2),
        (
            (19, 38),
            (34, 53, 110, 168, 69, 154, 37, 25, 72, 125),
            (29,),
            (31, 22),
        ),
        (1, 19, 171, 969, 3876),
        79_475,
    ),
    (
        (1, 3, 9, 10, 11, 12, 14),
        20,
        5,
        0,
        3,
        (39, 23, 16),
        (1, 2, 5),
        1,
        "3f74e3c972622579b4c1770ae9cbbecfa3bf58f33d0ff96f7dbaafcd506c5fe0",
        39,
        1,
        (39, 23, 16),
        None,
        (5, 3, 5, 4),
        (
            (19, 156, 182, 130, 132),
            (17, 21, 120),
            (26, 58, 31, 69, 15),
            (46, 78, 32, 29),
        ),
        (1, 20, 190, 1140),
        27_316,
    ),
    (
        (2, 5, 9, 11, 12, 13, 14),
        21,
        8,
        0,
        5,
        (16, 23, 19, 40, 46),
        (1, 2, 4, 5, 7),
        1,
        "ef1a8e8ffe8c32140eeccc972f4dead72ad75bb9e7b5c2723ebb8b5440d091d4",
        16,
        1,
        (16, 23, 19, 40, 46),
        None,
        (1, 3, 3, 3, 1, 1, 1, 1, 2),
        (
            (17,),
            (50, 69, 81),
            (37, 156, 78),
            (31, 80, 38),
            (75,),
            (54,),
            (24,),
            (20,),
            (25, 32),
        ),
        (1, 21, 210, 1330, 5985, 20349),
        463_622,
    ),
    (
        (2, 3, 4, 5, 6, 7, 8),
        13,
        6,
        0,
        3,
        (22, 26, 33),
        (1, 2, 9),
        4,
        "efd87dabb7c352d9946e8d92f7f4bc2130c7f0b345dcafb7b875a705f70b3587",
        None,
        0,
        None,
        None,
        (4, 2, 3, 1),
        (
            (18, 52, 65, 91),
            (27, 78),
            (36, 39, 20),
            (44,),
        ),
        (1, 13, 78, 286),
        4_527,
    ),
)
EXPECTED_RANK1_RESULTS: tuple[object, ...] | None = (
    (
        (2, 8, 9, 10, 11, 13, 14),
        19,
        1,
        0,
        18,
        (),
        (),
        (17, 23, 46, 24, 29, 31, 22, 38, 34, 53, 110, 168, 69, 154, 37, 25, 72, 125),
    ),
    (
        (1, 3, 9, 10, 11, 12, 14),
        39,
        1,
        0,
        19,
        (),
        (),
        (23, 17, 19, 16, 46, 26, 21, 78, 58, 120, 156, 182, 31, 32, 69, 15, 29, 130, 132),
    ),
    (
        (2, 5, 9, 11, 12, 13, 14),
        16,
        1,
        0,
        20,
        (),
        (),
        (23, 17, 19, 40, 50, 46, 20, 37, 25, 75, 31, 32, 54, 156, 69, 80, 24, 38, 81, 78),
    ),
)
EXPECTED_MINIMUM_SEED_FAMILIES: tuple[tuple[tuple[int, ...], ...], ...] | None = (
    (
        (
            (17, 23, 46, 24),
            (17, 23, 46, 29),
        ),
        ((39, 23, 16),),
        ((16, 23, 19, 40, 46),),
        (
            (22, 26, 33),
            (22, 18, 36),
            (22, 18, 44),
            (22, 18, 33),
        ),
    )
)
EXPECTED_TOTAL_COUNTS: tuple[int, ...] | None = (
    4,
    73,
    25,
    15,
    58,
    8,
    2,
    1,
    16,
    292,
    1298,
    73,
    34_661,
    574_940,
    3,
    0,
    57,
)
EXPECTED_EXTREMA: tuple[object, ...] | None = (
    (
        "11273502240/3720701",
        (2, 8, 9, 10, 11, 13, 14),
        168,
        3030,
    ),
    3030,
    5,
    3,
    4,
)
EXPECTED_SYSTEM_DIGEST: str | None = (
    "43edbe3acf8ec1fbc66301d8c52500afe083c212db2242f76fc8279cf38d7124"
)

CERTIFIED_RANK1_APEX = {
    (2, 8, 9, 10, 11, 13, 14): 19,
    (1, 3, 9, 10, 11, 12, 14): 39,
    (2, 5, 9, 11, 12, 13, 14): 16,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_scalar():
    require(file_sha256(SCALAR_PATH) == SCALAR_SHA256, "scalar battery changed")
    require(
        file_sha256(SCALAR_OUTPUT_PATH) == SCALAR_OUTPUT_SHA256,
        "scalar battery transcript changed",
    )
    spec = importlib.util.spec_from_file_location("j6_suffix_scalar_locked", SCALAR_PATH)
    require(spec is not None and spec.loader is not None, "cannot load scalar battery")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_scalar()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def ceiling(value: F) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


def interval_mass(intervals: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in intervals), F(0))


def profile_apex(
    root: dict[str, object],
    gate: tuple[int, ...],
    apex: int,
) -> dict[str, object]:
    """Build a complete exact activation profile for one gate apex."""

    body = root["body"]
    literal = M.R.subtract_local(root["good"], apex)
    direct, components, direct_mass = M.T.CORE.good_norm(
        tuple(sorted((*body, apex)))
    )
    mass = interval_mass(literal)
    require(
        literal == direct
        and len(literal) == components
        and mass == direct_mass > 0,
        f"literal/direct apex mismatch: {body}, {apex}",
    )

    outside_speeds = [
        speed
        for speed in range(FIRST_EXTERNAL, BASE_HORIZON + 1)
        if speed not in gate
    ]
    outside_rows = M.T.coverages_many(literal, outside_speeds)
    ranked = sorted(outside_rows, key=lambda item: (-item[0], item[1]))
    q5 = ranked[4][0]
    require(
        q5 > mass / 7,
        f"outside-gate q5 misses discrepancy limit: {body}, {apex}",
    )
    threshold = S2 * components / (7 * (q5 - mass / 7))
    tail_first = max(BASE_HORIZON + 1, ceiling(threshold))
    if tail_first > BASE_HORIZON + 1:
        outside_rows.extend(
            M.T.coverages_many(
                literal,
                [
                    speed
                    for speed in range(BASE_HORIZON + 1, tail_first)
                    if speed not in gate
                ],
            )
        )
    require(
        mass / 7 + S2 * components / (7 * tail_first) <= q5,
        f"outside-gate tail did not seal: {body}, {apex}",
    )
    ranked = sorted(outside_rows, key=lambda item: (-item[0], item[1]))
    outside_top5 = tuple(ranked[:5])
    require(
        all(speed not in gate for _, speed in outside_top5),
        f"gate speed survived outside ledger: {body}, {apex}",
    )

    gate_speeds = [speed for speed in gate if speed != apex]
    gate_rows = M.T.coverages_many(literal, gate_speeds)
    gate_coverages = {speed: value for value, speed in gate_rows}
    require(
        len(gate_coverages) == len(gate) - 1,
        f"gate coverage ledger changed: {body}, {apex}",
    )
    for speed, value in gate_coverages.items():
        require(
            value == M.T.coverage(literal, speed),
            f"gate vector/scalar mismatch: {body}, {apex}, {speed}",
        )
    ranked_candidates = tuple(
        sorted(
            (*outside_top5, *gate_rows),
            key=lambda item: (-item[0], item[1]),
        )
    )

    by_speed = {speed: value for value, speed in outside_rows}
    controls = M.select_controls(ranked, by_speed, tail_first)
    for speed in controls:
        require(
            by_speed[speed] == M.T.coverage(literal, speed),
            f"outside vector/scalar mismatch: {body}, {apex}, {speed}",
        )
    return {
        "body": body,
        "apex": apex,
        "m": mass,
        "r": components,
        "outside_top5": outside_top5,
        "gate_coverages": gate_coverages,
        "ranked_candidates": ranked_candidates,
        "threshold": threshold,
        "tail_first": tail_first,
        "outside_controls": len(controls),
        "gate_controls": len(gate_coverages),
        "direct_controls": 1,
    }


def activation_margin(
    profile: dict[str, object],
    gate: tuple[int, ...],
    active: set[int],
) -> F:
    """Exact scalar margin after excluding the active set and this apex."""

    del gate
    values: list[F] = []
    for value, speed in profile["ranked_candidates"]:
        if speed in active:
            continue
        values.append(value)
        if len(values) == 5:
            break
    require(len(values) == 5, "activation candidate set too small")
    return profile["m"] - sum(values, F(0))


def closure(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
    seed: tuple[int, ...],
) -> dict[str, object]:
    """Compute the order-independent least activation fixed point."""

    active = set(seed)
    rounds: list[tuple[int, ...]] = []
    round_margins: list[tuple[F, ...]] = []
    predicate_checks = 0
    while True:
        additions: list[int] = []
        margins: list[F] = []
        for apex in gate:
            if apex in active:
                continue
            margin = activation_margin(profiles[apex], gate, active)
            predicate_checks += 1
            if margin > 0:
                additions.append(apex)
                margins.append(margin)
        if not additions:
            break
        rounds.append(tuple(additions))
        round_margins.append(tuple(margins))
        active.update(additions)
    return {
        "active": frozenset(active),
        "rounds": tuple(rounds),
        "round_margins": tuple(round_margins),
        "checks": predicate_checks,
    }


def minimum_seed_family(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
) -> dict[str, object]:
    """Exhaust through the first successful size and return every minimum seed."""

    tested_by_size: list[int] = []
    total_closure_checks = 0
    for size in range(len(gate) + 1):
        tested = 0
        successes: list[tuple[int, ...]] = []
        first_result: dict[str, object] | None = None
        for seed in combinations(gate, size):
            tested += 1
            result = closure(gate, profiles, seed)
            total_closure_checks += result["checks"]
            if len(result["active"]) == len(gate):
                successes.append(seed)
                if first_result is None:
                    first_result = result
        tested_by_size.append(tested)
        require(
            tested == comb(len(gate), size),
            "seed layer was not exhausted",
        )
        if successes:
            require(first_result is not None, "missing minimum-seed witness")
            return {
                "seed": successes[0],
                "minimum_seeds": tuple(successes),
                "rounds": first_result["rounds"],
                "round_margins": first_result["round_margins"],
                "tested_by_size": tuple(tested_by_size),
                "closure_checks": total_closure_checks,
            }
    raise RuntimeError("full gate failed even as its own seed")


def first_successful_seed_containing(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
    size: int,
    required_apex: int,
) -> tuple[int, ...] | None:
    """Return the deterministic first successful seed of one size containing apex."""

    for seed in combinations(gate, size):
        if required_apex not in seed:
            continue
        if len(closure(gate, profiles, seed)["active"]) == len(gate):
            return seed
    return None


def profile_line(
    profile: dict[str, object],
    gate: tuple[int, ...],
) -> str:
    return (
        f"E={','.join(map(str, profile['body']))};a={profile['apex']};"
        f"m={ftext(profile['m'])};r={profile['r']};"
        f"threshold={ftext(profile['threshold'])};"
        f"tail_first={profile['tail_first']};outside_top5="
        + ",".join(
            f"{speed}:{ftext(value)}"
            for value, speed in profile["outside_top5"]
        )
        + ";gate="
        + ",".join(
            f"{speed}:{ftext(profile['gate_coverages'][speed])}"
            for speed in gate
            if speed != profile["apex"]
        )
        + "\n"
    )


def main() -> None:
    coverage_hard_by_body = {
        body: open_count
        for body, _, _, open_count in M.EXPECTED_ROOT_COUNTS
    }
    roots = [M.global_root(body) for body in M.BATTERY_ROOTS]
    all_profiles: list[dict[str, object]] = []
    root_results: list[tuple[object, ...]] = []
    rank1_results: list[tuple[object, ...]] = []
    search_rows: list[dict[str, object]] = []
    for root in roots:
        gate = tuple(speed for _, speed in root["top"][: root["K"]])
        profiles = {
            apex: profile_apex(root, gate, apex) for apex in gate
        }
        all_profiles.extend(profiles.values())
        zero = closure(gate, profiles, ())
        exact = minimum_seed_family(gate, profiles)
        seed_ranks = tuple(gate.index(apex) + 1 for apex in exact["seed"])
        minimum_seed_hash = hashlib.sha256()
        minimum_seed_hash.update(
            (
                "LRC14/j6/suffix-bootstrap/minimum-seeds/"
                + ",".join(map(str, root["body"]))
                + "/v1\n"
            ).encode()
        )
        for seed in exact["minimum_seeds"]:
            minimum_seed_hash.update(
                (",".join(map(str, seed)) + "\n").encode()
            )
        minimum_seed_digest = minimum_seed_hash.hexdigest()
        certified_apex = CERTIFIED_RANK1_APEX.get(root["body"])
        certified_minimum = tuple(
            seed
            for seed in exact["minimum_seeds"]
            if certified_apex is not None and certified_apex in seed
        )
        plus_one = None
        if certified_apex is not None and not certified_minimum:
            plus_one = first_successful_seed_containing(
                gate,
                profiles,
                len(exact["seed"]) + 1,
                certified_apex,
            )
            require(
                plus_one is not None,
                f"no minimum-plus-one seed contains certified apex: {root['body']}",
            )
        result = (
            root["body"],
            len(gate),
            coverage_hard_by_body[root["body"]],
            len(zero["active"]),
            len(exact["seed"]),
            exact["seed"],
            seed_ranks,
            len(exact["minimum_seeds"]),
            minimum_seed_digest,
            certified_apex,
            len(certified_minimum),
            None if not certified_minimum else certified_minimum[0],
            plus_one,
            tuple(len(items) for items in exact["rounds"]),
            exact["rounds"],
            exact["tested_by_size"],
            exact["closure_checks"],
            exact["minimum_seeds"],
        )
        root_results.append(result)
        search_rows.append(
            {
                "body": root["body"],
                "gate": gate,
                "profiles": profiles,
                "zero": zero,
                "exact": exact,
                "result": result,
            }
        )
        if certified_apex is not None:
            rank1 = closure(gate, profiles, (certified_apex,))
            residual_hard = tuple(
                apex for apex in gate if apex not in rank1["active"]
            )
            rank1_results.append(
                (
                    root["body"],
                    certified_apex,
                    len(rank1["active"]),
                    len(rank1["active"]) - 1,
                    len(residual_hard),
                    tuple(len(items) for items in rank1["rounds"]),
                    rank1["rounds"],
                    residual_hard,
                )
            )

    coverage_hard_total = sum(row[2] for row in root_results)
    minimum_seed_total = sum(row[4] for row in root_results)
    activated_total = sum(row[1] - row[4] for row in root_results)
    total_counts = (
        len(roots),
        sum(row[1] for row in root_results),
        coverage_hard_total,
        minimum_seed_total,
        activated_total,
        sum(row[7] for row in root_results),
        sum(row[10] for row in root_results),
        sum(row[12] is not None for row in root_results),
        sum(root["controls"] for root in roots),
        sum(row["outside_controls"] for row in all_profiles),
        sum(row["gate_controls"] for row in all_profiles),
        sum(row["direct_controls"] for row in all_profiles),
        sum(sum(result[15]) for result in root_results),
        sum(result[16] for result in root_results),
        sum(result[2] for result in rank1_results),
        sum(result[3] for result in rank1_results),
        sum(result[4] for result in rank1_results),
    )
    maximum_threshold = max(
        all_profiles,
        key=lambda row: (
            row["threshold"],
            tuple(-x for x in row["body"]),
            -row["apex"],
        ),
    )
    extrema = (
        (
            ftext(maximum_threshold["threshold"]),
            maximum_threshold["body"],
            maximum_threshold["apex"],
            maximum_threshold["tail_first"],
        ),
        max(row["tail_first"] for row in all_profiles),
        max(result[4] for result in root_results),
        min(result[4] for result in root_results),
        max(result[7] for result in root_results),
    )

    digest = hashlib.sha256()
    digest.update(b"LRC14/j6/suffix-bootstrap-minseed-battery/v1\n")
    for row in search_rows:
        gate = row["gate"]
        digest.update(
            (
                f"ROOT={','.join(map(str, row['body']))};"
                f"gate={','.join(map(str, gate))}\n"
            ).encode()
        )
        for apex in gate:
            digest.update(profile_line(row["profiles"][apex], gate).encode())
        exact = row["exact"]
        digest.update(
            (
                f"seed={','.join(map(str, exact['seed']))};"
                f"minimum_seeds={exact['minimum_seeds']};"
                f"rounds={exact['rounds']};"
                f"round_margins="
                + "|".join(
                    ",".join(ftext(value) for value in values)
                    for values in exact["round_margins"]
                )
                + f";tested={exact['tested_by_size']};"
                f"checks={exact['closure_checks']}\n"
            ).encode()
        )
    digest.update(f"rank1={tuple(rank1_results)}\n".encode())
    system_digest = digest.hexdigest()

    locked_root_results = tuple(root_results)
    locked_rank1_results = tuple(rank1_results)
    if EXPECTED_ROOT_RESULTS is not None:
        require(
            tuple(row[:-1] for row in locked_root_results)
            == EXPECTED_ROOT_RESULTS,
            "minimum-seed root results changed",
        )
    if EXPECTED_MINIMUM_SEED_FAMILIES is not None:
        require(
            tuple(row[-1] for row in locked_root_results)
            == EXPECTED_MINIMUM_SEED_FAMILIES,
            "minimum-seed families changed",
        )
    if EXPECTED_RANK1_RESULTS is not None:
        require(
            locked_rank1_results == EXPECTED_RANK1_RESULTS,
            "proved-rank-one cascades changed",
        )
    if EXPECTED_TOTAL_COUNTS is not None:
        require(total_counts == EXPECTED_TOTAL_COUNTS, "total counts changed")
    if EXPECTED_EXTREMA is not None:
        require(extrema == EXPECTED_EXTREMA, "extrema changed")
    if EXPECTED_SYSTEM_DIGEST is not None:
        require(system_digest == EXPECTED_SYSTEM_DIGEST, "system digest changed")

    print("LRC14 j6 suffix bootstrap minimum-seed battery")
    print(f"root_results={locked_root_results}")
    print(f"proved_rank1_cascades={locked_rank1_results}")
    print(f"total_counts={total_counts}")
    print(f"extrema={extrema}")
    print(f"activation_system_sha256={system_digest}")
    if EXPECTED_TOTAL_COUNTS is None:
        print("mode=DISCOVERY (expected constants are not yet locked)")
    else:
        print("mode=LOCKED")
    print("scope=4/3432 roots;exact minimum scalar-bootstrap seeds;not uniform;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
