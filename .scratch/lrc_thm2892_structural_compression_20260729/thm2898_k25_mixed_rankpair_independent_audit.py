#!/usr/bin/env python3
"""Independent audit of the mixed rank-pair/parity closure at the K=25 root.

The audited scalar/BFS reconstruction is hash-pinned and loaded only as a
foundation.  This program does not import the owner K25 program, atlas,
bootstrap code, or any owner scratch module.  It independently recomputes
every parent-carrier B2 cap used by the THM-2897 rank-pair schedule and every
H4 child cap on the remaining parity branches.
"""

from __future__ import annotations

import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FOUNDATION = (
    ROOT
    / ".scratch/lrc_thm2892_structural_compression_20260729"
    / "thm2898_k25_independent_audit.py"
)
FOUNDATION_SHA = (
    "559380aeb3d266b74ebc041a66ee0da9b533a0f34bf873361e41bcf8e6cfb071"
)

EXPECTED_PARITY = (23, 27, 19, 46, 18, 17)
EXPECTED_PARITY_STATUS = (
    "EXACT_FAIL",
    "LOWER_REFUTED",
    "LOWER_REFUTED",
    "LOWER_REFUTED",
    "EXACT_FAIL",
    "EXACT_FAIL",
)
EXPECTED_RANKPAIR_STEPS = (
    (34, 100, 156, 125),
    (32, 110),
)
EXPECTED_H4 = (13, 13, 10, 11, 7, 3)
EXPECTED_AGGREGATE = (6, 2, 6, 280, 279, 280, 280, 0, 2561)
EXPECTED_RANK_AGGREGATE = (8, 149, 103, 46, 6, 40, 699)
EXPECTED_ORDERED = (
    23,
    27,
    19,
    46,
    18,
    168,
    182,
    34,
    100,
    156,
    125,
    63,
    44,
    29,
    92,
    32,
    110,
    130,
    72,
    17,
    25,
    50,
    54,
    38,
    37,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_foundation():
    require(FOUNDATION.is_file(), "independent K25 foundation is missing")
    require(
        hashlib.sha256(FOUNDATION.read_bytes()).hexdigest() == FOUNDATION_SHA,
        "independent K25 foundation changed",
    )
    spec = importlib.util.spec_from_file_location(
        "thm2898_k25_independent_foundation",
        FOUNDATION,
    )
    require(spec is not None and spec.loader is not None, "cannot load foundation")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


A = load_foundation()
E = A.E


def ftext(value: F | None) -> str:
    if value is None:
        return "none"
    return f"{value.numerator}/{value.denominator}"


def rank_pair_audit(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
    prefix: tuple[int, ...],
) -> tuple[list[dict[str, object]], tuple[int, ...]]:
    """Audit THM-2897 at one scalar-closed state.

    A cheap actual pair is a lower bound for B2.  Thus a nonpositive
    ``m-q5-2*U(pair)`` rigorously refutes the certificate without a full cap
    computation.  Only positive lower margins trigger a sealed exact B2 call.
    """

    prior = set(prefix)
    rows: list[dict[str, object]] = []
    passing: list[int] = []
    for apex in gate:
        if apex in prior:
            continue
        excluded = prior | {apex}
        top5 = A.available_top(profiles[apex], excluded)
        q5 = top5[4][0]
        witness = (top5[0][1], top5[1][1])
        residual = E.subtract_many(profiles[apex]["carrier"], witness)
        witness_union = profiles[apex]["m"] - E.interval_mass(residual)
        direct = E.good_set(tuple(sorted((*A.BODY, apex, *witness))))
        require(residual == direct, f"rank-pair witness mismatch at {apex}")
        lower_margin = profiles[apex]["m"] - q5 - 2 * witness_union

        exact_cap = None
        exact_margin = None
        maximizer = None
        paid = 0
        if lower_margin <= 0:
            status = "LOWER_REFUTED"
        else:
            ranked, tail = E.sealed_ranking(
                profiles[apex]["carrier"],
                excluded,
                E.PAIR_HORIZON,
                5,
            )
            require(ranked[4] == top5[4], f"q5 ranking mismatch at {apex}")
            head, maximizer, paid = E.finite_pair_cap(
                profiles[apex]["carrier"],
                ranked,
            )
            exact_cap = max(head, ranked[0][0] + tail)
            exact_margin = profiles[apex]["m"] - q5 - 2 * exact_cap
            status = "PASS" if exact_margin > 0 else "EXACT_FAIL"
            if status == "PASS":
                passing.append(apex)

        rows.append(
            {
                "prefix": prefix,
                "apex": apex,
                "q5_speed": top5[4][1],
                "q5": q5,
                "witness": witness,
                "witness_union": witness_union,
                "lower_margin": lower_margin,
                "status": status,
                "exact_cap": exact_cap,
                "exact_margin": exact_margin,
                "maximizer": maximizer,
                "paid": paid,
            }
        )
    return rows, tuple(passing)


def rank_line(row: dict[str, object]) -> str:
    return (
        f"P={row['prefix']};a={row['apex']};"
        f"q5={row['q5_speed']}:{ftext(row['q5'])};"
        f"W={row['witness']}:{ftext(row['witness_union'])};"
        f"lower={ftext(row['lower_margin'])};status={row['status']};"
        f"B2={ftext(row['exact_cap'])};"
        f"margin={ftext(row['exact_margin'])};"
        f"max={row['maximizer']};paid={row['paid']}\n"
    )


def main() -> None:
    root = A.root_profile40(A.BODY)
    gate = tuple(speed for _, speed in root["top"][: root["K"]])
    require(gate == A.EXPECTED_GATE and root["K"] == 25, "K25 gate changed")
    profiles = {
        apex: A.apex_profile(root, gate, apex)
        for apex in gate
    }
    search = A.minimum_paid_path(gate, profiles)
    seed_bank = tuple(step[1] for step in search["steps"])
    require(
        seed_bank == EXPECTED_PARITY and search["minimum"] == 6,
        "scalar seed bank changed",
    )

    index = {speed: bit for bit, speed in enumerate(gate)}
    full = (1 << len(gate)) - 1
    mask, rounds0, margins0 = A.scalar_closure(gate, profiles, 0)
    require(not rounds0 and not margins0, "empty prefix scalar-activated")
    prefix: list[int] = []
    seed_cursor = 0
    all_rank_rows: list[dict[str, object]] = []
    state_rows: list[tuple[object, ...]] = []
    schedule: list[tuple[object, ...]] = []
    parity_branches: list[dict[str, object]] = []
    parity_pairs: list[tuple[dict[str, object], dict[str, object]]] = []

    while mask != full:
        require(
            set(prefix) == set(A.labels_from_mask(gate, mask)),
            "ordered prefix and scalar-closed mask disagree",
        )
        prior = tuple(prefix)
        audit, passing = rank_pair_audit(gate, profiles, prior)
        all_rank_rows.extend(audit)
        state_rows.append(
            (
                prior,
                len(audit),
                sum(row["status"] == "LOWER_REFUTED" for row in audit),
                sum(row["status"] != "LOWER_REFUTED" for row in audit),
                passing,
            )
        )

        if passing:
            pass_data = tuple(
                (
                    row["apex"],
                    ftext(row["exact_cap"]),
                    ftext(row["exact_margin"]),
                    row["maximizer"],
                    row["paid"],
                )
                for row in audit
                if row["status"] == "PASS"
            )
            for apex in passing:
                require(apex not in prefix, "rank-pair repeated an apex")
                prefix.append(apex)
                mask |= 1 << index[apex]
            target, rounds, margins = A.scalar_closure(gate, profiles, mask)
            for scalar_round in rounds:
                for apex in scalar_round:
                    require(apex not in prefix, "scalar round repeated an apex")
                    prefix.append(apex)
            schedule.append(
                (
                    "RANKPAIR",
                    prior,
                    pass_data,
                    rounds,
                    tuple(
                        tuple(ftext(value) for value in row)
                        for row in margins
                    ),
                )
            )
            mask = target
            continue

        while (
            seed_cursor < len(seed_bank)
            and mask & (1 << index[seed_bank[seed_cursor]])
        ):
            seed_cursor += 1
        require(seed_cursor < len(seed_bank), "scalar seed bank exhausted")
        apex = seed_bank[seed_cursor]
        seed_cursor += 1
        selected = next(row for row in audit if row["apex"] == apex)
        require(selected["status"] != "PASS", "paid apex has rank-pair proof")

        excluded = set(prior) | {apex}
        top5 = A.available_top(profiles[apex], excluded)
        scalar_margin = profiles[apex]["m"] - sum(
            (value for value, _ in top5),
            F(0),
        )
        require(scalar_margin <= 0, "paid apex was scalar-active")
        branch = {
            "body": A.BODY,
            "rank": len(parity_branches) + 1,
            "apex": apex,
            "prefix": (*prior, apex),
            "excluded": excluded,
            "carrier": profiles[apex]["carrier"],
            "m": profiles[apex]["m"],
            "r": profiles[apex]["r"],
        }
        h4 = E.h4_profile(branch)
        pairs = [
            E.pair_profile(h4, pair)
            for pair in combinations(h4["H"], 2)
        ]
        require(
            all(pair["adaptive_closed"] for pair in pairs),
            "parity branch left a hard H4 residual",
        )
        parity_branches.append(h4)
        parity_pairs.extend((h4, pair) for pair in pairs)

        prefix.append(apex)
        mask |= 1 << index[apex]
        target, rounds, margins = A.scalar_closure(gate, profiles, mask)
        for scalar_round in rounds:
            for activated in scalar_round:
                require(activated not in prefix, "scalar round repeated an apex")
                prefix.append(activated)
        schedule.append(
            (
                "PARITY",
                prior,
                apex,
                selected["status"],
                rounds,
                tuple(
                    tuple(ftext(value) for value in row)
                    for row in margins
                ),
            )
        )
        mask = target

    parity_apices = tuple(
        step[2] for step in schedule if step[0] == "PARITY"
    )
    parity_status = tuple(
        step[3] for step in schedule if step[0] == "PARITY"
    )
    rank_steps = tuple(
        tuple(row[0] for row in step[2])
        for step in schedule
        if step[0] == "RANKPAIR"
    )
    aggregate = (
        len(parity_branches),
        len(rank_steps),
        sum(len(step) for step in rank_steps),
        len(parity_pairs),
        sum(pair["direct_margin"] > 0 for _, pair in parity_pairs),
        sum(pair["pair_margin"] > 0 for _, pair in parity_pairs),
        sum(pair["adaptive_closed"] for _, pair in parity_pairs),
        sum(not pair["adaptive_closed"] for _, pair in parity_pairs),
        sum(pair["paid"] for _, pair in parity_pairs),
    )
    rank_aggregate = (
        len(state_rows),
        len(all_rank_rows),
        sum(row["status"] == "LOWER_REFUTED" for row in all_rank_rows),
        sum(row["status"] != "LOWER_REFUTED" for row in all_rank_rows),
        sum(row["status"] == "PASS" for row in all_rank_rows),
        sum(row["status"] == "EXACT_FAIL" for row in all_rank_rows),
        sum(row["paid"] for row in all_rank_rows),
    )

    require(parity_apices == EXPECTED_PARITY, "parity apex path changed")
    require(parity_status == EXPECTED_PARITY_STATUS, "parity statuses changed")
    require(rank_steps == EXPECTED_RANKPAIR_STEPS, "rank-pair steps changed")
    require(
        tuple(len(branch["H"]) for branch in parity_branches) == EXPECTED_H4,
        "mixed-schedule H4 sizes changed",
    )
    require(aggregate == EXPECTED_AGGREGATE, "mixed aggregate changed")
    require(
        rank_aggregate == EXPECTED_RANK_AGGREGATE,
        "rank-pair aggregate changed",
    )
    require(tuple(prefix) == EXPECTED_ORDERED, "ordered closure changed")

    digest = hashlib.sha256(b"LRC14/THM2898/independent-mixed/v1\n")
    digest.update(
        (
            f"E={A.BODY};gate={gate};seed_bank={seed_bank};"
            f"search={(search['minimum'], search['states'], search['edges'], search['closure_calls'])};"
            f"schedule={tuple(schedule)};states={tuple(state_rows)};"
            f"ordered={tuple(prefix)}\n"
        ).encode()
    )
    for row in all_rank_rows:
        digest.update(rank_line(row).encode())
    for branch, pair in parity_pairs:
        digest.update(
            (
                f"a={branch['apex']};P={branch['prefix']};"
                f"L={pair['hpair']};h={ftext(pair['m'])};"
                f"d3={ftext(pair['direct_margin'])};"
                f"d21={ftext(pair['pair_margin'])};paid={pair['paid']}\n"
            ).encode()
        )

    print("THM-2898 UNIQUE K25 MIXED RANK-PAIR INDEPENDENT AUDIT")
    print(
        f"root={A.BODY};m={ftext(root['m'])};r={len(root['carrier'])};"
        f"K={root['K']};gate={gate}"
    )
    print(
        f"search={(search['minimum'], search['states'], search['edges'], search['closure_calls'])};"
        f"seed_bank={seed_bank}"
    )
    print(f"schedule={tuple(schedule)}")
    print(f"state_rows={tuple(state_rows)}")
    print(f"ordered_prefix={tuple(prefix)}")
    print(f"H4_sizes={tuple(len(branch['H']) for branch in parity_branches)}")
    print(f"aggregate={aggregate}")
    print(f"rank_aggregate={rank_aggregate}")
    print(f"independent_mixed_ledger_sha256={digest.hexdigest()}")
    print("scope=unique K25 root only;not uniform j6;not LRC14")
    print("all_independent_mixed_controls=PASS")


if __name__ == "__main__":
    main()
