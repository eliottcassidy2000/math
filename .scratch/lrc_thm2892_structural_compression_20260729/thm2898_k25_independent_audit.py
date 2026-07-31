#!/usr/bin/env python3
"""Independent exact audit of the proposed THM-2898 unique K25 closure.

This audit imports only the already tracked standalone THM-2895 independent
rational-interval engine.  It does not import the THM-2898 owner program,
the all-root atlas, the bootstrap implementation, or any owner scratch
module.  It independently reconstructs the root, all activation profiles,
the closed-state BFS, and every paid H4-pair residual.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import deque
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENGINE = (
    ROOT
    / "04-computation/lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.py"
)
ENGINE_SHA = "9023a4042dc8def3f8781668e721549972fb1458d07f2fab89bf7ac3a6f745cc"

BODY = (1, 8, 10, 11, 12, 13, 14)
EXPECTED_GATE = (
    23,
    27,
    19,
    46,
    18,
    17,
    25,
    34,
    38,
    100,
    63,
    156,
    29,
    125,
    130,
    44,
    37,
    50,
    92,
    168,
    72,
    32,
    110,
    54,
    182,
)
EXPECTED_PAID = (23, 27, 19, 46, 18, 17)
EXPECTED_H4 = (13, 13, 10, 11, 7, 5)
EXPECTED_AGGREGATE = (6, 287, 286, 287, 287, 0)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_engine():
    require(
        hashlib.sha256(ENGINE.read_bytes()).hexdigest() == ENGINE_SHA,
        "independent interval engine changed",
    )
    spec = importlib.util.spec_from_file_location("thm2898_independent_engine", ENGINE)
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E = load_engine()


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def root_profile40(body: tuple[int, ...]) -> dict[str, object]:
    """Independently seal the top forty needed to detect a K=25 gate."""

    carrier = E.good_set(body)
    mass = E.interval_mass(carrier)
    rows = E.coverages(
        carrier,
        range(E.FIRST_EXTERNAL, E.ROOT_HORIZON + 1),
    )
    base = sorted(rows, key=lambda item: (-item[0], item[1]))
    q40 = base[39][0]
    require(q40 > mass / 7, "root rank forty misses limiting density")
    threshold = E.S2 * len(carrier) / (7 * (q40 - mass / 7))
    tail_first = max(E.ROOT_HORIZON + 1, E.ceiling(threshold))
    if tail_first > E.ROOT_HORIZON + 1:
        rows.extend(
            E.coverages(
                carrier,
                range(E.ROOT_HORIZON + 1, tail_first),
            )
        )
    require(
        E.discrepancy_tail(carrier, tail_first) <= q40,
        "root top-forty tail failed",
    )
    ranked = sorted(rows, key=lambda item: (-item[0], item[1]))
    top = tuple(ranked[:40])
    gate = None
    gate_margin = None
    for index in range(35):
        margin = mass - sum(
            (value for value, _ in top[index : index + 6]),
            F(0),
        )
        if gate is None and margin > 0:
            gate = index
            gate_margin = margin
    require(gate is not None and gate_margin is not None, "root gate absent")
    require(
        top[0][0] == E.scalar_coverage(carrier, top[0][1])
        and top[39][0] == E.scalar_coverage(carrier, top[39][1]),
        "root scalar controls failed",
    )
    return {
        "body": body,
        "carrier": carrier,
        "m": mass,
        "r": len(carrier),
        "top": top,
        "K": gate,
        "gate_margin": gate_margin,
        "tail_first": tail_first,
    }


def apex_profile(
    root: dict[str, object],
    gate: tuple[int, ...],
    apex: int,
) -> dict[str, object]:
    """Build exact top-five data for every gate-exclusion state."""

    carrier = E.subtract_one(root["carrier"], apex)
    direct = E.good_set(tuple(sorted((*BODY, apex))))
    require(carrier == direct, "apex literal/direct reconstruction failed")
    mass = E.interval_mass(carrier)
    outside_speeds = [
        speed
        for speed in range(E.FIRST_EXTERNAL, E.ROOT_HORIZON + 1)
        if speed not in gate
    ]
    outside_rows = E.coverages(carrier, outside_speeds)
    ranked = sorted(outside_rows, key=lambda item: (-item[0], item[1]))
    q5 = ranked[4][0]
    require(q5 > mass / 7, "outside-gate q5 misses limiting density")
    threshold = E.S2 * len(carrier) / (7 * (q5 - mass / 7))
    tail_first = max(E.ROOT_HORIZON + 1, E.ceiling(threshold))
    if tail_first > E.ROOT_HORIZON + 1:
        outside_rows.extend(
            E.coverages(
                carrier,
                [
                    speed
                    for speed in range(E.ROOT_HORIZON + 1, tail_first)
                    if speed not in gate
                ],
            )
        )
    require(
        E.discrepancy_tail(carrier, tail_first) <= q5,
        "outside-gate q5 tail failed",
    )
    outside_ranked = sorted(
        outside_rows,
        key=lambda item: (-item[0], item[1]),
    )
    outside_top5 = tuple(outside_ranked[:5])
    gate_rows = E.coverages(
        carrier,
        [speed for speed in gate if speed != apex],
    )
    ranked_candidates = tuple(
        sorted(
            (*outside_top5, *gate_rows),
            key=lambda item: (-item[0], item[1]),
        )
    )
    require(
        len(ranked_candidates) == len(gate) - 1 + 5,
        "apex candidate ledger has wrong size",
    )
    for value, speed in (
        outside_top5[0],
        outside_top5[4],
        gate_rows[0],
        gate_rows[-1],
    ):
        require(
            value == E.scalar_coverage(carrier, speed),
            "apex scalar/vector control failed",
        )
    return {
        "apex": apex,
        "carrier": carrier,
        "m": mass,
        "r": len(carrier),
        "outside_top5": outside_top5,
        "ranked_candidates": ranked_candidates,
        "tail_first": tail_first,
    }


def available_top(
    profile: dict[str, object],
    excluded: set[int],
    count: int = 5,
) -> tuple[tuple[F, int], ...]:
    rows = tuple(
        (value, speed)
        for value, speed in profile["ranked_candidates"]
        if speed not in excluded
    )[:count]
    require(len(rows) == count, "not enough activation candidates")
    return rows


def scalar_closure(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
    start: int,
) -> tuple[int, tuple[tuple[int, ...], ...], tuple[tuple[F, ...], ...]]:
    index = {speed: bit for bit, speed in enumerate(gate)}
    active = start
    rounds: list[tuple[int, ...]] = []
    margins_by_round: list[tuple[F, ...]] = []
    while True:
        additions: list[int] = []
        margins: list[F] = []
        for bit, apex in enumerate(gate):
            if active & (1 << bit):
                continue
            excluded = {
                speed
                for speed, other_bit in index.items()
                if active & (1 << other_bit)
            }
            excluded.add(apex)
            top = available_top(profiles[apex], excluded)
            margin = profiles[apex]["m"] - sum(
                (value for value, _ in top),
                F(0),
            )
            if margin > 0:
                additions.append(bit)
                margins.append(margin)
        if not additions:
            return active, tuple(rounds), tuple(margins_by_round)
        rounds.append(tuple(gate[bit] for bit in additions))
        margins_by_round.append(tuple(margins))
        for bit in additions:
            active |= 1 << bit


def minimum_paid_path(
    gate: tuple[int, ...],
    profiles: dict[int, dict[str, object]],
) -> dict[str, object]:
    """Unit-cost BFS on scalar-closed masks."""

    full = (1 << len(gate)) - 1
    start, rounds, margins = scalar_closure(gate, profiles, 0)
    require(not rounds and not margins, "empty state unexpectedly activates")
    queue = deque([start])
    distance = {start: 0}
    previous: dict[
        int,
        tuple[
            int,
            int,
            tuple[tuple[int, ...], ...],
            tuple[tuple[F, ...], ...],
        ],
    ] = {}
    closure_cache: dict[
        int,
        tuple[int, tuple[tuple[int, ...], ...], tuple[tuple[F, ...], ...]],
    ] = {start: (start, rounds, margins)}
    edges = 0
    while queue and full not in distance:
        state = queue.popleft()
        for bit, apex in enumerate(gate):
            if state & (1 << bit):
                continue
            raw = state | (1 << bit)
            if raw not in closure_cache:
                closure_cache[raw] = scalar_closure(gate, profiles, raw)
            target, scalar_rounds, scalar_margins = closure_cache[raw]
            edges += 1
            if target in distance:
                continue
            distance[target] = distance[state] + 1
            previous[target] = (
                state,
                apex,
                scalar_rounds,
                scalar_margins,
            )
            if target == full:
                break
            queue.append(target)
    require(full in distance, "full gate not reachable")
    steps = []
    state = full
    while state != start:
        old, apex, scalar_rounds, scalar_margins = previous[state]
        steps.append((old, apex, scalar_rounds, scalar_margins))
        state = old
    steps.reverse()
    return {
        "minimum": distance[full],
        "states": len(distance),
        "edges": edges,
        "closure_calls": len(closure_cache),
        "steps": tuple(steps),
    }


def labels_from_mask(gate: tuple[int, ...], mask: int) -> tuple[int, ...]:
    return tuple(
        speed for bit, speed in enumerate(gate) if mask & (1 << bit)
    )


def main() -> None:
    root = root_profile40(BODY)
    gate = tuple(speed for _, speed in root["top"][: root["K"]])
    require(gate == EXPECTED_GATE and root["K"] == 25, "K25 gate changed")
    profiles = {
        apex: apex_profile(root, gate, apex)
        for apex in gate
    }
    search = minimum_paid_path(gate, profiles)
    paid = tuple(step[1] for step in search["steps"])
    require(
        search["minimum"] == 6 and paid == EXPECTED_PAID,
        "minimum paid path changed",
    )

    ordered_prefix: list[int] = []
    branch_rows: list[dict[str, object]] = []
    pair_rows: list[tuple[dict[str, object], dict[str, object]]] = []
    path_rows = []
    for rank, (state, apex, rounds, margins) in enumerate(
        search["steps"],
        start=1,
    ):
        require(
            set(ordered_prefix) == set(labels_from_mask(gate, state)),
            "ordered prefix disagrees with closed BFS state",
        )
        prior = tuple(ordered_prefix)
        excluded = set(prior) | {apex}
        top5 = available_top(profiles[apex], excluded)
        scalar_margin = profiles[apex]["m"] - sum(
            (value for value, _ in top5),
            F(0),
        )
        require(scalar_margin <= 0, "paid apex was already scalar-active")

        branch = {
            "body": BODY,
            "rank": rank,
            "apex": apex,
            "prefix": (*prior, apex),
            "excluded": excluded,
            "carrier": profiles[apex]["carrier"],
            "m": profiles[apex]["m"],
            "r": profiles[apex]["r"],
        }
        h4 = E.h4_profile(branch)
        local = [
            E.pair_profile(h4, pair)
            for pair in combinations(h4["H"], 2)
        ]
        branch_rows.append(h4)
        pair_rows.extend((h4, pair) for pair in local)
        path_rows.append(
            (
                prior,
                apex,
                rounds,
                tuple(tuple(ftext(value) for value in row) for row in margins),
                ftext(scalar_margin),
                len(h4["H"]),
                len(local),
                sum(pair["direct_margin"] > 0 for pair in local),
                sum(pair["pair_margin"] > 0 for pair in local),
            )
        )

        ordered_prefix.append(apex)
        for scalar_round in rounds:
            for activated in scalar_round:
                require(
                    activated not in ordered_prefix,
                    "scalar cascade repeated a gate label",
                )
                ordered_prefix.append(activated)

    require(
        len(ordered_prefix) == 25 and set(ordered_prefix) == set(gate),
        "paid path and scalar cascades do not exhaust the gate",
    )
    aggregate = (
        len(branch_rows),
        len(pair_rows),
        sum(pair["direct_margin"] > 0 for _, pair in pair_rows),
        sum(pair["pair_margin"] > 0 for _, pair in pair_rows),
        sum(pair["adaptive_closed"] for _, pair in pair_rows),
        sum(not pair["adaptive_closed"] for _, pair in pair_rows),
    )
    require(
        tuple(len(branch["H"]) for branch in branch_rows) == EXPECTED_H4,
        "paid H4 core sizes changed",
    )
    require(aggregate == EXPECTED_AGGREGATE, "paid pair census changed")

    digest = hashlib.sha256(b"LRC14/THM2898/independent-K25/v1\n")
    digest.update(
        (
            f"E={BODY};gate={gate};"
            f"search={(search['minimum'],search['states'],search['edges'],search['closure_calls'])};"
            f"path={tuple(path_rows)};ordered={tuple(ordered_prefix)}\n"
        ).encode()
    )
    for branch, pair in pair_rows:
        digest.update(
            (
                f"a={branch['apex']};P={branch['prefix']};"
                f"L={pair['hpair']};h={ftext(pair['m'])};"
                f"d3={ftext(pair['direct_margin'])};"
                f"d21={ftext(pair['pair_margin'])}\n"
            ).encode()
        )

    print("THM-2898 UNIQUE K25 INDEPENDENT AUDIT")
    print(
        f"root={BODY};m={ftext(root['m'])};r={len(root['carrier'])};"
        f"K={root['K']};gate={gate}"
    )
    print(
        f"search={(search['minimum'],search['states'],search['edges'],search['closure_calls'])}"
    )
    print(f"paid_path={paid}")
    print(f"path_rows={tuple(path_rows)}")
    print(f"ordered_prefix={tuple(ordered_prefix)}")
    print(f"H4_sizes={tuple(len(branch['H']) for branch in branch_rows)}")
    print(f"aggregate={aggregate}")
    print(f"independent_ledger_sha256={digest.hexdigest()}")
    print("scope=unique K25 root only;not uniform j6;not LRC14")
    print("all_independent_controls=PASS")


if __name__ == "__main__":
    main()
