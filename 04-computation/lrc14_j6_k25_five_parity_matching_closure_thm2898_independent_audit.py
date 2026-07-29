#!/usr/bin/env python3
"""Standalone independent audit of the unique K=25 matching closure.

The audit reconstructs the root, all apex profiles, the scalar-closed BFS,
all paid parity branches, and the matching repair.  It imports only the
tracked standalone THM-2895 independent rational-interval engine.  It does
not import the owner K25 program, atlas, bootstrap code, matching probe, or
any scratch helper.  No rank-pair activation is used in this proof.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import deque
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENGINE = (
    ROOT
    / "04-computation"
    / "lrc14_j6_suffix_parity_flag_closure_thm2895_independent_audit.py"
)
ENGINE_SHA = (
    "9023a4042dc8def3f8781668e721549972fb1458d07f2fab89bf7ac3a6f745cc"
)

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
EXPECTED_LEDGER = (
    "4eb7bb8f0d4fdb02cb6c2dc26eacce23b2e41749aefdbc05eac5506b89efc970"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_engine():
    require(ENGINE.is_file(), "tracked THM-2895 independent engine is missing")
    require(
        hashlib.sha256(ENGINE.read_bytes()).hexdigest() == ENGINE_SHA,
        "tracked THM-2895 independent engine changed",
    )
    spec = importlib.util.spec_from_file_location(
        "thm2898_k25_independent_engine",
        ENGINE,
    )
    require(spec is not None and spec.loader is not None, "cannot load engine")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


E = load_engine()


def ftext(value: F | None) -> str:
    if value is None:
        return "none"
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
    """Run unit-cost BFS on scalar-closed masks."""

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


def matching_seed_audit(
    profile: dict[str, object],
    prefix: tuple[int, ...],
    apex: int,
) -> dict[str, object]:
    """Test q5+M_(2,2)<h through THM-2897's exact heavy-edge cutoff."""

    excluded = set(prefix) | {apex}
    top5 = available_top(profile, excluded)
    ranked, tail = E.sealed_ranking(
        profile["carrier"],
        excluded,
        E.PAIR_HORIZON,
        5,
    )
    require(tuple(ranked[:5]) == top5, "matching top-five ranking mismatch")
    h = profile["m"]
    q1 = ranked[0][0]
    q5 = ranked[4][0]
    head, maximizer, paid = E.finite_pair_cap(profile["carrier"], ranked)
    b2 = max(head, q1 + tail)
    target = h - q5
    edge_floor = target - b2
    gamma = E.S2 * profile["r"] / 7
    delta = edge_floor - q1 - h / 7
    base = {
        "prefix": prefix,
        "apex": apex,
        "h": h,
        "r": profile["r"],
        "q1_speed": ranked[0][1],
        "q1": q1,
        "q5_speed": ranked[4][1],
        "q5": q5,
        "b2": b2,
        "b2_maximizer": maximizer,
        "b2_paid": paid,
        "target": target,
        "edge_floor": edge_floor,
        "delta": delta,
        "cutoff": None,
        "candidate_pairs": 0,
        "heavy_edges": (),
        "disjoint_edge_pairs": (),
        "best_disjoint": None,
        "matching_margin": None,
        "hostile_pairs": 0,
        "hostile_matching": None,
        "status": "NONPOSITIVE",
    }
    if delta <= 0:
        return base

    cutoff = E.ceiling(gamma / delta) - 1
    require(
        E.FIRST_EXTERNAL <= cutoff <= E.PAIR_HORIZON,
        "matching cutoff escaped independently sealed horizon",
    )
    by_speed = {speed: value for value, speed in ranked}
    vertices = sorted(
        (
            (by_speed[speed], speed)
            for speed in range(E.FIRST_EXTERNAL, cutoff + 1)
            if speed not in excluded
        ),
        key=lambda item: (-item[0], item[1]),
    )
    edges: list[tuple[F, int, int]] = []
    candidate_pairs = 0
    for first_index, (first_value, first) in enumerate(vertices[:-1]):
        if first_value + vertices[first_index + 1][0] < edge_floor:
            break
        after_first = E.subtract_one(profile["carrier"], first)
        for second_value, second in vertices[first_index + 1 :]:
            if first_value + second_value < edge_floor:
                break
            candidate_pairs += 1
            survivor = E.subtract_one(after_first, second)
            weight = h - E.interval_mass(survivor)
            if weight >= edge_floor:
                x, y = sorted((first, second))
                edges.append((weight, x, y))
    edges.sort(key=lambda row: (-row[0], row[1], row[2]))

    hostile_pairs = 0
    hostile_matching = None
    disjoint_edge_pairs = []
    for first_index, first_edge in enumerate(edges[:-1]):
        for second_edge in edges[first_index + 1 :]:
            total = first_edge[0] + second_edge[0]
            if {first_edge[1], first_edge[2]}.isdisjoint(
                {second_edge[1], second_edge[2]}
            ):
                witness = (
                    total,
                    (first_edge[1], first_edge[2]),
                    (second_edge[1], second_edge[2]),
                )
                disjoint_edge_pairs.append(witness)
                if total >= target:
                    hostile_pairs += 1
                    if hostile_matching is None or witness > hostile_matching:
                        hostile_matching = witness
    disjoint_edge_pairs.sort(reverse=True)
    best_disjoint = (
        disjoint_edge_pairs[0] if disjoint_edge_pairs else None
    )
    matching_margin = (
        None if best_disjoint is None else target - best_disjoint[0]
    )
    base.update(
        {
            "cutoff": cutoff,
            "candidate_pairs": candidate_pairs,
            "heavy_edges": tuple(edges),
            "disjoint_edge_pairs": tuple(disjoint_edge_pairs),
            "best_disjoint": best_disjoint,
            "matching_margin": matching_margin,
            "hostile_pairs": hostile_pairs,
            "hostile_matching": hostile_matching,
            "status": "PASS" if hostile_matching is None else "OPEN",
        }
    )
    return base


def matching_line(row: dict[str, object]) -> str:
    edge_hash = hashlib.sha256(b"LRC14/THM2898/matching-heavy-edges/v1\n")
    for weight, x, y in row["heavy_edges"]:
        edge_hash.update(f"E={x},{y};W={ftext(weight)}\n".encode())
    return (
        f"P={row['prefix']};a={row['apex']};h={ftext(row['h'])};"
        f"r={row['r']};q1={row['q1_speed']}:{ftext(row['q1'])};"
        f"q5={row['q5_speed']}:{ftext(row['q5'])};"
        f"B2={ftext(row['b2'])}:{row['b2_maximizer']}:{row['b2_paid']};"
        f"T={ftext(row['target'])};L={ftext(row['edge_floor'])};"
        f"delta={ftext(row['delta'])};N={row['cutoff']};"
        f"candidate_pairs={row['candidate_pairs']};"
        f"edges={len(row['heavy_edges'])};"
        f"disjoint_pairs={len(row['disjoint_edge_pairs'])};"
        f"best_disjoint={row['best_disjoint']};"
        f"matching_margin={ftext(row['matching_margin'])};"
        f"hostile_pairs={row['hostile_pairs']};"
        f"hostile_matching={row['hostile_matching']};"
        f"status={row['status']};edge_digest={edge_hash.hexdigest()}\n"
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
    seed_bank = tuple(step[1] for step in search["steps"])
    search_result = (
        search["minimum"],
        search["states"],
        search["edges"],
        search["closure_calls"],
    )
    require(
        search_result == (6, 68_407, 323_776, 68_407)
        and seed_bank == EXPECTED_PAID,
        "scalar seed bank changed",
    )

    index = {speed: bit for bit, speed in enumerate(gate)}
    full = (1 << len(gate)) - 1
    mask, rounds0, margins0 = scalar_closure(gate, profiles, 0)
    require(not rounds0 and not margins0, "empty prefix scalar-activated")
    prefix: list[int] = []
    seed_cursor = 0
    matching_rows: list[dict[str, object]] = []
    schedule: list[tuple[object, ...]] = []
    branches: list[dict[str, object]] = []
    pairs: list[tuple[dict[str, object], dict[str, object]]] = []

    while mask != full:
        require(
            set(prefix) == set(labels_from_mask(gate, mask)),
            "ordered prefix and scalar-closed mask disagree",
        )
        prior = tuple(prefix)
        while (
            seed_cursor < len(seed_bank)
            and mask & (1 << index[seed_bank[seed_cursor]])
        ):
            seed_cursor += 1
        require(seed_cursor < len(seed_bank), "scalar seed bank exhausted")
        apex = seed_bank[seed_cursor]
        seed_cursor += 1

        matching = matching_seed_audit(profiles[apex], prior, apex)
        matching_rows.append(matching)
        if matching["status"] == "PASS":
            prefix.append(apex)
            mask |= 1 << index[apex]
            target, rounds, margins = scalar_closure(gate, profiles, mask)
            for scalar_round in rounds:
                prefix.extend(scalar_round)
            schedule.append(
                (
                    "MATCHING",
                    prior,
                    apex,
                    matching["cutoff"],
                    len(matching["heavy_edges"]),
                    matching["best_disjoint"],
                    ftext(matching["matching_margin"]),
                    rounds,
                    tuple(
                        tuple(ftext(value) for value in row)
                        for row in margins
                    ),
                )
            )
            mask = target
            continue

        excluded = set(prior) | {apex}
        top5 = available_top(profiles[apex], excluded)
        require(
            profiles[apex]["m"] - sum((v for v, _ in top5), F(0)) <= 0,
            "paid apex was scalar-active",
        )
        branch = {
            "body": BODY,
            "rank": len(branches) + 1,
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
        require(all(pair["adaptive_closed"] for pair in local), "hard H4 pair")
        branches.append(h4)
        pairs.extend((h4, pair) for pair in local)

        prefix.append(apex)
        mask |= 1 << index[apex]
        target, rounds, margins = scalar_closure(gate, profiles, mask)
        for scalar_round in rounds:
            prefix.extend(scalar_round)
        schedule.append(
            (
                "PARITY",
                prior,
                apex,
                matching["status"],
                rounds,
                tuple(
                    tuple(ftext(value) for value in row)
                    for row in margins
                ),
            )
        )
        mask = target

    parity_apices = tuple(step[2] for step in schedule if step[0] == "PARITY")
    matching_steps = tuple(
        step[2] for step in schedule if step[0] == "MATCHING"
    )
    aggregate = (
        len(branches),
        len(matching_steps),
        len(pairs),
        sum(pair["direct_margin"] > 0 for _, pair in pairs),
        sum(pair["pair_margin"] > 0 for _, pair in pairs),
        sum(pair["adaptive_closed"] for _, pair in pairs),
        sum(not pair["adaptive_closed"] for _, pair in pairs),
        sum(pair["paid"] for _, pair in pairs),
    )

    require(parity_apices == (23, 27, 19, 46, 18), "five-parity path changed")
    require(matching_steps == (17,), "matching seed changed")
    require(tuple(len(row["H"]) for row in branches) == (13, 13, 10, 11, 7),
            "five-parity H4 sizes changed")
    require(aggregate == (5, 1, 277, 276, 277, 277, 0, 2549),
            "five-parity aggregate changed")
    require(
        tuple(prefix)
        == (
            23,
            27,
            19,
            46,
            18,
            168,
            182,
            17,
            25,
            156,
            125,
            44,
            130,
            72,
            54,
            34,
            100,
            92,
            32,
            110,
            38,
            63,
            29,
            37,
            50,
        ),
        "ordered scalar closure changed",
    )
    final_matching = matching_rows[-1]
    require(final_matching["apex"] == 17, "final matching apex changed")
    require(final_matching["cutoff"] == 1844, "matching cutoff changed")
    require(final_matching["status"] == "PASS", "matching certificate failed")
    require(len(final_matching["heavy_edges"]) == 3, "heavy graph changed")
    require(
        len(final_matching["disjoint_edge_pairs"]) == 1,
        "disjoint matching count changed",
    )
    require(
        final_matching["matching_margin"] == F(69186919, 39666627000),
        "matching margin changed",
    )

    digest = hashlib.sha256(b"LRC14/THM2898/independent-matching/v1\n")
    digest.update(
        (
            f"E={BODY};gate={gate};seed_bank={seed_bank};"
            f"search={(search['minimum'], search['states'], search['edges'], search['closure_calls'])};"
            f"schedule={tuple(schedule)};ordered={tuple(prefix)};"
            f"aggregate={aggregate}\n"
        ).encode()
    )
    for row in matching_rows:
        digest.update(matching_line(row).encode())
    for branch, pair in pairs:
        digest.update(
            (
                f"a={branch['apex']};P={branch['prefix']};"
                f"L={pair['hpair']};h={ftext(pair['m'])};"
                f"d3={ftext(pair['direct_margin'])};"
                f"d21={ftext(pair['pair_margin'])};paid={pair['paid']}\n"
            ).encode()
        )
    ledger_digest = digest.hexdigest()
    require(ledger_digest == EXPECTED_LEDGER, "independent ledger changed")

    print("THM-2898 UNIQUE K25 FIVE-PARITY MATCHING INDEPENDENT AUDIT")
    print(
        f"root={BODY};m={ftext(root['m'])};r={len(root['carrier'])};"
        f"K={root['K']};gate={gate}"
    )
    print(
        f"search={(search['minimum'], search['states'], search['edges'], search['closure_calls'])};"
        f"seed_bank={seed_bank}"
    )
    print(f"schedule={tuple(schedule)}")
    print(f"ordered_prefix={tuple(prefix)}")
    print(f"H4_sizes={tuple(len(row['H']) for row in branches)}")
    print(f"aggregate={aggregate}")
    for row in matching_rows:
        print("matching=" + matching_line(row).strip())
    print(f"independent_matching_ledger_sha256={ledger_digest}")
    print("scope=unique K25 root only;not uniform j6;not LRC14")
    print("all_independent_matching_controls=PASS")


if __name__ == "__main__":
    main()
