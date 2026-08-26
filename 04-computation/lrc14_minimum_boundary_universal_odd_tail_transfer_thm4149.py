#!/usr/bin/env python3
"""Exact arithmetic audit for THM-4149's minimum-one/two boundary.

THM-4148 already handles first-window bodies with minimum at least three.
The width inequality leaves only two small boundary universes.  This audit
checks their explicit physical-sheet clocks, the complete 68-ratio residual
partition, the hostile first-window trap at minimum one, and the enlarged
eleven-body census.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction as Q
from math import comb, gcd


DELTA = Q(1, 14)
GATE = Q(2, 189)
EXPECTED_SEMANTIC_SHA256 = (
    "90de52f9075568cead74de79315daf0494ad7f1f49e5ae7303eff0eae7881f6d"
)
EXPECTED_SEMANTIC_FNV64 = "6784b7e0b01a759d"


def gap(speed: int, phase: Q) -> Q:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def clearance(speeds: tuple[int, ...], phase: Q) -> Q:
    return min(gap(speed, phase) for speed in speeds)


def physical_phase(body_phase: Q, sheet: int) -> Q:
    assert sheet in (0, 1)
    return (body_phase + sheet) / 2


def both_sheets_bad(p: int, q: int, body_phase: Q) -> bool:
    return all(
        min(gap(p, physical_phase(body_phase, sheet)),
            gap(q, physical_phase(body_phase, sheet))) < DELTA
        for sheet in (0, 1)
    )


def cross_walls(p: int, q: int) -> tuple[Q, ...]:
    walls = {Q(0), Q(1)}
    for sheet_shift in (Q(0), Q(1, 2)):
        for speed in (p, q):
            for integer in range(speed):
                for sign in (-1, 1):
                    phase_wall = (Q(integer) + sign * DELTA) / speed
                    walls.add((2 * (phase_wall - sheet_shift)) % 1)
    return tuple(sorted(walls))


def cross_components(p: int, q: int) -> tuple[tuple[Q, Q], ...]:
    walls = cross_walls(p, q)
    cells = [
        (left, right)
        for left, right in zip(walls, walls[1:])
        if both_sheets_bad(p, q, (left + right) / 2)
    ]
    components: list[tuple[Q, Q]] = []
    for left, right in cells:
        if components and components[-1][1] == left and both_sheets_bad(p, q, left):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    for left, right in components:
        assert not both_sheets_bad(p, q, left)
        assert not both_sheets_bad(p, q, right)
    return tuple(components)


def fnv64(data: bytes) -> str:
    value = 14695981039346656037
    for byte in data:
        value ^= byte
        value = value * 1099511628211 & ((1 << 64) - 1)
    return f"{value:016x}"


def main() -> None:
    assert 27 * (13 * 1 - 11) >= 4 * 1 * 11
    assert 27 * (13 * 1 - 12) < 4 * 1 * 12
    assert 27 * (13 * 2 - 20) >= 4 * 2 * 20
    assert 27 * (13 * 2 - 21) < 4 * 2 * 21

    pairs = tuple(
        (p, q)
        for q in range(3, 26, 2)
        for p in range(1, q, 2)
        if gcd(p, q) == 1
    )
    assert len(pairs) == 68

    # Minimum one: every admissible body is contained in {1,...,11}.
    full_one = tuple(range(1, 12))
    one_clocks = {
        "ordinary": (Q(1, 14), 0),
        "one_tail": (Q(6, 77), 1),
        "one_thirteen": (Q(3, 14), 0),
    }
    assert tuple(clearance(full_one, y) for y, _ in one_clocks.values()) == (
        Q(1, 14), Q(6, 77), Q(1, 14)
    )

    one_categories: dict[str, list[tuple[int, int]]] = {
        key: [] for key in one_clocks
    }
    one_tail_minima = {key: Q(1, 2) for key in one_clocks}
    for p, q in pairs:
        if p >= 3:
            category = "ordinary"
        elif q != 13:
            category = "one_tail"
        else:
            category = "one_thirteen"
        y, sheet = one_clocks[category]
        x = physical_phase(y, sheet)
        tail_gap = min(gap(p, x), gap(q, x))
        row_gap = clearance(tuple(2 * h for h in full_one) + (p, q), x)
        assert tail_gap >= DELTA and row_gap >= DELTA
        one_categories[category].append((p, q))
        one_tail_minima[category] = min(one_tail_minima[category], tail_gap)
    assert tuple(len(one_categories[key]) for key in one_clocks) == (56, 11, 1)
    assert tuple(one_tail_minima[key] for key in one_clocks) == (
        Q(3, 28), Q(1, 14), Q(3, 28)
    )

    # The first window really is trapped for (1,13); the third clock is an
    # isolated safe body phase owned from opposite sides by speeds 5 and 9.
    hostile_component = cross_components(1, 13)[0]
    first_window_one = (Q(1, 14), Q(13, 154))
    assert hostile_component == (Q(6, 91), Q(8, 91))
    assert hostile_component[0] < first_window_one[0]
    assert first_window_one[1] < hostile_component[1]
    isolated = Q(3, 14)
    epsilon = Q(1, 10000)
    assert clearance(full_one, isolated) == DELTA
    assert gap(5, isolated) == gap(9, isolated) == DELTA
    assert clearance(full_one, isolated - epsilon) < DELTA
    assert clearance(full_one, isolated + epsilon) < DELTA

    # Minimum two: every admissible body lies in {2,...,20}; its whole common
    # interval [1/28,13/280] is safe.  Two endpoints and both sheets suffice.
    full_two = tuple(range(2, 21))
    common_two = (Q(1, 28), Q(13, 280))
    for y in common_two:
        assert clearance(full_two, y) == DELTA
    two_clocks = {
        "q_at_most_23": (Q(1, 28), 1),
        "q_25_p_at_least_5": (Q(1, 28), 0),
        "q_25_p_small": (Q(13, 280), 1),
    }
    two_categories: dict[str, list[tuple[int, int]]] = {
        key: [] for key in two_clocks
    }
    two_tail_minima = {key: Q(1, 2) for key in two_clocks}
    for p, q in pairs:
        if q <= 23:
            category = "q_at_most_23"
        elif p >= 5:
            category = "q_25_p_at_least_5"
        else:
            category = "q_25_p_small"
        y, sheet = two_clocks[category]
        x = physical_phase(y, sheet)
        tail_gap = min(gap(p, x), gap(q, x))
        row_gap = clearance(tuple(2 * h for h in full_two) + (p, q), x)
        assert tail_gap >= DELTA and row_gap >= DELTA
        two_categories[category].append((p, q))
        two_tail_minima[category] = min(two_tail_minima[category], tail_gap)
    assert tuple(len(two_categories[key]) for key in two_clocks) == (58, 8, 2)
    assert two_categories["q_25_p_small"] == [(1, 25), (3, 25)]
    assert gap(25, physical_phase(Q(1, 28), 1)) == Q(3, 56) < DELTA
    assert gap(1, physical_phase(Q(1, 28), 0)) == Q(1, 56) < DELTA
    assert gap(3, physical_phase(Q(1, 28), 0)) == Q(3, 56) < DELTA
    assert two_tail_minima["q_25_p_small"] == Q(9, 112) > DELTA

    rows: list[tuple[int, int, int]] = []
    for minimum in range(1, 1001):
        last_cap = 351 * minimum // (4 * minimum + 27)
        if last_cap - minimum >= 10:
            rows.append((minimum, last_cap, comb(last_cap - minimum, 10)))
    assert [row[0] for row in rows] == list(range(1, 71))
    assert rows[0] == (1, 11, 1)
    assert rows[1] == (2, 20, 43_758)
    assert max(row[1] for row in rows) == 80
    all_eleven = sum(row[2] for row in rows)
    added = rows[0][2] + rows[1][2]
    assert added == 43_759
    assert all_eleven == 60_301_653_510

    ledger = (
        "gate=2/189;m1cap=11;m2cap=20;pairs=68;"
        "m1partition=56,11,1;m1clocks=1/14:0,6/77:1,3/14:0;"
        "m1bodygaps=1/14,6/77,1/14;m1hostile=6/91..8/91;m1isolated=3/14;"
        "m2partition=58,8,2;m2clocks=1/28:1,1/28:0,13/280:1;"
        "m2hostile=3/56;m2rescue=9/112;"
        "all11=60301653510;added=43759;mrange=1..70;maxlabel=80"
    )
    semantic = hashlib.sha256(ledger.encode()).hexdigest()
    semantic_fnv = fnv64(ledger.encode())
    assert semantic == EXPECTED_SEMANTIC_SHA256
    assert semantic_fnv == EXPECTED_SEMANTIC_FNV64

    print("THM4149_MINIMUM_BOUNDARY_ODD_TAIL_TRANSFER_20260825")
    print("status=PASS;scope=exact rational minimum-one/two boundary")
    print("width_gate=2/189;minimum_caps=(1:11,2:20)")
    print("primitive_residual_pairs=68")
    print("minimum_one_partition=(56,11,1)")
    print("minimum_one_clocks=((1/14,0),(6/77,1),(3/14,0))")
    print("minimum_one_body_gaps=(1/14,6/77,1/14)")
    print("minimum_one_first_window_hostile=(1,13;6/91,8/91)")
    print("minimum_one_isolated_clock=3/14")
    print("minimum_two_partition=(58,8,2)")
    print("minimum_two_clocks=((1/28,1),(1/28,0),(13/280,1))")
    print("minimum_two_hostile=3/56;minimum_two_rescue=9/112")
    print(f"added_eleven_bodies={added}")
    print("all_width_body_minima=1..70;all_width_max_label=80")
    print(f"all_width_eleven_bodies={all_eleven}")
    print(f"semantic_sha256={semantic}")
    print(f"semantic_fnv64={semantic_fnv}")


if __name__ == "__main__":
    main()
