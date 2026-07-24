#!/usr/bin/env python3
"""Exact finite ledger for THM-2146.

Universe: every subset C of {1,...,13} with |C| <= 6, at threshold h=3/41.
Consequence object: the exact Haar measure of

    G_h(C) = {t in R/Z : ||ct|| >= h for every c in C}.

Two independent exact interval implementations are compared on every one of
the 4096 inputs:

  1. endpoint cells, classified by a rational midpoint;
  2. direct construction and merging of the closed danger intervals.

The second table applies only the union bound for the extra (noncore) speeds.
It is therefore a proved base measure, not a decorrelation estimate.
"""

from fractions import Fraction as Q
from itertools import combinations
from math import comb

H = Q(3, 41)
CORE = tuple(range(1, 14))


def circle_distance(x: Q) -> Q:
    y = x % 1
    return min(y, 1 - y)


def endpoint_safe_data(speeds: tuple[int, ...]) -> tuple[Q, Q, int]:
    """Safe measure, longest interval, number of positive-length components."""
    if not speeds:
        return Q(1), Q(1), 1
    endpoints = {Q(0), Q(1)}
    for v in speeds:
        for k in range(v):
            endpoints.add(((Q(k) - H) / v) % 1)
            endpoints.add(((Q(k) + H) / v) % 1)
    points = sorted(endpoints)
    safe_lengths: list[Q] = []
    for left, right in zip(points, points[1:]):
        middle = (left + right) / 2
        if all(circle_distance(v * middle) >= H for v in speeds):
            safe_lengths.append(right - left)
    return sum(safe_lengths, Q(0)), max(safe_lengths, default=Q(0)), len(safe_lengths)


def merged_safe_data(speeds: tuple[int, ...]) -> tuple[Q, Q, int]:
    """Independent check by clipping and merging danger intervals in [0,1]."""
    if not speeds:
        return Q(1), Q(1), 1
    intervals: list[tuple[Q, Q]] = []
    for v in speeds:
        # k=-1 and k=v catch the two intervals which meet the endpoints.
        for k in range(-1, v + 1):
            left = max(Q(0), (Q(k) - H) / v)
            right = min(Q(1), (Q(k) + H) / v)
            if left < right:
                intervals.append((left, right))
    intervals.sort()
    merged: list[list[Q]] = []
    for left, right in intervals:
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    danger = sum((right - left for left, right in merged), Q(0))
    safe_lengths: list[Q] = []
    cursor = Q(0)
    for left, right in merged:
        if cursor < left:
            safe_lengths.append(left - cursor)
        cursor = max(cursor, right)
    if cursor < 1:
        safe_lengths.append(Q(1) - cursor)
    return Q(1) - danger, max(safe_lengths, default=Q(0)), len(safe_lengths)


def main() -> None:
    expected_minima = {
        0: Q(1),
        1: Q(35, 41),
        2: Q(59, 82),
        3: Q(1615, 2706),
        4: Q(239, 492),
        5: Q(2729, 7380),
        6: Q(153101, 568260),
    }
    expected_winners = {
        0: [()],
        1: [(v,) for v in CORE],
        2: [(1, 12)],
        3: [(1, 11, 12)],
        4: [(1, 7, 8, 9)],
        5: [(1, 5, 7, 8, 9)],
        6: [(1, 5, 7, 8, 9, 11)],
    }
    minima: dict[int, Q] = {}
    minimizers: dict[int, list[tuple[tuple[int, ...], Q, int]]] = {}
    checked = 0

    print("THM-2146 EXACT SMALL-CORE UNION LEDGER")
    print("threshold h=3/41; universe C subset {1,...,13}, |C|<=6")
    print("endpoint-cell method == merged-danger-interval method on every input")
    print()
    print("j  inputs  min measure       equality sets                 longest      positive-length components")

    for j in range(7):
        best: Q | None = None
        winners: list[tuple[tuple[int, ...], Q, int]] = []
        for speeds in combinations(CORE, j):
            data_a = endpoint_safe_data(speeds)
            data_b = merged_safe_data(speeds)
            assert data_a == data_b, (speeds, data_a, data_b)
            checked += 1
            measure, longest, components = data_a
            if best is None or measure < best:
                best = measure
                winners = [(speeds, longest, components)]
            elif measure == best:
                winners.append((speeds, longest, components))
        assert best is not None
        minima[j] = best
        minimizers[j] = winners
        assert best == expected_minima[j]
        assert [winner[0] for winner in winners] == expected_winners[j]
        if j == 1:
            equality = "all 13 singletons"
            longest = "varies"
            components = "varies"
        else:
            equality = ",".join(str(w[0]) for w in winners)
            longest = str(winners[0][1])
            components = str(winners[0][2])
        print(
            f"{j}  {comb(13, j):>6}  {str(best):<16}  "
            f"{equality:<29} {longest:<12} {components}"
        )

    assert checked == 4096
    print()
    print("DEFECT-SENSITIVE SIX-SPEED RESIDUAL BASE")
    print("Choose seven noncore speeds as the seven-speed body.")
    print("d  core slots  extra residual speeds  proved safe base  factor over 5/41")
    for defect in range(7, 14):
        core_slots = 13 - defect
        extra = defect - 7
        base = minima[core_slots] - extra * 2 * H
        factor = base / Q(5, 41)
        print(
            f"{defect}  {core_slots:>10}  {extra:>21}  "
            f"{str(base):<16} {factor}"
        )

    print()
    print("checked=4096; independent exact methods agree; all assertions passed")
    print("scope: these bases do NOT bound the aggregate covariance with the seven-speed body")


if __name__ == "__main__":
    main()
