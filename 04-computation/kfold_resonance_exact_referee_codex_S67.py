#!/usr/bin/env python3
"""Exact rational referee for the THM-1095 k=3 residue-template refutation.

The original experiment evaluated the resonance lattice with a finite cutoff.
This referee instead computes every danger intersection by exact rational
interval arithmetic, reconstructs centered moments by subset expansion, and
finds m_7 by exhaustive product-ordered enumeration.  A second endpoint-cell
engine independently checks every intersection measure.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, product


LAMBDA = F(1, 14)
P = 2 * LAMBDA


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def merge(intervals: list[tuple[F, F]]) -> tuple[tuple[F, F], ...]:
    answer: list[tuple[F, F]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if answer and left <= answer[-1][1]:
            answer[-1] = (answer[-1][0], max(answer[-1][1], right))
        else:
            answer.append((left, right))
    return tuple(answer)


def danger(speed: int) -> tuple[tuple[F, F], ...]:
    intervals: list[tuple[F, F]] = []
    radius = LAMBDA / speed
    for numerator in range(speed):
        centre = F(numerator, speed)
        left, right = centre - radius, centre + radius
        if left < 0:
            intervals.extend(((F(0), right), (left + 1, F(1))))
        elif right > 1:
            intervals.extend(((left, F(1)), (F(0), right - 1)))
        else:
            intervals.append((left, right))
    answer = merge(intervals)
    require(sum((b - a for a, b in answer), F(0)) == P,
            "one-runner danger mass")
    return answer


def intersect(
    left: tuple[tuple[F, F], ...], right: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    answer: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        stop = min(left[i][1], right[j][1])
        if start < stop:
            answer.append((start, stop))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(answer)


def intersection_measure_interval(speeds: tuple[int, ...]) -> F:
    intervals = ((F(0), F(1)),)
    for speed in speeds:
        intervals = intersect(intervals, danger(speed))
    return sum((right - left for left, right in intervals), F(0))


def in_danger(speed: int, time: F) -> bool:
    residue = (speed * time) % 1
    return min(residue, 1 - residue) < LAMBDA


def intersection_measure_cells(speeds: tuple[int, ...]) -> F:
    endpoints = {F(0), F(1)}
    for speed in speeds:
        for left, right in danger(speed):
            endpoints.add(left)
            endpoints.add(right)
    ordered = sorted(endpoints)
    total = F(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(in_danger(speed, midpoint) for speed in speeds):
            total += right - left
    return total


def intersection_measure(speeds: tuple[int, ...]) -> F:
    first = intersection_measure_interval(speeds)
    second = intersection_measure_cells(speeds)
    require(first == second, f"intersection engines disagree at {speeds}")
    return first


def centered_moment(speeds: tuple[int, ...]) -> F:
    total = F(0)
    for size in range(len(speeds) + 1):
        for indices in combinations(range(len(speeds)), size):
            chosen = tuple(speeds[index] for index in indices)
            mass = F(1) if not chosen else intersection_measure(chosen)
            total += (-P) ** (len(speeds) - size) * mass
    return total


def multiplicative_minimum(
    speeds: tuple[int, int, int], product_limit: int
) -> tuple[int, tuple[tuple[int, int, int], ...], int]:
    best: int | None = None
    witnesses: list[tuple[int, int, int]] = []
    checked = 0
    for absolute in product(range(1, product_limit + 1), repeat=3):
        size = absolute[0] * absolute[1] * absolute[2]
        if size > product_limit or any(value % 7 == 0 for value in absolute):
            continue
        for signs in product((-1, 1), repeat=3):
            vector = tuple(sign * value for sign, value in zip(signs, absolute))
            checked += 1
            if sum(speed * value for speed, value in zip(speeds, vector)) != 0:
                continue
            if best is None or size < best:
                best = size
                witnesses = [vector]
            elif size == best:
                witnesses.append(vector)
    require(best is not None, f"no resonance found for {speeds}")
    # Every vector of product below `best` occurs in the finite enumeration;
    # finding a witness therefore proves minimality, not just an upper bound.
    return best, tuple(sorted(set(witnesses))), checked


PAIR_CONTROLS = (
    ((1, 15), F(6, 49)),
    ((1, 2), F(5, 49)),
    ((1, 3), F(4, 49)),
    ((1, 4), F(3, 49)),
    ((1, 16), F(5, 49)),
    ((15, 16), F(5, 49)),
    ((3, 5), F(6, 49)),
    ((3, 19), F(6, 49)),
    ((17, 19), F(6, 49)),
)

TRIPLES = (
    (1, 2, 3), (1, 2, 17), (1, 16, 17), (15, 16, 17),
    (1, 2, 5), (1, 2, 19), (1, 16, 19), (15, 16, 19),
    (1, 2, 9), (1, 2, 23), (1, 16, 23), (15, 16, 23),
    (1, 2, 11), (1, 2, 25), (1, 16, 25),
)


def main() -> None:
    print("THM-1095 exact rational centered-moment referee")
    print("threshold 1/14; two independent exact intersection engines")
    print("pair controls: delta(a,b)*a*b")
    for pair, expected in PAIR_CONTROLS:
        value = centered_moment(pair) * pair[0] * pair[1]
        require(value == expected, f"pair control mismatch at {pair}")
        print(f"  {pair}: {value}")

    print("triple exact values: A residue delta m7 delta*m7 witnesses")
    rows = {}
    for speeds in TRIPLES:
        delta = centered_moment(speeds)
        minimum, witnesses, checked = multiplicative_minimum(speeds, 36)
        row = (tuple(value % 14 for value in speeds), delta, minimum,
               delta * minimum, witnesses, checked)
        rows[speeds] = row
        print(
            f"  {speeds} {row[0]} {delta} {minimum} {delta * minimum} "
            f"{witnesses} checked={checked}"
        )

    first = rows[(1, 2, 3)]
    second = rows[(1, 2, 17)]
    require(first[0] == second[0] == (1, 2, 3), "residue-class mismatch")
    require(first[1] == F(61, 2058) and first[2] == 1,
            "first decisive exact row")
    require(second[1] == -F(37, 11662) and second[2] == 8,
            "second decisive exact row")
    require(first[3] > 0 > second[3], "exact sign reversal missing")
    print("decisive same-residue sign reversal")
    print(f"  (1,2,3): delta*m7 = {first[3]}")
    print(f"  (1,2,17): delta*m7 = {second[3]}")
    print("exact one-parameter residue-pair law")
    for c in (3, 5, 9, 11):
        low = rows[(1, 2, c)]
        high = rows[(1, 2, c + 14)]
        require(low[1] == F(70 - 3 * c, 686 * c),
                f"low-c formula mismatch at {c}")
        require(high[1] == F(14 - 3 * (c + 14), 686 * (c + 14)),
                f"high-c formula mismatch at {c + 14}")
        require(low[2] == (c - 1) // 2, f"low-c m7 mismatch at {c}")
        require(high[2] == (c + 13) // 2,
                f"high-c m7 mismatch at {c + 14}")
        print(
            f"  c={c},{c + 14} residue={c}: "
            f"delta={low[1]},{high[1]} m7={low[2]},{high[2]}"
        )
    same_m7 = rows[(1, 16, 17)]
    require(first[0] == same_m7[0] and first[2] == same_m7[2] == 1,
            "same-residue-and-m7 comparison mismatch")
    require(first[1] != same_m7[1],
            "delta unexpectedly determined by residues and m7")
    print("stronger same-residue-and-m7 mismatch")
    print(f"  delta(1,2,3)={first[1]} != delta(1,16,17)={same_m7[1]}")
    print("verdict delta*m7 is not a function of the speed residues modulo 14")
    print("verdict delta is not a function of (speed residues modulo 14,m7)")


if __name__ == "__main__":
    main()
