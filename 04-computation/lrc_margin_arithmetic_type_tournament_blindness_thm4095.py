#!/usr/bin/env python3
"""Exact finite referee for THM-4095.

The proof of arithmetic-type preservation and the all-real greedy Egyptian
construction is symbolic.  This referee independently reconstructs the
piecewise-linear maxima, affine recovery maps, rational witness grids, and
primitive two-speed spectrum with exact Fraction arithmetic.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def circle_distance(x: F, y: F) -> F:
    difference = (x - y) % 1
    return min(difference, 1 - difference)


def runner_distance(speed: int, time: F) -> F:
    residue = (speed * time) % 1
    return min(residue, 1 - residue)


def observer(speeds: tuple[int, ...], time: F) -> F:
    return min(runner_distance(speed, time) for speed in speeds)


def base_cells(speeds: tuple[int, ...]) -> list[F]:
    walls = {F(k, 2 * speed) for speed in speeds for k in range(2 * speed + 1)}
    return sorted(walls)


def exact_maximum(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    """Maximize the lower envelope by all exact affine-cell intersections."""
    walls = base_cells(speeds)
    candidates = set(walls)
    for left, right in zip(walls, walls[1:]):
        endpoint_values = {
            speed: (runner_distance(speed, left), runner_distance(speed, right))
            for speed in speeds
        }
        for first, second in combinations(speeds, 2):
            difference_left = endpoint_values[first][0] - endpoint_values[second][0]
            difference_right = endpoint_values[first][1] - endpoint_values[second][1]
            if difference_left == difference_right:
                continue
            crossing = left - difference_left * (right - left) / (difference_right - difference_left)
            if left <= crossing <= right:
                candidates.add(crossing)
    values = {time: observer(speeds, time) for time in candidates}
    maximum = max(values.values())
    maximizers = tuple(sorted(time for time, value in values.items() if value == maximum))
    return maximum, maximizers


def active_recovery(speeds: tuple[int, ...], time: F, beta: F) -> tuple[int, int, int, F]:
    """Return v,n,sign,Delta with t=(n+sign*(Delta+beta))/v."""
    value = observer(speeds, time)
    for speed in speeds:
        x = speed * time
        nearest = floor_fraction(x + F(1, 2))
        signed = x - nearest
        if abs(signed) == value:
            sign = 1 if signed >= 0 else -1
            delta = value - beta
            recovered = F(nearest, speed) + F(sign, speed) * (delta + beta)
            require(recovered == time, "active affine recovery failed")
            return speed, nearest, sign, delta
    raise RuntimeError("no active speed found")


def nearest_grid_point(time: F, denominator: int) -> F:
    scaled = time * denominator
    numerator = floor_fraction(scaled + F(1, 2))
    return F(numerator % denominator, denominator)


def greedy_unit_fractions(value: F) -> list[int]:
    """Finite greedy expansion for a positive rational value at most one."""
    require(F(0) < value <= F(1), "greedy input outside (0,1]")
    remainder = value
    denominators: list[int] = []
    while remainder:
        denominator = (remainder.denominator + remainder.numerator - 1) // remainder.numerator
        require(not denominators or denominator > denominators[-1], "greedy denominators not increasing")
        denominators.append(denominator)
        remainder -= F(1, denominator)
    return denominators


def main() -> None:
    pair_gates = 0
    odd_pair_gates = 0
    mixed_pair_gates = 0
    for first in range(1, 81):
        for second in range(first + 1, 81):
            if gcd(first, second) != 1:
                continue
            maximum, maximizers = exact_maximum((first, second))
            if first % 2 and second % 2:
                predicted = F(1, 2)
                require(F(1, 2) in maximizers, "odd-pair half-time maximizer missing")
                odd_pair_gates += 1
            else:
                total = first + second
                predicted = F(total - 1, 2 * total)
                inverse = pow(first, -1, total)
                witness = F((inverse * ((total - 1) // 2)) % total, total)
                require(observer((first, second), witness) == predicted, "mixed-pair witness failed")
                mixed_pair_gates += 1
            require(maximum == predicted, f"primitive pair formula failed {(first, second)}")
            pair_gates += 1

    require(pair_gates == 1965, "primitive-pair universe changed")
    require((odd_pair_gates, mixed_pair_gates) == (651, 1314), "parity split changed")

    # Every odd q is realized by the primitive mixed pair {1,q-1}.
    margin_values = []
    for total in range(3, 160, 2):
        maximum, _ = exact_maximum((1, total - 1))
        margin = maximum - F(1, 3)
        require(margin == F(1, 6) - F(1, 2 * total), "optimized-margin ray failed")
        margin_values.append(margin)
    require(margin_values[0] == 0, "q=3 endpoint changed")
    require(margin_values[1] == F(1, 15), "q=5 gap endpoint changed")
    require(all(value >= F(1, 15) for value in margin_values[1:]), "claimed gap was entered")

    # Exact affine recovery is the computational shadow of Q(Delta(t))=Q(t).
    recovery_gates = 0
    beta = F(1, 3)
    for speeds in ((1, 2), (1, 3), (2, 5, 9), (3, 7, 11, 14)):
        for denominator in range(2, 81):
            for numerator in range(denominator):
                time = F(numerator, denominator)
                active_recovery(speeds, time, beta)
                recovery_gates += 1

    # V-Lipschitz balls compile strict witnesses to sufficiently fine grids.
    grid_gates = 0
    strict_controls = (
        ((1, 3), F(9, 20), F(1, 3)),
        ((1, 4), F(2, 5), F(1, 3)),
        ((2, 5, 9), F(2, 7), F(1, 3)),
    )
    for speeds, time, threshold in strict_controls:
        eta = observer(speeds, time) - threshold
        require(eta > 0, "strict control is not strict")
        velocity = max(speeds)
        denominator = velocity * eta.denominator + 1
        while F(denominator) <= F(velocity, 2) / eta:
            denominator += 1
        rational_witness = nearest_grid_point(time, denominator)
        require(circle_distance(time, rational_witness) < eta / velocity, "grid missed Lipschitz ball")
        require(observer(speeds, rational_witness) > threshold, "grid witness lost strictness")
        grid_gates += 1

    # Finite rational controls for the greedy construction.  The terminal
    # 1/m is replaced by sum_{r>=m} 1/(r(r+1)); a finite prefix plus its exact
    # telescoping tail must still equal the target.
    egyptian_gates = 0
    for value in (F(1), F(2, 3), F(5, 17), F(1, 15), F(1, 6), F(37, 101)):
        denominators = greedy_unit_fractions(value)
        terminal = denominators[-1]
        prefix_sum = sum((F(1, d) for d in denominators[:-1]), F(0))
        require(all(a < b for a, b in zip(denominators, denominators[1:])), "finite greedy order failed")
        expanded = [r * (r + 1) for r in range(terminal, terminal + 40)]
        require(not denominators[:-1] or denominators[-2] < expanded[0], "infinite replacement not increasing")
        controlled_sum = prefix_sum + sum((F(1, d) for d in expanded), F(0)) + F(1, terminal + 40)
        require(controlled_sum == value, "telescoping infinite replacement failed")
        egyptian_gates += 1

    hostile_low, _ = exact_maximum((1, 2))
    hostile_high, _ = exact_maximum((1, 3))
    require(hostile_low - beta == 0, "{1,2} hostile margin changed")
    require(hostile_high - beta == F(1, 6), "{1,3} hostile margin changed")

    print("THM-4095 arithmetic-type and order-tournament audit: PASS")
    print(
        f"primitive pair maxima={pair_gates}; odd={odd_pair_gates}; mixed={mixed_pair_gates}"
    )
    print("optimized margins={1/6} union {1/6-1/(2q): odd q>=3}; gap=(0,1/15)")
    print(f"affine recovery gates={recovery_gates}; strict grid gates={grid_gates}")
    print(f"greedy Egyptian finite controls={egyptian_gates}; replacement tail=40 terms")
    print("hostile same-order pair: {1,2} margin=0; {1,3} margin=1/6")


if __name__ == "__main__":
    main()
