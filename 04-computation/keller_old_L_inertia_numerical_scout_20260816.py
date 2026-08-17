"""Numerical monodromy scout for the old L divisor in the fixed Keller tower.

This is deliberately a numerical discovery tool, not a proof companion.  It
continues all 3^n inverse points around a small transverse loop about the
generic point (a,b,c)=(2/27,1,1) of L=0 and reports the resulting cycle types.

The inverse chart is the exact one used by THM-3533.  Matching uses the
chordal metric on each affine coordinate, so the two branches tending to
infinity remain numerically visible.  Repeated radii and step counts are
hostile controls against a continuation artefact.
"""

from __future__ import annotations

import cmath
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linear_sum_assignment


Point = tuple[complex, complex, complex]


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def l_value(point: Point) -> complex:
    a, b, c = point
    return 27 * a * a * c * c - 18 * a * b * c + 16 * a + b**3 * c - b * b


def forward(point: Point) -> Point:
    x, y, z = point
    unit = 1 + x * y
    return (
        unit**3 * z + y * y * unit * (4 + 3 * x * y),
        y + 3 * x * unit * unit * z + 3 * x * y * y * (4 + 3 * x * y),
        2 * x - 3 * x * x * y - x**3 * z,
    )


def relative_error(left: Point, right: Point) -> float:
    return max(
        abs(x - y) / max(1.0, abs(x), abs(y)) for x, y in zip(left, right)
    )


def inverse_points(target: Point) -> tuple[list[Point], float]:
    a, b, c = target
    lead = l_value(target)
    linear = 4 - 3 * b * c
    scale = max(abs(lead), abs(linear), abs(2 * c), 1.0)
    roots = np.roots(
        np.asarray([lead / scale, 0.0, linear / scale, -2 * c / scale], dtype=complex)
    )
    answer: list[Point] = []
    residual = 0.0
    for raw_x in roots:
        x = complex(raw_x)
        denominator = 12 * a * x * x - b * b * x * x + b * x + 2
        require(abs(denominator) > 1.0e-13, ("inverse denominator", target, x))
        y = b - 3 * a * x * (9 * a * c * x - b * x + 2) / denominator
        require(abs(x) > 1.0e-13, ("inverse x", target, x))
        z = (2 * x - 3 * x * x * y - c) / x**3
        point = (x, y, z)
        residual = max(residual, relative_error(forward(point), target))
        answer.append(point)
    require(len(answer) == 3, ("inverse degree", len(answer)))
    return answer, residual


def inverse_levels(target: Point, depth: int) -> tuple[list[list[Point]], float]:
    levels: list[list[Point]] = []
    parents = [target]
    residual = 0.0
    for _ in range(depth):
        children: list[Point] = []
        for parent in parents:
            new_children, new_residual = inverse_points(parent)
            children.extend(new_children)
            residual = max(residual, new_residual)
        levels.append(children)
        parents = children
    return levels, residual


def projective_pair(value: complex) -> tuple[complex, float]:
    size = abs(value)
    normalizer = np.hypot(1.0, size)
    return value / normalizer, 1.0 / normalizer


def chordal(left: complex, right: complex) -> float:
    left_num, left_den = projective_pair(left)
    right_num, right_den = projective_pair(right)
    return abs(left_num * right_den - right_num * left_den)


def point_distance(left: Point, right: Point) -> float:
    return float(sum(chordal(x, y) ** 2 for x, y in zip(left, right)) ** 0.5)


def distance_matrix(left: list[Point], right: list[Point]) -> np.ndarray:
    return np.asarray(
        [[point_distance(x, y) for y in right] for x in left], dtype=float
    )


def reorder_by_continuation(
    previous: list[Point], current: list[Point]
) -> tuple[list[Point], float]:
    costs = distance_matrix(previous, current)
    rows, columns = linear_sum_assignment(costs)
    require(list(rows) == list(range(len(previous))), ("assignment rows", rows))
    column_for_row = np.empty(len(previous), dtype=int)
    column_for_row[rows] = columns
    selected = [current[int(column_for_row[row])] for row in range(len(previous))]
    return selected, float(max(costs[row, column_for_row[row]] for row in rows))


def endpoint_permutation(initial: list[Point], final: list[Point]) -> tuple[int, ...]:
    costs = distance_matrix(final, initial)
    rows, columns = linear_sum_assignment(costs)
    require(list(rows) == list(range(len(final))), ("endpoint rows", rows))
    permutation = np.empty(len(final), dtype=int)
    permutation[rows] = columns
    require(float(max(costs[row, permutation[row]] for row in rows)) < 2.0e-6,
            ("endpoint mismatch", float(np.max(costs))))
    return tuple(int(value) for value in permutation)


def cycle_lengths(permutation: tuple[int, ...]) -> tuple[int, ...]:
    seen = [False] * len(permutation)
    lengths: list[int] = []
    for start in range(len(permutation)):
        if seen[start]:
            continue
        cursor = start
        length = 0
        while not seen[cursor]:
            seen[cursor] = True
            cursor = permutation[cursor]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


@dataclass(frozen=True)
class Run:
    radius: float
    steps: int
    cycle_rows: tuple[tuple[int, tuple[int, ...], int], ...]
    max_step_cost: float
    max_forward_residual: float


def run(radius: float, steps: int, depth: int = 4) -> Run:
    tracked: list[list[Point]] | None = None
    initial: list[list[Point]] | None = None
    max_step_cost = 0.0
    max_forward_residual = 0.0
    for step in range(steps + 1):
        parameter = radius * cmath.exp(2j * np.pi * step / steps)
        target = (complex(2 / 27) + parameter, 1.0 + 0j, 1.0 + 0j)
        levels, residual = inverse_levels(target, depth)
        max_forward_residual = max(max_forward_residual, residual)
        if tracked is None:
            tracked = [
                sorted(level, key=lambda point: tuple((abs(z), z.real, z.imag) for z in point))
                for level in levels
            ]
            initial = [list(level) for level in tracked]
            continue
        for level_index, level in enumerate(levels):
            tracked[level_index], step_cost = reorder_by_continuation(
                tracked[level_index], level
            )
            max_step_cost = max(max_step_cost, step_cost)

    require(tracked is not None and initial is not None, "empty continuation")
    rows = []
    for level_index, (first, last) in enumerate(zip(initial, tracked), start=1):
        permutation = endpoint_permutation(first, last)
        cycles = cycle_lengths(permutation)
        exponent = len(permutation) - len(cycles)
        rows.append((level_index, cycles, exponent))
    return Run(
        radius=radius,
        steps=steps,
        cycle_rows=tuple(rows),
        max_step_cost=max_step_cost,
        max_forward_residual=max_forward_residual,
    )


def main() -> None:
    runs = (
        run(1.0e-3, 240),
        run(3.0e-4, 480),
        run(1.0e-4, 720),
    )
    baseline = runs[0].cycle_rows
    for item in runs[1:]:
        require(item.cycle_rows == baseline, ("inconsistent cycle rows", runs))
    require(baseline[0][1] == (2, 1), ("level-one positive control", baseline[0]))
    require(all(row[2] % 2 == (1 if row[0] == 1 else 0) for row in baseline),
            ("THM-3531 old-L parity", baseline))
    print("== fixed Keller old-L numerical inertia scout ==")
    for item in runs:
        print(
            f"radius={item.radius:.1e};steps={item.steps};"
            f"max_step_cost={item.max_step_cost:.3e};"
            f"max_forward_residual={item.max_forward_residual:.3e}"
        )
        for depth, cycles, exponent in item.cycle_rows:
            histogram = {length: cycles.count(length) for length in sorted(set(cycles))}
            print(
                f"depth={depth};degree={3**depth};cycles={histogram};"
                f"tame_permutation_exponent={exponent};parity={exponent % 2}"
            )
    print("status=VERIFIED-NUMERICAL DISCOVERY ONLY;no exact inertia theorem or index claim")


if __name__ == "__main__":
    main()
