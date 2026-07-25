#!/usr/bin/env python3
"""Exact hostile controls for THM-2174.

The all-row statements in THM-2174 are elementary.  This companion checks
the displayed defect-six pump, the fixed-core endpoint currents, and the two
safe measures by independent exact interval algorithms.
"""

from fractions import Fraction
from functools import reduce
from math import gcd, lcm


RADIUS = Fraction(1, 14)
CORE = (1, 2, 3, 4, 5, 6, 8)
FAR = (3361, 3362, 6721, 6722, 10081, 10082)
PUMPED_FAR = (1681, 1682, 3361, 3362, 5041, 5042)
ROW = CORE + FAR
PUMPED_ROW = CORE + PUMPED_FAR
BASE = 2
LOW_LEVEL = 4
HIGH_LEVEL = 5

EXPECTED_OLD_MEASURE = Fraction(
    1561405750435498559, 10390707539702618590
)
EXPECTED_NEW_MEASURE = Fraction(
    317645844187362436, 2113439446871390435
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def circle_fraction(x: Fraction) -> Fraction:
    return x % 1


def is_safe(speed: int, time: Fraction) -> bool:
    phase = circle_fraction(speed * time)
    return min(phase, 1 - phase) >= RADIUS


def boundary_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    points = {Fraction(0), Fraction(1)}
    for speed in speeds:
        for index in range(speed):
            points.add(
                circle_fraction(
                    Fraction(14 * index - 1, 14 * speed)
                )
            )
            points.add(
                circle_fraction(
                    Fraction(14 * index + 1, 14 * speed)
                )
            )
    return tuple(sorted(points))


def positive_safe_components(
    speeds: tuple[int, ...],
) -> tuple[tuple[Fraction, Fraction], ...]:
    points = boundary_points(speeds)
    components: list[tuple[Fraction, Fraction]] = []
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if all(is_safe(speed, midpoint) for speed in speeds):
            if components and components[-1][1] == left:
                components[-1] = (components[-1][0], right)
            else:
                components.append((left, right))
    return tuple(components)


def primitive_h(x: Fraction) -> Fraction:
    """The THM-2162 primitive, normalized by H(0)=0."""
    x = circle_fraction(x)
    if x <= Fraction(1, 14):
        return Fraction(6, 7) * x
    if x <= Fraction(13, 14):
        return Fraction(1, 14) - x / 7
    return (
        -Fraction(3, 49)
        + Fraction(6, 7) * (x - Fraction(13, 14))
    )


def endpoint_numerator(
    components: tuple[tuple[Fraction, Fraction], ...],
    residue: int,
) -> Fraction:
    return sum(
        (
            primitive_h(residue * right)
            - primitive_h(residue * left)
            for left, right in components
        ),
        Fraction(0),
    )


def owner(row: tuple[int, ...], level: int) -> tuple[int, ...]:
    power = BASE**level
    return tuple(index for index, value in enumerate(row) if value >= power)


def tie_mask(row: tuple[int, ...], level: int) -> tuple[int, ...]:
    power = BASE**level
    quotient = tuple(value // power for value in row)
    return tuple(
        index
        for index in range(len(row) - 1)
        if quotient[index] < quotient[index + 1]
    )


def carry(
    row: tuple[int, ...], relation: tuple[int, ...], level: int
) -> int:
    power = BASE**level
    return -sum(
        coefficient * (value // power)
        for coefficient, value in zip(relation, row)
    )


def pump(row: tuple[int, ...]) -> tuple[int, ...]:
    low_power = BASE**LOW_LEVEL
    high_power = BASE**HIGH_LEVEL
    return tuple(
        value % low_power + low_power * (value // high_power)
        for value in row
    )


def crossing_relations() -> tuple[tuple[int, ...], ...]:
    relations = []
    for pair_index in range(3):
        relation = [0] * 13
        relation[0] = -1
        relation[7 + 2 * pair_index] = -1
        relation[8 + 2 * pair_index] = 1
        relations.append(tuple(relation))
    return tuple(relations)


def danger_arcs(
    speeds: tuple[int, ...],
) -> list[tuple[Fraction, Fraction]]:
    arcs: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        half_width = Fraction(1, 14 * speed)
        for index in range(speed):
            left = Fraction(index, speed) - half_width
            right = Fraction(index, speed) + half_width
            left %= 1
            width = 2 * half_width
            right = left + width
            if right <= 1:
                arcs.append((left, right))
            else:
                arcs.append((left, Fraction(1)))
                arcs.append((Fraction(0), right - 1))
    return arcs


def safe_measure_by_union(speeds: tuple[int, ...]) -> Fraction:
    arcs = sorted(danger_arcs(speeds))
    total = Fraction(0)
    current_left: Fraction | None = None
    current_right: Fraction | None = None
    for left, right in arcs:
        if current_left is None:
            current_left, current_right = left, right
        elif left <= current_right:
            current_right = max(current_right, right)
        else:
            total += current_right - current_left
            current_left, current_right = left, right
    if current_left is not None:
        total += current_right - current_left
    return 1 - total


def safe_measure_by_cells(speeds: tuple[int, ...]) -> Fraction:
    points = boundary_points(speeds)
    return sum(
        (
            right - left
            for left, right in zip(points, points[1:])
            if all(
                is_safe(speed, (left + right) / 2)
                for speed in speeds
            )
        ),
        Fraction(0),
    )


def vector_rank_over_q(
    vectors: tuple[tuple[int, ...], ...], q: int
) -> int:
    matrix = [[entry % q for entry in vector] for vector in vectors]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (
                row
                for row in range(rank, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, q)
        matrix[rank] = [(inverse * x) % q for x in matrix[rank]]
        for row in range(len(matrix)):
            if row != rank and matrix[row][column]:
                scale = matrix[row][column]
                matrix[row] = [
                    (x - scale * y) % q
                    for x, y in zip(matrix[row], matrix[rank])
                ]
        rank += 1
    return rank


def main() -> None:
    require(
        all(left < right for left, right in zip(ROW, ROW[1:])),
        "source row is not strictly ordered",
    )
    require(
        all(
            left < right
            for left, right in zip(PUMPED_ROW, PUMPED_ROW[1:])
        ),
        "pumped row is not strictly ordered",
    )
    require(reduce(gcd, ROW) == 1, "source row is not primitive")
    require(reduce(gcd, PUMPED_ROW) == 1, "pumped row is not primitive")
    require(pump(ROW) == PUMPED_ROW, "displayed pump is wrong")

    endpoint_modulus = reduce(lcm, (14 * speed for speed in CORE))
    require(endpoint_modulus == 1680, "core endpoint modulus changed")
    require(
        all(
            (old - new) % endpoint_modulus == 0
            for old, new in zip(ROW, PUMPED_ROW)
        ),
        "coordinatewise endpoint phases are not preserved",
    )
    require(
        owner(ROW, LOW_LEVEL) == owner(ROW, HIGH_LEVEL),
        "owner sidecar does not repeat",
    )
    require(
        tie_mask(ROW, LOW_LEVEL) == tie_mask(ROW, HIGH_LEVEL),
        "tie sidecar does not repeat",
    )

    relations = crossing_relations()
    for relation in relations:
        require(sum(a * v for a, v in zip(relation, ROW)) == 0, "old relation")
        require(
            sum(a * v for a, v in zip(relation, PUMPED_ROW)) == 0,
            "pumped relation",
        )
        require(
            carry(ROW, relation, LOW_LEVEL)
            == carry(ROW, relation, HIGH_LEVEL)
            == 0,
            "relation carry does not repeat",
        )
    require(vector_rank_over_q(relations, 2) == 3, "relations lost rank")

    components = positive_safe_components(CORE)
    core_mass = sum(
        (right - left for left, right in components), Fraction(0)
    )
    require(core_mass == Fraction(27, 70), "core mass changed")
    require(len(components) == 16, "core BV component count changed")
    require(
        all(
            endpoint.denominator != 0
            and endpoint_modulus % endpoint.denominator == 0
            for component in components
            for endpoint in component
        ),
        "a core endpoint escaped the common modulus",
    )

    current_one = endpoint_numerator(components, 1)
    current_two = endpoint_numerator(components, 2)
    require(current_one == -Fraction(27, 490), "residue-one current")
    require(current_two == -Fraction(27, 245), "residue-two current")
    for old, new in zip(FAR, PUMPED_FAR):
        residue = old % endpoint_modulus
        require(residue == new % endpoint_modulus, "far residue changed")
        numerator = endpoint_numerator(components, residue)
        require(numerator != 0, "selected endpoint current vanished")
        require(numerator / old != numerator / new, "current retained scale")

    old_union = safe_measure_by_union(ROW)
    old_cells = safe_measure_by_cells(ROW)
    new_union = safe_measure_by_union(PUMPED_ROW)
    new_cells = safe_measure_by_cells(PUMPED_ROW)
    require(old_union == old_cells == EXPECTED_OLD_MEASURE, "old measure")
    require(new_union == new_cells == EXPECTED_NEW_MEASURE, "new measure")
    require(new_union > old_union > 0, "safe-measure separation changed")

    print("THM-2174 exact endpoint-phase scale obstruction")
    print(f"core={CORE}")
    print(f"core_endpoint_modulus={endpoint_modulus}")
    print(f"core_mass={core_mass}, BV_components={len(components)}")
    print(f"C_1={current_one}, C_2={current_two}")
    print(f"source_far={FAR}")
    print(f"pumped_far={PUMPED_FAR}")
    print(
        "repeated_owner="
        + str(owner(ROW, LOW_LEVEL))
        + ", repeated_tie_mask="
        + str(tie_mask(ROW, LOW_LEVEL))
    )
    print("independent_height_one_crossing_relations=3")
    print(f"source_safe_measure={old_union}")
    print(f"pumped_safe_measure={new_union}")
    print(f"measure_difference={new_union-old_union}")
    print("independent_exact_evaluators=interval_union,boundary_cells")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
