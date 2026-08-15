#!/usr/bin/env python3
"""Exact controls for THM-3446's weighted-depth Hensel lattice action."""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import product
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "317fbb4d6c074e76a3886a8653f4b6bb6423f97c135cff44c2df212722e85bb2"
ROOT = Path(__file__).resolve().parents[1]
BASE_THEOREM = ROOT / "01-canon/theorems/THM-3444-commuting-smooth-hensel-vector-field-lattice-action.md"
EXPECTED_BASE_SHA256 = "d6128768b4591273577c351ddad1977506667204ea39008ff24f741d77a4fb32"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def compose(first, second):
    return {point: first[second[point]] for point in second}


def power_map(mapping, exponent: int):
    answer = {point: point for point in mapping}
    base = mapping
    while exponent:
        if exponent & 1:
            answer = compose(answer, base)
        base = compose(base, base)
        exponent //= 2
    return answer


def generated_orbits(points, generators):
    unseen = set(points)
    orbits = []
    while unseen:
        start = min(unseen)
        orbit = {start}
        frontier = [start]
        while frontier:
            point = frontier.pop()
            for mapping in generators:
                image = mapping[point]
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        unseen.difference_update(orbit)
        orbits.append(tuple(sorted(orbit)))
    return tuple(sorted(orbits))


def fibre_points(prime: int, dimension: int, base_depth: int, level: int,
                 base):
    modulus = prime ** level
    scale = prime ** base_depth
    width = prime ** (level - base_depth)
    return tuple(
        tuple((base[i] + scale * offset[i]) % modulus for i in range(dimension))
        for offset in product(range(width), repeat=dimension)
    )


def affine_case(prime: int, dimension: int, depths, level: int):
    base_depth = min(depths)
    modulus = prime ** level
    points = fibre_points(
        prime, dimension, base_depth, level, tuple(range(dimension))
    )
    generators = []
    orders = []
    for axis, depth in enumerate(depths):
        step = prime ** depth
        order = prime ** (level - depth)
        generators.append({
            point: point[:axis] + ((point[axis] + step) % modulus,) + point[axis + 1:]
            for point in points
        })
        orders.append(order)
    generators = tuple(generators)
    identity = {point: point for point in points}
    for mapping, order in zip(generators, orders):
        require(power_map(mapping, order) == identity,
                ("affine exponent", prime, dimension, depths, level, order))
        if order > 1:
            require(power_map(mapping, order // prime) != identity,
                    ("affine exact order", prime, depths, order))
    orbits = generated_orbits(points, generators)
    group_size = 1
    for order in orders:
        group_size *= order
    expected_points = prime ** (dimension * (level - base_depth))
    expected_banks = expected_points // group_size
    bank_exponent = (
        (dimension - len(depths)) * (level - base_depth)
        + sum(depth - base_depth for depth in depths)
    )
    require(expected_banks == prime ** bank_exponent,
            ("bank identity", prime, dimension, depths, level))
    require(len(points) == expected_points,
            ("affine points", prime, dimension, depths, level))
    require({len(orbit) for orbit in orbits} == {group_size},
            ("affine free", prime, dimension, depths, level))
    require(len(orbits) == expected_banks,
            ("affine banks", prime, dimension, depths, level))
    return (prime, dimension, tuple(depths), level, len(points),
            group_size, len(orbits), bank_exponent)


def affine_cases():
    return tuple(
        affine_case(*case)
        for case in (
            (3, 2, (1, 2), 3),
            (3, 3, (1, 2), 3),
            (3, 3, (1, 2, 3), 4),
            (5, 2, (2, 2), 3),
            (2, 2, (2, 3), 4),
        )
    )


def nonlinear_cases():
    prime = 3
    first_depth = 1
    second_depth = 2
    level = 3
    modulus = prime ** level
    first_step = prime ** first_depth
    second_step = prime ** second_depth
    rows = []
    for base_x, base_y in product(range(prime), repeat=2):
        points = fibre_points(prime, 2, first_depth, level, (base_x, base_y))

        def first(point):
            x, y = point
            return ((x + first_step) % modulus,
                    (y + 2 * first_step * x + first_step * first_step) % modulus)

        def second(point):
            x, y = point
            return (x, (y + second_step) % modulus)

        maps = tuple({point: action(point) for point in points}
                     for action in (first, second))
        require(compose(maps[0], maps[1]) == compose(maps[1], maps[0]),
                ("nonlinear commute", base_x, base_y))
        require(power_map(maps[0], 9) == {point: point for point in points},
                ("nonlinear shallow order", base_x, base_y))
        require(power_map(maps[1], 3) == {point: point for point in points},
                ("nonlinear deep order", base_x, base_y))
        orbits = generated_orbits(points, maps)
        require({len(orbit) for orbit in orbits} == {27},
                ("nonlinear weighted orbit", base_x, base_y))
        require(len(orbits) == 3, ("nonlinear weighted bank", base_x, base_y))
        rows.append((base_x, base_y, len(points), 27, len(orbits)))
    return tuple(rows)


def delayed_dependence():
    prime = 3
    rows = []
    for level in (2, 3):
        modulus = prime ** level
        points = fibre_points(prime, 1, 1, level, (0,))
        first = {point: ((point[0] + prime) % modulus,) for point in points}
        if level == 2:
            maps = (first,)
            relation = None
        else:
            second = {point: ((point[0] + prime * prime) % modulus,)
                      for point in points}
            maps = (first, second)
            relation = compose(power_map(first, prime), power_map(second, prime - 1))
            require(relation == {point: point for point in points},
                    "delayed dependence relation")
        orbits = generated_orbits(points, maps)
        require(len(orbits) == 1 and len(orbits[0]) == len(points),
                ("delayed image orbit", level, orbits))
        rows.append((level, len(points), len(orbits[0]), relation is not None))
    return tuple(rows)


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert node found")
    require(lf_sha256(BASE_THEOREM) == EXPECTED_BASE_SHA256,
            "THM-3444 dependency changed")

    affine = affine_cases()
    nonlinear = nonlinear_cases()
    delayed = delayed_dependence()
    semantic_payload = (affine, nonlinear, delayed)
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("weighted-depth commuting Hensel lattice controls")
    print(f"affine_rows={affine}")
    print(f"nonlinear_rows={nonlinear}")
    print(f"delayed_dependence={delayed}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
