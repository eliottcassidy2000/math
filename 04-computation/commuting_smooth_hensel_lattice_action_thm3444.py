#!/usr/bin/env python3
"""Exact controls for THM-3444's commuting smooth Hensel lattice action."""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import product
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "a3393fafed345ae96bea50a5dcc4dfb7dadb56320e1d7cefd57a800d3ed615c1"
ROOT = Path(__file__).resolve().parents[1]
BASE_THEOREM = ROOT / "01-canon/theorems/THM-3442-smooth-hensel-fibre-vector-field-orbit-law.md"
EXPECTED_BASE_SHA256 = "5ca39bf0693f7d314fd2907d113eeb8aaca613a5f6e02171cd5a1a4700b3a11a"


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


def fibre_points(prime: int, dimension: int, carry_depth: int, level: int,
                 base):
    modulus = prime ** level
    scale = prime ** carry_depth
    width = prime ** (level - carry_depth)
    return tuple(
        tuple((base[i] + scale * offset[i]) % modulus for i in range(dimension))
        for offset in product(range(width), repeat=dimension)
    )


def affine_translation_cases():
    rows = []
    cases = (
        (3, 2, 2, 1, 3),
        (3, 3, 2, 1, 3),
        (5, 4, 3, 2, 3),
        (2, 3, 2, 2, 4),
    )
    for prime, dimension, rank, carry_depth, level in cases:
        modulus = prime ** level
        scale = prime ** carry_depth
        points = fibre_points(
            prime, dimension, carry_depth, level, tuple(range(dimension))
        )
        generators = []
        for axis in range(rank):
            generators.append({
                point: point[:axis] + ((point[axis] + scale) % modulus,) + point[axis + 1:]
                for point in points
            })
        generators = tuple(generators)
        orbits = generated_orbits(points, generators)
        width = prime ** (level - carry_depth)
        expected_orbit = width ** rank
        expected_count = width ** (dimension - rank)
        require(len(points) == width ** dimension, ("affine points", cases))
        require({len(orbit) for orbit in orbits} == {expected_orbit},
                ("affine orbit", prime, dimension, rank))
        require(len(orbits) == expected_count,
                ("affine bank", prime, dimension, rank, len(orbits)))
        for mapping in generators:
            require(power_map(mapping, width) == {point: point for point in points},
                    ("affine exponent", prime, rank))
        rows.append((prime, dimension, rank, carry_depth, level,
                     len(points), expected_orbit, len(orbits)))
    return tuple(rows)


def nonlinear_conjugate_cases():
    prime = 3
    carry_depth = 1
    level = 3
    modulus = prime ** level
    scale = prime ** carry_depth
    rows = []
    for base_x, base_y in product(range(prime), repeat=2):
        points = fibre_points(prime, 2, carry_depth, level, (base_x, base_y))

        def first(point):
            x, y = point
            return ((x + scale) % modulus,
                    (y + 2 * scale * x + scale * scale) % modulus)

        def second(point):
            x, y = point
            return (x, (y + scale) % modulus)

        maps = tuple({point: action(point) for point in points}
                     for action in (first, second))
        require(compose(maps[0], maps[1]) == compose(maps[1], maps[0]),
                ("nonlinear commute", base_x, base_y))
        identity = {point: point for point in points}
        require(all(power_map(mapping, 9) == identity for mapping in maps),
                ("nonlinear exponent", base_x, base_y))
        orbits = generated_orbits(points, maps)
        require(len(points) == 81 and len(orbits) == 1 and len(orbits[0]) == 81,
                ("nonlinear regular", base_x, base_y, len(orbits)))

        # At the first lift the two translation vectors are (1,2*xbar),(0,1).
        determinant = 1
        require(determinant % prime, ("nonlinear carry", base_x))
        rows.append((base_x, base_y, len(points), len(orbits[0])))
    return tuple(rows)


def dependent_hostile():
    prime = 5
    carry_depth = 1
    level = 2
    modulus = prime ** level
    scale = prime ** carry_depth
    points = fibre_points(prime, 2, carry_depth, level, (0, 0))
    maps = tuple(
        {
            point: ((point[0] + multiple * scale) % modulus, point[1])
            for point in points
        }
        for multiple in (1, 2)
    )
    orbits = generated_orbits(points, maps)
    relation = compose(power_map(maps[0], 3), maps[1])
    identity = {point: point for point in points}
    require(relation == identity, "dependent relation")
    require(len(orbits) == 5 and {len(orbit) for orbit in orbits} == {5},
            ("dependent orbits", orbits))
    return len(points), tuple(sorted(len(orbit) for orbit in orbits)), (3, 1)


def mmul2(left, right, modulus: int):
    return tuple(
        tuple(sum(left[i][h] * right[h][j] for h in range(2)) % modulus
              for j in range(2))
        for i in range(2)
    )


def projective_map(carry_depth: int, level: int):
    modulus = 2 ** level
    scale = 2 ** carry_depth
    e = ((0, 1), (1, 1))
    u = tuple(
        tuple((int(i == j) + scale * e[i][j]) % modulus for j in range(2))
        for i in range(2)
    )
    slopes = tuple(scale * value for value in range(2 ** (level - carry_depth)))

    def move(slope):
        top = (u[0][0] + u[0][1] * slope) % modulus
        bottom = (u[1][0] + u[1][1] * slope) % modulus
        return bottom * pow(top, -1, modulus) % modulus

    return slopes, {slope: move(slope) for slope in slopes}


def two_adic_product_boundary():
    rows = []
    for carry_depth, level, expected_orbit, expected_count in (
        (1, 3, 8, 2),
        (2, 4, 16, 1),
    ):
        slopes, projective = projective_map(carry_depth, level)
        scale = 2 ** carry_depth
        modulus = 2 ** level
        width = 2 ** (level - carry_depth)
        affine = tuple(scale * value for value in range(width))
        points = tuple(product(slopes, affine))
        first = {(x, y): (projective[x], y) for x, y in points}
        second = {(x, y): (x, (y + scale) % modulus) for x, y in points}
        orbits = generated_orbits(points, (first, second))
        require({len(orbit) for orbit in orbits} == {expected_orbit},
                ("two-adic orbit", carry_depth, orbits))
        require(len(orbits) == expected_count,
                ("two-adic count", carry_depth, len(orbits)))
        if carry_depth == 1:
            require(power_map(first, 2) == {point: point for point in points},
                    "two-adic depth-one kernel")
        rows.append((carry_depth, level, len(points), expected_orbit, len(orbits)))
    return tuple(rows)


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert node found")
    require(lf_sha256(BASE_THEOREM) == EXPECTED_BASE_SHA256,
            "THM-3442 dependency changed")

    affine = affine_translation_cases()
    nonlinear = nonlinear_conjugate_cases()
    dependent = dependent_hostile()
    two_adic = two_adic_product_boundary()
    semantic_payload = (affine, nonlinear, dependent, two_adic)
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("commuting smooth Hensel lattice-action controls")
    print(f"affine_rows={affine}")
    print(f"nonlinear_rows={nonlinear}")
    print(f"dependent_hostile={dependent}")
    print(f"two_adic_product={two_adic}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
