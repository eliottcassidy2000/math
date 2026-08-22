#!/usr/bin/env python3
"""Exact controls for THM-3442's smooth Hensel-fibre orbit law."""

from __future__ import annotations

import ast
from hashlib import sha256
from itertools import product
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "c72801a6c4fcb534b9e71ab0bb23da9163fc1051c45895b9bcd65241fd0a1c49"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def cycles_of(mapping):
    unseen = set(mapping)
    cycles = []
    while unseen:
        start = min(unseen)
        point = start
        cycle = []
        while point not in cycle:
            require(point in mapping, ("escaped", point))
            cycle.append(point)
            point = mapping[point]
        require(point == start, ("tail", start, point))
        for item in cycle:
            unseen.discard(item)
        cycles.append(tuple(cycle))
    return tuple(sorted(cycles))


def affine_fibre(prime: int, dimension: int, carry_depth: int, level: int,
                 base, action):
    modulus = prime ** level
    scale = prime ** carry_depth
    width = prime ** (level - carry_depth)
    points = tuple(
        tuple((base[i] + scale * offsets[i]) % modulus for i in range(dimension))
        for offsets in product(range(width), repeat=dimension)
    )
    mapping = {point: action(point, modulus) for point in points}
    require(set(mapping.values()) == set(points), ("not fibre permutation", prime, base))
    return points, cycles_of(mapping)


def translation_cases():
    rows = []
    cases = (
        (3, 1, 1, 4),
        (3, 2, 1, 3),
        (3, 3, 1, 3),
        (5, 2, 2, 4),
        (2, 2, 2, 4),
    )
    for prime, dimension, carry_depth, level in cases:
        scale = prime ** carry_depth
        base = tuple(range(dimension))

        def action(point, modulus, step=scale):
            return ((point[0] + step) % modulus,) + point[1:]

        points, cycles = affine_fibre(
            prime, dimension, carry_depth, level, base, action
        )
        orbit_length = prime ** (level - carry_depth)
        orbit_count = prime ** ((dimension - 1) * (level - carry_depth))
        require(len(points) == prime ** (dimension * (level - carry_depth)),
                ("translation fibre", prime, dimension, carry_depth, level))
        require({len(cycle) for cycle in cycles} == {orbit_length},
                ("translation length", prime, dimension, cycles))
        require(len(cycles) == orbit_count,
                ("translation count", prime, dimension, len(cycles), orbit_count))
        rows.append((prime, dimension, carry_depth, level, len(points), orbit_length, len(cycles)))
    return tuple(rows)


def shear_cases():
    rows = []
    prime = 5
    carry_depth = 1
    level = 3
    scale = prime ** carry_depth
    for base_y in range(prime):
        base = (2, base_y)

        def action(point, modulus, step=scale):
            x, y = point
            return ((x + step * y * y) % modulus, y)

        points, cycles = affine_fibre(prime, 2, carry_depth, level, base, action)
        if base_y == 0:
            require(all(len(cycle) == 1 for cycle in cycles),
                    ("zero fibre moved", cycles[:3]))
            require(len(cycles) == len(points), ("zero fibre count", len(cycles)))
        else:
            expected_length = prime ** (level - carry_depth)
            expected_count = prime ** (level - carry_depth)
            require({len(cycle) for cycle in cycles} == {expected_length},
                    ("free shear length", base_y, cycles[:3]))
            require(len(cycles) == expected_count,
                    ("free shear count", base_y, len(cycles)))
        rows.append((base_y, len(points), tuple(sorted({len(cycle) for cycle in cycles})), len(cycles)))
    return tuple(rows)


def torus_cases():
    rows = []
    for prime, carry_depth, level in ((3, 1, 4), (5, 1, 3), (2, 2, 4)):
        modulus = prime ** level
        scale = prime ** carry_depth
        multiplier = 1 + scale
        for residue in range(1, prime):
            base = (residue,)

            def action(point, modulus, unit=multiplier):
                return (unit * point[0] % modulus,)

            points, cycles = affine_fibre(prime, 1, carry_depth, level, base, action)
            expected = prime ** (level - carry_depth)
            require(len(points) == expected and len(cycles) == 1
                    and len(cycles[0]) == expected,
                    ("torus torsor", prime, carry_depth, level, residue))
            rows.append((prime, carry_depth, level, residue, len(points), len(cycles[0])))

        # The same formula on A^1 has vector field x*d/dx and fixes the whole
        # first lifting fibre above its zero.
        first_level = carry_depth + 1
        zero_base = (0,)

        def zero_action(point, modulus, unit=multiplier):
            return (unit * point[0] % modulus,)

        points, cycles = affine_fibre(
            prime, 1, carry_depth, first_level, zero_base, zero_action
        )
        require(len(cycles) == prime and all(len(cycle) == 1 for cycle in cycles),
                ("dilation zero hostile", prime, cycles))
    return tuple(rows)


def inverse2(matrix, modulus: int):
    a, b = matrix[0]
    c, d = matrix[1]
    determinant = (a * d - b * c) % modulus
    inverse = pow(determinant, -1, modulus)
    return (
        (d * inverse % modulus, -b * inverse % modulus),
        (-c * inverse % modulus, a * inverse % modulus),
    )


def mmul2(left, right, modulus: int):
    return tuple(
        tuple(sum(left[i][h] * right[h][j] for h in range(2)) % modulus
              for j in range(2))
        for i in range(2)
    )


def mpow2(matrix, exponent: int, modulus: int):
    answer = ((1, 0), (0, 1))
    base = matrix
    while exponent:
        if exponent & 1:
            answer = mmul2(answer, base, modulus)
        base = mmul2(base, base, modulus)
        exponent //= 2
    return answer


def projective_fibre_cycles(carry_depth: int, level: int):
    prime = 2
    modulus = prime ** level
    scale = prime ** carry_depth
    e = ((0, 1), (1, 1))
    u = tuple(
        tuple((int(i == j) + scale * e[i][j]) % modulus for j in range(2))
        for i in range(2)
    )
    width = prime ** (level - carry_depth)
    slopes = tuple(scale * value for value in range(width))

    def move(slope):
        top = (u[0][0] + u[0][1] * slope) % modulus
        bottom = (u[1][0] + u[1][1] * slope) % modulus
        return bottom * pow(top, -1, modulus) % modulus

    mapping = {slope: move(slope) for slope in slopes}
    cycles = cycles_of(mapping)
    return u, mpow2(u, prime ** (level - carry_depth), modulus), tuple(sorted(len(c) for c in cycles))


def two_adic_boundary():
    bad = projective_fibre_cycles(1, 3)
    good = projective_fibre_cycles(2, 4)
    require(bad[2] == (2, 2), ("p2 c1", bad))
    require(good[2] == (4,), ("p2 c2", good))
    return bad, good


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert node found")

    translations = translation_cases()
    shears = shear_cases()
    torus = torus_cases()
    two_boundary = two_adic_boundary()
    semantic_payload = (translations, shears, torus, two_boundary)
    semantic_sha256 = sha256(repr(semantic_payload).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic_sha256))

    print("smooth Hensel-fibre vector-field orbit controls")
    print(f"translation_rows={translations}")
    print(f"shear_rows={shears}")
    print(f"torus_rows={torus}")
    print(f"two_adic_boundary={two_boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_lf_sha256={lf_sha256(source)}")


if __name__ == "__main__":
    main()
