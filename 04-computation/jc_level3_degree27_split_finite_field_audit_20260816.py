#!/usr/bin/env python3
"""Independent split-fibre audit of the fixed-F degree-27 x-eliminant gate.

This companion does not use the characteristic-zero nested-algebra
interpolation.  It searches directly over F_101 for a target with three
explicit first preimages and nine explicit second preimages, multiplies the
nine third-level cubic cores, and checks that the resulting degree-27
polynomial is squarefree.  Hence every cubic is separable and the three
degree-nine blocks are pairwise coprime at this specialization.
"""

from __future__ import annotations

import hashlib


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


P = 101
EXPECTED_TARGET = (93, 28, 83)
EXPECTED_LEVEL_ONE = [(10, 17, 82), (39, 90, 66), (52, 36, 24)]
EXPECTED_LEVEL_TWO_BLOCKS = [
    [(16, 70, 96), (90, 98, 7), (96, 9, 53)],
    [(50, 45, 99), (72, 58, 17), (80, 32, 9)],
    [(2, 57, 13), (21, 83, 49), (78, 15, 61)],
]
EXPECTED_PRODUCT_SHA256 = "fd4e337b6b7c5bb043741c72e39aaec2edd16e8cd9a23a216d68c099bf0f8c67"


def inv(value: int) -> int:
    value %= P
    require(value != 0, "attempted inversion of zero")
    return pow(value, P - 2, P)


def l_value(point: tuple[int, int, int]) -> int:
    a, b, c = point
    return (27 * a * a * c * c - 18 * a * b * c + 16 * a + b**3 * c - b * b) % P


def fmap(point: tuple[int, int, int]) -> tuple[int, int, int]:
    x, y, z = point
    unit = (1 + x * y) % P
    return (
        (unit**3 * z + y * y * unit * (4 + 3 * x * y)) % P,
        (y + 3 * x * unit * unit * z + 3 * x * y * y * (4 + 3 * x * y)) % P,
        (2 * x - 3 * x * x * y - x**3 * z) % P,
    )


def core_coefficients(target: tuple[int, int, int]) -> list[int]:
    a, b, c = target
    return [(-2 * c) % P, (4 - 3 * b * c) % P, 0, l_value(target)]


def evaluate(coefficients: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % P
    return result


def roots(coefficients: list[int]) -> list[int]:
    return [value for value in range(P) if evaluate(coefficients, value) == 0]


def inverse_point(target: tuple[int, int, int], x: int) -> tuple[int, int, int] | None:
    a, b, c = target
    denominator = ((12 * a - b * b) * x * x + b * x + 2) % P
    if denominator == 0 or x == 0:
        return None
    y = (
        b
        - 3 * a * x * ((9 * a * c - b) * x + 2) * inv(denominator)
    ) % P
    z = ((2 * x - c - 3 * x * x * y) * inv(x**3)) % P
    point = (x, y, z)
    require(fmap(point) == target, "finite-field inverse formula failed")
    return point


def trim(poly: list[int]) -> list[int]:
    result = [value % P for value in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] = (result[i + j] + x * y) % P
    return trim(result)


def derivative(poly: list[int]) -> list[int]:
    return trim([(index * poly[index]) % P for index in range(1, len(poly))])


def remainder(numerator: list[int], denominator: list[int]) -> list[int]:
    out = trim(numerator)
    divisor = trim(denominator)
    require(divisor != [0], "polynomial division by zero")
    inverse_lead = inv(divisor[-1])
    while len(out) >= len(divisor) and out != [0]:
        shift = len(out) - len(divisor)
        scale = out[-1] * inverse_lead % P
        for index, coefficient in enumerate(divisor):
            out[index + shift] = (out[index + shift] - scale * coefficient) % P
        out = trim(out)
    return out


def gcd(left: list[int], right: list[int]) -> list[int]:
    x, y = trim(left), trim(right)
    while y != [0]:
        x, y = y, remainder(x, y)
    scale = inv(x[-1])
    return trim([(scale * value) % P for value in x])


def block_product(points: list[tuple[int, int, int]]) -> list[int]:
    result = [1]
    for point in points:
        result = multiply(result, core_coefficients(point))
    return result


def split_level(target: tuple[int, int, int]) -> list[tuple[int, int, int]] | None:
    if l_value(target) == 0:
        return None
    x_roots = roots(core_coefficients(target))
    if len(x_roots) != 3:
        return None
    points = [inverse_point(target, x) for x in x_roots]
    if any(point is None for point in points):
        return None
    result = [point for point in points if point is not None]
    require(len(set(result)) == 3, "split fibre repeated a point")
    return result


def find_witness() -> tuple[
    tuple[int, int, int],
    list[tuple[int, int, int]],
    list[list[tuple[int, int, int]]],
]:
    # Deterministic permutation of F_101^3; expected identity-Frobenius density
    # is about 1/1296 at depth two.
    state = 1
    seen: set[int] = set()
    for _ in range(P**3):
        state = (742939 * state + 2311) % (P**3)
        if state in seen:
            break
        seen.add(state)
        target = (state % P, (state // P) % P, (state // (P * P)) % P)
        first = split_level(target)
        if first is None:
            continue
        second_blocks: list[list[tuple[int, int, int]]] = []
        for point in first:
            child = split_level(point)
            if child is None:
                break
            second_blocks.append(child)
        if len(second_blocks) == 3:
            second_points = [point for block in second_blocks for point in block]
            if all(l_value(point) != 0 for point in second_points):
                trial_blocks = [block_product(block) for block in second_blocks]
                trial_product = [1]
                for block in trial_blocks:
                    trial_product = multiply(trial_product, block)
                if len(trial_product) == 28 and gcd(trial_product, derivative(trial_product)) == [1]:
                    return target, first, second_blocks
    raise RuntimeError("no fully split depth-two target found")


target, level_one, level_two_blocks = find_witness()
require(target == EXPECTED_TARGET, "deterministic witness target changed")
require(level_one == EXPECTED_LEVEL_ONE, "level-one split fibre changed")
require(level_two_blocks == EXPECTED_LEVEL_TWO_BLOCKS, "level-two split blocks changed")
level_two = [point for block in level_two_blocks for point in block]
require(len(level_two) == 9 and len(set(level_two)) == 9, "depth-two fibre is not nine points")

blocks = [block_product(block) for block in level_two_blocks]
require(all(len(block) - 1 == 9 for block in blocks), "a degree-nine block lost degree")
product = [1]
for block in blocks:
    product = multiply(product, block)
require(len(product) - 1 == 27, "degree-27 product lost degree")
require(gcd(product, derivative(product)) == [1], "degree-27 product is not squarefree")
for i in range(3):
    require(gcd(blocks[i], derivative(blocks[i])) == [1], "a block is inseparable")
    for j in range(i):
        require(gcd(blocks[i], blocks[j]) == [1], "two degree-nine blocks meet")

ledger = "\n".join(f"{index}:{coefficient}" for index, coefficient in enumerate(product))
digest = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(digest == EXPECTED_PRODUCT_SHA256, "degree-27 finite-field ledger changed")

print("== independent F_101 split-fibre degree-27 audit ==")
print(f"target={target}; level-one points={level_one}")
print(f"level-two blocks={level_two_blocks}")
print("fibre census: 3 first preimages and 3+3+3 second preimages")
print("three degree-nine blocks: squarefree and pairwise coprime")
print("direct nine-cubic product: degree 27 and squarefree")
print(f"degree-27 ascending-coefficient sha256={digest}")
print("therefore generic x-coordinate separation/coprimality is nonempty")
print("scope: finite-field genericity certificate only; no discriminant expansion")
print("all independent audit checks passed")
