#!/usr/bin/env python3
"""Square/triangular tag sidecar for THM-3825's prime-colour decoder."""

from __future__ import annotations

import ast
from collections import Counter
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True:
        raise RuntimeError(message)


def factor(n: int) -> tuple[tuple[int, int], ...]:
    factors: list[tuple[int, int]] = []
    divisor = 2
    while divisor * divisor <= n:
        exponent = 0
        while n % divisor == 0:
            n //= divisor
            exponent += 1
        if exponent:
            factors.append((divisor, exponent))
        divisor += 1
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


def admissible_sum(total: int) -> bool:
    return total >= 3 and all(
        prime % 3 == 2 and exponent <= 2
        for prime, exponent in factor(total)
    )


def triangular_index(value: int) -> int | None:
    discriminant = 1 + 8 * value
    root = isqrt(discriminant)
    if root * root != discriminant or root % 2 == 0:
        return None
    return (root - 1) // 2


def fundamental_pell(discriminant: int) -> tuple[int, int]:
    """Least positive solution of x^2-discriminant*y^2=1."""
    root = isqrt(discriminant)
    if root * root == discriminant:
        raise RuntimeError("Pell discriminant must be nonsquare")
    remainder = 0
    denominator = 1
    partial = root
    x_previous, x = 1, partial
    y_previous, y = 0, 1
    while x * x - discriminant * y * y != 1:
        remainder = denominator * partial - remainder
        denominator = (discriminant - remainder * remainder) // denominator
        partial = (root + remainder) // denominator
        x_previous, x = x, partial * x + x_previous
        y_previous, y = y, partial * y + y_previous
    return x, y


def pure_power_of_three(value: int) -> bool:
    while value % 3 == 0:
        value //= 3
    return value == 1


EXPECTED_TRIANGULAR_EXPONENTS = [
    (1, 3, 0, 7),
    (2, 3, 1, 14),
    (1, 10, 1, 77),
    (2, 9, 1, 66),
    (3, 8, 2, 98),
    (9, 13, 0, 76),
    (9, 13, 1, 132),
    (5, 18, 5, 1701),
    (10, 13, 1, 138),
    (3, 43, 2, 1196),
    (13, 33, 2, 828),
    (8, 45, 1, 741),
    (7, 52, 4, 4778),
    (33, 38, 2, 1278),
    (37, 45, 0, 532),
    (39, 46, 1, 969),
    (41, 44, 2, 1665),
    (21, 73, 0, 892),
    (37, 63, 0, 775),
    (39, 68, 1, 1497),
    (6, 109, 0, 1609),
    (9, 136, 1, 3885),
    (66, 113, 3, 9666),
    (37, 150, 0, 2617),
    (46, 159, 0, 2869),
    (101, 113, 3, 11556),
    (25, 264, 2, 18206),
    (2, 293, 5, 110565),
    (106, 189, 3, 20709),
    (53, 258, 6, 158921),
    (5, 314, 2, 23606),
]


sums = tuple(total for total in range(3, 357) if admissible_sum(total))
atlas = tuple(
    (first, total - first)
    for total in sums
    for first in range(1, (total + 1) // 2)
    if first < total - first and gcd(first, total - first) == 1
)

gate(len(sums) == 94, "admissible-sum census")
gate(len(atlas) == 5855, "primitive-pair census")

bases: dict[int, tuple[int, int]] = {}
for first, second in atlas:
    base = first**3 + second**3
    gate(base % 3 != 0, "prime-colour base is 3-free")
    gate(base not in bases, "primitive prime-colour base is injective")
    bases[base] = (first, second)

gate(len(bases) == 5855, "primitive-value census")
gate(
    tuple(
        (first, second, isqrt(first**3 + second**3))
        for first, second in atlas
        if isqrt(first**3 + second**3) ** 2 == first**3 + second**3
    )
    == ((56, 65, 671),),
    "unique primitive square base",
)

powers = tuple(3**exponent for exponent in range(156))
unordered_square_hits: list[tuple[int, int, int]] = []
oriented_square_hits: list[tuple[int, int, int]] = []
unordered_triangular_hits: list[tuple[int, int, int, int]] = []
oriented_triangular_hits: list[tuple[int, int, int, int]] = []

for base, (first, second) in bases.items():
    for kappa in range(78):
        address = powers[kappa] * base
        square = isqrt(address) ** 2 == address
        expected_square = (first, second) == (56, 65) and kappa % 2 == 0
        gate(square == expected_square, "unordered square/parity equivalence")
        if square:
            unordered_square_hits.append((first, second, kappa))
        index = triangular_index(address)
        if index is not None:
            unordered_triangular_hits.append((first, second, kappa, index))

    for exponent in range(156):
        address = powers[exponent] * base
        square = isqrt(address) ** 2 == address
        expected_square = (first, second) == (56, 65) and exponent % 2 == 0
        gate(square == expected_square, "oriented square/orientation equivalence")
        if square:
            oriented_square_hits.append((first, second, exponent))
        index = triangular_index(address)
        if index is not None:
            oriented_triangular_hits.append((first, second, exponent, index))

gate(
    unordered_square_hits
    == [(56, 65, exponent) for exponent in range(0, 78, 2)],
    "complete unordered square-tag census",
)
gate(
    oriented_square_hits
    == [(56, 65, exponent) for exponent in range(0, 156, 2)],
    "complete oriented square-tag census",
)
gate(
    unordered_triangular_hits == EXPECTED_TRIANGULAR_EXPONENTS,
    "complete unordered triangular-tag census",
)
gate(
    oriented_triangular_hits == EXPECTED_TRIANGULAR_EXPONENTS,
    "complete oriented triangular-tag census",
)

for first, second, exponent, index in EXPECTED_TRIANGULAR_EXPONENTS:
    base = first**3 + second**3
    address = powers[exponent] * base
    gate(address == index * (index + 1) // 2,
         "triangular tagged-address reconstruction")
    gate(isqrt(address) ** 2 != address,
         "no square-triangular tagged address")
    pell_label, pell_sheet = divmod(exponent, 2)
    gate(
        (2 * index + 1) ** 2
        - 8 * 3**pell_sheet * base * (3**pell_label) ** 2
        == 1,
        "two-sheet Pell reconstruction",
    )
    discriminant = 8 * 3**pell_sheet * base
    pell_point = (2 * index + 1, 3**pell_label)
    gate(fundamental_pell(discriminant) == pell_point,
         "triangular tag is the fundamental Pell point")
    next_x = pell_point[0] ** 2 + discriminant * pell_point[1] ** 2
    next_y = 2 * pell_point[0] * pell_point[1]
    gate(next_x * next_x - discriminant * next_y * next_y == 1,
         "next Pell point reconstruction")
    gate(next_y == 2 * pell_point[0] * pell_point[1],
         "next Pell ordinate formula")
    gate(not pure_power_of_three(next_y),
         "immediate exit from the pure-3 Pell ordinate")

exponent_histogram = dict(sorted(Counter(
    exponent for _, _, exponent, _ in EXPECTED_TRIANGULAR_EXPONENTS
).items()))
gate(exponent_histogram == {0: 8, 1: 9, 2: 7, 3: 3, 4: 1, 5: 2, 6: 1},
     "triangular exponent histogram")
gate(sum(exponent % 2 == 0 for _, _, exponent, _ in EXPECTED_TRIANGULAR_EXPONENTS) == 17,
     "oriented triangular epsilon-zero count")
gate(sum(exponent % 2 == 1 for _, _, exponent, _ in EXPECTED_TRIANGULAR_EXPONENTS) == 14,
     "oriented triangular epsilon-one count")

triangular_lookup = {
    (first, second, exponent): index
    for first, second, exponent, index in EXPECTED_TRIANGULAR_EXPONENTS
}
triangular_successors = [
    (first, second, exponent, index,
     triangular_lookup[(first, second, exponent + 1)])
    for first, second, exponent, index in EXPECTED_TRIANGULAR_EXPONENTS
    if (first, second, exponent + 1) in triangular_lookup
]
gate(triangular_successors == [(9, 13, 0, 76, 132)],
     "unique triangular label-successor edge")
gate(3 * (76 * 77 // 2) == 132 * 133 // 2,
     "triangular tripling reconstruction")
gate((2 * 132 + 1) ** 2 - 3 * (2 * 76 + 1) ** 2 == -2,
     "negative-Pell triangular-tripling address")

semantic = {
    "universe": {
        "primitive_pairs": 5855,
        "unordered_channels": 78,
        "oriented_channels": 156,
    },
    "squares": {
        "primitive_pair": [56, 65],
        "unordered_count": 39,
        "unordered_rule": "kappa even",
        "oriented_count": 78,
        "oriented_rule": "epsilon=0",
    },
    "triangular_exponent_hits": EXPECTED_TRIANGULAR_EXPONENTS,
    "triangular_exponent_histogram": exponent_histogram,
    "pell_reframe": "X^2-8*3^epsilon*m*Y^2=1; Y=3^kappa",
    "pell_fundamental_hits": 31,
    "pell_first_step": "Y_2=2XY is not a power of three",
    "triangular_label_successor": [9, 13, 0, 1, 76, 132],
    "scope": "finite arithmetic tag sidecar only; no loneliness implication",
}
semantic_blob = json.dumps(
    semantic, sort_keys=True, separators=(",", ":")
).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3825-prime-colour-square-triangular-tag-sidecar")
print("universe=primitive_pairs:5855;unordered_channels:78;oriented_channels:156")
print("square_unordered=pair:(56,65);kappa_even;count=39")
print("square_oriented=pair:(56,65);epsilon=0;count=78")
print("triangular_shared_exponent_hits=" + ";".join(
    f"({first},{second};e={exponent};T_{index})"
    for first, second, exponent, index in EXPECTED_TRIANGULAR_EXPONENTS
))
print("triangular_counts=unordered:31;oriented:31")
print("triangular_exponent_histogram=0:8,1:9,2:7,3:3,4:1,5:2,6:1")
print("oriented_triangular_bits=epsilon0:17;epsilon1:14")
print("triangular_pell_sheets=X^2-8*3^epsilon*m*Y^2=1;Y=3^kappa")
print("triangular_pell_fundamental_hits=31;next_Y=2XY;pure3_returns_at_next_step=none")
print("triangular_label_successor=(9,13;e:0->1;T_76->T_132);count=1;negative_pell=265^2-3*153^2=-2")
print("square_triangular_tagged_addresses=none")
print("scope=finite_arithmetic_tag_sidecar_only;no_LRC14_row_excluded")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
