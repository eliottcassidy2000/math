#!/usr/bin/env python3
"""Square/triangular finite sidecar for THM-3818's cube-address atlas."""

from __future__ import annotations

import ast
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
    root = isqrt(1 + 8 * value)
    if root * root != 1 + 8 * value or root % 2 == 0:
        return None
    return (root - 1) // 2


sums = tuple(total for total in range(3, 357) if admissible_sum(total))
atlas = tuple(
    (p, total - p)
    for total in sums
    for p in range(1, (total + 1) // 2)
    if p < total - p and gcd(p, total - p) == 1
)

gate(len(sums) == 94, "admissible-sum census")
gate(len(atlas) == 5855, "primitive-pair census")
for p, q in atlas:
    gate(p < q and gcd(p, q) == 1 and p + q <= 356,
         "declared primitive-pair universe")
    gate(admissible_sum(p + q), "declared inert cube-free sum universe")
    cofactor = p * p - p * q + q * q
    gate(gcd(p + q, cofactor) == 1, "coprime cube-sum factors")

square_shells = tuple(total for total in sums if isqrt(total) ** 2 == total)
triangular_shells = tuple(
    total for total in sums if triangular_index(total) is not None
)
square_shell_rows = sum(p + q in square_shells for p, q in atlas)
triangular_shell_rows = sum(p + q in triangular_shells for p, q in atlas)

gate(square_shells == (4, 25, 100, 121, 289), "square sum shells")
gate(square_shell_rows == 222, "square-shell ratio count")
gate(triangular_shells == (10, 55, 253), "triangular sum shells")
gate(triangular_shell_rows == 132, "triangular-shell ratio count")

square_addresses: list[tuple[int, int, int]] = []
triangular_addresses: list[tuple[int, int, int]] = []
for p, q in atlas:
    address = p**3 + q**3
    square_root = isqrt(address)
    if square_root * square_root == address:
        square_addresses.append((p, q, square_root))
    index = triangular_index(address)
    if index is not None:
        triangular_addresses.append((p, q, index))

expected_triangular_addresses = [
    (1, 3, 7),
    (9, 13, 76),
    (37, 45, 532),
    (21, 73, 892),
    (37, 63, 775),
    (6, 109, 1609),
    (37, 150, 2617),
    (46, 159, 2869),
]
gate(square_addresses == [(56, 65, 671)], "unique square cube address")
gate(triangular_addresses == expected_triangular_addresses,
     "complete triangular cube-address census")
for p, q, index in triangular_addresses:
    gate(p**3 + q**3 == index * (index + 1) // 2,
         "triangular-address reconstruction")
    gate(isqrt(p**3 + q**3) ** 2 != p**3 + q**3,
         "no square-triangular cube address")
gate(56**3 + 65**3 == 671**2 == 450241,
     "square-address reconstruction")

# THM-3819's square-triangular/Pell pair recurrence, stopped at the exact
# THM-3743 l1 cap.  Only the middle pair belongs to the inert atlas.
pell_pairs: list[tuple[int, int]] = []
a, b = 1, 5
while a + b <= 356:
    pell_pairs.append((a, b))
    a, b = b, 6 * b - a
pell_intersection = tuple(pair for pair in pell_pairs if pair in set(atlas))
gate(pell_pairs == [(1, 5), (5, 29), (29, 169)], "capped Pell pairs")
gate(pell_intersection == ((5, 29),), "Pell/inert-atlas intersection")
gate(8 * 9 // 2 == 6**2 == 36, "square-triangular conductor identity")
gate(12 * 13 // 2 == 78, "triangular labelled-placement identity")
gate(5**3 + 29**3 == 24514, "merged Pell cube address")

semantic = {
    "universe": "THM-3818: 94 admissible sums, 5855 primitive pairs",
    "sum_shells": {
        "square": [4, 25, 100, 121, 289],
        "square_rows": 222,
        "triangular": [10, 55, 253],
        "triangular_rows": 132,
    },
    "square_address": [56, 65, 671],
    "triangular_addresses": expected_triangular_addresses,
    "square_triangular_addresses": [],
    "pell_intersection": [5, 29],
    "scope": "finite arithmetic sidecar only; no loneliness implication",
}
semantic_blob = json.dumps(
    semantic, sort_keys=True, separators=(",", ":")
).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3818-square-triangular-cube-address-sidecar")
print("universe=admissible_sums:94;primitive_pairs:5855")
print("square_sum_shells=4,25,100,121,289;rows=222")
print("triangular_sum_shells=10,55,253;rows=132")
print("square_addresses=(56,65;sqrt=671);count=1")
print("triangular_addresses=" + ";".join(
    f"({p},{q};T_{index})" for p, q, index in triangular_addresses
))
print("square_triangular_addresses=none")
print("pell_cap_pairs=(1,5),(5,29),(29,169);atlas_intersection=(5,29)")
print("classic_join=T_8=36=6^2;T_12=78;cube_address_5_29=24514")
print("scope=finite_arithmetic_sidecar_only;no_LRC14_row_excluded")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
