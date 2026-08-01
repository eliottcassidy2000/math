#!/usr/bin/env python3
"""Exact finite controls for THM-3042's common-quotient owner criterion."""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


N = 12
R = tuple(range(N))
PRODUCT = tuple(product(R, repeat=2))
DIVISORS = tuple(d for d in range(1, N + 1) if N % d == 0)


def add(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return ((x[0] + y[0]) % N, (x[1] + y[1]) % N)


def mul(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return ((x[0] * y[0]) % N, (x[1] * y[1]) % N)


def congruence_order(g: int) -> frozenset[tuple[int, int]]:
    require(N % g == 0, f"{g} does not divide {N}")
    return frozenset((r, c) for r, c in PRODUCT if (r - c) % g == 0)


def generated_by_one(d: int, b: tuple[int, int]) -> frozenset[tuple[int, int]]:
    base = congruence_order(d)
    return frozenset(add(a, ((k * b[0]) % N, (k * b[1]) % N)) for a in base for k in R)


def generated_by_two(
    d: int, b1: tuple[int, int], b2: tuple[int, int]
) -> frozenset[tuple[int, int]]:
    base = congruence_order(d)
    return frozenset(
        add(
            a,
            (
                (k * b1[0] + ell * b2[0]) % N,
                (k * b1[1] + ell * b2[1]) % N,
            ),
        )
        for a in base
        for k in R
        for ell in R
    )


def check_ring(B: frozenset[tuple[int, int]], label: str) -> None:
    require((0, 0) in B and (1, 1) in B, f"{label}: zero/unit missing")
    require(all((r, r) in B for r in R), f"{label}: diagonal missing")
    for x in B:
        require(((-x[0]) % N, (-x[1]) % N) in B, f"{label}: inverse missing")
        for y in B:
            require(add(x, y) in B, f"{label}: addition not closed")
            require(mul(x, y) in B, f"{label}: multiplication not closed")


def kernel_first(B: frozenset[tuple[int, int]]) -> frozenset[int]:
    return frozenset(r for r in R if (r, 0) in B)


def kernel_second(B: frozenset[tuple[int, int]]) -> frozenset[int]:
    return frozenset(c for c in R if (0, c) in B)


def ideal_multiples(g: int) -> frozenset[int]:
    return frozenset((g * k) % N for k in R)


def conductor_singleton(B: frozenset[tuple[int, int]]) -> frozenset[int]:
    return frozenset(
        s
        for s in R
        if (s, 0) in B
        and all(mul((s, 0), z) in B for z in PRODUCT)
    )


# Every congruence order is a literal subdirect fibre product.  Its two
# kernels and conductor singleton slice are the same ideal (g).
for g in DIVISORS:
    B = congruence_order(g)
    check_ring(B, f"B_{g}")
    require({r for r, _ in B} == set(R), f"B_{g}: first projection")
    require({c for _, c in B} == set(R), f"B_{g}: second projection")
    ideal = ideal_multiples(g)
    require(kernel_first(B) == ideal, f"B_{g}: first kernel")
    require(kernel_second(B) == ideal, f"B_{g}: second kernel")
    require(conductor_singleton(B) == ideal, f"B_{g}: conductor slice")
    require(((1, 0) in B) == (g == 1), f"B_{g}: singleton criterion")


# Exhaust every one-generator enlargement of every congruence order.
single_cases = 0
for d in DIVISORS:
    for b in PRODUCT:
        delta = b[0] - b[1]
        g = gcd(N, d, delta)
        actual = generated_by_one(d, b)
        predicted = congruence_order(g)
        require(actual == predicted, f"one-generator mismatch d={d}, b={b}, g={g}")
        check_ring(actual, f"generated d={d}, b={b}")
        single_cases += 1


# Selected two-generator tests include a jointly Bezout but individually
# nonsplitting pair.
two_controls = (
    (12, (0, 4), (0, 9)),
    (12, (0, 4), (0, 6)),
    (6, (2, 5), (7, 1)),
    (4, (0, 2), (1, 1)),
)
for d, b1, b2 in two_controls:
    g = gcd(N, d, b1[0] - b1[1], b2[0] - b2[1])
    actual = generated_by_two(d, b1, b2)
    predicted = congruence_order(g)
    require(actual == predicted, f"two-generator mismatch d={d}, {b1}, {b2}")
    check_ring(actual, f"two-generator d={d}, {b1}, {b2}")


# Sharp preserve / shrink / split controls for d=4.
sharp = []
for label, b, expected_g in (
    ("preserve", (1, 1), 4),
    ("shrink", (0, 2), 2),
    ("split", (0, 1), 1),
):
    B = generated_by_one(4, b)
    require(B == congruence_order(expected_g), f"{label}: wrong quotient")
    require(((1, 0) in B) == (expected_g == 1), f"{label}: wrong owner bit")
    sharp.append((label, expected_g, len(B), int((1, 0) in B)))


semantic = "|".join(
    [f"N={N}", f"div={','.join(map(str, DIVISORS))}", f"single={single_cases}"]
    + [f"{label}:{g}:{size}:{owner}" for label, g, size, owner in sharp]
)

print("THM3042_SUBDIRECT_GRAPH_ORDER_COMMON_QUOTIENT_EXACT")
print(f"finite_ring_modulus {N}")
print(f"congruence_subdirect_products {len(DIVISORS)}")
print(f"one_generator_enlargements {single_cases}")
print(f"two_generator_controls {len(two_controls)}")
for label, g, size, owner in sharp:
    print(f"sharp_{label} quotient_size {g} order_size {size} singleton {owner}")
print(f"semantic_sha256 {sha256(semantic.encode()).hexdigest()}")
print("all_runtime_checks_explicit 1")
