#!/usr/bin/env python3
"""Exact two-jet locality and Catalan sharpness controls for THM-3214."""

from __future__ import annotations

import ast
import hashlib
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = (
    ROOT
    / "01-canon"
    / "theorems"
    / "THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return.md"
)
EXPECTED_DEPENDENCY_SHA256 = (
    "69ae0b5e4c9526d643bc976220a4752e33632248c80394ca4cf9c9620105878f"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def add(left: list[int], right: list[int]) -> list[int]:
    length = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(length)
    ]


def scale(value: int, coefficients: list[int]) -> list[int]:
    return [value * coefficient for coefficient in coefficients]


def multiply(
    left: list[int], right: list[int], length: int | None = None
) -> list[int]:
    if length is None:
        length = len(left) + len(right) - 1
    answer = [0] * length
    for i, a_value in enumerate(left):
        for j, b_value in enumerate(right):
            if i + j >= length:
                break
            answer[i + j] += a_value * b_value
    return answer


def power(coefficients: list[int], exponent: int, length: int) -> list[int]:
    answer = [1] + [0] * (length - 1)
    factor = coefficients[:length]
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply(answer, factor, length)
        factor = multiply(factor, factor, length)
        remaining //= 2
    return answer


def p1(f: list[int], g: list[int]) -> int:
    return f[0] * g[1] - f[1] * g[0]


def pseudo_by_numerator(f: list[int], g: list[int]) -> list[int]:
    """Compute E(f,g) by forming the numerator and deleting two zeros."""
    length = min(len(f), len(g))
    f = f[:length]
    g = g[:length]
    connection = p1(f, g)
    numerator = add(scale(f[0] ** 2, g), scale(-f[0] * g[0], f))
    shifted = [0] + scale(-connection, f[:-1])
    numerator = add(numerator, shifted)[:length]
    require(numerator[:2] == [0, 0], ("pseudo numerator", numerator[:2]))
    return numerator[2:]


def closed_coefficient(f: list[int], g: list[int], index: int) -> int:
    return (
        f[0] ** 2 * g[index + 2]
        - f[0] * g[0] * f[index + 2]
        - p1(f, g) * f[index + 1]
    )


def chain(f: list[int], g: list[int], steps: int) -> list[list[int]]:
    """Return [R_-1,R_0,...,R_steps] for R_-1=g and R_0=f."""
    values = [g[:], f[:]]
    for _ in range(steps):
        values.append(pseudo_by_numerator(values[-1], values[-2]))
    return values


def epsilon(index: int) -> int:
    values = { -1: 1, 0: 1 }
    for current in range(index):
        values[current + 1] = -values[current - 1]
    return values[index]


dependency_hash = lf_sha256(DEPENDENCY)
require(
    dependency_hash == EXPECTED_DEPENDENCY_SHA256,
    ("dependency hash", dependency_hash),
)


# Formula (4), calculated independently from numerator formation.
COEFFICIENT_TABLES = 11
coefficient_checks = 0
for seed in range(COEFFICIENT_TABLES):
    f = [(-1) ** (seed + index) * ((seed + 2) * (index + 1) + index**2)
         for index in range(13)]
    g = [(-1) ** index * ((seed + 3) * (index + 2) - 2 * index + 1)
         for index in range(13)]
    result = pseudo_by_numerator(f, g)
    for index, value in enumerate(result):
        require(value == closed_coefficient(f, g, index),
                ("coefficient formula", seed, index))
        coefficient_checks += 1


# Formula (9): alter every coefficient strictly beyond the required jet.
locality_pairs = 0
for steps in range(7):
    for height in range(9):
        cutoff = height + 2 * steps
        last_index = cutoff + 6
        f = [(-1) ** index * (index + 2) for index in range(last_index + 1)]
        g = [(-1) ** (index + 1) * (2 * index + 3)
             for index in range(last_index + 1)]
        f_hostile = f[:]
        g_hostile = g[:]
        for index in range(cutoff + 1, last_index + 1):
            f_hostile[index] += (index + 5) ** 3
            g_hostile[index] -= (index + 7) ** 2
        left = chain(f, g, steps)[steps + 1]
        right = chain(f_hostile, g_hostile, steps)[steps + 1]
        require(left[: height + 1] == right[: height + 1],
                ("jet locality", steps, height))
        locality_pairs += 1
require(locality_pairs == 63, ("locality count", locality_pairs))


# Catalan coefficients and the closed regular orbit (14).
MAX_SERIES_INDEX = 36
catalan = [comb(2 * index, index) // (index + 1)
           for index in range(MAX_SERIES_INDEX + 1)]
one = [1] + [0] * MAX_SERIES_INDEX
catalan_chain = chain(catalan, one, 8)
catalan_orbit_coefficients = 0
for index in range(9):
    actual = catalan_chain[index + 1]
    expected = scale(
        epsilon(index), power(catalan, 2 * index + 1, len(actual))
    )
    require(actual == expected, ("Catalan orbit", index))
    require(actual[0] == epsilon(index), ("Catalan pivot", index))
    catalan_orbit_coefficients += len(actual)


# Pivot sharpness (18)--(21), including every displayed leading perturbation.
pivot_leading_checks = 0
for depth in range(1, 11):
    last_index = 2 * depth + 8
    base = catalan[: last_index + 1]
    perturbed = base[:]
    perturbed[2 * depth] += 1
    base_chain = chain(base, one[: last_index + 1], depth)
    perturbed_chain = chain(perturbed, one[: last_index + 1], depth)
    for index in range(depth + 1):
        difference = [
            perturbed_chain[index + 1][place]
            - base_chain[index + 1][place]
            for place in range(len(base_chain[index + 1]))
        ]
        first = 2 * (depth - index)
        require(difference[:first] == [0] * first,
                ("pivot perturbation order", depth, index))
        require(difference[first] == epsilon(index),
                ("pivot perturbation lead", depth, index))
        pivot_leading_checks += 1
    require(
        [row[0] for row in perturbed_chain[1:depth + 1]]
        == [row[0] for row in base_chain[1:depth + 1]],
        ("earlier pivots", depth),
    )
    require(base_chain[depth + 1][0] == epsilon(depth),
            ("base kth pivot", depth))
    require(perturbed_chain[depth + 1][0] == 2 * epsilon(depth),
            ("sharp kth pivot", depth))


# Connection sharpness (22)--(24), plus the k=0 boundary.
connection_leading_checks = 0
connection_steps = 0
base_zero = catalan[:10]
perturbed_zero = base_zero[:]
perturbed_zero[1] += 1
require(p1(perturbed_zero, one[:10]) - p1(base_zero, one[:10]) == -1,
        "connection k=0")
connection_steps += 1
for depth in range(1, 11):
    last_index = 2 * depth + 9
    base = catalan[: last_index + 1]
    perturbed = base[:]
    perturbed[2 * depth + 1] += 1
    base_chain = chain(base, one[: last_index + 1], depth)
    perturbed_chain = chain(perturbed, one[: last_index + 1], depth)
    for index in range(depth + 1):
        difference = [
            perturbed_chain[index + 1][place]
            - base_chain[index + 1][place]
            for place in range(len(base_chain[index + 1]))
        ]
        first = 2 * (depth - index) + 1
        require(difference[:first] == [0] * first,
                ("connection perturbation order", depth, index))
        require(difference[first] == epsilon(index),
                ("connection perturbation lead", depth, index))
        connection_leading_checks += 1
    require(
        [row[0] for row in perturbed_chain[1:depth + 2]]
        == [row[0] for row in base_chain[1:depth + 2]],
        ("connection pivots unchanged", depth),
    )
    base_connection = p1(base_chain[depth + 1], base_chain[depth])
    hostile_connection = p1(
        perturbed_chain[depth + 1], perturbed_chain[depth]
    )
    expected_difference = -epsilon(depth) * epsilon(depth - 1)
    require(hostile_connection - base_connection == expected_difference,
            ("sharp connection", depth))
    connection_steps += 1


source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.parse(source).body),
        "assert-sensitive top-level test")

print("dependency_thm3192_sha256=" + dependency_hash)
print("universal_coefficient_tables=" + str(COEFFICIENT_TABLES))
print("universal_coefficient_checks=" + str(coefficient_checks))
print("locality_pairs=" + str(locality_pairs))
print("catalan_orbit_steps=8")
print("catalan_orbit_coefficients=" + str(catalan_orbit_coefficients))
print("catalan_signs=" + ",".join(str(epsilon(index)) for index in range(11)))
print("pivot_sharpness_steps=10")
print("pivot_leading_checks=" + str(pivot_leading_checks))
print("connection_sharpness_steps=" + str(connection_steps))
print("connection_leading_checks=" + str(connection_leading_checks))
print("jet_budget=(pivot_k:2k,connection_k:2k+1)")
print("regular_hostile=(C,C+z^(2k),C+z^(2k+1))")
print("scope=sharp_universal_information_budget_not_simultaneous_chart_nonvanishing")
