#!/usr/bin/env python3
"""Exact companion for proved THM-3484.

The script treats the three explicit degree-seven residue polynomials as a
self-contained ternary sequence.  It derives the compressed order-22
annihilator, proves finite linear complexity at least 22 with an exact Hankel
determinant, and audits the cubic Fourier degree drop.  Identifying this word
with the private-gradient determinant is the proved THM-3482 specialization;
no LRC current or LRC(14) assertion is made here.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from json import dumps
from math import comb
from pathlib import Path
import sys


STATE_POLYNOMIALS = {
    0: (0, 0, 0, 0, -2048, 4096, 8192, -16384),
    1: (0, 0, 0, 0, 2048, 4096, -24576, -16384),
    2: (0, -256, 3072, -15360, 40960, -61440, 49152, -16384),
}
EXPECTED_DENOMINATOR = (
    1, -1, 0, -7, 7, 0, 21, -21, 0, -35, 35, 0,
    35, -35, 0, -21, 21, 0, 7, -7, 0, -1, 1,
)
EXPECTED_NUMERATOR = (
    0, -34816, -338432, -28657152, -335106048, -313495296,
    -3294224640, -9788725248, -4813545216, -26592198912,
    -36296777728, -11030363648, -39322132992, -27167668224,
    -5002493952, -11956875264, -3768176640, -338608896,
    -502246656, -40216576, -239360, -186624,
)
EXPECTED_HANKEL_DETERMINANT = -(
    2 ** 382 * 3 ** 191 * 5 ** 22 * 7 ** 8 * 61 ** 7
)
EXPECTED_SEMANTIC_SHA256 = "0c9cac55b66f8a8aa4241ad728e3336dd000cff25544f3df7eef9af120012518"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def poly_eval(poly: tuple[int, ...], value: int) -> int:
    total = 0
    for coefficient in reversed(poly):
        total = total * value + coefficient
    return total


def determinant(matrix: list[list[int]]) -> int:
    size = len(matrix)
    require(size and all(len(row) == size for row in matrix), "square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (row for row in range(pivot_index, size) if work[row][pivot_index]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = work[row][column] * pivot - work[row][pivot_index] * work[pivot_index][column]
                require(numerator % previous == 0, (pivot_index, numerator, previous))
                work[row][column] = numerator // previous
        previous = pivot
        for row in range(pivot_index + 1, size):
            work[row][pivot_index] = 0
    return sign * work[-1][-1]


def sequence_value(index: int) -> int:
    require(index >= 0, index)
    return poly_eval(STATE_POLYNOMIALS[index % 3], index)


def poly_multiply(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    values = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            values[i + j] += a * b
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def poly_power(base: tuple[int, ...], exponent: int) -> tuple[int, ...]:
    result = (1,)
    for _unused in range(exponent):
        result = poly_multiply(result, base)
    return result


def recurrence_residual(coefficients: tuple[int, ...], values: tuple[int, ...], index: int) -> int:
    return sum(coefficients[j] * values[index - j] for j in range(len(coefficients)))


def first_failure(coefficients: tuple[int, ...], values: tuple[int, ...]) -> tuple[int, int] | None:
    for index in range(len(coefficients) - 1, len(values)):
        residual = recurrence_residual(coefficients, values, index)
        if residual:
            return index, residual
    return None


def berlekamp_massey(values: tuple[int, ...]) -> tuple[int, tuple[Fraction, ...]]:
    connection = [Fraction(1)]
    previous_connection = [Fraction(1)]
    length = 0
    displacement = 1
    previous_discrepancy = Fraction(1)
    for index, raw_value in enumerate(values):
        discrepancy = Fraction(raw_value)
        for j in range(1, length + 1):
            discrepancy += connection[j] * values[index - j]
        if discrepancy == 0:
            displacement += 1
            continue
        old_connection = connection[:]
        scale = -discrepancy / previous_discrepancy
        needed = len(previous_connection) + displacement
        if len(connection) < needed:
            connection.extend(Fraction(0) for _unused in range(needed - len(connection)))
        for j, coefficient in enumerate(previous_connection):
            connection[j + displacement] += scale * coefficient
        if 2 * length <= index:
            length = index + 1 - length
            previous_connection = old_connection
            previous_discrepancy = discrepancy
            displacement = 1
        else:
            displacement += 1
    return length, tuple(connection)


def omega_add(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return left[0] + right[0], left[1] + right[1]


def omega_scale(value: tuple[Fraction, Fraction], scalar: int | Fraction) -> tuple[Fraction, Fraction]:
    return value[0] * scalar, value[1] * scalar


def fourier_components() -> tuple[tuple[tuple[Fraction, Fraction], ...], ...]:
    one = (Fraction(1), Fraction(0))
    omega = (Fraction(0), Fraction(1))
    omega_squared = (Fraction(-1), Fraction(-1))
    weights = (
        (one, one, one),
        (one, omega_squared, omega),
        (one, omega, omega_squared),
    )
    components = []
    for colour in range(3):
        coefficients = []
        for degree in range(8):
            total = (Fraction(0), Fraction(0))
            for state in range(3):
                total = omega_add(
                    total,
                    omega_scale(weights[colour][state], STATE_POLYNOMIALS[state][degree]),
                )
            coefficients.append(omega_scale(total, Fraction(1, 3)))
        components.append(tuple(coefficients))
    return tuple(components)


def component_degree(component: tuple[tuple[Fraction, Fraction], ...]) -> int:
    return max(index for index, value in enumerate(component) if value != (0, 0))


def security_report(source: Path) -> tuple[str, ...]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert present")
    forbidden = {
        "eval", "exec", "compile", "open", "system", "popen", "run", "Popen",
        "write_text", "write_bytes", "unlink", "remove", "rename",
    }
    called = {
        node.func.id for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    called.update(
        node.func.attr for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not called & forbidden, ("forbidden calls", sorted(called & forbidden)))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {"__future__", "ast", "fractions", "hashlib", "json", "math", "pathlib", "sys"}
    require(imports <= allowed, ("unexpected imports", sorted(imports - allowed)))
    return tuple(sorted(imports))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    denominator = [0] * 23
    for j in range(8):
        coefficient = (-1) ** j * comb(7, j)
        denominator[3 * j] += coefficient
        denominator[3 * j + 1] -= coefficient
    denominator = tuple(denominator)
    require(denominator == EXPECTED_DENOMINATOR, denominator)

    values = tuple(sequence_value(index) for index in range(5001))
    require(all(recurrence_residual(denominator, values, index) == 0
                for index in range(22, len(values))), "order-22 recurrence")

    numerator = tuple(
        sum(denominator[j] * values[index - j]
            for j in range(min(index, len(denominator) - 1) + 1))
        for index in range(22)
    )
    require(numerator == EXPECTED_NUMERATOR, numerator)

    naive = poly_power((1, 0, 0, -1), 8)
    require(len(naive) - 1 == 24, len(naive) - 1)
    require(first_failure(naive, values) is None, "naive order-24 annihilator")
    trivial_short = poly_multiply(poly_power((1, -1), 7), poly_power((1, 1, 1), 7))
    cubic_short = poly_multiply(poly_power((1, -1), 8), poly_power((1, 1, 1), 6))
    short_failures = (first_failure(trivial_short, values), first_failure(cubic_short, values))
    require(all(failure is not None for failure in short_failures), short_failures)

    bm_length, bm_connection = berlekamp_massey(values[:100])
    require(bm_length == 22, bm_length)
    require(bm_connection == tuple(Fraction(value) for value in denominator), bm_connection)

    hankel = [[values[row + column] for column in range(22)] for row in range(22)]
    hankel_determinant = determinant(hankel)
    require(hankel_determinant == EXPECTED_HANKEL_DETERMINANT, hankel_determinant)
    hankel_hash = sha256(str(hankel_determinant).encode("ascii")).hexdigest()

    components = fourier_components()
    degrees = tuple(component_degree(component) for component in components)
    leads = tuple(component[degree] for component, degree in zip(components, degrees))
    require(degrees == (7, 6, 6), degrees)
    require(leads[0] == (Fraction(-16384), Fraction(0)), leads[0])
    require(leads[1] == (Fraction(32768, 3), Fraction(24576)), leads[1])
    require(leads[2] == (Fraction(-40960, 3), Fraction(-24576)), leads[2])

    value_digest = sha256(dumps(values, separators=(",", ":")).encode("ascii")).hexdigest()
    semantic_payload = {
        "state_polynomials": STATE_POLYNOMIALS,
        "denominator": denominator,
        "numerator": numerator,
        "naive_degree": 24,
        "minimal_degree": 22,
        "short_failures": short_failures,
        "bm_length": bm_length,
        "hankel_factorization": ((2, 382), (3, 191), (5, 22), (7, 8), (61, 7)),
        "hankel_sign": -1,
        "hankel_sha256": hankel_hash,
        "fourier_degrees": degrees,
        "fourier_leads": leads,
        "value_window": (0, 5000),
        "value_digest": value_digest,
        "security": security_report(Path(__file__)),
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"), default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, semantic_hash)

    print("THM-3484 TERNARY WEIGHTED-DETERMINANT MINIMAL RECURRENCE EXACT COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print(f"STATE_POLYNOMIALS_CONSTANT_FIRST: {STATE_POLYNOMIALS}")
    print("GENERATING_DENOMINATOR: (1-z)*(1-z^3)^7 = (1-z)^8*(1+z+z^2)^7")
    print(f"ORDER22_RECURRENCE_COEFFICIENTS: {denominator}")
    print(f"GENERATING_NUMERATOR_COEFFICIENTS: {numerator}")
    print(f"CUBIC_FOURIER_DEGREES_AND_LEADS_IN_Q(omega): {(degrees, leads)}")
    print(f"SHORTENED_FACTOR_HOSTILES_FIRST_FAILURES: {short_failures}")
    print(f"BERLEKAMP_MASSEY_FIRST100: order={bm_length} coefficients={bm_connection}")
    print("HANKEL22_DETERMINANT: -2^382*3^191*5^22*7^8*61^7")
    print(f"HANKEL22_SHA256: {hankel_hash}")
    print(f"VALUES_K0_TO5000_SHA256: {value_digest}")
    print("HARMONIC_LANES: each residue class mod3 has natural and harmonic coefficient 1/3")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: the explicit ternary degree-seven word has minimal recurrence order 22 because its common leading term lowers both nontrivial cubic Fourier colours to degree six; no LRC current or LRC(14) conclusion follows")


if __name__ == "__main__":
    main()
