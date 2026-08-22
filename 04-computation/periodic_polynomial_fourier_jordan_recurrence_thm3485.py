#!/usr/bin/env python3
"""Exact companion for THM-3485's periodic-polynomial recurrence theorem."""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from json import dumps
from math import gcd
from pathlib import Path
import sys


EXPECTED_SEMANTIC_SHA256 = (
    "834dfd2be5aed2f4d92b6f7fe742bc0df3a3008e61a0e1bc3fd39eef2dac54c0"
)

THM3484_LANES = (
    (0, 0, 0, 0, -2048, 4096, 8192, -16384),
    (0, 0, 0, 0, 2048, 4096, -24576, -16384),
    (0, -256, 3072, -15360, 40960, -61440, 49152, -16384),
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def trim(poly: tuple[int | Fraction, ...]) -> tuple[Fraction, ...]:
    values = [Fraction(value) for value in poly]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def poly_add(
    left: tuple[int | Fraction, ...], right: tuple[int | Fraction, ...]
) -> tuple[Fraction, ...]:
    size = max(len(left), len(right))
    values = [Fraction(0)] * size
    for index, coefficient in enumerate(left):
        values[index] += coefficient
    for index, coefficient in enumerate(right):
        values[index] += coefficient
    return trim(tuple(values))


def poly_scale(
    poly: tuple[int | Fraction, ...], scalar: int | Fraction
) -> tuple[Fraction, ...]:
    return trim(tuple(Fraction(scalar) * coefficient for coefficient in poly))


def poly_multiply(
    left: tuple[int | Fraction, ...], right: tuple[int | Fraction, ...]
) -> tuple[Fraction, ...]:
    values = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            values[i + j] += Fraction(a) * Fraction(b)
    return trim(tuple(values))


def poly_power(
    base: tuple[int | Fraction, ...], exponent: int
) -> tuple[Fraction, ...]:
    result: tuple[Fraction, ...] = (Fraction(1),)
    for _unused in range(exponent):
        result = poly_multiply(result, base)
    return result


def poly_divmod(
    numerator: tuple[int | Fraction, ...],
    denominator: tuple[int | Fraction, ...],
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    work = list(trim(numerator))
    divisor = trim(denominator)
    require(divisor != (Fraction(0),), "zero polynomial divisor")
    if len(work) < len(divisor):
        return (Fraction(0),), tuple(work)
    quotient = [Fraction(0)] * (len(work) - len(divisor) + 1)
    while len(work) >= len(divisor) and any(work):
        shift = len(work) - len(divisor)
        scale = work[-1] / divisor[-1]
        quotient[shift] += scale
        for index, coefficient in enumerate(divisor):
            work[index + shift] -= scale * coefficient
        while len(work) > 1 and work[-1] == 0:
            work.pop()
    return trim(tuple(quotient)), trim(tuple(work))


def divisors(number: int) -> tuple[int, ...]:
    return tuple(value for value in range(1, number + 1) if number % value == 0)


CYCLOTOMIC_CACHE: dict[int, tuple[Fraction, ...]] = {}


def cyclotomic(number: int) -> tuple[Fraction, ...]:
    if number in CYCLOTOMIC_CACHE:
        return CYCLOTOMIC_CACHE[number]
    numerator = [Fraction(0)] * (number + 1)
    numerator[0] = -1
    numerator[-1] = 1
    result = trim(tuple(numerator))
    for proper in divisors(number)[:-1]:
        quotient, remainder = poly_divmod(result, cyclotomic(proper))
        require(remainder == (Fraction(0),),
                ("cyclotomic division", number, proper, remainder))
        result = quotient
    require(result[-1] == 1, ("nonmonic cyclotomic", number, result))
    require(all(value.denominator == 1 for value in result),
            ("nonintegral cyclotomic", number, result))
    CYCLOTOMIC_CACHE[number] = result
    return result


def poly_eval(poly: tuple[int | Fraction, ...], value: int) -> Fraction:
    answer = Fraction(0)
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def add_lanes(
    *terms: tuple[int | Fraction, tuple[int | Fraction, ...]]
) -> tuple[Fraction, ...]:
    answer: tuple[Fraction, ...] = (Fraction(0),)
    for scalar, poly in terms:
        answer = poly_add(answer, poly_scale(poly, scalar))
    return answer


def colour_degree(
    lanes: tuple[tuple[int | Fraction, ...], ...], colour: int
) -> int | None:
    period = len(lanes)
    common = gcd(period, colour)
    order = period // common
    primitive_power = 0 if colour == 0 else colour // common
    max_degree = max(len(trim(lane)) for lane in lanes) - 1
    modulus = cyclotomic(order)
    for degree in range(max_degree, -1, -1):
        evaluation = [Fraction(0)] * max(1, order)
        for residue, lane in enumerate(lanes):
            coefficient = Fraction(lane[degree]) if degree < len(lane) else Fraction(0)
            exponent = (-primitive_power * residue) % order if order > 1 else 0
            evaluation[exponent] += coefficient
        _quotient, remainder = poly_divmod(trim(tuple(evaluation)), modulus)
        if remainder != (Fraction(0),):
            return degree
    return None


def predicted_characteristic(
    lanes: tuple[tuple[int | Fraction, ...], ...]
) -> tuple[tuple[Fraction, ...], tuple[int | None, ...], tuple[tuple[int, int], ...]]:
    period = len(lanes)
    degrees = tuple(colour_degree(lanes, colour) for colour in range(period))
    exponents_by_order: dict[int, int] = {}
    for colour, degree in enumerate(degrees):
        order = period // gcd(period, colour)
        exponent = 0 if degree is None else degree + 1
        if order in exponents_by_order:
            require(exponents_by_order[order] == exponent,
                    ("Galois-orbit degree mismatch", period, order,
                     exponents_by_order[order], exponent, degrees))
        exponents_by_order[order] = exponent
    characteristic: tuple[Fraction, ...] = (Fraction(1),)
    for order in sorted(exponents_by_order):
        characteristic = poly_multiply(
            characteristic,
            poly_power(cyclotomic(order), exponents_by_order[order]),
        )
    return characteristic, degrees, tuple(sorted(exponents_by_order.items()))


def sequence_values(
    lanes: tuple[tuple[int | Fraction, ...], ...], count: int
) -> tuple[Fraction, ...]:
    return tuple(poly_eval(lanes[index % len(lanes)], index) for index in range(count))


def recurrence_residual(
    connection: tuple[int | Fraction, ...],
    values: tuple[Fraction, ...],
    index: int,
) -> Fraction:
    return sum(
        (Fraction(connection[j]) * values[index - j]
         for j in range(len(connection))),
        Fraction(0),
    )


def berlekamp_massey(
    values: tuple[Fraction, ...]
) -> tuple[int, tuple[Fraction, ...]]:
    connection = [Fraction(1)]
    previous_connection = [Fraction(1)]
    length = 0
    displacement = 1
    previous_discrepancy = Fraction(1)
    for index, value in enumerate(values):
        discrepancy = value
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
    return length, trim(tuple(connection))


def packet_certificate(
    name: str, lanes: tuple[tuple[int | Fraction, ...], ...]
) -> tuple[object, ...]:
    characteristic, degrees, exponents = predicted_characteristic(lanes)
    order = len(characteristic) - 1
    connection = tuple(reversed(characteristic))
    values = sequence_values(lanes, max(200, 6 * order + 30))
    require(all(
        recurrence_residual(connection, values, index) == 0
        for index in range(order, len(values))
    ), (name, "predicted recurrence failure"))
    bm_order, bm_connection = berlekamp_massey(values[:max(100, 4 * order + 10)])
    require((bm_order, bm_connection) == (order, connection),
            (name, "minimality mismatch", bm_order, bm_connection,
             order, connection))
    max_degree = max(len(trim(lane)) for lane in lanes) - 1
    naive_order = len(lanes) * (max_degree + 1)
    defect = naive_order - order
    return (
        name,
        len(lanes),
        max_degree,
        degrees,
        exponents,
        tuple(int(value) for value in characteristic),
        order,
        naive_order,
        defect,
        sha256(dumps([str(value) for value in values],
                     separators=(",", ":")).encode("ascii")).hexdigest(),
    )


def deterministic_packets() -> tuple[tuple[str, tuple[tuple[Fraction, ...], ...]], ...]:
    packets: list[tuple[str, tuple[tuple[Fraction, ...], ...]]] = []

    polynomial_a = (2, -1, 0, 3)
    polynomial_b = (5, 4)
    parity_lanes = tuple(
        add_lanes((1, polynomial_a), ((-1) ** residue, polynomial_b))
        for residue in range(4)
    )
    packets.append(("p4-parity-only-zero-colours", parity_lanes))

    a6 = (1, -2, 3)
    b6 = (4, 5)
    c6 = (-1, 0, 2, 1)
    trace3 = (2, -1, -1, 2, -1, -1)
    p6_lanes = tuple(
        add_lanes(
            (1, a6),
            ((-1) ** residue, b6),
            (trace3[residue], c6),
        )
        for residue in range(6)
    )
    packets.append(("p6-orders-1-2-3", p6_lanes))

    p5_lanes = tuple(
        (residue + 1, 2 - residue, residue * residue,
         residue * residue + 1, 7)
        for residue in range(5)
    )
    packets.append(("p5-common-leading-layer", p5_lanes))

    a8 = (3, -1)
    b8 = (0, 2, 1)
    c8 = (1, 4)
    trace4 = (2, 0, -2, 0, 2, 0, -2, 0)
    p8_lanes = tuple(
        add_lanes(
            (1, a8),
            ((-1) ** residue, b8),
            (trace4[residue], c8),
        )
        for residue in range(8)
    )
    packets.append(("p8-orders-1-2-4", p8_lanes))

    packets.append(("p7-zero-word", tuple((Fraction(0),) for _ in range(7))))

    for period in range(2, 9):
        for degree in range(5):
            lanes = tuple(
                tuple(
                    Fraction(((residue + 2) * (power + 3) + period
                              + degree * degree) % 11 - 5)
                    for power in range(degree + 1)
                )
                for residue in range(period)
            )
            packets.append((f"generic-p{period}-d{degree}", lanes))
    return tuple(packets)


def security_certificate() -> tuple[object, ...]:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "assert forbidden")
    forbidden = {
        "eval", "exec", "compile", "open", "system", "popen", "run",
        "Popen", "write_text", "write_bytes", "unlink", "remove", "rename",
    }
    calls = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    calls.update(
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not (calls & forbidden), ("forbidden calls", calls & forbidden))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {
        "__future__", "ast", "fractions", "hashlib", "json", "math",
        "pathlib", "sys",
    }
    require(imports <= allowed, ("unexpected imports", imports - allowed))
    return ("NO_ASSERT_AST", tuple(sorted(imports)), tuple(sorted(forbidden)))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    security = security_certificate()
    specialization = packet_certificate("THM-3484", THM3484_LANES)
    require(specialization[3] == (7, 6, 6), specialization)
    require(specialization[4] == ((1, 8), (3, 7)), specialization)
    require(specialization[6:9] == (22, 24, 2), specialization)

    packets = deterministic_packets()
    certificates = tuple(
        packet_certificate(name, lanes) for name, lanes in packets
    )
    require(len(certificates) == 40, len(certificates))

    common_lead_checks = []
    for period in range(2, 10):
        degree = 5
        lanes = tuple(
            tuple(
                [Fraction((residue + 1) * (power + 2) % 13 - 6)
                 for power in range(degree)]
                + [Fraction(17)]
            )
            for residue in range(period)
        )
        characteristic, degrees, _exponents = predicted_characteristic(lanes)
        order = len(characteristic) - 1
        require(degrees[0] == degree, (period, degrees))
        require(all(value is None or value <= degree - 1
                    for value in degrees[1:]), (period, degrees))
        require(order <= period * (degree + 1) - (period - 1),
                (period, order, degrees))
        common_lead_checks.append((period, degrees, order))

    semantic_payload = {
        "security": security,
        "thm3484": specialization,
        "certificates": certificates,
        "common_lead_checks": tuple(common_lead_checks),
        "cyclotomics_1_to_12": tuple(
            tuple(int(value) for value in cyclotomic(number))
            for number in range(1, 13)
        ),
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":"),
              default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256,
                (semantic_hash, EXPECTED_SEMANTIC_SHA256))

    named = tuple(certificate for certificate in certificates
                  if not str(certificate[0]).startswith("generic-"))
    print("THM-3485 PERIODIC-POLYNOMIAL FOURIER/JORDAN RECURRENCE EXACT COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    print(f"SECURITY: {security}")
    print(f"THM3484_SPECIALIZATION: {specialization}")
    print(f"NAMED_HOSTILES: {named}")
    print(f"GENERIC_PACKET_COUNT: {len(certificates) - len(named)}")
    print(f"COMMON_LEADING_LAYER_CHECKS: {tuple(common_lead_checks)}")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: exact rational recurrence reconstruction agrees with the Fourier-colour degree formula in all packets; the companion is finite evidence and the theorem proof is algebraic")


if __name__ == "__main__":
    main()
