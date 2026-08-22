#!/usr/bin/env python3
"""Exact finite-ring controls for THM-3679's p-adic scalar seam."""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


CHECKS = 0
FINITE_CASES = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def row_difference(left: int, right: int) -> list[int]:
    row = [0] * 9
    row[left] = 1
    row[right] = -1
    return row


UNITS = tuple(range(6))
BLOCKERS = (6, 7, 8)
CHARTS = ((0, 1), (2, 3), (4, 5), (1, 0))


def packet_matrix(source: int) -> Matrix:
    target_a, target_b = tuple(label for label in BLOCKERS if label != source)
    rows: list[list[int]] = []
    for graft_a, graft_b in CHARTS:
        rows.append(row_difference(target_a, graft_a))
        rows.append(row_difference(target_b, graft_b))
    return Matrix(rows)


def smith_profile(matrix: Matrix) -> tuple[int, tuple[int, ...], int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    entries = tuple(
        abs(int(diagonal[index, index]))
        for index in range(min(diagonal.rows, diagonal.cols))
        if diagonal[index, index] != 0
    )
    return len(entries), entries, matrix.cols - len(entries)


def valuation_capped(value: int, prime: int, depth: int) -> int:
    modulus = prime**depth
    value %= modulus
    if value == 0:
        return depth
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def finite_ring_controls() -> tuple[tuple[int, int], ...]:
    global FINITE_CASES
    prime = 3
    ledger: list[tuple[int, int]] = []
    for depth in range(1, 5):
        modulus = prime**depth
        units = tuple(value for value in range(modulus) if value % prime)
        cases = 0
        for source_speed in range(modulus):
            source_valuation = valuation_capped(source_speed, prime, depth)
            for other_sum in range(modulus):
                other_valuation = valuation_capped(other_sum, prime, depth)
                observed = any(
                    (other_sum + unit * source_speed) % modulus == 0
                    for unit in units
                )
                expected = source_valuation == other_valuation
                require(observed == expected,
                        ("one-source unit equation", depth, source_speed,
                         other_sum, observed, expected))
                cases += 1

        for total_speed in range(modulus):
            for scalar in units:
                observed = scalar * total_speed % modulus == 0
                require(observed == (total_speed == 0),
                        ("two-source scalar checksum", depth, total_speed,
                         scalar, observed))
                cases += 1

        ledger.append((depth, cases))
        FINITE_CASES += cases

    for depth in range(1, 4):
        lower = prime**depth
        upper = prime ** (depth + 1)
        for left in range(upper):
            for digit in range(prime):
                right = (left + digit * lower) % upper
                difference = (right - left) % upper
                require(difference % lower == 0,
                        ("Bockstein divisibility", depth, left, right))
                next_digit = (difference // lower) % prime
                require((next_digit == 0) == (difference == 0),
                        ("Bockstein zero upgrade", depth, left, right,
                         next_digit))
                FINITE_CASES += 1
    return tuple(ledger)


def canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def main() -> None:
    single_profiles = tuple(
        (source, smith_profile(packet_matrix(source))) for source in BLOCKERS
    )
    require(all(profile == (7, (1,) * 7, 2)
                for _source, profile in single_profiles),
            ("single-source Smith profiles", single_profiles))

    pair_profiles = []
    for first_index, first in enumerate(BLOCKERS):
        for second in BLOCKERS[first_index + 1:]:
            combined = packet_matrix(first).col_join(packet_matrix(second))
            pair_profiles.append(((first, second), smith_profile(combined)))
    pair_profiles = tuple(pair_profiles)
    require(all(profile == (8, (1,) * 8, 1)
                for _pair, profile in pair_profiles),
            ("two-source Smith profiles", pair_profiles))

    all_profile = smith_profile(
        packet_matrix(6).col_join(packet_matrix(7)).col_join(packet_matrix(8))
    )
    require(all_profile == (8, (1,) * 8, 1),
            ("three-source Smith profile", all_profile))

    finite_ledger = finite_ring_controls()

    source = Path(__file__).read_bytes()
    require(b"\r\n" not in source, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.decode("utf-8")))),
            "Python assert node present")

    semantic = (single_profiles, pair_profiles, all_profile,
                finite_ledger, FINITE_CASES)
    print("== THM-3679 p-adic scalar-seam exact controls ==")
    print(f"single_source_smith={single_profiles}")
    print(f"two_source_smith={pair_profiles}")
    print(f"three_source_smith={all_profile}")
    print(f"finite_ring_ledger_p3={finite_ledger};cases={FINITE_CASES}")
    print(f"semantic_sha256={canonical_digest(semantic)}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS;scope=finite-ring controls, not an LRC current lift")


if __name__ == "__main__":
    main()
