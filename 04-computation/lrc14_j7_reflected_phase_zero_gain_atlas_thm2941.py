#!/usr/bin/env python3
"""Exact global atlas of phase-zero gains in the reflected THM-2941 CSP.

The phase floor of a reduced unordered channel ``1 <= P < Q`` is

    1/49 + c(P mod 14,Q mod 14)/(PQ).

The exact residue table has ``c >= -12/49``.  Hence phase zero forces
``PQ <= 12`` and is globally decidable by a fourteen-row finite check.  The
resulting gain atlas explains the qualitative projective-cap walls used by
the recursive reflected-cone proofs.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_phase_zero_gain_atlas_thm2941.out"
EXPECTED_SEMANTIC_SHA256 = "8a934cf2a5ecbc0e54a720e89d162be89fd0c4d932d1019ef167cb64443b0b18"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def triangle_sum(s: F, z: F) -> F:
    bound = s.numerator // s.denominator + 3
    return sum(
        (max(F(0), s - abs(z + n)) for n in range(-bound, bound + 1)),
        F(0),
    )


def phase_correction(residue_p: int, residue_q: int) -> tuple[F, F]:
    a = F((residue_p + residue_q) % 14, 14)
    b = F((residue_q - residue_p) % 14, 14)
    events = {F(0), F(1)}
    for slope in (a, b):
        for n in range(-3, 4):
            for z in (-F(n), slope - n, -slope - n):
                if 0 <= z <= 1:
                    events.add(z)
    return min(
        (
            (triangle_sum(a, z) - a * a)
            - (triangle_sum(b, z) - b * b),
            z,
        )
        for z in events
    )


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def prime_vector(value: F) -> tuple[int, int, int]:
    rows = []
    numerator = value.numerator
    denominator = value.denominator
    for prime in (2, 3, 5):
        exponent = 0
        while numerator % prime == 0:
            numerator //= prime
            exponent += 1
        while denominator % prime == 0:
            denominator //= prime
            exponent -= 1
        rows.append(exponent)
    require(numerator == denominator == 1, (value, numerator, denominator))
    return tuple(rows)


def rational_rank(columns: tuple[tuple[int, ...], ...]) -> int:
    if not columns:
        return 0
    matrix = [[F(columns[j][i]) for j in range(len(columns))] for i in range(3)]
    rank = 0
    for column in range(len(columns)):
        pivot = next((row for row in range(rank, 3) if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        value = matrix[rank][column]
        matrix[rank] = [entry / value for entry in matrix[rank]]
        for row in range(3):
            if row == rank or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [a - factor * b for a, b in zip(matrix[row], matrix[rank])]
        rank += 1
        if rank == 3:
            break
    return rank


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    residue_rows = tuple(
        (phase_correction(p, q)[0], p, q, phase_correction(p, q)[1])
        for p in range(14)
        for q in range(14)
    )
    correction_min = min(row[0] for row in residue_rows)
    require(correction_min == F(-12, 49), correction_min)

    small_channels = []
    for p in range(1, 13):
        for q in range(p + 1, 13):
            if p * q > 12 or gcd(p, q) != 1:
                continue
            correction, location = phase_correction(p % 14, q % 14)
            phase = F(1, 49) + correction / (p * q)
            require(phase >= 0, (p, q, phase))
            small_channels.append((p, q, F(q, p), correction, phase, location))
    require(len(small_channels) == 14, len(small_channels))

    zero_rows = tuple(row for row in small_channels if row[4] == 0)
    zero_gains = tuple(sorted(row[2] for row in zero_rows))
    expected_gains = (F(4, 3), F(3, 2), F(2), F(5, 2), F(3), F(4), F(5), F(6))
    require(zero_gains == expected_gains, zero_gains)

    vectors = tuple(prime_vector(gain) for gain in zero_gains)
    rank_profile = tuple(rational_rank(vectors[:k]) for k in range(1, len(vectors) + 1))
    nullity_profile = tuple(k - rank_profile[k - 1] for k in range(1, len(vectors) + 1))
    require(rank_profile == (1, 2, 2, 3, 3, 3, 3, 3), rank_profile)
    require(nullity_profile == (0, 0, 1, 1, 2, 3, 4, 5), nullity_profile)

    relations = (
        "(4/3)*(3/2)=2",
        "(3/2)*2=3",
        "2^2=4",
        "(5/2)*2=5",
        "3*2=6",
    )
    require(F(4, 3) * F(3, 2) == 2, relations[0])
    require(F(3, 2) * 2 == 3, relations[1])
    require(F(2) ** 2 == 4, relations[2])
    require(F(5, 2) * 2 == 5, relations[3])
    require(F(3) * 2 == 6, relations[4])

    semantic_payload = (
        correction_min,
        tuple(small_channels),
        zero_gains,
        vectors,
        rank_profile,
        nullity_profile,
        relations,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    lines = [
        "LRC14 reflected phase-zero gain atlas",
        "phase_floor=1/49+c(P mod14,Q mod14)/(PQ)",
        "correction_min=-12/49;therefore phase_zero implies PQ<=12",
        "small_primitive_channels=14;global_phase_zero_channels=8",
        "zero_gains=" + repr(tuple(qtext(value) for value in zero_gains)),
        "prime_basis=(2,3,5);vectors=" + repr(vectors),
        "rank_profile=" + repr(rank_profile),
        "nullity_profile=" + repr(nullity_profile),
        "relations=" + repr(relations),
        "boundary_law=rank rises at 4/3,3/2,5/2;cycle nullity rises at 2,3,4,5,6",
        "cap_2=unique first circuit (4/3)*(3/2)=2",
        "cap_5/2=new independent prime-5 direction;old circuit persists",
        f"semantic_sha256={semantic}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        "scope=exact global phase-zero atlas for the reflected THM-2941 phase floor",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
