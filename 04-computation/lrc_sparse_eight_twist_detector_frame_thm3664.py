#!/usr/bin/env python3
"""Exact split-group certificate for the THM-3664 eight-twist frame."""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENTS = (
    ROOT / "04-computation/lrc14_relation_residue_pushforward_thm2334.py",
    ROOT / "05-knowledge/results/lrc14_relation_residue_pushforward_thm2334.out",
    ROOT / "04-computation/lrc_exceptional_detector_simple_spectrum_thm3661.py",
    ROOT / "05-knowledge/results/lrc_exceptional_detector_simple_spectrum_thm3661.out",
)
EXPECTED_PARENT_HASHES = (
    "ce220ede12175b6851810782c880f95048fe9f4643cc4f52f47a7f4d8dcb7b0c",
    "d2d9b49db9ef3eabf7e3ae17cea247da554a9f8df2abfc3907243d317b21fec1",
    "7c8ebacaff5318b58ed1588d2e38edf1bc3c14ea81beec069d781450db986251",
    "61893d6e628232e877c73d6102f042654327fddfe7d51d1cf79c1aa6bc65419f",
)
EXPECTED_SEMANTIC_SHA256 = "24c09ef5a3681ccbf9878276e1e7fda2efabc947b0e86c06a6bf594e088561a6"

P = 13
N = P * P
MOD = 755373809845391722745761
ZETA = 123453826432109539554819
X_PLUS = frozenset(((12, 0), (0, 11), (6, 5), (9, 3)))
X_MINUS = frozenset(((0, 12), (12, 1), (6, 7), (3, 9)))
PROJECTIVE_REPRESENTATIVES = tuple((1, slope) for slope in range(P)) + ((0, 1),)
EXPECTED_LINE_NORMS = (
    13, 53248, 13, 223093, 13, 9477, 137917,
    320437, 13, 36517, 9477, 2197, 13, 13,
)
EXPECTED_INTEGER_DETERMINANT = 18258988969361052598805681298174465892061184
EXPECTED_INTEGER_SQRT = 4273053822427357700928
EXPECTED_MODULAR_DETERMINANT = 669422013050837354847410


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def subtract(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return ((left[0] - right[0]) % P, (left[1] - right[1]) % P)


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant."""
    work = [row[:] for row in matrix]
    size = len(work)
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        require(pivot_row is not None, ("singular norm matrix", column))
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = (
                    work[row][target] * pivot
                    - work[row][column] * work[column][target]
                )
                require(numerator % previous == 0,
                        ("Bareiss nonexact division", column, row, target))
                work[row][target] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def reduce_cyclotomic(coefficients: list[int]) -> tuple[int, ...]:
    """Reduce modulo Phi_13=1+x+...+x^12 into the basis 1,...,x^11."""
    values = coefficients[:] + [0] * max(0, 12 - len(coefficients))
    for degree in range(len(values) - 1, 11, -1):
        coefficient = values[degree]
        if coefficient:
            values[degree] = 0
            shift = degree - 12
            for offset in range(12):
                values[shift + offset] -= coefficient
    return tuple(values[:12])


def cyclotomic_norm(coefficients: tuple[int, ...]) -> int:
    """Norm of P(zeta_13), as a multiplication determinant."""
    columns = []
    for shift in range(12):
        product = [0] * (len(coefficients) + shift)
        for degree, coefficient in enumerate(coefficients):
            product[degree + shift] = coefficient
        columns.append(reduce_cyclotomic(product))
    matrix = [[columns[column][row] for column in range(12)]
              for row in range(12)]
    return bareiss_determinant(matrix)


def main() -> None:
    parent_hashes = tuple(lf_sha256(path) for path in PARENTS)
    require(parent_hashes == EXPECTED_PARENT_HASHES,
            ("parent hashes", parent_hashes))
    require(pow(ZETA, P, MOD) == 1 and ZETA != 1, "zeta order")

    labels = tuple((r0, r1) for r0 in range(P) for r1 in range(P))
    signed_detector = {
        label: int(label in X_PLUS) - int(label in X_MINUS)
        for label in labels
    }
    detector = {label: value % MOD for label, value in signed_detector.items()}
    require(sum(abs(value) for value in signed_detector.values()) == 8,
            "eight-point support")
    require(sum(signed_detector.values()) == 0, "detector mean")

    fourier = {}
    phase = {}
    for frequency in labels:
        for label in labels:
            phase[(frequency, label)] = pow(
                ZETA,
                -(frequency[0] * label[0] + frequency[1] * label[1]) % P,
                MOD,
            )
        fourier[frequency] = sum(
            detector[label] * phase[(frequency, label)] for label in labels
        ) % MOD
    zero_set = tuple(frequency for frequency in labels if fourier[frequency] == 0)
    require(zero_set == ((0, 0),), ("Fourier zero set", zero_set))

    # For every profile basis vector delta_y, Fourier(C_g delta_y) equals
    # g_hat times Fourier(delta_y).  Linearity makes this the universal
    # sparse-detector transfer identity for all 169-coordinate profiles.
    transfer_digest_rows = []
    for frequency in labels:
        row = []
        for source in labels:
            direct = sum(
                detector[subtract(target, source)] * phase[(frequency, target)]
                for target in labels
            ) % MOD
            expected = fourier[frequency] * phase[(frequency, source)] % MOD
            require(direct == expected,
                    ("universal transfer", frequency, source, direct, expected))
            row.append(direct)
        transfer_digest_rows.append(tuple(row))
    transfer_digest = digest_json(tuple(transfer_digest_rows))

    # Exact characteristic-zero norms, one Galois orbit per projective line.
    line_polynomials = []
    line_norms = []
    for frequency in PROJECTIVE_REPRESENTATIVES:
        coefficients = [0] * P
        for label, value in signed_detector.items():
            exponent = -(
                frequency[0] * label[0] + frequency[1] * label[1]
            ) % P
            coefficients[exponent] += value
        require(sum(coefficients) == 0, ("line polynomial mean", frequency))
        line_polynomials.append(tuple(coefficients))
        line_norms.append(cyclotomic_norm(tuple(coefficients)))
    line_norms = tuple(line_norms)
    require(line_norms == EXPECTED_LINE_NORMS,
            ("projective line norms", line_norms))
    require(min(line_norms) == 13 and all(value > 0 for value in line_norms),
            ("line norm positivity/minimum", line_norms))

    integer_determinant = 1
    for value in line_norms:
        integer_determinant *= value
    require(integer_determinant == EXPECTED_INTEGER_DETERMINANT,
            ("integer determinant", integer_determinant))
    integer_sqrt = isqrt(integer_determinant)
    require(integer_sqrt == EXPECTED_INTEGER_SQRT
            and integer_sqrt * integer_sqrt == integer_determinant,
            ("square determinant", integer_sqrt))
    require(integer_determinant % MOD == EXPECTED_MODULAR_DETERMINANT,
            ("modular determinant lift", integer_determinant % MOD))

    # The basis gauge is genuine: every GL(2,13) transform produces a
    # different signed mask.  Thus the sparse criterion is basis-covariant,
    # not a canonical identification with the two-current address chart.
    transformed_masks = set()
    linear_count = 0
    for a in range(P):
        for b in range(P):
            for c in range(P):
                for d in range(P):
                    if (a * d - b * c) % P == 0:
                        continue
                    linear_count += 1
                    plus = tuple(sorted(
                        ((a * x + b * y) % P, (c * x + d * y) % P)
                        for x, y in X_PLUS
                    ))
                    minus = tuple(sorted(
                        ((a * x + b * y) % P, (c * x + d * y) % P)
                        for x, y in X_MINUS
                    ))
                    transformed_masks.add((plus, minus))
    require(linear_count == 26208 and len(transformed_masks) == linear_count,
            ("GL2 orbit/stabilizer", linear_count, len(transformed_masks)))
    mask_orbit_digest = digest_json(tuple(sorted(transformed_masks)))

    semantic = digest_json((
        parent_hashes,
        P, N, MOD, ZETA,
        tuple(sorted(X_PLUS)), tuple(sorted(X_MINUS)),
        tuple((frequency, fourier[frequency]) for frequency in labels),
        zero_set, transfer_digest,
        PROJECTIVE_REPRESENTATIVES,
        tuple(line_polynomials), line_norms,
        integer_determinant, integer_sqrt,
        linear_count, mask_orbit_digest,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3664 LRC sparse eight-twist detector frame ==")
    print(f"parents_sha256_lf={parent_hashes}")
    print(f"twist_group=F13^2;size:{N};detector_support:8;mean:0")
    print(f"fourier_zero_set={zero_set};nontrivial_support:{N - 1}")
    print(f"universal_transfer_pairs={N * N};sha256={transfer_digest}")
    print(f"projective_representatives={PROJECTIVE_REPRESENTATIVES}")
    print(f"projective_line_norms={line_norms};minimum:13")
    print(f"integer_augmentation_determinant={integer_determinant};square:{integer_sqrt}^2")
    print(f"integer_determinant_mod_p={integer_determinant % MOD}")
    print("pointwise_multiplier_bounds=13/8^11<=abs(g_hat)<=8")
    print("l2_frame_bounds=(169/8^22)||H-Hbar||^2<=||g*H||^2<=64||H-Hbar||^2")
    print(f"GL2_basis_masks={linear_count};distinct:{len(transformed_masks)};sha256={mask_orbit_digest}")
    print("criterion=nonzero target iff some translated signed eight-twist imbalance is nonzero")
    print("reconstruction=H-Hbar=h*(g*H), h_hat(0)=0, h_hat(q)=1/g_hat(q)")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=typed target-twist group after a chosen basis;not a physical address identification or LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
