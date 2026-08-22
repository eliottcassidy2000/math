#!/usr/bin/env python3
"""Exact finite certificate for the support-minimal THM-3665 detector."""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARENTS = (
    ROOT / "04-computation/lrc14_relation_residue_pushforward_thm2334.py",
    ROOT / "05-knowledge/results/lrc14_relation_residue_pushforward_thm2334.out",
)
EXPECTED_PARENT_HASHES = (
    "ce220ede12175b6851810782c880f95048fe9f4643cc4f52f47a7f4d8dcb7b0c",
    "d2d9b49db9ef3eabf7e3ae17cea247da554a9f8df2abfc3907243d317b21fec1",
)
EXPECTED_SEMANTIC_SHA256 = "7130760feb9c252c6f7d1b2dc5fdd2b224f3469718aa4363380511de88c586d3"

P = 13
N = P * P
MOD = 755373809845391722745761
ZETA = 123453826432109539554819
MASK = {(0, 0): 1, (1, 0): 1, (0, 1): -2}
PROJECTIVE_REPRESENTATIVES = tuple((1, slope) for slope in range(P)) + ((0, 1),)
EXPECTED_LINE_NORMS = (
    13, 13, 35503, 20969, 18603, 6773, 26039,
    169, 26039, 6773, 18603, 20969, 35503, 53248,
)
EXPECTED_INTEGER_DETERMINANT = 9072750758202804116713220630638041462433831593234432
EXPECTED_SQUARE_FACTOR = 26417870930056774319865792
EXPECTED_MODULAR_DETERMINANT = 408889968250631621108841


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
    detector = {label: MASK.get(label, 0) % MOD for label in labels}
    require(sum(MASK.values()) == 0 and len(MASK) == 3
            and sum(abs(value) for value in MASK.values()) == 4,
            "three-site mask")

    phase = {}
    fourier = {}
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

    transfer_rows = []
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
        transfer_rows.append(tuple(row))
    transfer_digest = digest_json(tuple(transfer_rows))

    # Every two-site mean-zero mask is a scalar translate of delta_0-delta_d.
    # Its zero characters are the annihilator of d, which has 13 elements.
    two_site_profiles = []
    for difference in labels[1:]:
        annihilator = tuple(
            frequency for frequency in labels
            if (frequency[0] * difference[0]
                + frequency[1] * difference[1]) % P == 0
        )
        require(len(annihilator) == P and annihilator[0] == (0, 0),
                ("two-site annihilator", difference, annihilator))
        two_site_profiles.append((difference, annihilator))
    require(len(two_site_profiles) == N - 1
            and {len(profile[1]) - 1 for profile in two_site_profiles} == {12},
            "two-site obstruction census")
    two_site_digest = digest_json(tuple(two_site_profiles))

    line_polynomials = []
    line_norms = []
    for frequency in PROJECTIVE_REPRESENTATIVES:
        coefficients = [0] * P
        for label, value in MASK.items():
            exponent = -(
                frequency[0] * label[0] + frequency[1] * label[1]
            ) % P
            coefficients[exponent] += value
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
    require(integer_determinant == 13 * EXPECTED_SQUARE_FACTOR ** 2,
            ("13 times square", integer_determinant))
    require(integer_determinant % MOD == EXPECTED_MODULAR_DETERMINANT,
            ("modular determinant", integer_determinant % MOD))

    # Exact combinatorial classification used in the trigonometric extremum
    # proof.  Write u in [-6,6] and delta=pi*k/13 with |k|<=13.
    spectral_cases = []
    min_candidates = []
    max_candidates = []
    for u, v in labels:
        if (u, v) == (0, 0):
            continue
        centred_u = u if u <= 6 else u - P
        k = (2 * v - centred_u) % (2 * P)
        if k > P:
            k -= 2 * P
        require(-12 <= centred_u <= 12 and -12 <= k <= 13,
                ("centred spectral coordinates", u, v, centred_u, k))
        spectral_cases.append(((u, v), centred_u, k))
        if k == 0 and abs(centred_u) == 2:
            min_candidates.append((u, v))
        if centred_u == 0 and abs(k) == 12:
            max_candidates.append((u, v))
    require(tuple(min_candidates) == ((2, 1), (11, 12)), min_candidates)
    require(tuple(max_candidates) == ((0, 6), (0, 7)), max_candidates)
    spectral_case_digest = digest_json(tuple(spectral_cases))

    # Every ordered basis gives a distinct support-minimal mask.
    basis_masks = set()
    basis_count = 0
    for a0 in range(P):
        for a1 in range(P):
            for b0 in range(P):
                for b1 in range(P):
                    if (a0 * b1 - a1 * b0) % P == 0:
                        continue
                    basis_count += 1
                    basis_masks.add(tuple(sorted((
                        ((0, 0), 1),
                        ((a0, a1), 1),
                        ((b0, b1), -2),
                    ))))
    require(basis_count == 26208 and len(basis_masks) == basis_count,
            ("basis mask count", basis_count, len(basis_masks)))

    semantic = digest_json((
        parent_hashes, P, N, MOD, ZETA,
        tuple(sorted(MASK.items())),
        tuple((frequency, fourier[frequency]) for frequency in labels),
        zero_set, transfer_digest,
        two_site_digest,
        PROJECTIVE_REPRESENTATIVES,
        tuple(line_polynomials), line_norms,
        integer_determinant, EXPECTED_SQUARE_FACTOR,
        tuple(spectral_cases), tuple(min_candidates), tuple(max_candidates),
        basis_count,
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3665 LRC support-minimal three-twist target detector ==")
    print(f"parents_sha256_lf={parent_hashes}")
    print("mask=delta_(0,0)+delta_(1,0)-2delta_(0,1);support:3;l1:4;mean:0")
    print(f"fourier_zero_set={zero_set};nontrivial_support:{N - 1}")
    print(f"universal_transfer_pairs={N * N};sha256={transfer_digest}")
    print(f"two_site_differences={N - 1};nontrivial_zeros_each:12;sha256={two_site_digest}")
    print("general_minimum_support_on_Fp^d=d+1")
    print(f"projective_line_norms={line_norms};minimum:13")
    print(f"integer_augmentation_determinant={integer_determinant}=13*{EXPECTED_SQUARE_FACTOR}^2")
    print(f"integer_determinant_mod_p={integer_determinant % MOD}")
    print(f"spectral_case_sha256={spectral_case_digest}")
    print(f"exact_min_multiplier=4*sin(pi/13)^2;frequencies={tuple(min_candidates)}")
    print(f"exact_max_multiplier=4*cos(pi/26);frequencies={tuple(max_candidates)}")
    print("sharp_l2_frame_bounds=16*sin(pi/13)^4 and 16*cos(pi/26)^2")
    print(f"ordered_basis_detectors={basis_count};distinct:{len(basis_masks)}")
    print("criterion=nonzero target iff H(s)+H(s-a)-2H(s-b) is nonzero for some s")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=abstract target-twist group;not inherited physical exceptional support or LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
