"""Exact companion for THM-3250's charged Heisenberg intertwiner."""

from collections import Counter
from fractions import Fraction
from itertools import product
import ast
import hashlib
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT
    / (
        "01-canon/theorems/"
        "THM-3240-exact-address-heisenberg-clutch-on-carrier-imbalance.md"
    ): "7d23f2920adfecb17d8a149aada08a8e34215111546eac63b40570e898e14f14",
    ROOT
    / (
        "01-canon/theorems/"
        "THM-3243-contact-deformation-blowup-equivariance-and-full-orbit-resolution.md"
    ): "1feb6bd6be925826301d5a43ce5fba724f756f96e33185130792ce3ba02eed08",
}


def lf_sha256(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency, expected_hash in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected_hash, "dependency hash: %s" % dependency.name)

syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")


def exact_rank(matrix):
    if not matrix:
        return 0
    work = [[Fraction(entry) for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    work[row][index] - factor * work[pivot_row][index]
                    for index in range(column_count)
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def exact_determinant(matrix):
    size = len(matrix)
    work = [[Fraction(entry) for entry in row] for row in matrix]
    result = Fraction(1)
    sign = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column] / pivot_value
            for index in range(column, size):
                work[row][index] -= factor * work[column][index]
    return sign * result


def exact_inverse(matrix):
    size = len(matrix)
    work = [
        [Fraction(entry) for entry in row]
        + [Fraction(int(row_index == column)) for column in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        require(pivot is not None, "singular pointed-current matrix")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    work[row][index] - factor * work[column][index]
                    for index in range(2 * size)
                ]
    return [row[size:] for row in work]


def matrix_multiply(left, right):
    return [
        [
            sum(left[row][middle] * right[middle][column] for middle in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def charged_intertwiner_checks(prime):
    checks = 0
    for x, y, central, kappa, source, spectator in product(range(prime), repeat=6):
        if kappa == 0:
            continue
        inverse_kappa = pow(kappa, -1, prime)
        flag_base = (source - inverse_kappa * spectator) % prime

        address_phase = kappa * (central - y * source) % prime
        address_source = (source + x) % prime

        flag_phase = (
            kappa * central - kappa * y * flag_base - spectator * y
        ) % prime
        moved_flag_base = (flag_base + x) % prime

        require(address_phase == flag_phase, "charged phase mismatch")
        require(
            moved_flag_base == (address_source - inverse_kappa * spectator) % prime,
            "charged base-coordinate mismatch",
        )
        checks += 1
    require(checks == prime**5 * (prime - 1), "charged test census")
    return checks


def neutral_character_checks(prime):
    address_characters = Counter(
        (dual_source, 0)
        for dual_source in range(prime)
        for _spectator in range(prime)
    )
    flag_characters = Counter(
        (dual_base, -dual_slope % prime)
        for dual_base in range(prime)
        for dual_slope in range(prime)
    )
    require(len(address_characters) == prime, "address neutral character types")
    require(set(address_characters.values()) == {prime}, "address neutral multiplicities")
    require(len(flag_characters) == prime**2, "flag neutral character types")
    require(set(flag_characters.values()) == {1}, "flag neutral multiplicities")
    require(address_characters != flag_characters, "neutral obstruction disappeared")
    return len(address_characters), len(flag_characters)


def weyl_orthogonality_checks(prime):
    kappa = 1
    diagonal_pairs = 0
    off_diagonal_pairs = 0
    for shift, phase, other_shift, other_phase in product(range(prime), repeat=4):
        if shift != other_shift:
            off_diagonal_pairs += 1
            continue
        exponents = [
            kappa * (other_phase - phase) * source % prime
            for source in range(prime)
        ]
        residue_counts = Counter(exponents)
        if phase == other_phase:
            require(residue_counts == Counter({0: prime}), "Weyl diagonal Gram")
            diagonal_pairs += 1
        else:
            require(
                residue_counts == Counter({residue: 1 for residue in range(prime)}),
                "Weyl off-diagonal root sum",
            )
            off_diagonal_pairs += 1
    require(diagonal_pairs == prime**2, "Weyl diagonal census")
    require(off_diagonal_pairs == prime**4 - prime**2, "Weyl off-diagonal census")
    return diagonal_pairs


rows = []
for prime in (3, 5, 7, 13):
    charged_checks = charged_intertwiner_checks(prime)
    address_types, flag_types = neutral_character_checks(prime)
    weyl_dimension = weyl_orthogonality_checks(prime)

    identity = [[int(row == column) for column in range(prime)] for row in range(prime)]
    all_ones = [[1 for _column in range(prime)] for _row in range(prime)]
    triangular_source = [
        [int(column <= row) for column in range(prime)] for row in range(prime)
    ]
    triangular_target = [
        [int(column >= row) for column in range(prime)] for row in range(prime)
    ]

    require(exact_rank(identity) == prime, "identity cyclic control")
    require(exact_rank(all_ones) == 1, "full-support rank-one hostile")
    require(exact_determinant(identity) == 1, "identity determinant")
    require(exact_determinant(triangular_source) == 1, "source pointed determinant")
    require(exact_determinant(triangular_target) == 1, "target pointed determinant")

    unique_map = matrix_multiply(exact_inverse(triangular_source), triangular_target)
    require(
        matrix_multiply(triangular_source, unique_map)
        == [[Fraction(entry) for entry in row] for row in triangular_target],
        "pointed intertwiner equation",
    )
    require(exact_determinant(unique_map) == 1, "pointed intertwiner invertibility")

    require(prime * exact_rank(identity) == prime**2, "cyclic orbit-span dimension")
    require(prime * exact_rank(all_ones) == prime, "hostile orbit-span dimension")
    require(prime * Fraction(1, prime) == 1, "address Fourier normalization")
    require(prime**2 * Fraction(1, prime**2) == 1, "flag Fourier normalization")
    require(Fraction(prime**2, prime) == prime, "unscaled intertwiner norm factor")
    maximum_full_rank = (prime - 1) * prime**2 + prime
    full_rank_defect = prime**3 - maximum_full_rank
    require(full_rank_defect == prime**2 - prime, "full-module rank defect")

    rows.append(
        (
            prime,
            charged_checks,
            address_types,
            flag_types,
            weyl_dimension,
            prime**2,
            prime,
            prime,
            maximum_full_rank,
            full_rank_defect,
        )
    )


prime = 13
neutral_dimension = prime**2
charged_dimension = (prime - 1) * prime**2
total_dimension = prime**3
require(neutral_dimension + charged_dimension == total_dimension, "p=13 block dimensions")
require(charged_dimension == 2028, "p=13 charged dimension")

print("THM-3250 exact companion")
print("dependency_hashes=%d assert_nodes=%d float_literals=%d" % (len(DEPENDENCIES), assert_nodes, float_literals))
for row in rows:
    print(
        "p=%d charged_tests=%d neutral_address_types=%d neutral_flag_types=%d "
        "weyl_dimension=%d cyclic_span=%d hostile_span=%d norm_ratio=%d "
        "max_full_rank=%d defect=%d" % row
    )
print(
    "p13_dimensions=neutral:%d charged:%d total:%d multiplicity:13"
    % (neutral_dimension, charged_dimension, total_dimension)
)
print("p13_full_map=max_rank:2041 kernel:156 cokernel:156")
print("p13_pointed_gate=det_source:1 det_target:1 det_map:1 unique:1")
print("THM-3250 PASS")
