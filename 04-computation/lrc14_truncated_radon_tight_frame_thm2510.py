#!/usr/bin/env python3
"""Exact referee for THM-2510's truncated-Radon tight-frame identity.

The proof is symbolic.  This dependency-free companion checks the complete
13 x 7 operator Gram matrix, its rank and vertical kernel, the sharp integral
dipole boundary, the full 42-cut energy multiplier, and the THM-2506 control.
It uses explicit RuntimeError gates, so normal and optimized Python execute
the same checks.
"""

from itertools import combinations, product


P = 13
Q = 7
FREE_ROWS = Q - 1
SLOPES = tuple(range(1, P))
CUTS = tuple(product(range(1, Q), range(Q)))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def zero_array() -> list[list[int]]:
    return [[0] * Q for _ in range(P)]


def from_free_coordinates(coordinates: list[int]) -> list[list[int]]:
    require(len(coordinates) == P * FREE_ROWS, "coordinate length drifted")
    defect = zero_array()
    index = 0
    for h in range(P):
        row_sum = 0
        for r in range(FREE_ROWS):
            value = coordinates[index]
            index += 1
            defect[h][r] = value
            row_sum += value
        defect[h][Q - 1] = -row_sum
    require(all(sum(row) == 0 for row in defect), "row-zero chart failed")
    return defect


def basis_array(h: int, r: int) -> list[list[int]]:
    defect = zero_array()
    defect[h][r] = 1
    defect[h][Q - 1] = -1
    return defect


def radon(
    defect: list[list[int]], tau: int, a: int = 1, c: int = 0
) -> list[int]:
    return [
        sum(
            defect[(v - tau * ((a * r + c) % Q)) % P][r]
            for r in range(Q)
        )
        for v in range(P)
    ]


def image(defect: list[list[int]]) -> tuple[int, ...]:
    return tuple(value for tau in SLOPES for value in radon(defect, tau))


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right))


def rhs_bilinear(
    left: list[list[int]], right: list[list[int]]
) -> int:
    ambient = P * sum(
        left[h][r] * right[h][r] for h in range(P) for r in range(Q)
    )
    left_sums = [sum(left[h][r] for h in range(P)) for r in range(Q)]
    right_sums = [sum(right[h][r] for h in range(P)) for r in range(Q)]
    return ambient - sum(a * b for a, b in zip(left_sums, right_sums))


def rhs_energy(defect: list[list[int]]) -> int:
    return rhs_bilinear(defect, defect)


def cut_bank_energy(defect: list[list[int]]) -> int:
    return sum(
        value * value
        for a, c in CUTS
        for tau in SLOPES
        for value in radon(defect, tau, a, c)
    )


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(value * inverse) % prime for value in work[rank]]
        for row in range(row_count):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                (a - factor * b) % prime
                for a, b in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def lcg_arrays(count: int) -> list[list[list[int]]]:
    state = 0x2510
    arrays = []
    for _ in range(count):
        coordinates = []
        for _ in range(P * FREE_ROWS):
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            coordinates.append(state % 7 - 3)
        arrays.append(from_free_coordinates(coordinates))
    return arrays


def main() -> None:
    bases = [basis_array(h, r) for h in range(P) for r in range(FREE_ROWS)]
    images = [image(defect) for defect in bases]
    require(len(bases) == 78, "row-zero basis dimension drifted")
    require(all(len(vector) == 156 for vector in images), "image size drifted")

    gram_checks = 0
    for first in range(len(bases)):
        for second in range(first, len(bases)):
            require(
                dot(images[first], images[second])
                == rhs_bilinear(bases[first], bases[second]),
                "tight-frame Gram identity failed",
            )
            gram_checks += 1
    require(gram_checks == 3_081, "Gram universe drifted")

    operator_rows = [
        [images[column][row] for column in range(len(images))]
        for row in range(len(images[0]))
    ]
    operator_rank = rank_mod(operator_rows, 1_000_003)
    require(operator_rank == 72, "Radon operator rank drifted")
    for r in range(FREE_ROWS):
        vertical = zero_array()
        for h in range(P):
            vertical[h][r] = 1
            vertical[h][Q - 1] = -1
        require(image(vertical) == (0,) * 156, "vertical kernel vector survived")

    single_row_controls = 0
    minimum = None
    equality_count = 0
    for free in product((-1, 0, 1), repeat=FREE_ROWS):
        if not any(free):
            continue
        coordinates = [0] * (P * FREE_ROWS)
        coordinates[:FREE_ROWS] = free
        defect = from_free_coordinates(coordinates)
        energy = dot(image(defect), image(defect))
        require(energy == rhs_energy(defect), "single-row identity failed")
        minimum = energy if minimum is None else min(minimum, energy)
        equality_count += energy == 24
        single_row_controls += 1
    require(
        (single_row_controls, minimum, equality_count) == (728, 24, 42),
        "sharp single-row boundary drifted",
    )

    baseline = tuple(range(-3, 4))
    require(sum(baseline) == 0, "vertical baseline is not row-zero")
    equality_controls = 0
    for h, (r0, r1), epsilon in product(
        range(P), combinations(range(Q), 2), (-1, 1)
    ):
        defect = [list(baseline) for _ in range(P)]
        defect[h][r0] += epsilon
        defect[h][r1] -= epsilon
        require(rhs_energy(defect) == 24, "dipole equality form lost sharpness")
        require(dot(image(defect), image(defect)) == 24, "dipole image energy failed")
        equality_controls += 1
    require(equality_controls == 546, "dipole equality universe drifted")

    random_arrays = lcg_arrays(200)
    for defect in random_arrays:
        require(
            dot(image(defect), image(defect)) == rhs_energy(defect),
            "deterministic row-zero control failed",
        )
    for defect in random_arrays[:20]:
        require(
            cut_bank_energy(defect) == len(CUTS) * rhs_energy(defect),
            "42-cut energy multiplier failed",
        )

    control = zero_array()
    control[0][5] = 1
    control[0][3] = -1
    control[1][5] = 1
    control[1][4] = -1
    require(rhs_energy(control) == 46, "THM-2506 control energy drifted")
    per_cut = {
        sum(value * value for tau in SLOPES for value in radon(control, tau, a, c))
        for a, c in CUTS
    }
    require(per_cut == {46}, "THM-2506 cut energies are not uniform")
    require(cut_bank_energy(control) == 1_932, "stored quadratic energy mismatch")

    print("THM-2510 truncated-Radon tight-frame exact companion: PASS")
    print(f"row_zero_dimension={len(bases)}; operator_rank={operator_rank}; kernel_dimension=6")
    print(f"exact_symmetric_Gram_checks={gram_checks}; image_coordinates=156")
    print(
        "single_row_integral_controls="
        f"{single_row_controls}; minimum={minimum}; equality_count={equality_count}"
    )
    print(f"sharp_dipole_plus_vertical_baseline_controls={equality_controls}")
    print("deterministic_row_zero_controls=200; full_42_cut_controls=20")
    print("THM-2506_two_row_energy=46; full_42_cut_energy=1932")
    print("counting_norm_identity=sum_tau||R_tau d||^2=13||d-E_h d||^2")
    print("integral_nonvertical_floor=24; equality=one horizontal dipole plus vertical kernel")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
