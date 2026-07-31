#!/usr/bin/env python3
"""Exact companion for THM-2810 factorial-Hankel faithfulness.

All checks are integer or finite-field exact.  The universal statements in
the theorem are proved algebraically; the finite controls here pin signs,
indices, characteristic-p boundaries, and the explicit cyclic hostile.
"""

from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matmul(left, right, modulus=None):
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    require(all(len(row) == middle for row in left), "left matrix shape")
    require(all(len(row) == columns for row in right), "right matrix shape")
    result = []
    for i in range(rows):
        row = []
        for j in range(columns):
            value = sum(left[i][k] * right[k][j] for k in range(middle))
            row.append(value if modulus is None else value % modulus)
        result.append(tuple(row))
    return tuple(result)


def transpose(matrix):
    return tuple(zip(*matrix))


def determinant_bareiss(matrix):
    work = [list(row) for row in matrix]
    size = len(work)
    require(all(len(row) == size for row in work), "determinant shape")
    if size == 0:
        return 1
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row
                    for row in range(pivot_index + 1, size)
                    if work[row][pivot_index] != 0
                ),
                None,
            )
            require(swap is not None, "singular Bareiss pivot")
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for i in range(pivot_index + 1, size):
            for j in range(pivot_index + 1, size):
                numerator = (
                    work[i][j] * pivot
                    - work[i][pivot_index] * work[pivot_index][j]
                )
                require(
                    numerator % previous == 0,
                    "Bareiss division stopped being exact",
                )
                work[i][j] = numerator // previous
        previous = pivot
        for i in range(pivot_index + 1, size):
            work[i][pivot_index] = 0
    return sign * work[-1][-1]


def rank_mod_prime(matrix, prime):
    work = [[value % prime for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    rank = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(rank, rows)
                if work[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [value * inverse % prime for value in work[rank]]
        for row in range(rows):
            if row == rank:
                continue
            multiple = work[row][column]
            if multiple:
                work[row] = [
                    (work[row][j] - multiple * work[rank][j]) % prime
                    for j in range(columns)
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def hankel(size, modulus=None):
    rows = []
    for i in range(size):
        row = []
        for j in range(size):
            value = factorial(i + j)
            row.append(value if modulus is None else value % modulus)
        rows.append(tuple(row))
    return tuple(rows)


def factorization(size):
    diagonal = tuple(
        tuple(factorial(i) if i == j else 0 for j in range(size))
        for i in range(size)
    )
    pascal = tuple(
        tuple(comb(i, j) if j <= i else 0 for j in range(size))
        for i in range(size)
    )
    return matmul(matmul(matmul(diagonal, pascal), transpose(pascal)), diagonal)


def char_zero_controls():
    determinants = []
    for degree in range(11):
        size = degree + 1
        matrix = hankel(size)
        require(
            matrix == factorization(size),
            f"D P P^T D factorization failed at degree {degree}",
        )
        determinant = determinant_bareiss(matrix)
        expected = 1
        for index in range(size):
            expected *= factorial(index) ** 2
        require(
            determinant == expected,
            f"factorial Hankel determinant failed at degree {degree}",
        )
        determinants.append(determinant)
    return tuple(determinants)


def characteristic_p_controls():
    rows = []
    for prime in (2, 3, 5, 7, 11, 13):
        full = hankel(prime, prime)
        next_block = hankel(prime + 1, prime)
        determinant = determinant_bareiss(full) % prime
        expected = 1
        for index in range(prime):
            expected = expected * factorial(index) ** 2 % prime
        require(
            determinant == expected != 0
            and rank_mod_prime(full, prime) == prime
            and rank_mod_prime(next_block, prime) == prime,
            f"characteristic-{prime} truncation boundary changed",
        )
        require(
            all(next_block[prime][column] == 0 for column in range(prime + 1))
            and all(
                next_block[row][prime] == 0 for row in range(prime + 1)
            ),
            f"characteristic-{prime} s^p annihilator row changed",
        )
        rows.append((prime, determinant, prime, prime))
    return tuple(rows)


def cyclic_alias_controls():
    checks = 0
    for exponent in range(1, 7):
        for period in range(2, 18):
            value = factorial(exponent) - factorial(exponent + period)
            require(value != 0, "cyclic alias accidentally preserved factorial height")
            checks += 1
    hostile = factorial(1) - factorial(14)
    require(hostile == -87_178_291_199, "C13 cyclic hostile changed")
    return checks, hostile


def sparse_support_control():
    minus_first = 6 * factorial(1) - factorial(3)
    minus_square = (
        36 * factorial(2)
        - 12 * factorial(4)
        + factorial(6)
    )
    plus_first = 6 * factorial(1) + factorial(3)
    require(
        (minus_first, minus_square, plus_first) == (0, 504, 12),
        "two-slot factorial hostile changed",
    )
    gaussian_minus_second = 2 * minus_first
    gaussian_minus_fourth = comb(4, 2) * minus_square
    gaussian_plus_second = 2 * plus_first
    require(
        (
            gaussian_minus_second,
            gaussian_minus_fourth,
            gaussian_plus_second,
        )
        == (0, 3024, 24),
        "two-charge lift hostile changed",
    )
    return (
        minus_first,
        minus_square,
        plus_first,
        gaussian_minus_second,
        gaussian_minus_fourth,
        gaussian_plus_second,
    )


def main():
    determinants = char_zero_controls()
    prime_rows = characteristic_p_controls()
    alias_checks, cyclic_hostile = cyclic_alias_controls()
    sparse = sparse_support_control()

    print("THM-2810 FACTORIAL-HANKEL FAITHFULNESS")
    print("status=VERIFIED-EXACT candidate; NC2/GMC2 already proved elsewhere")
    print(
        "identity=H_d=D*Pascal*Pascal^T*D; "
        "det(H_d)=product_(j=0)^d(j!)^2"
    )
    print(
        f"char0_degrees=0..10 determinants_nonzero={len(determinants)}/"
        f"{len(determinants)} det_d10={determinants[-1]}"
    )
    print(
        "consequence=Ann(L)=0; no proper algebra quotient or fixed "
        "finite-dimensional linear/Prony carrier"
    )
    print(
        "finite_state_invoice=matching moments 0..2d forces "
        "state_dimension>=d+1"
    )
    print(
        f"cyclic_alias_checks={alias_checks}; "
        f"L(s*(1-s^13))={cyclic_hostile}; two_charge={2 * cyclic_hostile}"
    )
    print(
        "characteristic_p_boundary="
        + ";".join(
            f"p{prime}:det={determinant},rank_p={rank_p},rank_p1={rank_p1}"
            for prime, determinant, rank_p, rank_p1 in prime_rows
        )
        + "; Ann(L_p)=(s^p)"
    )
    print(
        "same_support_control="
        f"L(H-)={sparse[0]},L(H-^2)={sparse[1]},L(H+)={sparse[2]}; "
        f"Gaussian=(M2-,M4-,M2+)=({sparse[3]},{sparse[4]},{sparse[5]})"
    )
    print(
        "scope=no fixed tower-wide bounded carrier; growing-height, "
        "clock-dependent, nonlinear-resultant, and Frobenius whole-face "
        "mechanisms remain available; HYP-8765 not refuted"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
