#!/usr/bin/env python3
"""Exact finite controls for THM-2870.

The script works over prime fields using dependency-free row reduction.
Every truth-bearing check uses ``require`` so optimized execution retains
the complete audit.
"""

from collections import Counter
from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank_mod(matrix, prime):
    """Return the row rank over F_prime."""
    if not matrix:
        return 0
    work = [[entry % prime for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    require(
        all(len(row) == column_count for row in work),
        "ragged matrix",
    )
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if work[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            entry * inverse % prime for entry in work[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            factor = work[row][column] % prime
            if factor:
                work[row] = [
                    (left - factor * right) % prime
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def matrix_multiply(left, right, prime):
    require(left and right, "empty matrix multiplication")
    inner = len(right)
    require(all(len(row) == inner for row in left), "left shape")
    width = len(right[0])
    require(all(len(row) == width for row in right), "right shape")
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(inner))
            % prime
            for j in range(width)
        ]
        for i in range(len(left))
    ]


def identity(size):
    return [
        [1 if row == column else 0 for column in range(size)]
        for row in range(size)
    ]


def subtract(left, right, prime):
    require(len(left) == len(right), "subtraction row mismatch")
    return [
        [
            (a - b) % prime
            for a, b in zip(left_row, right_row)
        ]
        for left_row, right_row in zip(left, right)
    ]


def cyclic_convolution_matrix(mask):
    """C[i,j]=mask[i-j], so (Cf)(i)=sum_g mask[g]f(i-g)."""
    size = len(mask)
    return [
        [mask[(row - column) % size] for column in range(size)]
        for row in range(size)
    ]


def diagonal_matrix(mask):
    size = len(mask)
    return [
        [mask[row] if row == column else 0 for column in range(size)]
        for row in range(size)
    ]


def equation_matrix(c_matrix, mask, orientation, prime):
    """Linear system for CT=TD or DT=TC in row-major T variables."""
    size = len(mask)
    require(
        len(c_matrix) == size
        and all(len(row) == size for row in c_matrix),
        "operator shape",
    )
    equations = []
    for i in range(size):
        for j in range(size):
            row = [0] * (size * size)
            if orientation == "forward":
                # (CT-TD)_(i,j)
                for k in range(size):
                    row[k * size + j] += c_matrix[i][k]
                row[i * size + j] -= mask[j]
            elif orientation == "reverse":
                # (DT-TC)_(i,j)
                row[i * size + j] += mask[i]
                for k in range(size):
                    row[i * size + k] -= c_matrix[k][j]
            else:
                raise RuntimeError("unknown intertwiner orientation")
            equations.append([entry % prime for entry in row])
    return equations


def check_projection_unit(c_matrix, mask, prime, label):
    """Check both Hom-space dimensions against r*dim ker(C-I)."""
    size = len(mask)
    require(rank_mod(c_matrix, prime) == size, f"{label}: C not a unit")
    c_minus_one = subtract(c_matrix, identity(size), prime)
    eigen_dimension = size - rank_mod(c_minus_one, prime)
    physical_rank = sum(mask)
    expected = physical_rank * eigen_dimension
    dimensions = []
    for orientation in ("forward", "reverse"):
        equations = equation_matrix(
            c_matrix, mask, orientation, prime
        )
        dimension = size * size - rank_mod(equations, prime)
        require(
            dimension == expected,
            (
                f"{label}: {orientation} dimension {dimension} "
                f"!= {expected}"
            ),
        )
        dimensions.append(dimension)
    return eigen_dimension, tuple(dimensions)


def all_matrices(size, prime):
    for entries in product(range(prime), repeat=size * size):
        yield [
            list(entries[row * size : (row + 1) * size])
            for row in range(size)
        ]


def abstract_two_by_two_controls():
    prime = 3
    checked = 0
    histogram = Counter()
    for c_matrix in all_matrices(2, prime):
        if rank_mod(c_matrix, prime) != 2:
            continue
        eigen_dimension = 2 - rank_mod(
            subtract(c_matrix, identity(2), prime),
            prime,
        )
        for physical_rank in range(3):
            mask = [1 if index < physical_rank else 0 for index in range(2)]
            expected = physical_rank * eigen_dimension
            for orientation in ("forward", "reverse"):
                equations = equation_matrix(
                    c_matrix, mask, orientation, prime
                )
                dimension = 4 - rank_mod(equations, prime)
                require(
                    dimension == expected,
                    "abstract two-by-two classification failed",
                )
            checked += 1
            histogram[(physical_rank, eigen_dimension)] += 1
    require(checked == 144, "wrong abstract control count")
    return checked, histogram


def c5_boolean_controls():
    prime = 5
    size = 5
    checked = 0
    histogram = Counter()
    for bits in product((0, 1), repeat=size):
        mass = sum(bits)
        if mass % prime == 0:
            continue
        c_matrix = cyclic_convolution_matrix(bits)
        eigen_dimension, dimensions = check_projection_unit(
            c_matrix, list(bits), prime, f"C5 mask {bits}"
        )
        require(0 < mass < size, "C5 unit mask should be proper")
        require(
            dimensions[0] == dimensions[1],
            "orientation dimensions differ",
        )
        histogram[(mass, eigen_dimension, dimensions[0])] += 1
        checked += 1
    require(checked == 30, "wrong C5 Boolean-mask count")
    return checked, histogram


def cayley_tournaments(size, prime):
    require(size % 2 == 1, "tournament order must be odd")
    inverse_pairs = [
        (value, (-value) % size)
        for value in range(1, (size + 1) // 2)
    ]
    histogram = Counter()
    count = 0
    for choices in product((0, 1), repeat=len(inverse_pairs)):
        mask = [0] * size
        for choice, pair in zip(choices, inverse_pairs):
            mask[pair[choice]] = 1
        require(
            all(
                mask[value] + mask[(-value) % size] == 1
                for value in range(1, size)
            ),
            "invalid Cayley tournament",
        )
        c_matrix = cyclic_convolution_matrix(mask)
        eigen_dimension, dimensions = check_projection_unit(
            c_matrix, mask, prime, f"C{size} tournament {choices}"
        )
        require(
            rank_mod(diagonal_matrix(mask), prime)
            == (size - 1) // 2,
            "wrong physical tournament rank",
        )
        histogram[(eigen_dimension, dimensions[0])] += 1
        count += 1
    require(
        count == 2 ** ((size - 1) // 2),
        "wrong tournament count",
    )
    return count, histogram


def sharp_c3_maps():
    prime = 3
    size = 3
    mask = [0, 1, 0]
    c_matrix = cyclic_convolution_matrix(mask)
    d_matrix = diagonal_matrix(mask)

    # T_+(f)=f(1)1 and T_-(f)=(sum f)delta_1.
    forward = [[1 if column == 1 else 0 for column in range(size)]
               for _row in range(size)]
    reverse = [[1 for _column in range(size)] if row == 1
               else [0 for _column in range(size)]
               for row in range(size)]

    require(
        matrix_multiply(c_matrix, forward, prime)
        == matrix_multiply(forward, d_matrix, prime),
        "sharp forward map failed",
    )
    require(
        matrix_multiply(d_matrix, reverse, prime)
        == matrix_multiply(reverse, c_matrix, prime),
        "sharp reverse map failed",
    )
    require(
        rank_mod(forward, prime) == rank_mod(reverse, prime) == 1,
        "sharp maps are not rank one",
    )
    eigen_dimension, dimensions = check_projection_unit(
        c_matrix, mask, prime, "C3 singleton"
    )
    require(
        eigen_dimension == 1 and dimensions == (1, 1),
        "C3 sharp dimension failed",
    )
    return eigen_dimension, dimensions


def mixed_prime_hostile():
    mask = [1, 0, 1, 0, 1, 0]
    c_matrix = cyclic_convolution_matrix(mask)
    require(sum(mask) % 2 == 1, "hostile mass is not a 2-unit")
    rank_two = rank_mod(c_matrix, 2)
    rank_large = rank_mod(c_matrix, 101)
    require(
        rank_two == rank_large == 2,
        "C6 mixed-prime hostile rank changed",
    )
    return rank_two


def format_histogram(histogram):
    return ";".join(
        f"{key}:{histogram[key]}" for key in sorted(histogram)
    )


def main():
    abstract_count, abstract_histogram = abstract_two_by_two_controls()
    c5_count, c5_histogram = c5_boolean_controls()
    c5_tournament_count, c5_tournament_histogram = cayley_tournaments(5, 5)
    c9_tournament_count, c9_tournament_histogram = cayley_tournaments(9, 3)
    c3_eigen_dimension, c3_dimensions = sharp_c3_maps()
    c6_rank = mixed_prime_hostile()

    print(
        "PRIME-POWER CONVOLUTION / PHYSICAL-DIAGONAL "
        "INTERTWINER AUDIT - exact"
    )
    print(
        "abstract F3^2 unit/projection controls: "
        f"{abstract_count}; hist={format_histogram(abstract_histogram)}"
    )
    print(
        "classification: CT=TD factors im(D)->ker(C-1); "
        "DT=TC factors coker(C-1)->im(D)"
    )
    print(
        "dimension in both orientations: "
        "rank(D)*dim ker(C-1); no semisimplicity used"
    )
    print(
        f"C5 proper Boolean p-unit masks: {c5_count}; "
        f"hist={format_histogram(c5_histogram)}"
    )
    print(
        f"C5 Cayley tournaments: {c5_tournament_count}; "
        f"hist={format_histogram(c5_tournament_histogram)}"
    )
    print(
        f"C9 Cayley tournaments: {c9_tournament_count}; "
        f"hist={format_histogram(c9_tournament_histogram)}"
    )
    print(
        "C3 singleton sharpness: "
        f"e={c3_eigen_dimension}; Hom dimensions={c3_dimensions}; "
        "both explicit maps have rank 1"
    )
    print(
        "C6 mixed-prime hostile: mass=3 is a 2-unit but "
        f"convolution rank={c6_rank}<6"
    )
    print(
        "Fourier boundary: convolution diagonalizes to hhat, not h; "
        "a prime-power unit hhat has no zero while a proper Boolean h does"
    )
    print(
        "scope: same-mask linear obstruction only; "
        "no physical translation, positive inverse, endpoint current, "
        "target/carry bridge, or LRC decrement"
    )
    print("PASS")


if __name__ == "__main__":
    main()
