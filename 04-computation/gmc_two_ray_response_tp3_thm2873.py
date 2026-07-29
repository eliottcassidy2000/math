#!/usr/bin/env python3
"""Exact companion for THM-2873.

For

    L(s^k) = k!,  f_j = s^j/j!,  d_j = f_(j+1)-f_j,

the theorem concerns the consecutive response determinant

    det [ L(d_(n+r) W^m) ]_(r=0,1,2; m=1,2,3)

for W = x*d_p + y*d_(p+h), h in {1,2}.  A common-factorial
normalization leaves seven integer polynomials.  After n=z+1 every
coefficient of every residual polynomial is strictly positive.

All truth-bearing gates raise explicitly, so normal and optimized Python
execute the same audit.
"""

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
from math import comb, factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def convolution(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            out[i + j] += first * second
    return out


def determinant_coefficients(entries):
    """Determinant of a 3x3 matrix whose entries are coefficient lists."""

    out = [0] * 7
    for column_permutation in permutations(range(3)):
        inversions = sum(
            column_permutation[i] > column_permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        sign = -1 if inversions % 2 else 1
        term = [1]
        for row, column in enumerate(column_permutation):
            term = convolution(term, entries[row][column])
        for exponent, coefficient in enumerate(term):
            out[exponent] += sign * coefficient
    return out


@lru_cache(maxsize=None)
def tensor(indices):
    """Closed inclusion-exclusion formula for L(prod_i d_(indices[i]))."""

    indices = tuple(indices)
    block_sizes = [index + 1 for index in indices]
    denominator = 1
    for block_size in block_sizes:
        denominator *= factorial(block_size)

    elementary = [1]
    for block_size in block_sizes:
        elementary.append(0)
        for degree in range(len(elementary) - 1, 0, -1):
            elementary[degree] += block_size * elementary[degree - 1]

    total_blocks = sum(block_sizes)
    numerator = sum(
        (-1) ** degree
        * elementary[degree]
        * factorial(total_blocks - degree)
        for degree in range(len(indices) + 1)
    )
    return Fraction(numerator, denominator)


def direct_two_ray_coefficients(row, pivot, gap):
    entries = []
    for shift in range(3):
        matrix_row = []
        for power in (1, 2, 3):
            column = []
            for high_count in range(power + 1):
                labels = (
                    (pivot,) * (power - high_count)
                    + (pivot + gap,) * high_count
                )
                column.append(
                    comb(power, high_count)
                    * tensor(tuple(sorted((row + shift,) + labels)))
                )
            matrix_row.append(column)
        entries.append(matrix_row)
    return determinant_coefficients(entries)


n_symbol, p_symbol, z_symbol = sp.symbols(
    "n p z", integer=True, nonnegative=True
)


def rising(base, length):
    out = sp.Integer(1)
    for shift in range(length):
        out *= base + shift
    return sp.expand(out)


def elementary_symmetric(values):
    out = [sp.Integer(1)]
    for value in values:
        out.append(sp.Integer(0))
        for degree in range(len(out) - 1, 0, -1):
            out[degree] = sp.expand(
                out[degree] + value * out[degree - 1]
            )
    return out


def normalized_entry(power, high_count, row_shift, gap):
    """The polynomial left after removing the common factorial blocks."""

    total_exponent = (
        n_symbol
        + row_shift
        + power * p_symbol
        + high_count * gap
    )
    block_sizes = (
        [n_symbol + row_shift + 1]
        + [p_symbol + 1] * (power - high_count)
        + [p_symbol + gap + 1] * high_count
    )
    elementary = elementary_symmetric(block_sizes)
    avoider_residual = sum(
        (-1) ** degree
        * elementary[degree]
        * rising(
            total_exponent + 1,
            power + 1 - degree,
        )
        for degree in range(power + 2)
    )
    factorial_shift = rising(
        n_symbol + power * p_symbol + 1,
        row_shift + high_count * gap,
    )
    return sp.expand(factorial_shift * avoider_residual)


def normalized_residuals(gap):
    entries = []
    for row_shift in range(3):
        matrix_row = []
        for power in (1, 2, 3):
            matrix_row.append(
                [
                    comb(power, high_count)
                    * normalized_entry(
                        power,
                        high_count,
                        row_shift,
                        gap,
                    )
                    for high_count in range(power + 1)
                ]
            )
        entries.append(matrix_row)
    return [sp.expand(value) for value in determinant_coefficients(entries)]


def common_factor(row, pivot, gap, high_count):
    return Fraction(
        factorial(row + pivot)
        * factorial(row + 2 * pivot)
        * factorial(row + 3 * pivot),
        factorial(row + 1)
        * factorial(row + 2)
        * factorial(row + 3)
        * factorial(pivot + 1) ** (6 - high_count)
        * factorial(pivot + gap + 1) ** high_count,
    )


def canonical_polynomial_digest(polynomials):
    records = []
    for residual_index, polynomial in enumerate(polynomials):
        shifted = sp.Poly(
            sp.expand(polynomial.subs(n_symbol, z_symbol + 1)),
            p_symbol,
            z_symbol,
        )
        for monomial, coefficient in shifted.terms():
            records.append(
                f"{residual_index}:{monomial[0]}:{monomial[1]}:{coefficient}"
            )
    return sha256(("\n".join(records) + "\n").encode()).hexdigest()


def singleton_curvature_polynomial(index, row):
    return (
        8 * index**5 * row
        - 4 * index**5
        + 28 * index**4 * row**2
        + 50 * index**4 * row
        - 5 * index**4
        + 38 * index**3 * row**3
        + 150 * index**3 * row**2
        + 173 * index**3 * row
        + 53 * index**3
        + 25 * index**2 * row**4
        + 152 * index**2 * row**3
        + 330 * index**2 * row**2
        + 302 * index**2 * row
        + 99 * index**2
        + 8 * index * row**5
        + 65 * index * row**4
        + 203 * index * row**3
        + 305 * index * row**2
        + 221 * index * row
        + 62 * index
        + row**6
        + 10 * row**5
        + 40 * row**4
        + 82 * row**3
        + 91 * row**2
        + 52 * row
        + 12
    )


def singleton_curvature(index, row):
    polynomial = singleton_curvature_polynomial(index, row)
    return Fraction(
        12
        * polynomial
        * factorial(2 * index + row)
        * factorial(3 * index + row),
        (row + 1)
        * (row + 2)
        * (row + 3)
        * factorial(index) ** 3
        * factorial(index + row + 1)
        * factorial(index + row + 2),
    )


def determinant(matrix):
    (a, b, c), (d, e, f), (g, h, i) = matrix
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def direct_singleton_curvature(index, row):
    matrix = [
        [
            tensor(tuple(sorted((row + shift,) + (index,) * power)))
            for power in (1, 2, 3)
        ]
        for shift in range(3)
    ]
    response_determinant = determinant(matrix)
    denominator = matrix[0][0] * matrix[1][0] * matrix[2][0]
    return response_determinant / denominator


def polarized_kernel(rows, labels):
    total = Fraction(0)
    positions = set(range(6))
    for singleton in range(6):
        remaining = positions - {singleton}
        for pair in combinations(sorted(remaining), 2):
            triple = tuple(sorted(remaining - set(pair)))
            matrix = [
                [
                    tensor(tuple(sorted((row, labels[singleton])))),
                    tensor(
                        tuple(
                            sorted(
                                (
                                    row,
                                    labels[pair[0]],
                                    labels[pair[1]],
                                )
                            )
                        )
                    ),
                    tensor(
                        tuple(
                            sorted(
                                (
                                    row,
                                    labels[triple[0]],
                                    labels[triple[1]],
                                    labels[triple[2]],
                                )
                            )
                        )
                    ),
                ]
                for row in rows
            ]
            total += determinant(matrix)
    return total


def universal_newton_scout():
    multiset_count = 0
    coefficient_count = 0
    minimum = None
    for labels in combinations_with_replacement(range(7), 6):
        degree = sum(labels)
        values = [
            polarized_kernel((row, row + 1, row + 2), labels)
            for row in range(1, degree + 2)
        ]
        for order in range(degree + 1):
            coefficient = values[0]
            require(
                coefficient > 0,
                f"universal Newton scout failed: {labels}, order={order}",
            )
            if minimum is None or coefficient < minimum:
                minimum = coefficient
            coefficient_count += 1
            values = [
                values[index + 1] - values[index]
                for index in range(len(values) - 1)
            ]
        multiset_count += 1
    return multiset_count, coefficient_count, minimum


def main():
    print("THM-2873 exact two-ray factorial-response TP3 certificate")

    residuals_by_gap = {}
    normalization_cells = 0
    for gap in (1, 2):
        residuals = normalized_residuals(gap)
        residuals_by_gap[gap] = residuals
        term_counts = []
        minima = []
        for high_count, residual in enumerate(residuals):
            shifted = sp.Poly(
                sp.expand(residual.subs(n_symbol, z_symbol + 1)),
                p_symbol,
                z_symbol,
            )
            coefficients = shifted.coeffs()
            require(coefficients, "empty residual certificate")
            require(
                all(coefficient > 0 for coefficient in coefficients),
                f"nonpositive residual coefficient: gap={gap}, "
                f"high_count={high_count}",
            )
            term_counts.append(len(coefficients))
            minima.append(min(coefficients))

        for row in range(1, 4):
            for pivot in range(4):
                direct = direct_two_ray_coefficients(row, pivot, gap)
                for high_count, direct_coefficient in enumerate(direct):
                    residual_value = residuals[high_count].subs(
                        {n_symbol: row, p_symbol: pivot}
                    )
                    reconstructed = (
                        Fraction(int(residual_value))
                        * common_factor(row, pivot, gap, high_count)
                    )
                    require(
                        reconstructed == direct_coefficient,
                        f"normalization mismatch: gap={gap}, row={row}, "
                        f"pivot={pivot}, high_count={high_count}",
                    )
                    normalization_cells += 1

        print(f"gap={gap}")
        print("shifted_term_counts=" + ",".join(map(str, term_counts)))
        print("shifted_min_coefficients=" + ",".join(map(str, minima)))
        print(
            "certificate_sha256="
            + canonical_polynomial_digest(residuals)
        )

    require(
        sp.expand(
            singleton_curvature_polynomial(p_symbol, z_symbol + 1)
            - singleton_curvature_polynomial(
                p_symbol,
                n_symbol,
            ).subs(n_symbol, z_symbol + 1)
        )
        == 0,
        "singleton shift substitution",
    )
    singleton_shifted = sp.Poly(
        sp.expand(
            singleton_curvature_polynomial(p_symbol, n_symbol).subs(
                n_symbol, z_symbol + 1
            )
        ),
        p_symbol,
        z_symbol,
    )
    require(
        all(coefficient > 0 for coefficient in singleton_shifted.coeffs()),
        "singleton n>=1 positivity",
    )
    for index in range(9):
        for row in range(6):
            require(
                singleton_curvature(index, row)
                == direct_singleton_curvature(index, row),
                f"singleton formula mismatch: index={index}, row={row}",
            )

    boundary_values = [
        singleton_curvature_polynomial(index, 0)
        for index in range(4)
    ]
    require(boundary_values == [12, 217, 748, 1143], "n=0 positive prefix")
    boundary_symbol = sp.symbols("k", integer=True, nonnegative=True)
    boundary_tail = sp.Poly(
        -sp.expand(
            singleton_curvature_polynomial(
                boundary_symbol + 4,
                0,
            )
        ),
        boundary_symbol,
    )
    require(
        boundary_tail.all_coeffs() == [4, 85, 667, 2305, 3002, 140],
        "n=0 negative tail certificate",
    )
    hostile_curvature = direct_singleton_curvature(4, 0)
    require(hostile_curvature == -4527600, "d4 n=0 hostile")

    # The semiring-output measures are not themselves TP2.  For W=d_1,
    # the power-one measure has mass 1 at index 1, while the power-two
    # prefix has mass 2 at indices 0 and 1.  The rows 0,1 minor is -2.
    raw_output_minor = 0 * 2 - 1 * 2
    require(raw_output_minor == -2, "raw output-measure hostile")

    scout_multisets, scout_coefficients, scout_minimum = (
        universal_newton_scout()
    )

    print(f"normalization_cells={normalization_cells}")
    print(
        "singleton_shifted_terms="
        f"{len(singleton_shifted.terms())};"
        f"minimum={min(singleton_shifted.coeffs())}"
    )
    print("singleton_n0_polynomial_values_j0_to_j3=12,217,748,1143")
    print("singleton_n0_negative_for_j_ge_4=true")
    print(f"hostile_d4_n0_curvature={hostile_curvature}")
    print(f"raw_output_measure_minor={raw_output_minor}")
    print(
        "universal_newton_scout="
        f"{scout_multisets}_label_multisets;"
        f"{scout_coefficients}_positive_coefficients;"
        f"minimum={scout_minimum}"
    )
    print("scope=two-ray gaps 1,2 proved; universal six-label TP3 remains open")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
