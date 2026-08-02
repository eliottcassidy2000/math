#!/usr/bin/env python3
"""Exact controls for THM-3122's labelled-deletion kernel ghost."""

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def partitions(total, largest=None):
    if total == 0:
        yield ()
        return
    if largest is None or largest > total:
        largest = total
    for first in range(largest, 0, -1):
        for tail in partitions(total - first, first):
            yield (first,) + tail


def lower(shape, part):
    answer = list(shape)
    answer.remove(part)
    if part > 1:
        answer.append(part - 1)
    return tuple(sorted(answer, reverse=True))


def hook_kernel(n):
    answer = {}
    for j in range(n - 1):
        answer[(n - j,) + (1,) * j] = (-1) ** j * comb(n, j)
    answer[(1,) * n] = (-1) ** (n - 1) * (n - 1)
    return answer


def labelled_delete(vector):
    answer = defaultdict(int)
    for shape, coefficient in vector.items():
        for part in set(shape):
            answer[lower(shape, part)] += (
                coefficient * part * shape.count(part)
            )
    return {shape: value for shape, value in answer.items() if value}


def derangements(n):
    if n == 0:
        return 1
    if n == 1:
        return 0
    a, b = 1, 0
    for k in range(2, n + 1):
        a, b = b, (k - 1) * (a + b)
    return b


def group_coefficient(n, support):
    """Coefficient of a permutation with this support size in O_n."""
    value = Fraction(0)
    for j in range(n - 1):
        k = n - j
        weight = (-1) ** j * factorial(n) // factorial(j)
        projection = Fraction(0)
        if support <= k:
            projection = Fraction(
                factorial(n - support),
                factorial(n) * factorial(k - support),
            )
        value += weight * ((1 if support == 0 else 0) - projection)
    require(value.denominator == 1, "group coefficient lost integrality")
    return value.numerator


def compose(left, right):
    return tuple(left[right[x]] for x in range(len(left)))


def inverse(permutation):
    answer = [0] * len(permutation)
    for x, y in enumerate(permutation):
        answer[y] = x
    return tuple(answer)


def cycle_permutation(n, cycle):
    answer = list(range(n))
    for x, y in zip(cycle, cycle[1:] + cycle[:1]):
        answer[x] = y
    return tuple(answer)


def is_derangement(permutation):
    return all(x != y for x, y in enumerate(permutation))


def corners(shape):
    for row, length in enumerate(shape):
        if row == len(shape) - 1 or shape[row + 1] < length:
            answer = list(shape)
            answer[row] -= 1
            if answer[row] == 0:
                answer.pop(row)
            yield tuple(answer)


@lru_cache(maxsize=None)
def standard_tableaux(shape):
    if not shape:
        return 1
    return sum(standard_tableaux(smaller) for smaller in corners(shape))


@lru_cache(maxsize=None)
def skew_tableaux_after_first_row(shape, row_length):
    """f^(shape/(row_length)); equals K_(shape,(row_length,1^j))."""
    if sum(shape) == row_length:
        return int(shape == (row_length,))
    if not shape or shape[0] < row_length:
        return 0
    return sum(
        skew_tableaux_after_first_row(smaller, row_length)
        for smaller in corners(shape)
        if smaller and smaller[0] >= row_length
    )


def carrier_scalar(n, shape):
    dimension = standard_tableaux(shape)
    value = Fraction(0)
    for j in range(n - 1):
        k = n - j
        coefficient = (-1) ** j * comb(n, j) * factorial(k)
        kostka = skew_tableaux_after_first_row(shape, k)
        value += coefficient * (1 - Fraction(kostka, dimension))
    require(value.denominator == 1, "irreducible scalar lost integrality")
    return value.numerator


def main():
    print("labelled_deletion_hook_kernel_and_derangement_laplacian")
    print("n derangements irreps positive zero min_positive max_positive")
    for n in range(2, 11):
        require(labelled_delete(hook_kernel(n)) == {},
                f"hook current escaped labelled deletion at n={n}")
        count = derangements(n)
        require(group_coefficient(n, 0) == count,
                f"identity coefficient is not D_n at n={n}")
        for support in range(2, n + 1):
            expected = -1 if support == n else 0
            require(group_coefficient(n, support) == expected,
                    f"support coefficient failed at {(n, support)}")

        if n >= 4:
            alpha = cycle_permutation(n, (0, 2, 1) + tuple(range(3, n)))
            transposition = tuple(1 if x == 0 else 0 if x == 1 else x
                                  for x in range(n))
            beta = compose(inverse(alpha), transposition)
            require(is_derangement(alpha) and is_derangement(beta),
                    f"transposition factors are not derangements at n={n}")
            require(compose(alpha, beta) == transposition,
                    f"derangement factorization lost (01) at n={n}")

        scalars = [carrier_scalar(n, shape) for shape in partitions(n)]
        require(all(value >= 0 for value in scalars),
                f"negative finite-control scalar at n={n}")
        positive = [value for value in scalars if value > 0]
        expected_zero = 2 if n == 3 else 1
        require(len(scalars) - len(positive) == expected_zero,
                f"unexpected equality space at n={n}")
        print(n, count, len(scalars), len(positive),
              len(scalars) - len(positive), min(positive), max(positive))
    print("group_coefficients=identity_Dn;partial_support_0;full_support_-1")
    print("all_exact_checks_pass")


if __name__ == "__main__":
    main()
