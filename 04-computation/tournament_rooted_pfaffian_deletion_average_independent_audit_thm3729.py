#!/usr/bin/env python3
"""Independent exact audit of rooted-Pfaffian averages (12)--(13).

This does not import the candidate probe.  It uses Fraction Gaussian
elimination, checks every tournament and every incident sign word for
2 <= n <= 4, and also checks the candidate's five-vertex hostile control.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import product


def identity(n):
    return [[Fraction(i == j) for j in range(n)] for i in range(n)]


def add(left, right):
    return [[a + b for a, b in zip(row_a, row_b)]
            for row_a, row_b in zip(left, right)]


def multiply(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))]
            for i in range(len(left))]


def determinant(matrix):
    n = len(matrix)
    if n == 0:
        return Fraction(1)
    work = [[Fraction(value) for value in row] for row in matrix]
    answer = Fraction(1)
    for column in range(n):
        pivot = next((row for row in range(column, n)
                      if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_value = work[column][column]
        answer *= pivot_value
        for entry in range(column, n):
            work[column][entry] /= pivot_value
        for row in range(column + 1, n):
            factor = work[row][column]
            for entry in range(column, n):
                work[row][entry] -= factor * work[column][entry]
    return answer


def inverse(matrix):
    n = len(matrix)
    work = [
        [Fraction(value) for value in row] + identity(n)[index]
        for index, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next((row for row in range(column, n)
                      if work[row][column]), None)
        if pivot is None:
            raise RuntimeError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [value / pivot_value for value in work[column]]
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b
                         for a, b in zip(work[row], work[column])]
    return [row[n:] for row in work]


def trace(matrix):
    return sum(matrix[i][i] for i in range(len(matrix)))


def tournament_matrix(n, mask):
    matrix = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            value = Fraction(1 if (mask >> bit) & 1 else -1)
            matrix[i][j] = value
            matrix[j][i] = -value
            bit += 1
    return matrix


def principal(matrix, kept):
    return [[matrix[i][j] for j in kept] for i in kept]


def outer(word):
    return [[Fraction(a * b) for b in word] for a in word]


def augment_skew(K, word):
    n = len(K)
    out = [[Fraction(0) for _ in range(n + 1)] for _ in range(n + 1)]
    for i in range(n):
        for j in range(n):
            out[i][j] = K[i][j]
        out[i][n] = Fraction(word[i])
        out[n][i] = Fraction(-word[i])
    return out


def check_case(K):
    n = len(K)
    if n < 2:
        raise RuntimeError("deletion-disc normalization requires n >= 2")
    M = add(identity(n), K)
    det_M = determinant(M)
    disc = det_M / (1 << (n - 1))

    inv_M = inverse(M)
    even_inverse = inverse(add(identity(n),
                               [[-value for value in row]
                                for row in multiply(K, K)]))
    if trace(inv_M) != trace(even_inverse):
        raise RuntimeError("trace (I-K^2)^-1 mismatch")

    deletion_determinants = []
    deletion_discs = []
    for vertex in range(n):
        kept = tuple(i for i in range(n) if i != vertex)
        deletion_det = determinant(add(identity(n - 1), principal(K, kept)))
        deletion_determinants.append(deletion_det)
        deletion_discs.append(deletion_det / (1 << (n - 2)))
        if inv_M[vertex][vertex] != deletion_det / det_M:
            raise RuntimeError("diagonal cofactor/deletion mismatch")

    odd_energies = []
    extension_discs = []
    for word in product((-1, 1), repeat=n):
        numerator = determinant(add(M, outer(word))) - det_M
        odd_energy = numerator / (1 << (n - 1))
        odd_energies.append(odd_energy)

        Khat = augment_skew(K, word)
        extension_disc = determinant(add(identity(n + 1), Khat)) / (1 << n)
        extension_discs.append(extension_disc)
        if odd_energy != 2 * extension_disc - disc:
            raise RuntimeError("one-root extension normalization mismatch")

    average_odd = sum(odd_energies) / (1 << n)
    if average_odd != disc * trace(inv_M):
        raise RuntimeError("Rademacher trace average mismatch")
    if average_odd != sum(deletion_discs) / 2:
        raise RuntimeError("deletion average (12) mismatch")

    average_extension = sum(extension_discs) / (1 << n)
    if average_extension != disc / 2 + sum(deletion_discs) / 4:
        raise RuntimeError("extension average (13) mismatch")

    return (
        n, det_M, disc, trace(inv_M), tuple(deletion_determinants),
        average_odd, average_extension,
    )


def main():
    records = []
    counts = []
    for n in range(2, 5):
        total = 1 << (n * (n - 1) // 2)
        for mask in range(total):
            records.append(check_case(tournament_matrix(n, mask)))
        counts.append((n, total))

    hostile = check_case(tournament_matrix(5, 0b1011010011))
    expected_hostile = (
        5, Fraction(32), Fraction(2), Fraction(7, 4),
        (Fraction(8), Fraction(8), Fraction(16), Fraction(16), Fraction(8)),
        Fraction(7, 2), Fraction(11, 4),
    )
    if hostile != expected_hostile:
        raise RuntimeError(("five-vertex hostile mismatch", hostile))

    semantic = (tuple(counts), tuple(records), hostile)
    print("exhaustive_counts=" + repr(tuple(counts)))
    print("n2_boundary_records=" + repr(tuple(
        record for record in records if record[0] == 2
    )))
    print("five_vertex_hostile=" + repr(hostile))
    print("semantic_sha256=" + sha256(repr(semantic).encode("ascii")).hexdigest())
    print("all_independent_deletion_checks=PASS")


if __name__ == "__main__":
    main()
