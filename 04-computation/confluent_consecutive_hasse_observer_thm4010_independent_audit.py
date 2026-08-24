#!/usr/bin/env python3
"""Exact independent audit for THM-4010's consecutive Hasse-jet sampler.

No THM-4000 producer module is imported.  The Hasse matrices, Bareiss
determinants, Euclidean Smith reductions, determinantal-divisor hostile,
polynomial division, and semantic transcript are implemented directly.
"""

from hashlib import sha256
from itertools import combinations
import json
from math import comb, factorial, gcd
import random
import sys


sys.stdout.reconfigure(newline="\n")

EXPECTED_SEMANTIC_SHA256 = "23c5b3fbbaeb75be24ffb5d2221abaff12697bade93a5b0d199c9f11c8d97af4"
GATES = 0


def gate(condition, label):
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(label)


def egcd(a, b):
    """Return nonnegative gcd and Bezout coefficients for signed inputs."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    gate(old_s * a + old_t * b == old_r, "Bezout identity")
    return old_r, old_s, old_t


def bareiss_det(matrix):
    """Fraction-free determinant for a square integer matrix."""
    n = len(matrix)
    gate(n > 0 and all(len(row) == n for row in matrix), "square determinant")
    work = [list(row) for row in matrix]
    if n == 1:
        return work[0][0]
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        pivot_row = next((row for row in range(pivot_index, n)
                          if work[row][pivot_index] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                numerator = (work[row][column] * pivot
                             - work[row][pivot_index] * work[pivot_index][column])
                gate(numerator % previous == 0, "Bareiss exact division")
                work[row][column] = numerator // previous
        for row in range(pivot_index + 1, n):
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def smith_diagonal(matrix):
    """Smith invariant factors via explicit unimodular Euclidean moves."""
    work = [list(row) for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    gate(rows > 0 and all(len(row) == columns for row in work), "rectangular Smith input")
    bound = min(rows, columns)

    for pivot_index in range(bound):
        location = next(((row, column)
                         for row in range(pivot_index, rows)
                         for column in range(pivot_index, columns)
                         if work[row][column] != 0), None)
        if location is None:
            break
        row, column = location
        work[pivot_index], work[row] = work[row], work[pivot_index]
        for current_row in work:
            current_row[pivot_index], current_row[column] = (
                current_row[column], current_row[pivot_index]
            )

        while True:
            changed = True
            while changed:
                changed = False
                for row in range(pivot_index + 1, rows):
                    b = work[row][pivot_index]
                    if b == 0:
                        continue
                    a = work[pivot_index][pivot_index]
                    if b % a == 0:
                        quotient = b // a
                        work[row] = [v - quotient * u
                                     for u, v in zip(work[pivot_index], work[row])]
                        gate(work[row][pivot_index] == 0, "Smith divisible row clear")
                        changed = True
                        continue
                    divisor, x, y = egcd(a, b)
                    old_pivot = work[pivot_index][:]
                    old_row = work[row][:]
                    work[pivot_index] = [x * u + y * v
                                         for u, v in zip(old_pivot, old_row)]
                    work[row] = [(-b // divisor) * u + (a // divisor) * v
                                 for u, v in zip(old_pivot, old_row)]
                    gate(work[pivot_index][pivot_index] == divisor
                         and work[row][pivot_index] == 0,
                         "Smith row gcd move")
                    changed = True

                for column in range(pivot_index + 1, columns):
                    b = work[pivot_index][column]
                    if b == 0:
                        continue
                    a = work[pivot_index][pivot_index]
                    if b % a == 0:
                        quotient = b // a
                        for row in range(rows):
                            work[row][column] -= quotient * work[row][pivot_index]
                        gate(work[pivot_index][column] == 0,
                             "Smith divisible column clear")
                        changed = True
                        continue
                    divisor, x, y = egcd(a, b)
                    old_pivot = [work[row][pivot_index] for row in range(rows)]
                    old_column = [work[row][column] for row in range(rows)]
                    for row in range(rows):
                        work[row][pivot_index] = x * old_pivot[row] + y * old_column[row]
                        work[row][column] = ((-b // divisor) * old_pivot[row]
                                             + (a // divisor) * old_column[row])
                    gate(work[pivot_index][pivot_index] == divisor
                         and work[pivot_index][column] == 0,
                         "Smith column gcd move")
                    changed = True

            pivot = work[pivot_index][pivot_index]
            gate(pivot != 0, "nonzero Smith pivot")
            offender = next(((row, column)
                             for row in range(pivot_index + 1, rows)
                             for column in range(pivot_index + 1, columns)
                             if work[row][column] % pivot != 0), None)
            if offender is None:
                break
            row, _column = offender
            work[pivot_index] = [u + v for u, v in zip(work[pivot_index], work[row])]

        if work[pivot_index][pivot_index] < 0:
            work[pivot_index] = [-entry for entry in work[pivot_index]]

    diagonal = tuple(abs(work[index][index]) for index in range(bound))
    gate(all(work[row][column] == 0
             for row in range(rows) for column in range(columns)
             if row != column), "Smith result diagonal")
    nonzero = tuple(entry for entry in diagonal if entry)
    gate(all(right % left == 0 for left, right in zip(nonzero, nonzero[1:])),
         "Smith divisibility chain")
    return diagonal


def hasse_matrix(a, m, k):
    size = (m + 1) * k
    return [
        [comb(degree, order) * (node ** (degree - order))
         if degree >= order else 0
         for degree in range(size)]
        for node in range(a, a + m + 1)
        for order in range(k)
    ]


def arbitrary_node_hasse_matrix(nodes, k):
    size = len(nodes) * k
    return [
        [comb(degree, order) * (node ** (degree - order))
         if degree >= order else 0
         for degree in range(size)]
        for node in nodes
        for order in range(k)
    ]


def ordinary_derivative_matrix(a, m, k):
    size = (m + 1) * k
    return [
        [factorial(order) * comb(degree, order) * (node ** (degree - order))
         if degree >= order else 0
         for degree in range(size)]
        for node in range(a, a + m + 1)
        for order in range(k)
    ]


def expected_determinant(m, k):
    return product(factorial(j) ** (k * k) for j in range(1, m + 1))


def product(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def factorization(value):
    value = abs(value)
    factors = []
    prime = 2
    while prime * prime <= value:
        exponent = 0
        while value % prime == 0:
            value //= prime
            exponent += 1
        if exponent:
            factors.append((prime, exponent))
        prime += 1 if prime == 2 else 2
    if value > 1:
        factors.append((value, 1))
    return tuple(factors)


def valuation(value, prime):
    exponent = 0
    while value and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def rank_mod_prime(matrix, prime):
    work = [[entry % prime for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows)
                      if work[row][column] != 0), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(inverse * entry) % prime for entry in work[rank]]
        for row in range(rows):
            if row == rank or work[row][column] == 0:
                continue
            multiplier = work[row][column]
            work[row] = [(u - multiplier * v) % prime
                         for u, v in zip(work[row], work[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def naive_smith(m, k):
    return tuple((factorial(j) ** k)
                 for j in range(m + 1) for _ in range(k))


def submatrix(matrix, row_indices, column_indices):
    return [[matrix[row][column] for column in column_indices]
            for row in row_indices]


def determinantal_divisors(matrix):
    """Exhaustive exact gcds of minors; used only on the rank-eight hostile."""
    n = len(matrix)
    answer = [1]
    for rank in range(1, n + 1):
        divisor = 0
        for rows in combinations(range(n), rank):
            for columns in combinations(range(n), rank):
                divisor = gcd(divisor, abs(bareiss_det(submatrix(matrix, rows, columns))))
                if divisor == 1:
                    break
            if divisor == 1:
                break
        answer.append(divisor)
    return tuple(answer)


def smith_from_divisors(divisors):
    result = []
    for left, right in zip(divisors, divisors[1:]):
        gate(left != 0 and right % left == 0, "determinantal divisor chain")
        result.append(right // left)
    return tuple(result)


def poly_trim(poly):
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add(left, right):
    size = max(len(left), len(right))
    return poly_trim([(left[i] if i < len(left) else 0)
                      + (right[i] if i < len(right) else 0)
                      for i in range(size)])


def poly_mul(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, u in enumerate(left):
        for j, v in enumerate(right):
            result[i + j] += u * v
    return poly_trim(result)


def poly_pow(poly, exponent):
    result = [1]
    base = list(poly)
    while exponent:
        if exponent & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        exponent >>= 1
    return result


def poly_eval(poly, value):
    answer = 0
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def hasse_eval(poly, value, order):
    return sum(coefficient * comb(degree, order) * value ** (degree - order)
               for degree, coefficient in enumerate(poly) if degree >= order)


def monic_divmod(dividend, divisor):
    gate(divisor and divisor[-1] == 1, "monic divisor")
    remainder = poly_trim(dividend)
    quotient = [0] * max(1, len(remainder) - len(divisor) + 1)
    while len(remainder) >= len(divisor) and remainder != [0]:
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1]
        quotient[shift] += coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        remainder = poly_trim(remainder)
    return poly_trim(quotient), poly_trim(remainder)


def node_polynomial(a, m):
    result = [1]
    for node in range(a, a + m + 1):
        result = poly_mul(result, [-node, 1])
    return result


def mutate_unimodular(matrix, seed):
    rng = random.Random(seed)
    result = [row[:] for row in matrix]
    n = len(result)
    for _ in range(5 * n):
        if rng.randrange(2) == 0:
            i, j = rng.sample(range(n), 2)
            if rng.randrange(3) == 0:
                result[i], result[j] = result[j], result[i]
            else:
                multiplier = rng.choice((-2, -1, 1, 2))
                result[i] = [u + multiplier * v for u, v in zip(result[i], result[j])]
        else:
            i, j = rng.sample(range(n), 2)
            if rng.randrange(3) == 0:
                for row in result:
                    row[i], row[j] = row[j], row[i]
            else:
                multiplier = rng.choice((-2, -1, 1, 2))
                for row in result:
                    row[i] += multiplier * row[j]
    return result


def atlas_row(m, k):
    matrix = hasse_matrix(0, m, k)
    determinant = abs(bareiss_det(matrix))
    diagonal = smith_diagonal(matrix)
    gate(determinant == expected_determinant(m, k), ("determinant formula", m, k))
    gate(product(diagonal) == determinant, ("Smith determinant", m, k))
    gate(all(prime <= m for entry in diagonal for prime, _exponent in factorization(entry)),
         ("prime support", m, k))
    return (m, k, len(matrix), determinant, diagonal, diagonal == naive_smith(m, k))


def main():
    # A substantial independent Smith atlas: 28 square matrices, up to rank 28.
    atlas = tuple(atlas_row(m, k) for m in range(0, 7) for k in range(1, 5))

    # A second exact path (gcds of every minor) covers every atlas matrix
    # through rank nine, including both first false corners.
    exhaustive_crosschecks = []
    for m, k, size, _determinant, diagonal, _naive in atlas:
        if size > 9:
            continue
        divisors = determinantal_divisors(hasse_matrix(0, m, k))
        recovered = smith_from_divisors(divisors)
        gate(recovered == diagonal, ("exhaustive minor crosscheck", m, k))
        exhaustive_crosschecks.append((m, k, divisors))

    # Correct general primewise information: exact mod-p rank and total valuation.
    primewise_controls = []
    for m, k, size, _determinant, diagonal, _naive in atlas:
        matrix = hasse_matrix(0, m, k)
        for prime in (2, 3, 5, 7):
            actual_rank = rank_mod_prime(matrix, prime)
            expected_rank = k * min(m + 1, prime)
            unit_count = sum(entry % prime != 0 for entry in diagonal)
            total_exponent = sum(valuation(entry, prime) for entry in diagonal)
            expected_exponent = k * k * sum(
                valuation(factorial(j), prime) for j in range(1, m + 1)
            )
            gate(actual_rank == expected_rank == unit_count,
                 ("primewise rank formula", m, k, prime))
            gate(total_exponent == expected_exponent,
                 ("primewise determinant valuation", m, k, prime))
            primewise_controls.append(
                (m, k, prime, actual_rank, total_exponent)
            )
        expected_units = size if m <= 1 else 2 * k
        gate(diagonal[:expected_units] == (1,) * expected_units,
             ("unit Smith prefix", m, k))
        if expected_units < size:
            gate(diagonal[expected_units] > 1, ("unit Smith prefix exact", m, k))

    # Translation is an integral unitriangular column change.
    translation_controls = []
    for a, m, k in ((1, 3, 2), (7, 4, 3), (23, 6, 4)):
        base = next(row for row in atlas if row[0:2] == (m, k))
        shifted = hasse_matrix(a, m, k)
        shifted_det = abs(bareiss_det(shifted))
        shifted_snf = smith_diagonal(shifted)
        gate(shifted_det == base[3] and shifted_snf == base[4],
             ("translation invariance", a, m, k))
        translation_controls.append((a, m, k, shifted_det, shifted_snf))

    # Minimal determinant-preserving hostile to the naive power/repeat Smith law.
    hostile = hasse_matrix(0, 3, 2)
    hostile_snf = smith_diagonal(hostile)
    hostile_naive = naive_smith(3, 2)
    gate(hostile_snf == (1, 1, 1, 1, 4, 4, 12, 108), "minimal hostile SNF")
    gate(hostile_naive == (1, 1, 1, 1, 4, 4, 36, 36), "minimal naive proposal")
    gate(product(hostile_snf) == product(hostile_naive), "hostile determinants agree")
    hostile_divisors = determinantal_divisors(hostile)
    gate(hostile_divisors == (1, 1, 1, 1, 1, 4, 16, 192, 20736),
         "hostile determinantal divisors")
    witness_rows = (0, 1, 2, 3, 4, 5, 7)
    witness_columns = (0, 1, 2, 3, 4, 5, 6)
    witness_minor = bareiss_det(submatrix(hostile, witness_rows, witness_columns))
    gate(witness_minor == 2112 and witness_minor % 576 != 0,
         "7-minor directly kills naive Smith law")

    # Minimality by rank: the only nontrivial lower-rank corner (m,k)=(2,2) passes.
    gate(next(row for row in atlas if row[0:2] == (2, 2))[4]
         == naive_smith(2, 2) == (1, 1, 1, 1, 4, 4),
         "lower-rank naive positive control")
    gate(all(row[5] for row in atlas if row[0] <= 1 or row[1] == 1),
         "unimodular and ordinary positive controls")

    # An exact partial Smith theorem: two arbitrary nodes, two Hasse jets.
    two_node_controls = []
    for distance in range(1, 13):
        matrix = arbitrary_node_hasse_matrix((0, distance), 2)
        divisor = distance * gcd(distance, 2)
        expected = (1, 1, divisor, distance ** 4 // divisor)
        actual = smith_diagonal(matrix)
        gate(actual == expected, ("two-node k=2 formula", distance))
        two_node_controls.append((distance, actual))

    # Hasse normalization is load-bearing beginning at order two.
    hasse = hasse_matrix(0, 2, 3)
    ordinary = ordinary_derivative_matrix(0, 2, 3)
    hasse_det = abs(bareiss_det(hasse))
    ordinary_det = abs(bareiss_det(ordinary))
    expected_ordinary_ratio = product(factorial(r) ** 3 for r in range(3))
    gate(ordinary_det == hasse_det * expected_ordinary_ratio,
         "ordinary derivative factorial mutation")
    gate(expected_ordinary_ratio == 8, "first Hasse normalization hostile")

    # Smith algorithm hostile: independent unimodular mutations preserve invariants.
    mutation_controls = []
    for seed, m, k in ((11, 2, 3), (29, 3, 2), (47, 4, 3)):
        matrix = hasse_matrix(0, m, k)
        expected = smith_diagonal(matrix)
        mutated = mutate_unimodular(matrix, seed)
        actual = smith_diagonal(mutated)
        gate(actual == expected, ("unimodular hostile", seed, m, k))
        mutation_controls.append((seed, m, k, actual))

    # Exact polynomial kernel, remainder, endpoint, and optimal-modulus controls.
    rng = random.Random(4000)
    sampler_controls = []
    for a, m, k in ((0, 0, 4), (0, 3, 2), (2, 2, 3), (7, 4, 2)):
        f = node_polynomial(a, m)
        kernel = poly_pow(f, k)
        rank = (m + 1) * k
        for trial in range(8):
            degree = rank + rng.randrange(1, 9)
            polynomial = [rng.randrange(-20, 21) for _ in range(degree + 1)]
            polynomial[-1] = rng.choice((-3, -2, -1, 1, 2, 3))
            quotient, remainder = monic_divmod(polynomial, kernel)
            gate(len(remainder) <= rank, "confluent remainder degree")
            gate(poly_add(poly_mul(kernel, quotient), remainder) == poly_trim(polynomial),
                 "confluent division identity")
            for node in range(a, a + m + 1):
                for order in range(k):
                    gate(hasse_eval(polynomial, node, order)
                         == hasse_eval(remainder, node, order),
                         "remainder preserves every Hasse jet")
            for evaluation_point in (a - 2, a + m + 2):
                modulus = abs(poly_eval(f, evaluation_point)) ** k
                difference = poly_eval(polynomial, evaluation_point) - poly_eval(remainder, evaluation_point)
                gate(modulus > 0 and difference % modulus == 0,
                     "off-node optimal congruence")
                gate(abs(poly_eval(kernel, evaluation_point)) == modulus,
                     "kernel generator attains optimal modulus")
        # At a sample node, order-zero data give exact equality; modulus zero is not used.
        node = a + m // 2
        gate(poly_eval(f, node) == 0, "sample endpoint has zero F value")
        gate(hasse_eval(kernel, node, 0) == 0, "sample value exact")
        # F itself is the minimal hostile to retaining only values when k>1.
        if k > 1:
            gate(all(poly_eval(f, node_) == 0 for node_ in range(a, a + m + 1)),
                 "ordinary values vanish on F")
            gate(any(hasse_eval(f, node_, 1) != 0 for node_ in range(a, a + m + 1)),
                 "first jets detect F")
            _q, f_remainder = monic_divmod(f, kernel)
            gate(f_remainder == f, "F is not divisible by F^k")
        sampler_controls.append((a, m, k, rank, tuple(kernel)))

    false_atlas = tuple((row[0], row[1]) for row in atlas if not row[5])
    gate(false_atlas and false_atlas[0] == (2, 3), "atlas-order first naive failure")
    # Rank-order minimum is (3,2), rank 8; (2,3) has rank 9.
    false_by_rank = min(((row[2], row[0], row[1]) for row in atlas if not row[5]))
    gate(false_by_rank == (8, 3, 2), "rank-minimal naive failure")

    semantic = {
        "atlas": [[row[0], row[1], row[2], str(row[3]), list(row[4]), row[5]]
                  for row in atlas],
        "exhaustive_crosschecks": [[m, k, list(divisors)]
                                   for m, k, divisors in exhaustive_crosschecks],
        "translation_controls": [[a, m, k, str(det), list(snf)]
                                 for a, m, k, det, snf in translation_controls],
        "primewise_controls": [list(row) for row in primewise_controls],
        "hostile_snf": list(hostile_snf),
        "hostile_naive": list(hostile_naive),
        "hostile_determinantal_divisors": list(hostile_divisors),
        "hostile_7minor": [list(witness_rows), list(witness_columns), witness_minor],
        "two_node_k2": [[distance, list(diagonal)]
                        for distance, diagonal in two_node_controls],
        "hasse_ordinary_determinant_ratio": expected_ordinary_ratio,
        "mutations": [[seed, m, k, list(diagonal)]
                      for seed, m, k, diagonal in mutation_controls],
        "sampler_controls": [[a, m, k, rank, list(kernel)]
                             for a, m, k, rank, kernel in sampler_controls],
        "first_false_by_rank": list(false_by_rank),
    }
    semantic_digest = sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":")
    ).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        gate(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("CONFLUENT_CONSECUTIVE_HASSE_OBSERVER_THM4010_INDEPENDENT_AUDIT")
    print("status=FINITE-EXACT_CANONICAL_PASS;KERNEL_AND_DETERMINANT_PROVED;GENERAL_SNF_OPEN")
    print("kernel=equal_jets_iff_difference_in_F_a_m_power_k_ZX")
    print("off_node_modulus=abs(F_a_m(B))^k_optimal;sample_node=value_exact_modulus_zero_not_used")
    print("determinant=product_1_le_j_le_m(j!)^(k^2);translation_independent=True")
    print("primewise=rank_mod_p:k*min(m+1,p);sum_smith_vp:k^2*sum_j_vp(j!)")
    print("unit_invariant_factors=N_if_m_le_1_else_exactly_2k")
    print(f"atlas_cases={len(atlas)};max_rank={max(row[2] for row in atlas)}")
    print(f"exhaustive_all_minor_crosschecks={len(exhaustive_crosschecks)};max_rank=9")
    for m, k, rank, determinant, diagonal, naive_match in atlas:
        print(f"atlas(m={m},k={k},N={rank})=det:{determinant};snf:{diagonal};naive:{naive_match}")
    print(f"rank_minimal_naive_hostile=(N,m,k)={false_by_rank}")
    print(f"hostile_actual_snf={hostile_snf}")
    print(f"hostile_naive_snf={hostile_naive}")
    print(f"hostile_determinantal_divisors={hostile_divisors}")
    print(f"hostile_7minor=rows{witness_rows};cols{witness_columns};det={witness_minor};not_divisible_by=576")
    print("partial_smith=two_nodes_distance_d_k2:diag(1,1,d*gcd(d,2),d^3/gcd(d,2))")
    print(f"two_node_controls={tuple(two_node_controls)}")
    print(f"hasse_vs_ordinary_first_ratio={expected_ordinary_ratio}")
    print(f"semantic_sha256={semantic_digest}")
    print(f"gates={GATES}")
    print("RESULT=PASS;GENERAL_CONSECUTIVE_SMITH_FORMULA=OPEN")


if __name__ == "__main__":
    main()
