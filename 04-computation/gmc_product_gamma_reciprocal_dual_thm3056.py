#!/usr/bin/env python3
"""Exact companion for THM-3056 (product-Gamma reciprocal duality).

Proof-level asymptotics are recorded, not inferred from floats.  Exact
rational arithmetic independently checks the hypergeometric coefficients,
the negative order-two Hankel wall, the single-factor determinant formula,
the integer width specialization, and its prime-incidence formula.  The
multi-factor all-order sign experiment is printed explicitly as a finite scout
and is not used as a theorem.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rising(a: Fraction, n: int) -> Fraction:
    value = Fraction(1)
    for j in range(n):
        value *= a + j
    return value


def reciprocal_weight(c: Fraction, shapes: tuple[Fraction, ...], n: int) -> Fraction:
    value = c**n
    for alpha in shapes:
        value *= rising(alpha, n)
    return 1 / value


def determinant_elimination(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    size = len(work)
    sign = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        for row in range(column + 1, size):
            multiplier = work[row][column] / pivot_value
            for target in range(column, size):
                work[row][target] -= multiplier * work[column][target]
    value = Fraction(sign)
    for diagonal in range(size):
        value *= work[diagonal][diagonal]
    return value


def permutation_sign(perm: tuple[int, ...]) -> int:
    inversions = sum(
        1
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
        if perm[i] > perm[j]
    )
    return -1 if inversions % 2 else 1


def determinant_leibniz(matrix: list[list[Fraction]]) -> Fraction:
    size = len(matrix)
    value = Fraction(0)
    for perm in permutations(range(size)):
        term = Fraction(permutation_sign(perm))
        for row, column in enumerate(perm):
            term *= matrix[row][column]
        value += term
    return value


def sign(value: Fraction) -> int:
    return (value > 0) - (value < 0)


def product(values) -> Fraction:
    answer = Fraction(1)
    for value in values:
        answer *= value
    return answer


def character(k: int) -> tuple[int, int, int]:
    harmonic = sum((Fraction(1, j) for j in range(1, k + 1)), Fraction(0))
    scale = factorial(k)
    A = scale * (harmonic - 1)
    B = scale * (k + 1 - 2 * harmonic)
    require(A.denominator == B.denominator == 1, f"nonintegral character at k={k}")
    return int(A), int(B), int(A + B)


def width_flag(k: int, t: int, n: int) -> int:
    A, B, I = character(k)
    value = 1
    for s in range(1, n):
        value *= (1 + s * t) ** I
    value *= (1 + n * t) ** B
    return value


hypergeometric_cells = 0
hankel_two_cells = 0
shape_families = (
    (Fraction(1),),
    (Fraction(1, 2),),
    (Fraction(3, 2),),
    (Fraction(1), Fraction(1)),
    (Fraction(1, 2), Fraction(3, 2)),
    (Fraction(1, 3), Fraction(2, 3), Fraction(4, 3)),
)
for c in (Fraction(1, 3), Fraction(1), Fraction(5, 2)):
    for shapes in shape_families:
        values = tuple(reciprocal_weight(c, shapes, n) for n in range(10))
        for n in range(9):
            predicted_ratio = 1 / (c * product(n + alpha for alpha in shapes))
            require(values[n + 1] / values[n] == predicted_ratio, "hypergeometric ratio mismatch")
            hypergeometric_cells += 1
        for n in range(1, 9):
            minor = values[n - 1] * values[n + 1] - values[n] ** 2
            require(minor < 0, "reciprocal order-two Hankel wall changed")
            hankel_two_cells += 1
# Closed single-factor contiguous determinant formula, plus a generalized
# sign audit.  The formula itself is proved by a falling-factorial alternant.
single_formula_cells = 0
single_generalized_cells = 0
for c in (Fraction(1, 2), Fraction(1), Fraction(7, 3)):
    for alpha in (Fraction(1, 3), Fraction(1), Fraction(5, 2)):
        for offset in range(5):
            for size in range(1, 7):
                matrix = [
                    [reciprocal_weight(c, (alpha,), offset + i + j) for j in range(size)]
                    for i in range(size)
                ]
                actual = determinant_elimination(matrix)
                exponent = size * offset + size * (size - 1)
                numerator = product(Fraction(factorial(h)) for h in range(1, size))
                denominator = product(
                    rising(alpha, offset + i + size - 1) for i in range(size)
                )
                expected = Fraction((-1) ** (size * (size - 1) // 2))
                expected *= c ** (-exponent) * numerator / denominator
                require(actual == expected, "single-factor closed determinant mismatch")
                single_formula_cells += 1

        indices = range(7)
        values = tuple(reciprocal_weight(c, (alpha,), n) for n in range(13))
        for size in range(1, 5):
            target = (-1) ** (size * (size - 1) // 2)
            for rows in combinations(indices, size):
                for columns in combinations(indices, size):
                    matrix = [[values[row + column] for column in columns] for row in rows]
                    actual = determinant_elimination(matrix)
                    if size <= 4:
                        require(actual == determinant_leibniz(matrix), "determinant paths disagree")
                    require(sign(actual) == target, "single-factor generalized sign mismatch")
                    single_generalized_cells += 1


# THM-3047 integer-width specialization.
width_cells = 0
collision_rows = 0
for k in range(2, 7):
    A, B, I = character(k)
    for t in (1, 2, 3):
        a = Fraction(1, t)
        shapes = (a,) * A + (a + 1,) * B
        values = tuple(width_flag(k, t, n) for n in range(9))
        reciprocal_values = tuple(reciprocal_weight(Fraction(t**I), shapes, n) for n in range(9))
        require(tuple(Fraction(1, value) for value in values) == reciprocal_values, "width dual mismatch")
        repeated = sum(1 for n in range(1, len(values)) if values[n] == values[n - 1])
        expected_repeated = int(k == 2)
        require(repeated == expected_repeated, f"width collision mismatch k={k}, t={t}")
        collision_rows += repeated
        require(all(values[n + 1] > values[n] for n in range(1, 8)), "width tail not increasing")
        width_cells += len(values)


def primes_through(limit: int) -> tuple[int, ...]:
    return tuple(
        n
        for n in range(2, limit + 1)
        if all(n % d for d in range(2, int(n**0.5) + 1))
    )


def valuation(value: int, prime: int) -> int:
    count = 0
    while value % prime == 0:
        value //= prime
        count += 1
    return count


def first_residue(t: int, modulus: int) -> int:
    for residue in range(1, modulus):
        if (1 + residue * t) % modulus == 0:
            return residue
    raise RuntimeError(f"no residue for t={t}, modulus={modulus}")


prime_cells = 0
for k in range(2, 6):
    A, B, I = character(k)
    for t in (1, 2, 3, 4, 6):
        for n in range(0, 18):
            value = width_flag(k, t, n)
            for prime in primes_through(43):
                if t % prime == 0:
                    expected = 0
                else:
                    expected = B * valuation(1 + n * t, prime)
                    power = prime
                    while power <= 1 + max(0, n - 1) * t:
                        residue = first_residue(t, power)
                        if n - 1 >= residue:
                            expected += I * (1 + (n - 1 - residue) // power)
                        power *= prime
                require(valuation(value, prime) == expected, "prime valuation mismatch")
                prime_cells += 1


# Finite-exact multi-factor generalized-Hankel scout.  No all-order theorem is
# inferred: generic sign-regular Hadamard closure is unavailable.
multi_scout_cells = 0
multi_shape_families = shape_families[3:] + (
    (Fraction(1),) * 5 + (Fraction(2),) * 2,
    (Fraction(1, 2),) * 5 + (Fraction(3, 2),) * 2,
)
for shapes in multi_shape_families:
    values = tuple(reciprocal_weight(Fraction(1), shapes, n) for n in range(13))
    for size in range(2, 5):
        target = (-1) ** (size * (size - 1) // 2)
        for rows in combinations(range(6), size):
            for columns in combinations(range(6), size):
                matrix = [[values[row + column] for column in columns] for row in rows]
                actual = determinant_elimination(matrix)
                require(actual == determinant_leibniz(matrix), "multi-factor determinant paths disagree")
                require(sign(actual) == target, "multi-factor finite scout sign mismatch")
                multi_scout_cells += 1


print("theorem=reciprocal product-Gamma series is _1F_J with order 1/J and type J*c^(-1/J)")
print(f"hypergeometric_ratio_cells={hypergeometric_cells};negative_H2_cells={hankel_two_cells}")
print(
    f"single_gamma_closed_formula_cells={single_formula_cells};"
    f"generalized_SSR_cells={single_generalized_cells};two determinant paths"
)
print(f"integer_width_cells={width_cells};duplicate_one_rows={collision_rows};only k=2")
print(f"prime_valuation_cells={prime_cells};prime shadow=all primes not dividing t")
print("support=Dirichlet abscissa 0;counting ~ log(x)/(I loglog(x));harmonic mass hypergeometric")
print("divisor_scar=THM2438 mean loss equals support harmonic mass")
print(f"multi_factor_sign_cells={multi_scout_cells};status=FINITE-EXACT_SCOUT_NOT_THEOREM")
