#!/usr/bin/env python3
"""Exact companion for THM-2412.

Checks the derivative/finite-difference basis dictionary, the exponential
eigenfunction laws, the Pascal half-row split, and the Catalan convolution
ladder using integers and Fractions only.
"""

from fractions import Fraction as F
from math import comb, factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def falling(x: F, k: int) -> F:
    value = F(1)
    for j in range(k):
        value *= x - j
    return value


def scaled_falling(x: F, k: int, h: F) -> F:
    value = F(1)
    for j in range(k):
        value *= x - j * h
    return value


def poly_value(coefficients: tuple[F, ...], x: F) -> F:
    value = F(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def derivative_coefficients(coefficients: tuple[F, ...]) -> tuple[F, ...]:
    if len(coefficients) <= 1:
        return (F(0),)
    return tuple(j * coefficients[j] for j in range(1, len(coefficients)))


def delta_value(coefficients: tuple[F, ...], x: F, order: int = 1) -> F:
    return sum(
        (-1) ** (order - j) * comb(order, j) * poly_value(coefficients, x + j)
        for j in range(order + 1)
    )


def stirling_second(size: int) -> list[list[int]]:
    table = [[0] * (size + 1) for _ in range(size + 1)]
    table[0][0] = 1
    for n in range(1, size + 1):
        for k in range(1, n + 1):
            table[n][k] = table[n - 1][k - 1] + k * table[n - 1][k]
    return table


def stirling_first_signed(size: int) -> list[list[int]]:
    table = [[0] * (size + 1) for _ in range(size + 1)]
    table[0][0] = 1
    for n in range(1, size + 1):
        for k in range(1, n + 1):
            table[n][k] = table[n - 1][k - 1] - (n - 1) * table[n - 1][k]
    return table


def convolution(left: list[int], right: list[int], n: int) -> int:
    return sum(left[k] * right[n - k] for k in range(n + 1))


# ---------------------------------------------------------------------------
# 1. The lowering bases and the Stirling change of coordinates.
# ---------------------------------------------------------------------------

degree = 12
S2 = stirling_second(degree)
s1 = stirling_first_signed(degree)
lowering_checks = 0

for k in range(1, degree + 1):
    for numerator in range(-7, 12):
        x = F(numerator, 3)
        require(
            falling(x + 1, k) - falling(x, k) == k * falling(x, k - 1),
            "Delta did not lower the falling-factorial basis exactly",
        )
        lowering_checks += 1

scaled_lowering_checks = 0
for h in (F(1), F(2, 3), F(-3, 5)):
    for k in range(1, degree + 1):
        for numerator in range(-5, 8):
            x = F(numerator, 4)
            require(
                (
                    scaled_falling(x + h, k, h)
                    - scaled_falling(x, k, h)
                )
                / h
                == k * scaled_falling(x, k - 1, h),
                "scaled difference quotient did not lower its basic sequence",
            )
            scaled_lowering_checks += 1

for j in range(degree + 1):
    for numerator in range(-7, 12):
        x = F(numerator, 3)
        require(
            x**j == sum(S2[j][k] * falling(x, k) for k in range(j + 1)),
            "second-kind Stirling basis change failed",
        )
        require(
            falling(x, j) == sum(s1[j][k] * x**k for k in range(j + 1)),
            "first-kind Stirling basis change failed",
        )

# A deterministic rational polynomial supplies a nontrivial all-coordinate
# Maclaurin/Gregory--Newton comparison.
coefficients = tuple(F((-1) ** j * (j + 2), j + 1) for j in range(9))
newton_coefficients = tuple(
    delta_value(coefficients, F(0), k) / factorial(k)
    for k in range(len(coefficients))
)

for k, newton_coefficient in enumerate(newton_coefficients):
    transformed = sum(
        coefficients[j] * S2[j][k] for j in range(k, len(coefficients))
    )
    require(newton_coefficient == transformed, "Maclaurin-to-Newton transform")

for j, maclaurin_coefficient in enumerate(coefficients):
    transformed = sum(
        newton_coefficients[k] * s1[k][j]
        for k in range(j, len(coefficients))
    )
    require(maclaurin_coefficient == transformed, "Newton-to-Maclaurin transform")

for numerator in range(-11, 18):
    x = F(numerator, 4)
    reconstructed = sum(
        newton_coefficients[k] * falling(x, k)
        for k in range(len(coefficients))
    )
    require(reconstructed == poly_value(coefficients, x), "Newton reconstruction")

# E=exp(D), Delta=E-I, and D=log(I+Delta), all finite on polynomials.
operator_checks = 0
for numerator in range(-9, 13):
    x = F(numerator, 5)
    derivative = derivative_coefficients(coefficients)
    exp_D = F(0)
    current = coefficients
    for order in range(len(coefficients)):
        exp_D += poly_value(current, x) / factorial(order)
        current = derivative_coefficients(current)
    require(exp_D == poly_value(coefficients, x + 1), "E != exp(D)")

    log_shift = sum(
        F((-1) ** (order + 1), order)
        * delta_value(coefficients, x, order)
        for order in range(1, len(coefficients))
    )
    require(log_shift == poly_value(derivative, x), "D != log(I+Delta)")
    operator_checks += 1


# ---------------------------------------------------------------------------
# 2. Continuous and discrete exponentials.
# ---------------------------------------------------------------------------

exponential_checks = 0
for base in range(-3, 8):
    for n in range(21):
        require(
            base**n
            == sum(comb(n, k) * (base - 1) ** k for k in range(n + 1)),
            "finite Newton exponential failed",
        )
        exponential_checks += 1

for n in range(31):
    require(
        2**n == sum(falling(F(n), k) / factorial(k) for k in range(n + 1)),
        "unit Delta exponential failed",
    )
    require(
        4**n == sum(3**k * comb(n, k) for k in range(n + 1)),
        "4^n Newton expansion failed",
    )
    if n < 30:
        require(2 ** (n + 1) - 2**n == 2**n, "Delta 2^n eigenvalue")
        require(4 ** (n + 1) - 4**n == 3 * 4**n, "Delta 4^n eigenvalue")

# A labelled tournament has one binary choice for every unordered pair.
tournament_checks = 0
for vertices in range(13):
    pairs = vertices * (vertices - 1) // 2
    require(sum(comb(pairs, k) for k in range(pairs + 1)) == 2**pairs,
            "tournament binary-coordinate count failed")
    tournament_checks += 1

switching_gauge_checks = 0
for vertices in range(1, 8):
    signs = {
        (i, j): 1 if (i * 5 + j * 3 + vertices) % 2 else -1
        for i in range(vertices)
        for j in range(i + 1, vertices)
    }
    for word in range(1 << vertices):
        x = tuple(1 if word & (1 << i) else -1 for i in range(vertices))
        energy = sum(signs[i, j] * x[i] * x[j] for i, j in signs)
        negated = tuple(-entry for entry in x)
        negated_energy = sum(
            signs[i, j] * negated[i] * negated[j] for i, j in signs
        )
        require(energy == negated_energy, "global switching gauge failed")
        switching_gauge_checks += 1


# ---------------------------------------------------------------------------
# 3. Complementary Pascal half-rows and the Catalan ladder.
# ---------------------------------------------------------------------------

limit = 50
A = [
    sum(comb(2 * n, k) for k in range(n + 1))
    for n in range(limit + 2)
]
B = [
    sum(comb(2 * n + 2, k) for k in range(n + 1))
    for n in range(limit + 1)
]
P = [4**n for n in range(limit + 2)]
Catalan = [comb(2 * n, n) // (n + 1) for n in range(limit + 2)]

require(A[:6] == [1, 3, 11, 42, 163, 638], "A032443 prefix mismatch")
require(B[:5] == [1, 5, 22, 93, 386], "A000346 prefix mismatch")
require(P[:5] == [1, 4, 16, 64, 256], "power-of-four prefix mismatch")

split_checks = 0
for n in range(limit + 1):
    central = comb(2 * n + 2, n + 1)
    require(A[n + 1] + B[n] == 4 ** (n + 1), "half-row sum failed")
    require(A[n + 1] - B[n] == central, "central tie-layer difference failed")
    require(A[n + 1] == (4 ** (n + 1) + central) // 2, "weak half")
    require(B[n] == (4 ** (n + 1) - central) // 2, "strict half")
    split_checks += 1

catalan_checks = 0
for n in range(limit + 1):
    require(convolution(Catalan, A, n) == P[n], "Catalan*A != 4^n")
    require(convolution(Catalan, P, n) == B[n], "Catalan*4^n != B")
    catalan_squared_A = sum(
        Catalan[i] * Catalan[j] * A[n - i - j]
        for i in range(n + 1)
        for j in range(n - i + 1)
    )
    require(catalan_squared_A == B[n], "Catalan^2*A != B")
    catalan_checks += 1

catalan_leakage_checks = 0
strict_half = [0] + B[:limit]
for n in range(1, limit + 1):
    cat = Catalan[n - 1]
    require(A[n] == 4 * A[n - 1] - cat, "weak-half Catalan leakage")
    require(
        strict_half[n] == 4 * strict_half[n - 1] + cat,
        "strict-half Catalan leakage",
    )
    catalan_leakage_checks += 1

# Hostile boundaries: ordinary powers are not the exact Delta-lowering basis,
# bases other than 2 have a different Delta eigenvalue, and deleting the
# central layer changes the half-row.
require((F(3) + 1) ** 2 - F(3) ** 2 == 7, "monomial hostile setup")
require(7 != 2 * F(3), "monomial basis accidentally lowered exactly")
require(3**8 - 3**7 == 2 * 3**7, "base-three hostile setup")
require(3**8 - 3**7 != 3**7, "base 2 was not unique at eigenvalue one")
require(A[5] == 638 and A[5] != 639, "sixth-term prefix guard failed")


print("theorem=THM-2412")
print("status=PROVED+VERIFIED-EXACT")
print(f"falling_factorial_lowering_checks={lowering_checks}")
print(f"scaled_falling_lowering_checks={scaled_lowering_checks}")
print(f"operator_dictionary_checks={operator_checks}")
print(f"finite_exponential_checks={exponential_checks}")
print("delta_unit=2^n;delta_eigenvalue=1")
print("four_power_newton=4^n=sum_k(3^k binom(n,k))")
print(f"tournament_binary_coordinate_checks={tournament_checks}")
print(f"switching_global_gauge_checks={switching_gauge_checks}")
print("A032443_prefix=1,3,11,42,163,638")
print("A000346_prefix=1,5,22,93,386")
print(f"central_half_split_checks={split_checks}")
print("half_sum=4^(n+1);half_difference=binom(2n+2,n+1)")
print(f"catalan_convolution_ladder_checks={catalan_checks}")
print(f"catalan_leakage_recurrence_checks={catalan_leakage_checks}")
print("catalan*A032443=4^n;catalan*4^n=A000346;catalan^2*A032443=A000346")
print("operation_collision=A_(n+1)+B_n=4^(n+1)=4*4^n")
print("hostiles=monomial-lowering;base-not-two;central-tie-deletion;sixth-term-639")
print("all_checks=PASS")
