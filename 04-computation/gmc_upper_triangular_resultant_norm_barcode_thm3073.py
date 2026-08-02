#!/usr/bin/env python3
"""Exact arithmetic companion for THM-3073's norm/phase barcode."""

from itertools import permutations, product
from math import factorial, gcd, prod


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def degree_tuples():
    return [
        values
        for length in range(1, 4)
        for values in product(range(1, 5), repeat=length)
    ]


def exact_order(modulus, exponent):
    return modulus // gcd(modulus, exponent)


def first_stage(modulus, exponent_at_stage):
    for stage in range(3, 1000):
        if exponent_at_stage(stage) % modulus == 0:
            return stage
    raise RuntimeError("first-stage search range exhausted")


def death_duration(modulus, width):
    running = 1
    for duration in range(1, 1000):
        running *= width + duration
        if running % modulus == 0:
            return duration
    raise RuntimeError("death-duration search range exhausted")


def poly_add(left, right, scale=1):
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, 0) + scale * coefficient
        if answer[monomial] == 0:
            del answer[monomial]
    return answer


def poly_mul(left, right):
    answer = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            monomial = (i + k, j + ell)
            answer[monomial] = answer.get(monomial, 0) + a * b
    return {monomial: value for monomial, value in answer.items() if value}


ZERO = {}
ONE = {(0, 0): 1}
C = {(1, 0): 1}
DVAR = {(0, 1): 1}


def permutation_sign(values):
    inversions = sum(
        values[i] > values[j]
        for i in range(len(values))
        for j in range(i + 1, len(values))
    )
    return -1 if inversions % 2 else 1


def polynomial_determinant(matrix):
    answer = ZERO
    for sigma in permutations(range(len(matrix))):
        term = ONE
        for row, column in enumerate(sigma):
            term = poly_mul(term, matrix[row][column])
        answer = poly_add(answer, term, permutation_sign(sigma))
    return answer


# The proposed factor S^N E^D has exactly the full resultant multidegree in
# every lower and upper coefficient block.
degree_cells = 0
degree_entries = 0
for lower_degrees in degree_tuples():
    for upper_degrees in degree_tuples():
        lower_product = prod(lower_degrees)
        upper_product = prod(upper_degrees)
        full_product = lower_product * upper_product
        for degree in lower_degrees:
            require(
                upper_product * (lower_product // degree)
                == full_product // degree,
                "lower multidegree mismatch",
            )
            degree_entries += 1
        for degree in upper_degrees:
            require(
                lower_product * (upper_product // degree)
                == full_product // degree,
                "upper multidegree mismatch",
            )
            degree_entries += 1
        degree_cells += 1


# The phase-torus map (u,v) -> u^N v^D has image size q/g and kernel qg,
# where g=gcd(q,N,D).  Check both GMC exponent families exactly.
phase_cells = 0
for modulus in (2, 3, 4, 5, 7, 13, 14, 26, 91):
    for stage in range(3, 21):
        exponent_pairs = (
            (stage, factorial(stage - 1)),
            (stage * (stage - 1), factorial(stage - 2)),
        )
        for lower_exponent, normal_exponent in exponent_pairs:
            image = {
                (lower_exponent * u + normal_exponent * v) % modulus
                for u in range(modulus)
                for v in range(modulus)
            }
            kernel = sum(
                (lower_exponent * u + normal_exponent * v) % modulus == 0
                for u in range(modulus)
                for v in range(modulus)
            )
            divisor = gcd(modulus, lower_exponent, normal_exponent)
            require(len(image) == modulus // divisor, "phase image mismatch")
            require(kernel == modulus * divisor, "phase kernel mismatch")
            phase_cells += 1


# Direct first-death stages for a phase placed in either input block.
direct_rows = {}
for modulus in (7, 13, 91):
    direct_rows[modulus] = (
        first_stage(modulus, lambda k: k),
        first_stage(modulus, lambda k: factorial(k - 1)),
        first_stage(modulus, lambda k: k * (k - 1)),
        first_stage(modulus, lambda k: factorial(k - 2)),
    )
require(direct_rows[7] == (7, 8, 7, 9), "C7 direct barcode changed")
require(direct_rows[13] == (13, 14, 13, 15), "C13 direct barcode changed")
require(direct_rows[91] == (91, 14, 14, 15), "C91 direct barcode changed")


# The iterated one-normal tower multiplies an old phase by k!/m!.
tower_cells = 0
for width in range(2, 21):
    for modulus in (2, 7, 13, 91):
        duration = death_duration(modulus, width)
        exponent = prod(range(width + 1, width + duration + 1))
        require(exponent % modulus == 0, "tower phase did not die")
        if duration > 1:
            previous = prod(range(width + 1, width + duration))
            require(previous % modulus != 0, "tower death was not first")
        if modulus in (2, 7, 13):
            require(
                duration == modulus - (width % modulus),
                "prime death-duration formula changed",
            )
        if modulus == 91:
            require(
                duration
                == max(7 - (width % 7), 13 - (width % 13)),
                "C91 death-duration formula changed",
            )
        tower_cells += 1


# Phase placement and cancellation hostiles.
require(exact_order(7, 8 * 7) == 1, "THM-3063 k=8 base should die")
require(exact_order(7, factorial(6)) == 7, "THM-3063 k=8 normal should survive")
require(exact_order(7, 7) == 1, "THM-3069 k=7 base should die")
require(exact_order(7, factorial(6)) == 7, "THM-3069 k=7 normal should survive")
require(exact_order(7, 8) == 7, "THM-3069 k=8 base should survive")
require(exact_order(7, factorial(7)) == 1, "THM-3069 k=8 normal should die")
require((4 * 3 + factorial(2)) % 7 == 0, "two-input cancellation hostile changed")
require(4 * 3 % 7 != 0 and factorial(2) % 7 != 0, "hostile factors died separately")


# If the lower equations are allowed to mix z, the nominal factorization is
# false.  The Sylvester determinant of y(y+cz), z(z+dy) is exactly 1-cd.
mixed_matrix = (
    (ONE, C, ZERO, ZERO),
    (ZERO, ONE, C, ZERO),
    (ZERO, DVAR, ONE, ZERO),
    (ZERO, ZERO, DVAR, ONE),
)
mixed_resultant = polynomial_determinant(mixed_matrix)
require(mixed_resultant == {(0, 0): 1, (1, 1): -1}, "mixed hostile changed")


print("THM-3073 UPPER-TRIANGULAR RESULTANT NORM BARCODE")
print(
    f"degree_cells={degree_cells} multidegree_entries={degree_entries} "
    "coordinate_power_sign=+1"
)
print(
    f"phase_kernel_cells={phase_cells} image=q/gcd(q,N,D) "
    "kernel=q*gcd(q,N,D)"
)
print("one_normal_phase=(u,v)->u^k*v^((k-1)!)")
print("two_normal_phase=(u,v)->u^(k(k-1))*v^((k-2)!)")
print("direct_death_C7=one:7/8,two:7/9")
print("direct_death_C13=one:13/14,two:13/15")
print("direct_death_C91=one:91/14,two:14/15")
print("tower_from_m2=durations:5,11,11;terminal_widths:7,13,13")
print("tower_from_m3=durations:4,10,10;terminal_widths:7,13,13")
print(f"tower_cells={tower_cells} order_after=N/gcd(N,k!/m!)")
print("placement_hostile=k8_two_normal:C7_base_dies,C7_normal_survives")
print("cancellation_hostile=k4_two_normal:C7_u=v_nontrivial,scalar_trivial")
print("mixed_lower_hostile=resultant(y(y+cz),z(z+dy))=1-cd")
print("scope=upper_triangular;fixed_composable_tower;pullback_is_not_norm")
print("all_exact_checks=PASS")
