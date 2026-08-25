#!/usr/bin/env python3
"""Exact arithmetic audit for the mixed-binomial second-moment exponent.

This is the exact verifier for THM-4036. It checks the signed
falling-factorial coefficients, evaluates the energy-peeling recurrence over
every order of (2,4,6,8), and records the exact leastness gap and target CRT
coordinates.  It is not a search for counterexamples.
"""

from fractions import Fraction
from itertools import permutations


DEGREES = (2, 4, 6, 8)
TARGET = 896_315_812_331_399
PUBLISHED_PREFIX = 2_000_000_000_000


def require(condition: bool, message: object) -> None:
    if not condition:
        raise AssertionError(message)


def multiply_by_linear(coefficients: list[int], root: int) -> list[int]:
    """Ascending coefficients of p(t) * (t-root)."""
    out = [0] * (len(coefficients) + 1)
    for exponent, coefficient in enumerate(coefficients):
        out[exponent] -= root * coefficient
        out[exponent + 1] += coefficient
    return out


def falling_coefficients(degree: int) -> tuple[int, ...]:
    coefficients = [1]
    for root in range(degree):
        coefficients = multiply_by_linear(coefficients, root)
    return tuple(coefficients)


EXPECTED = {
    2: (0, -1, 1),
    4: (0, -6, 11, -6, 1),
    6: (0, -120, 274, -225, 85, -15, 1),
    8: (0, -5040, 13068, -13132, 6769, -1960, 322, -28, 1),
}


def reciprocal_sum(degrees: tuple[int, ...]) -> Fraction:
    return sum((Fraction(1, degree) for degree in degrees), Fraction(0))


def energy_exponent(peeling_order: tuple[int, ...]) -> Fraction:
    """Exponent from E_(d,D) <= X^(1/d) E_D + X^eps B_D^2."""
    if len(peeling_order) == 1:
        return Fraction(1, peeling_order[0])
    degree = peeling_order[0]
    tail = peeling_order[1:]
    diagonal = Fraction(1, degree) + energy_exponent(tail)
    off_diagonal = 2 * reciprocal_sum(tail)
    return max(diagonal, off_diagonal)


def main() -> None:
    for degree in DEGREES:
        actual = falling_coefficients(degree)
        require(actual == EXPECTED[degree], (degree, actual))
        print(f"falling_coefficients_{degree}={actual}")

    best = None
    best_orders: list[tuple[int, ...]] = []
    for order in permutations(DEGREES):
        exponent = energy_exponent(order)
        print(f"order={order} energy_exponent={exponent}")
        if best is None or exponent < best:
            best = exponent
            best_orders = [order]
        elif exponent == best:
            best_orders.append(order)

    require(best == Fraction(13, 12), "wrong best energy exponent")
    require((2, 4, 6, 8) in best_orders, "canonical peeling order is not optimal")
    print(f"best_energy_exponent={best} best_orders={best_orders}")

    summatory_exponent = reciprocal_sum(DEGREES)
    support_exponent = 2 * summatory_exponent - best
    require(summatory_exponent == Fraction(25, 24), "wrong summatory exponent")
    require(support_exponent == 1, "wrong support exponent")
    print(f"summatory_exponent={summatory_exponent}")
    print(f"cauchy_schwarz_support_exponent={support_exponent}")

    gap_start = PUBLISHED_PREFIX + 1
    gap_end = TARGET - 1
    gap_count = gap_end - gap_start + 1
    require(gap_count == 894_315_812_331_398, "wrong unresolved leastness gap")
    print(f"leastness_gap=[{gap_start},{gap_end}] count={gap_count}")

    moduli = (3, 9, 11, 33, 99, 5, 7, 13, 17, 19, 23)
    print("target_residues=" + ",".join(f"{q}:{TARGET % q}" for q in moduli))
    print("PASS")


if __name__ == "__main__":
    main()
