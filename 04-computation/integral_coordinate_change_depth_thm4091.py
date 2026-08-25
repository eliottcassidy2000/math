#!/usr/bin/env python3
"""Exact finite controls for THM-4091.

The theorem itself is symbolic.  This companion exhausts integral coordinate
changes h(f)=h_1 f+...+h_4 f^4 with h_i in [-2,2], through output degree six.
It checks the derivative identity responsible for the exceptional depth-one
law, the finite LCM clearing matrix, the absence of a failure before degree
three, and the sharp hostile h=f+f^2 for every tested depth e>=2.

Reproduction:
  python -B 04-computation/integral_coordinate_change_depth_thm4091.py
  python -B -O 04-computation/integral_coordinate_change_depth_thm4091.py
"""

from fractions import Fraction
from itertools import product
from math import gcd


MAX_DEGREE = 6


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def multiply(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (MAX_DEGREE + 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i + j <= MAX_DEGREE:
                out[i + j] += a * b
    return tuple(out)


def powers(h: tuple[int, ...]) -> list[tuple[int, ...]]:
    result = [(1,) + (0,) * MAX_DEGREE]
    for _ in range(MAX_DEGREE):
        result.append(multiply(result[-1], h))
    return result


def lcm_up_to(n: int) -> int:
    value = 1
    for k in range(1, n + 1):
        value = value * k // gcd(value, k)
    return value


def main() -> None:
    coordinate_changes = 0
    derivative_gates = 0
    lcm_matrix_gates = 0
    low_degree_gates = 0
    depth_two_failures = 0
    first_depth_two_failure: tuple[int, int, tuple[int, ...], Fraction] | None = None

    for coefficients in product(range(-2, 3), repeat=4):
        h = (0,) + coefficients + (0,) * (MAX_DEGREE - 4)
        hp = tuple((i + 1) * h[i + 1] for i in range(MAX_DEGREE)) + (0,)
        hpowers = powers(h)
        coordinate_changes += 1

        for n in range(1, MAX_DEGREE + 1):
            for k in range(1, n + 1):
                coefficient = hpowers[k][n]

                # n[f^n]h^k = k[f^(n-1)]h^(k-1)h'.
                rhs_series = multiply(hpowers[k - 1], hp)
                require(
                    n * coefficient == k * rhs_series[n - 1],
                    "derivative identity failed",
                )
                require(
                    Fraction(n * coefficient, k).denominator == 1,
                    "depth-one composition entry is not integral",
                )
                derivative_gates += 1

                for exponent in range(1, 5):
                    clearing = lcm_up_to(MAX_DEGREE) ** exponent
                    entry = Fraction(clearing * coefficient, k**exponent)
                    require(entry.denominator == 1, "LCM matrix clearing failed")
                    lcm_matrix_gates += 1

                    if n <= 2:
                        scaled = Fraction(n**exponent * coefficient, k**exponent)
                        require(scaled.denominator == 1, "failure below degree three")
                        low_degree_gates += 1

                scaled_depth_two = Fraction(n * n * coefficient, k * k)
                if scaled_depth_two.denominator != 1:
                    depth_two_failures += 1
                    candidate = (n, k, coefficients, scaled_depth_two)
                    if first_depth_two_failure is None or candidate[:2] < first_depth_two_failure[:2]:
                        first_depth_two_failure = candidate

    require(first_depth_two_failure is not None, "depth-two hostile was not found")
    require(first_depth_two_failure[0:2] == (3, 2), "first possible failure is not (n,k)=(3,2)")

    # Sharp universal hostile from the proof: R_e=q^2/2^e and h=f+f^2.
    hostile_h = (0, 1, 1) + (0,) * (MAX_DEGREE - 2)
    hostile_powers = powers(hostile_h)
    hostile_rows: list[tuple[int, Fraction]] = []
    for exponent in range(2, 9):
        coefficient = Fraction(hostile_powers[2][3], 2**exponent)
        scaled = 3**exponent * coefficient
        require(scaled.denominator == 2 ** (exponent - 1), "unexpected hostile denominator")
        hostile_rows.append((exponent, scaled))

    print("THM-4091 INTEGRAL COORDINATE-CHANGE AUDIT")
    print("h(f)=h1*f+...+h4*f^4, each hi in [-2,2]; output degree <= 6")
    print()
    print(f"coordinate changes exhausted = {coordinate_changes}")
    print(f"derivative/depth-one gates   = {derivative_gates}")
    print(f"LCM matrix gates             = {lcm_matrix_gates}")
    print(f"degree-one/two gates         = {low_degree_gates}")
    print(f"depth-two failing entries    = {depth_two_failures}")
    print("first possible failure       = (n,k)=(3,2)")
    print()
    print("hostile R_e(q)=q^2/2^e, h(f)=f+f^2")
    for exponent, scaled in hostile_rows:
        print(f"e={exponent}: 3^e*[f^3]R_e(h) = {scaled}")
    print("EXACT CONTROLS: PASS")


if __name__ == "__main__":
    main()
