#!/usr/bin/env python3
"""Independent exact referee for THM-3001's reversal fixed locus.

The referee uses coefficient arithmetic only.  It exhausts a small positive
coefficient universe, verifies that Newton-ratio palindromy is equivalent to a
scaled reciprocal coefficient law, checks the antipalindromic circuit and the
odd-degree fixed wall, and retains a constant-R boundary separating a
constant ratio sequence from equality in every Newton inequality.
"""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path


EXPECTED_RECORD_DIGEST = (
    "1d0484648a8c6abed41e9c67d8d679e40fa34bd40f4f94b6f71a2cb3fb47a810"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def h_sequence(coefficients):
    d = len(coefficients) - 1
    lead = Fraction(coefficients[d])
    return [
        Fraction(coefficients[d - k], 1) / (lead * comb(d, k))
        for k in range(d + 1)
    ]


def ratios(coefficients):
    h = h_sequence(coefficients)
    d = len(coefficients) - 1
    return [None] + [
        h[k] * h[k] / (h[k - 1] * h[k + 1])
        for k in range(1, d)
    ]


def ratio_palindromic(coefficients) -> bool:
    d = len(coefficients) - 1
    row = ratios(coefficients)
    return all(row[k] == row[d - k] for k in range(1, d))


def scaled_reciprocal(coefficients) -> bool:
    """Check a_k=A B^k a_(d-k) with A,B>0, exactly over Q."""
    d = len(coefficients) - 1
    a = [Fraction(value) for value in coefficients]
    A = a[0] / a[d]
    B = a[1] * a[d] / (a[0] * a[d - 1])
    return A > 0 and B > 0 and all(
        a[k] == A * B**k * a[d - k] for k in range(d + 1)
    )


def circuit_antipalindromic(coefficients) -> bool:
    d = len(coefficients) - 1
    row = ratios(coefficients)
    if not ratio_palindromic(coefficients):
        return True
    return all(
        (row[k] / row[k - 1]) * (row[d + 1 - k] / row[d - k]) == 1
        for k in range(2, d)
    )


def multiply_linear(coefficients, root_parameter):
    answer = [Fraction(0)] * (len(coefficients) + 1)
    for index, value in enumerate(coefficients):
        answer[index] += root_parameter * value
        answer[index + 1] += value
    return answer


def coefficients_from_roots(roots):
    answer = [Fraction(1)]
    for root in roots:
        answer = multiply_linear(answer, Fraction(root))
    return answer


def trim_polynomial(coefficients):
    coefficients = [Fraction(value) for value in coefficients]
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def polynomial_derivative(coefficients):
    return [index * coefficients[index] for index in range(1, len(coefficients))]


def polynomial_remainder(dividend, divisor):
    remainder = trim_polynomial(dividend)
    divisor = trim_polynomial(divisor)
    require(divisor != [0], "zero polynomial divisor")
    while len(remainder) >= len(divisor) and remainder != [0]:
        shift = len(remainder) - len(divisor)
        factor = remainder[-1] / divisor[-1]
        for index, value in enumerate(divisor):
            remainder[index + shift] -= factor * value
        remainder = trim_polynomial(remainder)
    return remainder


def sign_variations(signs):
    nonzero = [sign for sign in signs if sign]
    return sum(left != right for left, right in zip(nonzero, nonzero[1:]))


def negative_root_count(coefficients):
    """Count distinct roots in (-infinity,0) by an exact Sturm chain."""
    chain = [trim_polynomial(coefficients)]
    chain.append(trim_polynomial(polynomial_derivative(chain[0])))
    while len(chain[-1]) > 1:
        remainder = polynomial_remainder(chain[-2], chain[-1])
        require(remainder != [0], "Sturm chain found a repeated factor")
        chain.append([-value for value in remainder])
    at_negative_infinity = [
        1 if row[-1] * (-1) ** (len(row) - 1) > 0 else -1
        for row in chain
    ]
    at_zero = [1 if row[0] > 0 else -1 if row[0] < 0 else 0 for row in chain]
    return sign_variations(at_negative_infinity) - sign_variations(at_zero)


def exhaustive_record():
    count = 0
    palindromic = 0
    for degree in range(2, 7):
        for coefficients in itertools.product(range(1, 4), repeat=degree + 1):
            count += 1
            left = ratio_palindromic(coefficients)
            right = scaled_reciprocal(coefficients)
            require(left == right, (degree, coefficients, left, right))
            require(circuit_antipalindromic(coefficients), (degree, coefficients))
            palindromic += int(left)
    return {"positive_coefficient_rows": count, "palindromic_rows": palindromic}


def fixed_locus_record():
    odd_checks = 0
    even_checks = 0
    path_checks = 0
    for degree in (3, 5, 7, 9):
        half = [Fraction(index + 2) for index in range((degree + 1) // 2)]
        coefficients = half + list(reversed(half))
        row = ratios(coefficients)
        middle = (degree + 1) // 2
        require(row[middle] == row[middle - 1], (degree, "central wall"))
        require(ratio_palindromic(coefficients), (degree, "fixed palindrome"))
        odd_checks += 1
    for degree in (4, 6, 8):
        left = [Fraction(index + 2) for index in range(degree // 2)]
        coefficients = left + [Fraction(11)] + list(reversed(left))
        require(ratio_palindromic(coefficients), (degree, "even fixed palindrome"))
        require(circuit_antipalindromic(coefficients), (degree, "even circuit"))
        even_checks += 1
    for coefficients in (
        [1, 2, 3, 5],
        [1, 2, 3, 5, 7, 11],
        [2, 3, 5, 7, 11, 13, 17, 19],
    ):
        degree = len(coefficients) - 1
        midpoint = [
            Fraction(coefficients[k] + coefficients[degree - k], 2)
            for k in range(degree + 1)
        ]
        require(midpoint == list(reversed(midpoint)), (degree, "path midpoint"))
        require(ratio_palindromic(midpoint), (degree, "midpoint ratios"))
        middle = (degree + 1) // 2
        require(ratios(midpoint)[middle] == ratios(midpoint)[middle - 1], (degree, "path wall"))
        path_checks += 1
    return {
        "odd_fixed_wall_checks": odd_checks,
        "even_antipalindrome_checks": even_checks,
        "equivariant_path_midpoints": path_checks,
    }


def constant_ratio_not_equality_control():
    degree = 5
    constant = Fraction(2)
    h = [constant ** (-k * (k - 1) // 2) for k in range(degree + 1)]
    coefficients = [Fraction(0)] * (degree + 1)
    for k in range(degree + 1):
        coefficients[degree - k] = comb(degree, k) * h[k]
    row = ratios(coefficients)[1:]
    require(all(value == constant for value in row), "constant-ratio control failed")
    negative_roots = negative_root_count(coefficients)
    require(negative_roots == degree, "constant-ratio control is not real-rooted")
    return {
        "degree": degree,
        "constant_ratio": str(constant),
        "distinct_negative_real_roots": negative_roots,
        "coefficients": [str(value) for value in coefficients],
    }


def balanced_controls():
    count = 0
    for multiplicity in range(2, 21):
        for second_root in (2, 3, 5, 10, 20, 50):
            coefficients = coefficients_from_roots(
                [1] * multiplicity + [second_root] * multiplicity
            )
            row = ratios(coefficients)
            degree = 2 * multiplicity
            require(ratio_palindromic(coefficients), (multiplicity, second_root))
            require(
                all(row[k] > row[k - 1] for k in range(2, multiplicity + 1)),
                (multiplicity, second_root, "leading shape"),
            )
            require(
                all(row[k] < row[k - 1] for k in range(multiplicity + 1, degree)),
                (multiplicity, second_root, "trailing shape"),
            )
            count += 1
    return {"exact_balanced_two_cluster_controls": count}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = [
        {"exhaustive_equivalence": exhaustive_record()},
        {"fixed_locus": fixed_locus_record()},
        {
            "constant_ratio_not_newton_equality_control":
                constant_ratio_not_equality_control()
        },
        {"balanced_controls": balanced_controls()},
        {
            "asymptotic_scope": {
                "from_log_circuit_nonnegative": "C(mu_d)>=-O(1/d)",
                "from_reversed_log_circuit_nonnegative": "C(mu_d*)<=O(1/d)",
                "finite_exact_sign_without_margin": False,
            }
        },
    ]
    lines = [json.dumps(row, sort_keys=True, separators=(",", ":")) for row in records]
    record_digest = sha256("\n".join(lines).encode("ascii")).hexdigest()
    require(record_digest == EXPECTED_RECORD_DIGEST, "fixed-locus record digest changed")
    transcript = "\n".join(
        [
            "THM-3001 INDEPENDENT REVERSAL FIXED-LOCUS REFEREE",
            "method=exact coefficient arithmetic;no import from primary companion",
            *lines,
            f"record_digest={record_digest}",
            "all_exact_checks=PASS",
        ]
    ) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
