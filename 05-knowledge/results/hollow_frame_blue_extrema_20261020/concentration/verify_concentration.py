#!/usr/bin/env python3
"""Exact checks of a formula-derived concentration theorem, not a frame census."""

import argparse
import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def coefficient_table(max_r, max_j):
    """Rows are coefficients of ((1+t)/(1-t))**r, via multiplication."""
    table = [[1] + [0] * max_j]
    for _ in range(max_r):
        previous = table[-1]
        row = []
        prefix = 0
        for j in range(max_j + 1):
            row.append(previous[j] + 2 * prefix)
            prefix += previous[j]
        table.append(row)
    return table


def composition_coefficient(r, j):
    if j == 0:
        return 1
    return sum(
        2**s * comb(r, s) * comb(j - 1, s - 1)
        for s in range(1, min(r, j) + 1)
    )


def next_odd_strictly_above(q):
    """Least odd integer > the rational number q, including odd-integer q."""
    half = (q + 1) / 2
    return 2 * (half.numerator // half.denominator) + 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-k", type=int, default=128)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.max_k >= 7, "--max-k must be at least 7")

    table = coefficient_table(args.max_k + 1, args.max_k // 2)
    inherited_variances = {
        1: Fraction(0),
        2: Fraction(1),
        3: Fraction(6),
        4: Fraction(75),
        5: Fraction(1500),
        6: Fraction(89685, 2),
        7: Fraction(1875825),
    }
    inherited_maxima = {1: 1, 2: 3, 3: 9, 4: 39, 5: 195}
    rows = []
    term_count = 0

    for k in range(1, args.max_k + 1):
        mu = factorial(k)
        normalized_terms = []
        majorant_terms = []
        for j in range(k // 2 + 1):
            r = k - 2 * j + 1
            c = table[r][j]
            require(
                c == composition_coefficient(r, j),
                f"Coefficient constructions disagree at k={k}, j={j}",
            )
            coefficient_bound = 2**j * comb(k - j, j)
            a = Fraction(factorial(k - 2 * j) * c, 4**j * mu)
            b = Fraction(factorial(k - j), 2**j * factorial(j) * mu)
            d = Fraction(1, factorial(j) * k**j)
            require(
                0 <= c <= coefficient_bound and 0 <= a <= b <= d,
                f"Uniform coefficient/term bounds fail at k={k}, j={j}",
            )
            normalized_terms.append(a)
            majorant_terms.append(d)
            term_count += 1

        relative_variance = sum(normalized_terms, Fraction(0)) - 1
        variance = relative_variance * mu**2
        q = mu * (1 + relative_variance)
        if k in inherited_variances:
            require(
                variance == inherited_variances[k],
                f"Inherited variance disagrees at k={k}",
            )

        if k == 1:
            require(relative_variance == 0 and q == 1, "k=1 exception failed")
            odd_lower_bound = 1
        else:
            require(
                normalized_terms[1] == Fraction(1, 2 * k),
                f"Exact j=1 term failed at k={k}",
            )
            tail = relative_variance - Fraction(1, 2 * k)
            finite_majorant = sum(majorant_terms[2:], Fraction(0))
            rational_majorant = Fraction(3, 2 * k * (3 * k - 1))
            require(
                0 <= tail <= finite_majorant <= rational_majorant,
                f"Uniform tail bounds fail at k={k}",
            )
            odd_lower_bound = next_odd_strictly_above(q)
            require(
                odd_lower_bound % 2 == 1
                and odd_lower_bound - 2 <= q < odd_lower_bound,
                f"Strict odd rounding failed at k={k}",
            )

        if k in inherited_maxima:
            require(
                odd_lower_bound <= inherited_maxima[k],
                f"Lower bound contradicts inherited maximum at k={k}",
            )

        if k <= 10 or k in {16, 32, 64, 128, args.max_k}:
            rows.append(
                {
                    "k": k,
                    "relative_variance_exact": str(relative_variance),
                    "k_times_relative_variance_decimal": float(
                        k * relative_variance
                    ),
                    "Q_exact": str(q),
                    "strict_odd_maximum_lower_bound": odd_lower_bound,
                    "provenance": "second-moment formula; not a frame census",
                }
            )

    require(next_odd_strictly_above(Fraction(7)) == 9, "Odd equality regression")
    require(next_odd_strictly_above(Fraction(5, 2)) == 3, "Fraction regression")
    require(next_odd_strictly_above(Fraction(8)) == 9, "Even equality regression")
    summary = {
        "status": "PASS",
        "checked_k_range": [1, args.max_k],
        "checked_terms": term_count,
        "arithmetic": "exact integer/rational; decimals only for display",
        "checks": [
            "independent coefficient constructions",
            "all allowed j coefficient and normalized-term bounds",
            "exact j=1 contribution",
            "finite sum and rational tail majorants",
            "inherited exact variances k=1,...,7",
            "strict odd-rounding rule",
            "consistency with inherited finite-exact maxima k=1,...,5",
        ],
        "frame_enumeration_performed": False,
        "rows": rows,
    }
    rendered = json.dumps(summary, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(
        f"PASS: k=1,...,{args.max_k}; {term_count} exact term checks; "
        "no frame enumeration."
    )
    for row in rows:
        print(
            f"k={row['k']}: k*Var(h)/(k!)^2="
            f"{row['k_times_relative_variance_decimal']:.12g}"
        )


if __name__ == "__main__":
    main()
