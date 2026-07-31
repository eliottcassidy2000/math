#!/usr/bin/env python3
"""Exact controls for THM-2991's arbitrarily delayed Newton-ratio return."""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def b_values(n: int, B: int) -> tuple[Fraction, ...]:
    values = []
    for k in range(2 * n + 1):
        elementary = sum(
            comb(n, j) * comb(n, k - j) * B**j
            for j in range(max(0, k - n), min(n, k) + 1)
        )
        values.append(Fraction(elementary, comb(2 * n, k)))
    return tuple(values)


def ratios(values: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(
        values[k] ** 2 / (values[k - 1] * values[k + 1])
        for k in range(1, len(values) - 1)
    )


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    handle = None
    if args.output is not None:
        handle = args.output.open("w", encoding="utf-8", newline="\n")
        sys.stdout = handle

    print("THM-2991 PF-INFINITY ARBITRARILY DELAYED NEWTON-RATIO RETURN")
    print("family=(x+1)^n*(x+B)^n;degree=2n;B_positive_integer")

    asymptotic_records = []
    for n in range(2, 41):
        limits = []
        for k in range(1, n):
            c_prev = Fraction(comb(n, k - 1), comb(2 * n, k - 1))
            c_now = Fraction(comb(n, k), comb(2 * n, k))
            c_next = Fraction(comb(n, k + 1), comb(2 * n, k + 1))
            direct = c_now**2 / (c_prev * c_next)
            closed = Fraction(
                (n - k + 1) * (2 * n - k),
                (2 * n - k + 1) * (n - k),
            )
            require(direct == closed, (n, k, direct, closed))
            limits.append(closed)
        require(
            all(limits[k] > limits[k - 1] for k in range(1, len(limits))),
            f"limit ladder failed at n={n}",
        )
        center_constant = Fraction(1, (n + 1) ** 2)
        asymptotic_records.append(
            f"n={n};L1={fraction_text(limits[0])};"
            f"Llast={fraction_text(limits[-1])};"
            f"center_over_B={fraction_text(center_constant)}"
        )

    witness_records = []
    B = 10**6
    for n in range(2, 41):
        values = b_values(n, B)
        rs = ratios(values)
        require(
            all(rs[k] > rs[k - 1] for k in range(1, n)),
            f"concrete prefix failed at n={n}",
        )
        require(
            all(
                values[2 * n - k]
                == Fraction(B) ** (n - k) * values[k]
                for k in range(2 * n + 1)
            ),
            f"coefficient reciprocity failed at n={n}",
        )
        require(
            all(rs[2 * n - k - 1] == rs[k - 1] for k in range(1, 2 * n)),
            f"ratio reciprocity failed at n={n}",
        )
        require(rs[n] == rs[n - 2] < rs[n - 1], f"return failed at n={n}")
        require(all(value > 1 for value in rs), f"strict ULC failed at n={n}")
        witness_records.append(
            f"n={n};B={B};first_return={n + 1};"
            f"R1={fraction_text(rs[0])};Rn={fraction_text(rs[n - 1])}"
        )

    small = b_values(2, 2)
    small_ratios = ratios(small)
    require(
        small == (Fraction(1), Fraction(3, 2), Fraction(13, 6), Fraction(3), Fraction(4))
        and small_ratios
        == (Fraction(27, 26), Fraction(169, 162), Fraction(27, 26)),
        "degree-four hostile changed",
    )

    asymptotic_digest = sha256("\n".join(asymptotic_records).encode()).hexdigest()
    witness_digest = sha256("\n".join(witness_records).encode()).hexdigest()
    print("asymptotic_n=2..40;records=39;digest=" + asymptotic_digest)
    print("concrete_n=2..40;B=1000000;records=39;digest=" + witness_digest)
    print("limit_ladder=STRICT;center_ratio_over_B=1/(n+1)^2")
    print("reciprocal_identity=b_(2n-k)=B^(n-k)*b_k;all_1677_coefficients=PASS")
    print("ratio_symmetry=R_(2n-k)=R_k;all_1599_circuits=PASS")
    print("first_return=n+1;concrete_witnesses=39;strict_ULC=PASS")
    print("degree4_b=1,3/2,13/6,3,4")
    print("degree4_R=27/26,169/162,27/26")
    print("scope=generic_PF_infinity_no_return_boundary;not_first_gap_core_counterexample")
    print("all_exact_checks=PASS")

    if handle is not None:
        handle.flush()
        handle.close()


if __name__ == "__main__":
    main()
