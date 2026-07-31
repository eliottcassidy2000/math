#!/usr/bin/env python3
"""Exact controls for THM-2991's delayed global Newton-ratio return."""

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


def normalized_coefficients(roots: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    elementary = [Fraction(1)]
    for root in roots:
        updated = [Fraction(0)] * (len(elementary) + 1)
        for k, value in enumerate(elementary):
            updated[k] += value
            updated[k + 1] += root * value
        elementary = updated
    degree = len(roots)
    return tuple(elementary[k] / comb(degree, k) for k in range(degree + 1))


def three_cluster_values(n: int, C: int) -> tuple[Fraction, ...]:
    roots = (Fraction(1),) * n + (Fraction(3),) + (Fraction(C),) * n
    return normalized_coefficients(roots)


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

    print("THM-2991 PF-INFINITY ARBITRARILY DELAYED GLOBAL NEWTON-RATIO RETURN")
    print("family=(x+1)^n*(x+3)*(x+C)^n;degree=2n+1;C_positive_integer")

    asymptotic_records = []
    t = Fraction(1, 3)
    for n in range(2, 61):
        degree = 2 * n + 1
        limits = []
        for k in range(1, n):
            c_prev = Fraction(comb(n, k - 1), comb(degree, k - 1))
            c_now = Fraction(comb(n, k), comb(degree, k))
            c_next = Fraction(comb(n, k + 1), comb(degree, k + 1))
            direct = c_now**2 / (c_prev * c_next)
            closed = Fraction(
                (n - k + 1) * (degree - k),
                (degree - k + 1) * (n - k),
            )
            require(direct == closed, (n, k, direct, closed))
            limits.append(closed)
        require(
            all(limits[k] > limits[k - 1] for k in range(1, len(limits))),
            f"limit ladder failed at n={n}",
        )
        center_constant = Fraction(1, (n + 3) * (n + 2))
        leading_limit = Fraction(2 * n * n, degree * (n - 1))
        final_limit = Fraction(2, degree) * (n + t) ** 2 / (n - 1 + 2 * t)
        gap = leading_limit - final_limit
        closed_gap = (
            2
            * t
            * (2 * n - t * (n - 1))
            / (degree * (n - 1) * (n - 1 + 2 * t))
        )
        require(gap == closed_gap > 0, (n, gap, closed_gap))
        asymptotic_records.append(
            f"n={n};L1={fraction_text(limits[0])};"
            f"Llast={fraction_text(limits[-1])};"
            f"center_over_C={fraction_text(center_constant)};"
            f"edge_gap={fraction_text(gap)}"
        )

    witness_records = []
    C = 10**6
    coefficient_checks = 0
    circuit_checks = 0
    for n in range(2, 61):
        values = three_cluster_values(n, C)
        rs = ratios(values)
        degree = 2 * n + 1
        require(
            all(rs[k] > rs[k - 1] for k in range(1, n)),
            f"concrete prefix failed at n={n}",
        )
        require(rs[-1] < rs[0], f"global edge return failed at n={n}")
        require(all(value > 1 for value in rs), f"strict ULC failed at n={n}")
        reciprocal_values = normalized_coefficients(
            (Fraction(1),) * n
            + (Fraction(1, 3),)
            + (Fraction(1, C),) * n
        )
        reciprocal_rs = ratios(reciprocal_values)
        require(rs[-1] == reciprocal_rs[0], f"reciprocal edge failed at n={n}")
        coefficient_checks += degree + 1
        circuit_checks += degree - 1
        witness_records.append(
            f"n={n};C={C};prefix={n};return={2*n};"
            f"R1={fraction_text(rs[0])};Rlast={fraction_text(rs[-1])}"
        )

    small = three_cluster_values(2, 20)
    small_ratios = ratios(small)
    require(
        small
        == (
            Fraction(1),
            Fraction(9),
            Fraction(607, 10),
            Fraction(2283, 10),
            Fraction(584),
            Fraction(1200),
        )
        and small_ratios
        == (
            Fraction(810, 607),
            Fraction(368449, 205470),
            Fraction(5212089, 3544880),
            Fraction(42632, 34245),
        )
        and small_ratios[0] < small_ratios[1]
        and small_ratios[-1] < small_ratios[0],
        "degree-five global-return control changed",
    )

    directional_only = normalized_coefficients(
        (Fraction(1),) * 3 + (Fraction(10**6),) * 3
    )
    directional_ratios = ratios(directional_only)
    require(
        directional_ratios[3] == directional_ratios[1]
        and directional_ratios[3] < directional_ratios[2]
        and directional_ratios[3] > directional_ratios[0],
        "two-cluster directional-only boundary changed",
    )

    asymptotic_digest = sha256("\n".join(asymptotic_records).encode()).hexdigest()
    witness_digest = sha256("\n".join(witness_records).encode()).hexdigest()
    print("asymptotic_n=2..60;records=59;digest=" + asymptotic_digest)
    print("concrete_n=2..60;C=1000000;records=59;digest=" + witness_digest)
    print("limit_ladder=STRICT;center_ratio_over_C=1/((n+3)(n+2))")
    print("last_edge_limit_gap=2t(2n-t(n-1))/(d(n-1)(n-1+2t));t=1/3;POSITIVE")
    print(f"reciprocal_last_edge=PASS;coefficient_checks={coefficient_checks}")
    print(f"strict_ULC_and_global_return=PASS;circuit_checks={circuit_checks}")
    print("degree5_b=1,9,607/10,2283/10,584,1200")
    print("degree5_R=810/607,368449/205470,5212089/3544880,42632/34245")
    print("two_cluster_boundary=directional_turn_only;opposite_inner_ratio_above_R1")
    print("scope=generic_PF_infinity_global_no_return_boundary;not_first_gap_core_counterexample")
    print("all_exact_checks=PASS")

    if handle is not None:
        handle.flush()
        handle.close()


if __name__ == "__main__":
    main()
