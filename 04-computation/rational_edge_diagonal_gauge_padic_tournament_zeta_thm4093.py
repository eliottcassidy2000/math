#!/usr/bin/env python3
"""Exact audit for THM-4093.

The symbolic theorem is proved in the companion markdown file.  This program
hostile-tests its three finite interfaces:

* vertex-ratio edge weights are diagonal similarity over Q;
* the cubic coefficient of det(I-uA) is minus the directed-triangle count;
* the resulting p-adic tangent has the claimed exact valuation.

All calculations are exact.  Assertions are expressed through ``require`` so
the optimized interpreter executes every gate as well.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from itertools import combinations
from math import gcd
from random import Random

import sympy as sp

sys.dont_write_bytecode = True


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def directed_triangle_count(adjacency: sp.Matrix) -> int:
    total = 0
    for i, j, k in combinations(range(adjacency.rows), 3):
        total += int(
            adjacency[i, j] * adjacency[j, k] * adjacency[k, i]
            + adjacency[i, k] * adjacency[k, j] * adjacency[j, i]
        )
    return total


def principal_cubic_sum(adjacency: sp.Matrix) -> int:
    total = 0
    for i, j, k in combinations(range(adjacency.rows), 3):
        # Six-term Leibniz determinant, deliberately separate from the two
        # directed-cycle products used by ``directed_triangle_count``.
        total += int(
            adjacency[i, i] * adjacency[j, j] * adjacency[k, k]
            + adjacency[i, j] * adjacency[j, k] * adjacency[k, i]
            + adjacency[i, k] * adjacency[j, i] * adjacency[k, j]
            - adjacency[i, k] * adjacency[j, j] * adjacency[k, i]
            - adjacency[i, j] * adjacency[j, i] * adjacency[k, k]
            - adjacency[i, i] * adjacency[j, k] * adjacency[k, j]
        )
    return total


def tournament_from_bits(n: int, bits: int) -> sp.Matrix:
    adjacency = sp.zeros(n)
    cursor = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> cursor) & 1:
                adjacency[i, j] = 1
            else:
                adjacency[j, i] = 1
            cursor += 1
    return adjacency


def stern_depth(a: int, b: int) -> int:
    divisor = gcd(a, b)
    a //= divisor
    b //= divisor
    quotient_sum = 0
    while b:
        quotient_sum += a // b
        a, b = b, a % b
    return quotient_sum - 1


def stern_tournament(n: int) -> sp.Matrix:
    adjacency = sp.zeros(n)
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if (-1) ** stern_depth(a, b) == 1:
                adjacency[a - 1, b - 1] = 1
            else:
                adjacency[b - 1, a - 1] = 1
    return adjacency


def apex_imbalance(adjacency: sp.Matrix) -> int:
    apex = adjacency.rows - 1
    incoming = sum(int(adjacency[i, apex]) for i in range(apex))
    outgoing = sum(int(adjacency[apex, i]) for i in range(apex))
    return incoming - outgoing


def valuation(integer: int, prime: int) -> int:
    require(integer != 0, "valuation of zero is not used")
    integer = abs(integer)
    answer = 0
    while integer % prime == 0:
        integer //= prime
        answer += 1
    return answer


def rational_edge_gauge_audit() -> int:
    rng = Random(4093)
    cases = 0
    x = sp.symbols("x")
    for n in range(2, 9):
        for _ in range(18):
            adjacency = sp.zeros(n)
            for i in range(n):
                for j in range(n):
                    if i != j and rng.randrange(3) == 0:
                        adjacency[i, j] = rng.randrange(-3, 4)
            labels = []
            for _i in range(n):
                numerator = 0
                while numerator == 0:
                    numerator = rng.randrange(-11, 12)
                labels.append(sp.Rational(numerator, rng.randrange(1, 12)))
            diagonal = sp.diag(*labels)
            weighted = sp.Matrix(
                n,
                n,
                lambda i, j: adjacency[i, j] * labels[i] / labels[j],
            )
            require(weighted == diagonal * adjacency * diagonal.inv(), "similarity")
            require(
                weighted.charpoly(x).as_expr() == adjacency.charpoly(x).as_expr(),
                "characteristic polynomial changed",
            )
            for power in range(1, n + 3):
                require(
                    sp.trace(weighted**power) == sp.trace(adjacency**power),
                    "closed-walk trace changed",
                )
            cases += 1
    return cases


def tournament_cubic_audit() -> tuple[int, int]:
    tournaments = 0
    cubic_checks = 0
    for n in range(1, 7):
        for bits in range(1 << (n * (n - 1) // 2)):
            adjacency = tournament_from_bits(n, bits)
            c3 = directed_triangle_count(adjacency)
            # The coefficient of u^3 in det(I-uA) is the negative sum of
            # the principal 3-by-3 minors.  Each such minor is one precisely
            # for a cyclic triple and zero for a transitive triple.
            require(principal_cubic_sum(adjacency) == c3, "cubic principal minor")
            tournaments += 1
            cubic_checks += max(0, n * (n - 1) * (n - 2) // 6)
    return tournaments, cubic_checks


def padic_tangent_audit(polynomial: sp.Poly, c3: int) -> int:
    checks = 0
    for prime in (2, 3, 5, 7, 11, 17, 19):
        if c3 % prime == 0:
            continue
        for m in range(1, 6):
            point = prime**m
            p_value = int(polynomial.eval(point))
            numerator = 1 - p_value
            require(p_value % prime != 0, "zeta denominator is not a p-adic unit")
            require(valuation(numerator, prime) == 3 * m, "wrong p-adic tangent")
            checks += 1
            unit_point = point * (prime + 1)
            unit_value = int(polynomial.eval(unit_point))
            require(unit_value % prime != 0, "unit-shifted denominator is not a unit")
            require(
                valuation(1 - unit_value, prime) == 3 * m,
                "unit-shifted p-adic tangent",
            )
            checks += 1
    return checks


def main() -> None:
    u = sp.symbols("u")
    gauge_cases = rational_edge_gauge_audit()
    tournament_count, cubic_checks = tournament_cubic_audit()

    tangent_checks = 0
    controls = []
    for n in (5, 9, 13):
        adjacency = stern_tournament(n)
        polynomial = sp.Poly(sp.expand((sp.eye(n) - u * adjacency).det()), u)
        c3 = directed_triangle_count(adjacency)
        require(polynomial.coeff_monomial(u**3) == -c3, "full determinant cubic")
        tangent_checks += padic_tangent_audit(polynomial, c3)
        controls.append(
            (n, apex_imbalance(adjacency), c3, sp.factor(polynomial.as_expr()))
        )

    require(
        sp.expand(controls[0][3] - (1 - u**3)) == 0,
        "Stern T5 is not the C3 control",
    )
    require(controls[1][1] == 0 and controls[1][2] == 17, "Stern T9 hostile drift")

    divisible_adjacency = tournament_from_bits(4, 4)
    divisible_polynomial = sp.Poly(
        sp.expand((sp.eye(4) - u * divisible_adjacency).det()), u
    )
    require(
        divisible_polynomial.as_expr() == 1 - 2 * u**3 - u**4,
        "p-divides-c3 hostile polynomial drift",
    )
    require(
        valuation(1 - int(divisible_polynomial.eval(2)), 2) == 5,
        "p-divides-c3 hostile valuation drift",
    )

    print("THM-4093 exact audit")
    print(f"rational diagonal-gauge cases: {gauge_cases}")
    print(f"tournaments exhausted through n=6: {tournament_count}")
    print(f"cyclic/transitive triple slots audited: {cubic_checks}")
    print(f"exact p-adic tangent checks: {tangent_checks}")
    for n, imbalance, c3, polynomial in controls:
        print(f"Stern T_{n}: apex B={imbalance}, c3={c3}, det(I-uA)={polynomial}")
    print("hostile: B(T_9)=0 but c3(T_9)=17 and zeta is nontrivial")
    print("sharpness hostile: P=1-2u^3-u^4, p=x=2 gives v_p(zeta(x)-1)=5")
    print("all exact gates PASS")


if __name__ == "__main__":
    main()
