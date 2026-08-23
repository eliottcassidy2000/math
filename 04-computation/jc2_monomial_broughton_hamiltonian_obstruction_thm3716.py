#!/usr/bin/env python3
"""Exact companion for THM-3716's monomial Broughton obstruction."""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


x, y = sp.symbols("x y")


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(left, x) * sp.diff(right, y)
        - sp.diff(left, y) * sp.diff(right, x)
    )


semantic_rows: list[str] = []

for m in range(2, 10):
    for n in range(1, 9):
        q = x + x**m * y**n
        qx = sp.diff(q, x)
        qy = sp.diff(q, y)
        z = x ** (m - 1) * y**n
        bezout_x = sp.Rational(m * m, n) * x ** (m - 2) * y ** (n + 1)
        bezout_y = m * z - 1
        tau_coefficient = sp.Rational(m * (m + n), n)
        tau = tau_coefficient * x ** (m - 2) * y**n

        gate(sp.expand(q - x * (1 + z)) == 0, "reducible zero fibre")
        gate(
            sp.expand(bezout_x * qy - bezout_y * qx) == 1,
            "gradient Bezout identity",
        )
        gate(
            sp.expand(sp.diff(bezout_x, y) - sp.diff(bezout_y, x) - tau) == 0,
            "complementary curl",
        )

        for i in range(0, 7):
            for j in range(0, 7):
                monomial = x**i * y**j
                expected = (
                    -j * x**i * y ** (j - 1)
                    + (n * i - m * j)
                    * x ** (i + m - 1)
                    * y ** (j + n - 1)
                )
                gate(
                    sp.expand(bracket(monomial, q) - expected) == 0,
                    "Hamiltonian two-edge formula",
                )

        coefficients: list[Fraction] = [
            -Fraction(m * (m + n), n * (n + 1))
        ]
        for k in range(1, 13):
            coefficients.append(
                -Fraction(m + (k + 1) * n, (k + 1) * n + 1)
                * coefficients[-1]
            )
        for k, coefficient in enumerate(coefficients):
            i_k = m - 2 + k * (m - 1)
            j_k = (k + 1) * n + 1
            gate(coefficient != 0, "tail coefficient nonzero")
            gate(
                n * i_k - m * j_k == -(m + (k + 2) * n),
                "tail transport coefficient",
            )

            partial = sum(
                sp.Rational(value.numerator, value.denominator)
                * x ** (m - 2 + index * (m - 1))
                * y ** ((index + 1) * n + 1)
                for index, value in enumerate(coefficients[: k + 1])
            )
            residual_coefficient = (
                -(m + (k + 2) * n)
                * sp.Rational(coefficient.numerator, coefficient.denominator)
            )
            residual = (
                residual_coefficient
                * x ** (m - 2 + (k + 1) * (m - 1))
                * y ** ((k + 2) * n)
            )
            gate(
                sp.expand(bracket(partial, q) - tau - residual) == 0,
                "finite tail residual",
            )

        semantic_rows.append(
            f"{m},{n}:" + hashlib.sha256(
                "|".join(
                    (
                        sp.srepr(q),
                        sp.srepr(bezout_x),
                        sp.srepr(bezout_y),
                        sp.srepr(tau),
                        ",".join(str(value) for value in coefficients),
                    )
                ).encode()
            ).hexdigest()
        )

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

print("theorem=THM-3716-monomial-Broughton-Hamiltonian-obstruction-family")
print("family=Q=x+x^m*y^n;m>=2;n>=1;characteristic=0")
print("gradient=unimodular;zero_fibre=reducible;coordinate=NO")
print("obstruction=tau=m(m+n)/n*x^(m-2)*y^n;hamiltonian_image=NO")
print("mechanism=isolated_two_edge_lattice_chain_with_nonterminating_coefficients")
print("hostile_grid=m:2..9;n:1..8;monomials:0..6x0..6;tail_depth:12=PASS")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
