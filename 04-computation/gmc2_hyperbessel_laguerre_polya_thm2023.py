#!/usr/bin/env python3
"""Exact coefficient audit for THM-2023.

Gauss multiplication identifies Phi_(p,q) with a scaled positive-parameter
0F_(p+q-1).  The real-negative zero theorem is external (Baricz--Singh), not
reproved numerically here.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


def denominator_parameters(p: int, q: int) -> list[sp.Rational | sp.Integer]:
    return (
        [sp.Integer(1)]
        + [sp.Rational(j, p) for j in range(1, p)]
        + [sp.Rational(j, q) for j in range(1, q)]
    )


def phi_coefficient(k: int, p: int, q: int) -> sp.Rational:
    return sp.Rational(1, factorial(p * k) * factorial(q * k))


def hyper_coefficient(k: int, p: int, q: int) -> sp.Expr:
    scale = sp.Integer(p) ** p * sp.Integer(q) ** q
    denominator = sp.Integer(factorial(k)) * scale**k
    for parameter in denominator_parameters(p, q):
        denominator *= sp.rf(parameter, k)
    return sp.cancel(1 / denominator)


def verify_pair(p: int, q: int, cutoff: int = 30) -> str:
    params = denominator_parameters(p, q)
    assert len(params) == p + q - 1
    assert all(parameter > 0 for parameter in params)
    for k in range(cutoff + 1):
        assert sp.cancel(phi_coefficient(k, p, q) - hyper_coefficient(k, p, q)) == 0
    return (
        f"Phi_({p},{q}) = 0F_{p+q-1}(-;{params}; "
        f"x/{p**p * q**q}), coefficients k=0..{cutoff} exact"
    )


def main() -> None:
    print("THM-2023 HYPER-BESSEL LAGUERRE-POLYA PARAMETER AUDIT")
    pairs = [(1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (1, 4), (3, 4), (3, 5)]
    for p, q in pairs:
        print("PASS", verify_pair(p, q))
    print("PASS all denominator parameters are positive, the hypothesis of the Baricz-Singh 0Fq zero theorem")
    print("SCOPE only Phi_(p,q), the rd-e=r boundary; Psi_r=sum y^j/(rj)! is a different Mittag-Leffler family")
    print("TOURNAMENT ANALYSIS not used: zero location is unary; candidate vertices (zeros, parameters, boundary sides, proof obligations) admit no pairwise orientation preserving the Laguerre-Polya predicate")
    print("CHALLENGED ASSUMPTION the numerical L-P pattern was open; it is a direct classical 0Fq corollary after the correct scaling")


if __name__ == "__main__":
    main()
