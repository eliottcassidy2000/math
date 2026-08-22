#!/usr/bin/env python3
"""Finite-exact image-algebra probe for the explicit THM-3556 packet.

The packet is rescaled to integral coordinates ``(L,T,W,R)=(L,T,2U,2S)``.
We verify its known cubic relation, exhibit a new quartic relation, and compare
the pullback Hilbert function through target degree seven with that of a
complete intersection of degrees three and four.  The agreement is finite
evidence only: this script does not assert that the two relations generate
the full image ideal.
"""

from __future__ import annotations

from itertools import combinations_with_replacement
from math import comb

import sympy as sp
from flint import fmpz_mat


def require(condition: bool, label: str) -> None:
    """Keep truth-bearing checks active under ``python -O``."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def multiply(first: dict[tuple[int, int], int],
             second: dict[tuple[int, int], int]) -> dict[tuple[int, int], int]:
    out: dict[tuple[int, int], int] = {}
    for (a, b), coefficient in first.items():
        for (c, d), other in second.items():
            exponent = (a + c, b + d)
            out[exponent] = out.get(exponent, 0) + coefficient * other
    return {exponent: coefficient for exponent, coefficient in out.items()
            if coefficient}


def product(polynomials: list[dict[tuple[int, int], int]]) -> dict[tuple[int, int], int]:
    out = {(0, 0): 1}
    for polynomial in polynomials:
        out = multiply(out, polynomial)
    return out


def dictionary(poly: sp.Expr, v: sp.Symbol,
               y: sp.Symbol) -> dict[tuple[int, int], int]:
    return {monomial: int(coefficient)
            for monomial, coefficient in sp.Poly(poly, v, y).terms()}


def cumulative_monomials(maximum_degree: int) -> int:
    if maximum_degree < 0:
        return 0
    return comb(maximum_degree + 4, 4)


def main() -> None:
    v, y = sp.symbols("v y")
    L0, T0, W0, R0 = sp.symbols("L T W R")

    W = sp.expand(2 + 2 * y - y**2 - 3 * v * y**2 + 9 * v * y)
    T = sp.expand(y**2 - 3 * v * W)
    R = sp.expand(2 * y**3 - 9 * v * W * y)
    L = sp.expand(v**2 * (4 * v * W - y**2))
    packet = [L, T, W, R]

    cubic = sp.expand(R0**2 - 4 * T0**3 - 27 * L0 * W0**2)
    quartic = sp.expand(
        27 * L0 * R0**2 - 243 * L0 * R0 * T0 - 81 * L0 * R0 * W0
        + 729 * L0 * R0 - 243 * L0 * T0 * W0**2
        + 486 * L0 * T0 * W0 - 243 * L0 * W0**3
        + 1215 * L0 * W0**2 - 1458 * L0 * W0
        + 3 * R0**2 * T0 + 7 * R0**2 * W0 - 18 * R0**2
        - 15 * R0 * T0**2 - 12 * R0 * T0 * W0 + 48 * R0 * T0
        + 3 * R0 * W0**2 - 16 * R0 * W0 + 36 * R0
        + 21 * T0**2 * W0**2 - 138 * T0**2 * W0 + 192 * T0**2
        + 6 * T0 * W0**3 - 36 * T0 * W0**2 + 48 * T0 * W0
        + W0**4 - 6 * W0**3 + 12 * W0**2 - 8 * W0
    )
    substitution = {L0: L, T0: T, W0: W, R0: R}
    require(sp.expand(cubic.subs(substitution)) == 0,
            "cubic packet relation")
    require(sp.expand(quartic.subs(substitution)) == 0,
            "quartic packet relation")
    require(sp.Poly(quartic, L0, T0, W0, R0).total_degree() == 4,
            "quartic target degree")
    require(len(sp.Poly(quartic, L0, T0, W0, R0).terms()) == 28,
            "quartic term count")

    packet_dicts = [dictionary(entry, v, y) for entry in packet]
    hilbert_rows = []
    for degree_cap in range(1, 8):
        pullbacks = []
        for degree in range(1, degree_cap + 1):
            for factors in combinations_with_replacement(range(4), degree):
                pullbacks.append(product([packet_dicts[index]
                                          for index in factors]))
        support = sorted(set().union(*[set(poly) for poly in pullbacks]),
                         key=lambda exponent: (sum(exponent), exponent))
        matrix_rows = [[poly.get(exponent, 0) for poly in pullbacks]
                       for exponent in support]
        rank = fmpz_mat(matrix_rows).rank()
        kernel = len(pullbacks) - rank
        complete_intersection_kernel = (
            cumulative_monomials(degree_cap - 3)
            + cumulative_monomials(degree_cap - 4)
            - cumulative_monomials(degree_cap - 7)
        )
        require(kernel == complete_intersection_kernel,
                f"degree-{degree_cap} (3,4)-complete-intersection count")
        hilbert_rows.append((degree_cap, len(pullbacks), rank, kernel))

    collision = {
        L0: -sp.Rational(6724, 3645),
        T0: sp.Rational(57, 20),
        W0: sp.Rational(27, 20),
        R0: sp.Rational(27, 20),
    }
    require(cubic.subs(collision) == 0 and quartic.subs(collision) == 0,
            "collision lies on both displayed relations")
    variables = (L0, T0, W0, R0)
    cubic_gradient = sp.Matrix([sp.diff(cubic, entry).subs(collision)
                                for entry in variables])
    quartic_gradient = sp.Matrix([sp.diff(quartic, entry).subs(collision)
                                  for entry in variables])
    require(cubic_gradient != sp.zeros(4, 1), "cubic gradient is nonzero")
    require(quartic_gradient == -sp.Rational(21, 5) * cubic_gradient,
            "relation gradients become dependent at the double point")

    print("THM-3556 PACKET IMAGE-RELATION PROBE")
    print("coordinates=(L,T,W=2U,R=2S)")
    print("relations_verified=degree_3_cusp,degree_4_28_term_relation")
    print("hilbert_rows_degree,total_monomials,pullback_rank,kernel="
          f"{hilbert_rows}")
    print("complete_intersection_3_4_counts_match_through_degree_7=True")
    print("collision_relation_gradient_ratio=quartic/cubic=-21/5")
    print("VERDICT=finite exact evidence for a (3,4) packet-image complete intersection")
    print("SCOPE=full image-ideal generation and primeness remain open")


if __name__ == "__main__":
    main()
