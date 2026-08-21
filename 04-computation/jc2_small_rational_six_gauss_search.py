#!/usr/bin/env python3
"""Exact first atlas for the genuinely coupled six-Gauss JC(2) chart.

Universe
--------
For

    M6=E_-(S)E_+(R)E_-(Z)E_+(V)E_-(U)E_+(W),

take U,V,W,Z independently from

    +/-{x,y,x^2,xy,y^2}.

For every one of the 10,000 bases, solve the *entire* coefficient space
R in C[x,y] of total degree at most six.  A second, frontier-facing box takes
U,V,Z from the same signed monomials, W=y^2+/-x, and the entire degree-eight
R space: 2,000 bases and 45 unknown coefficients.  The script proves
inconsistency by maximal-rank certificates modulo 1009: respectively
(28,29) and (45,46) for coefficient/augmented rank.  Maximal coefficient
rank makes these sound characteristic-zero certificates, not bad-prime
heuristics.  Since the top row cannot close, S need not be searched.

Run with both

    python3 04-computation/jc2_small_rational_six_gauss_search.py
    python3 -O 04-computation/jc2_small_rational_six_gauss_search.py
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from pathlib import Path

from flint import nmod_mat

import jc2_small_rational_four_gauss_search as g


Poly = g.Poly
ONE: Poly = {(0, 0): Q(1)}
ZERO: Poly = {}
PRIME = 1009


def require(label: str, condition: bool, checks: list[str]) -> None:
    if not condition:
        raise RuntimeError(f"FAIL: {label}")
    checks.append(label)


def modular_rank(rows: list[list[Q]], columns: int) -> int:
    data = [
        (entry.numerator * pow(entry.denominator, -1, PRIME)) % PRIME
        for row in rows
        for entry in row[:columns]
    ]
    return nmod_mat(len(rows), columns, data, PRIME).rank()


def transport_columns(
    first: Poly,
    second: Poly,
    zeroth: Poly,
    basis: list[g.Monomial],
) -> list[Poly]:
    """Columns of first*d_y-second*d_x+zeroth on a monomial basis."""

    answer = []
    for monomial in basis:
        column = g.add(
            g.multiply(first, g.derivative_monomial(monomial, 1)),
            g.multiply(second, g.derivative_monomial(monomial, 0)),
            Q(-1),
        )
        answer.append(g.add(column, g.multiply(zeroth, {monomial: Q(1)})))
    return answer


def six_word_data(U: Poly, V: Poly, W: Poly, Z: Poly) -> tuple[Poly, ...]:
    """Return A,B,epsilon,delta,C,D,kappa for the staged six-word."""

    A = g.add(ONE, g.multiply(U, V))
    B = g.add(V, g.multiply(W, A))
    epsilon = g.add(g.derivative(A, 1), g.derivative(B, 0), Q(-1))
    delta = g.add(
        g.derivative(U, 1),
        g.derivative(g.multiply(U, W), 0),
        Q(-1),
    )
    C = g.add(U, g.multiply(Z, A))
    D = g.add(g.add(ONE, g.multiply(U, W)), g.multiply(Z, B))
    kappa = g.add(
        delta,
        g.add(
            g.multiply(A, g.derivative(Z, 1)),
            g.multiply(B, g.derivative(Z, 0)),
            Q(-1),
        ),
    )
    kappa = g.add(kappa, g.multiply(epsilon, Z))
    return A, B, epsilon, delta, C, D, kappa


def main() -> None:
    checks: list[str] = []
    signed_monomials = tuple(
        {monomial: Q(sign)}
        for monomial in ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
        for sign in (-1, 1)
    )
    r_basis = g.monomial_basis(6)
    rank_histogram: Counter[tuple[int, int]] = Counter()
    systems = 0

    for U, V, W, Z in product(signed_monomials, repeat=4):
        systems += 1
        _, _, epsilon, _, C, D, kappa = six_word_data(U, V, W, Z)
        columns = transport_columns(C, D, kappa, r_basis)
        rows = g.rational_system_rows(columns, g.scalar(epsilon, Q(-1)))
        coefficient_rank = modular_rank(rows, len(r_basis))
        augmented_rank = modular_rank(rows, len(r_basis) + 1)
        rank_histogram[(coefficient_rank, augmented_rank)] += 1

    require("ten signed monomial parameters", len(signed_monomials) == 10, checks)
    require("R space has dimension 28", len(r_basis) == 28, checks)
    require("ten thousand six-word bases", systems == 10_000, checks)
    require(
        "every top-row system has maximal inconsistent rank",
        rank_histogram == Counter({(28, 29): 10_000}),
        checks,
    )

    # First two-scale boundary left open by the all-degree leading-form gate:
    # W has pure-y quadratic top but a lower x term.  U,V,Z remain signed
    # monomials, while R ranges over all 45 coefficients through degree eight.
    r_two_scale_basis = g.monomial_basis(8)
    two_scale_W = (
        {(0, 2): Q(1), (1, 0): Q(-1)},
        {(0, 2): Q(1), (1, 0): Q(1)},
    )
    two_scale_histogram: Counter[tuple[int, int]] = Counter()
    two_scale_systems = 0
    for U, V, W, Z in product(
        signed_monomials, signed_monomials, two_scale_W, signed_monomials
    ):
        two_scale_systems += 1
        _, _, epsilon, _, C, D, kappa = six_word_data(U, V, W, Z)
        columns = transport_columns(C, D, kappa, r_two_scale_basis)
        rows = g.rational_system_rows(columns, g.scalar(epsilon, Q(-1)))
        coefficient_rank = modular_rank(rows, len(r_two_scale_basis))
        augmented_rank = modular_rank(rows, len(r_two_scale_basis) + 1)
        two_scale_histogram[(coefficient_rank, augmented_rank)] += 1
    require("two-scale R space has dimension 45",
            len(r_two_scale_basis) == 45, checks)
    require("two thousand two-scale bases", two_scale_systems == 2_000, checks)
    require(
        "every two-scale top-row system has maximal inconsistent rank",
        two_scale_histogram == Counter({(45, 46): 2_000}),
        checks,
    )

    # Positive sign/control for the genuinely coupled equations, with all
    # six factors nonzero.  U=1,V=y,W=-1,Z=-1,R=1 gives the fifth-prefix
    # rows (1,0),(-y,1), with epsilon=1 and kappa=-1.  S=y repairs the
    # bottom row and the complete six-word is the identity.
    U = ONE
    V = {(0, 1): Q(1)}
    W = {(0, 0): Q(-1)}
    Z = {(0, 0): Q(-1)}
    R = ONE
    A, B, epsilon, _, C, D, kappa = six_word_data(U, V, W, Z)
    require("coupled control has epsilon one", epsilon == ONE, checks)
    require("coupled control has kappa minus one", kappa == {(0, 0): Q(-1)}, checks)
    E = g.add(A, g.multiply(R, C))
    F = g.add(B, g.multiply(R, D))
    require("R repairs the top row to identity", E == ONE and not F, checks)
    S = {(0, 1): Q(1)}
    final_bottom_first = g.add(C, g.multiply(S, E))
    final_bottom_second = g.add(D, g.multiply(S, F))
    require("coupled control repairs to identity bottom row",
            not final_bottom_first and final_bottom_second == ONE, checks)
    hostile_bottom_curl = g.add(g.derivative(C, 1), g.derivative(D, 0), Q(-1))
    require("omitting S leaves the advertised nonzero curl",
            hostile_bottom_curl == {(0, 0): Q(-1)}, checks)

    # A richer height-two tame control has nonconstant U,V,W and two equal,
    # nonzero polynomial defects.  It exercises the full staged mechanism
    # outside the signed-monomial base universe rather than collapsing by a
    # constant elementary relation.
    x = {(1, 0): Q(1)}
    y = {(0, 1): Q(1)}
    s = g.add(x, {(0, 2): Q(1, 2)})
    U, V, W, Z, R = s, y, y, ONE, {(0, 0): Q(-1)}
    A, B, epsilon, _, C, D, kappa = six_word_data(U, V, W, Z)
    require("height-two control has equal nonconstant defects",
            epsilon == s and kappa == s, checks)
    E = g.add(A, g.multiply(R, C))
    F = g.add(B, g.multiply(R, D))
    q = g.add(y, g.scalar(g.power(s, 2), Q(1, 2)))
    require("height-two R stage is minus dq",
            E == g.scalar(g.derivative(q, 0), Q(-1))
            and F == g.scalar(g.derivative(q, 1), Q(-1)), checks)
    S = g.scalar(g.power(s, 2), Q(-1, 2))
    final_C = g.add(C, g.multiply(S, E))
    final_D = g.add(D, g.multiply(S, F))
    target_Q = g.add(g.add(s, q), g.scalar(g.power(q, 2), Q(1, 2)))
    require("height-two S stage is d(s+q+q^2/2)",
            final_C == g.derivative(target_Q, 0)
            and final_D == g.derivative(target_Q, 1), checks)
    require("height-two final determinant is one",
            g.add(g.multiply(E, final_D), g.multiply(F, final_C), Q(-1)) == ONE,
            checks)

    semantic_lines = [
        "status=FINITE-EXACT; no JC(2) conclusion outside the stated universe",
        "word=E_-(S)E_+(R)E_-(Z)E_+(V)E_-(U)E_+(W)",
        "base=U,V,W,Z independently in +/-{x,y,x^2,xy,y^2}",
        "R=entire C[x,y] total-degree<=6; constant included; dim=28",
        f"systems={systems};prime={PRIME};rank_histogram={dict(rank_histogram)}",
        ("certificate=maximal coefficient rank 28 and augmented rank 29 modulo "
         "1009 implies characteristic-zero inconsistency"),
        "S_not_searched=top row never closes",
        ("two_scale_base=U,V,Z in +/-{x,y,x^2,xy,y^2};W=y^2+/-x;"
         "R=entire total-degree<=8,dim=45"),
        (f"two_scale_systems={two_scale_systems};prime={PRIME};"
         f"rank_histogram={dict(two_scale_histogram)}"),
        ("positive_control=U=1,V=y,W=-1,Z=-1,R=1,S=y has "
         "epsilon=1,kappa=-1 and gives identity"),
        ("hostile=omit R or S: top defect remains +1 or bottom curl remains "
         "-1, respectively"),
        ("height2_control=s=x+y^2/2;W=V=y,U=s,Z=1,R=-1,S=-s^2/2;"
         "epsilon=kappa=s;F=(-q,s+q+q^2/2),q=y+s^2/2"),
        "counterexample_survivors=0",
    ]
    semantic_digest = sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()
    source_digest = sha256(Path(__file__).read_bytes()).hexdigest()
    dependency_digest = sha256(Path(g.__file__).read_bytes()).hexdigest()
    for line in semantic_lines:
        print(line)
    print(f"checks={len(checks)}")
    print(f"semantic_sha256={semantic_digest}")
    print(f"dependency_sha256={dependency_digest}")
    print(f"source_sha256={source_digest}")


if __name__ == "__main__":
    main()
