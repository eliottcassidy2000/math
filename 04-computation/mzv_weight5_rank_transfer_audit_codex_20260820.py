"""Exact weight-five/type audit for Zagier-dimension claims.

Status: PROVED elementary identities + FINITE-EXACT rank audit.

This does not prove the classical Zagier dimension conjecture.  It identifies
the exact weight-five carrier, its change-of-basis determinant, and the
information lost when a formal/motivic or mod-p rank statement is pushed to
real periods without an injectivity theorem.
"""

from __future__ import annotations

from fractions import Fraction

from flint import arb, ctx
from sympy import Matrix, isprime
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


EXCEPTIONAL_CANDIDATE = 543_606_522_303_979
ANNOUNCED_NUMERATOR = 285_075_330_345_178


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def zagier_dimensions(top: int) -> list[int]:
    d = [0] * (top + 1)
    d[0] = 1
    if top >= 2:
        d[2] = 1
    for weight in range(3, top + 1):
        d[weight] = d[weight - 2] + d[weight - 3]
    return d


def hoffman_words(weight: int) -> list[tuple[int, ...]]:
    if weight == 0:
        return [()]
    out: list[tuple[int, ...]] = []
    for first in (2, 3):
        if first <= weight:
            out.extend((first,) + tail for tail in hoffman_words(weight - first))
    return out


def rank_mod_prime(matrix: Matrix, p: int) -> int:
    rows = [[int(x) % p for x in matrix.row(i)] for i in range(matrix.rows)]
    rank = 0
    col = 0
    while rank < len(rows) and col < matrix.cols:
        pivot = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = pow(rows[rank][col], -1, p)
        rows[rank] = [(inv * x) % p for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col]:
                q = rows[r][col]
                rows[r] = [(x - q * y) % p for x, y in zip(rows[r], rows[rank])]
        rank += 1
        col += 1
    return rank


def main() -> None:
    dims = zagier_dimensions(24)
    for weight in range(25):
        require(len(hoffman_words(weight)) == dims[weight], f"Padovan count failed at {weight}")
    print("Zagier recurrence / Hoffman {2,3}-word counts PASS through weight 24")
    print("d_0,...,d_12 =", dims[:13])
    print("weight-five Hoffman words =", hoffman_words(5))
    quotient_dims = [dims[w] - (dims[w - 2] if w >= 2 else 0) for w in range(25)]
    require(quotient_dims[5] == 1, "finite/symmetric quotient dimension drift")
    print("predicted dimensions after quotienting by zeta(2), weights 0..12 =", quotient_dims[:13])
    print("in particular, the weight-five finite/symmetric quotient has dimension 1, not 2")

    # Euler's double-zeta reductions, in the ordered target basis
    # (zeta(5), zeta(2)zeta(3)):
    #   zeta(2,3) =  9/2 zeta(5) - 2 zeta(2)zeta(3)
    #   zeta(3,2) = -11/2 zeta(5) + 3 zeta(2)zeta(3).
    rational_change = [
        [Fraction(9, 2), Fraction(-2)],
        [Fraction(-11, 2), Fraction(3)],
    ]
    determinant = (
        rational_change[0][0] * rational_change[1][1]
        - rational_change[0][1] * rational_change[1][0]
    )
    require(determinant == Fraction(5, 2), "weight-five determinant drift")
    integral_change = Matrix([[9, -4], [-11, 6]])
    snf = smith_normal_form(integral_change, domain=ZZ)
    require(abs(int(integral_change.det())) == 10, "cleared determinant drift")
    require([abs(int(snf[i, i])) for i in range(2)] == [1, 10], "SNF drift")
    print("weight-five change determinant = 5/2")
    print("cleared integral Smith factors = [1, 10]; only lattice primes 2 and 5 occur")

    p0 = EXCEPTIONAL_CANDIDATE
    require(isprime(p0), "reported exceptional integer is not prime")
    require(p0 not in (2, 5), "candidate unexpectedly belongs to weight-five determinant")
    print(f"reported exceptional candidate p0={p0} is prime")
    print("p0 cannot come from the canonical weight-five 2x2 basis change")

    # Reproduce the public numerical sidecar with a certified real ball.  It
    # excludes one rational having denominator exactly p0; it does not supply
    # the missing theorem that every rational value would have that form.
    ctx.dps = 80
    announced_gap = (
        arb(ANNOUNCED_NUMERATOR) / p0
        - arb(5).zeta() / (arb(2).zeta() * arb(3).zeta())
    )
    require(announced_gap.lower() > 0, "announced rational gap lost its sign")
    print("certified announced denominator-p0 gap =", announced_gap)
    print("this excludes exactly one denominator-p0 rational, not rationality")

    # Minimal rank-transfer hostile: full characteristic-zero rank with one
    # exceptional reduction.  The exception detects torsion in the chosen
    # integral presentation, not a kernel of an unrelated real period map.
    presentation = Matrix([[1, 0], [0, p0]])
    require(presentation.rank() == 2, "Q-rank hostile failed")
    require(rank_mod_prime(presentation, p0) == 1, "exceptional reduction hostile failed")
    for prime in (2, 3, 5, 7, 11, 13, 17, 19):
        require(rank_mod_prime(presentation, prime) == 2, f"good-prime hostile failed at {prime}")
    print("single-exception Smith hostile PASS: Q-rank 2, mod-p0 rank 1")

    # The dimension statement at weight five is exactly a two-column period
    # injectivity statement because zeta(2)zeta(3) is nonzero.
    print()
    print("TYPE CONCLUSION")
    print("classical d_5=2 iff zeta(5)/(zeta(2)zeta(3)) is irrational")
    print("formal/motivic rank 2 supplies the source; real period injectivity is the missing map")
    print("finite/symmetric and p-adic realizations kill the zeta(2)zeta(3) column")
    print("a mod-p exceptional-prime theorem alone supplies neither the real period map nor its injectivity")


if __name__ == "__main__":
    main()
