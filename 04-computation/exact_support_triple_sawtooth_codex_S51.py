#!/usr/bin/env python3
"""Exact 28-residue law for the THM-935 support mass M({1,2,N}).

At lambda = 1/14, put p = 2*lambda = 1/7 and let B_v be the bad arc
``||v*t|| < lambda``.  The subset-Moebius definition gives

    M({1,2,N})
      = mu(B_1 cap B_2 cap B_N)
        - p * sum_{pairs A} mu(cap_{v in A} B_v) + 2*p**3.

The proof implemented below is not a fitted interpolation.  It rescales the
two fixed windows B_1 and B_1 cap B_2 and evaluates two explicit periodic
tooth-discrepancy functions D and E.  This gives

    M({1,2,N}) = k[N mod 28] / (686*N),

where k is printed below.  Two independent exact intersection engines--an
endpoint sweep and a midpoint-cell count--then referee the formula.

Tournament analysis / challenged assumption
--------------------------------------------
The faithful finite quotient has vertices given by the 28 positions of the
moving tooth boundary relative to the fixed windows.  It preserves the exact
Moebius mass, its sign, endpoint chronology, and the antipodal involution
``r -> 28-r``.  A runner tournament does not: pair mass is symmetric, and the
first sign-bearing datum is a signed triple hyperedge.  Orienting pairs by a
label tie-break creates a transitive tournament whose score histogram, cycle
count, SCCs, edge flips, and Hamiltonian paths are artifacts.  The natural
object here is the signed support hypergraph decorated by a cyclic 28-residue
boundary walk, not a tournament on runners.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, groupby
from math import gcd


LAMBDA = F(1, 14)
P = 2 * LAMBDA


def lcm(values: tuple[int, ...]) -> int:
    result = 1
    for value in values:
        result = result * abs(value) // gcd(result, abs(value))
    return result


@lru_cache(maxsize=None)
def intersection_sweep(support: tuple[int, ...]) -> F:
    """Exact endpoint sweep for ``mu(cap B_v)``."""
    support = tuple(sorted(support))
    scale = lcm(support)
    denominator = 14 * scale
    events: list[tuple[int, int]] = []
    for speed in support:
        width = scale // speed
        for center in range(speed):
            low = (14 * center - 1) * width
            high = (14 * center + 1) * width
            if low < 0:
                events.extend(
                    [(0, 1), (high, -1), (denominator + low, 1), (denominator, -1)]
                )
            else:
                events.extend([(low, 1), (high, -1)])

    depth = 0
    previous = 0
    total = 0
    target = len(support)
    for coordinate, same_coordinate in groupby(sorted(events), key=lambda event: event[0]):
        if depth == target:
            total += coordinate - previous
        depth += sum(delta for _, delta in same_coordinate)
        previous = coordinate
    return F(total, denominator)


@lru_cache(maxsize=None)
def intersection_midpoints(support: tuple[int, ...]) -> F:
    """Independent exact midpoint count on the common endpoint grid."""
    support = tuple(sorted(support))
    scale = lcm(support)
    denominator = 14 * scale
    doubled = 2 * denominator
    count = 0
    for cell in range(denominator):
        midpoint = 2 * cell + 1
        if all(
            14 * min((speed * midpoint) % doubled, (-speed * midpoint) % doubled)
            < doubled
            for speed in support
        ):
            count += 1
    return F(count, denominator)


def direct_mass(intersection, n: int) -> F:
    """The support-three subset-Moebius mass, using either exact engine."""
    triple = (1, 2, n)
    pair_sum = sum(intersection(tuple(pair)) for pair in combinations(triple, 2))
    return intersection(triple) - P * pair_sum + 2 * P**3


def unit_fractional_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def integer_tooth_discrepancy(x: F) -> F:
    """D(x) = integral_[-x,x] 1_{||u||<1/14} du - 2x/7.

    D is one-periodic.  The three pieces below follow by intersecting
    ``[-x,x]`` with the teeth centered at the integers.
    """
    x = unit_fractional_part(x)
    if x <= LAMBDA:
        return F(12, 7) * x
    if x <= 1 - LAMBDA:
        return (1 - 2 * x) / 7
    return -F(12, 7) * (1 - x)


def half_integer_tooth_discrepancy(x: F) -> F:
    """E(x) for teeth centered at the half-integers, also one-periodic."""
    x = unit_fractional_part(x)
    if x <= F(3, 7):
        return -2 * x / 7
    if x <= F(4, 7):
        return F(12, 7) * x - F(6, 7)
    return 2 * (1 - x) / 7


def discrepancy_mass(n: int) -> F:
    """Mass obtained from the two proved tooth-discrepancy formulas.

    Here is the complete window reduction.  On the circle, ``B_1`` is the
    interval centered at zero of radius ``1/14``, while

        B_1 cap B_2 = (-1/28, 1/28).

    Put ``alpha=N/28``.  Rescaling either fixed window by ``u=N*t`` gives

        mu(B_1 cap B_N)             = p^2   + D(2*alpha)/N,
        mu(B_1 cap B_2 cap B_N)     = p^2/2 + D(alpha)/N,

    and ``mu(B_1 cap B_2)=p/2``.  The set ``B_2`` has two radius-``1/28``
    components, centered at 0 and 1/2.  Translation by 1/2 changes ``N*t``
    by ``N/2``.  Thus the second component sees another integer tooth when N
    is even and a half-integer tooth when N is odd:

        mu(B_2 cap B_N) = p^2 + 2*D(alpha)/N                 (N even),
        mu(B_2 cap B_N) = p^2 + (D(alpha)+E(alpha))/N        (N odd).

    Substitution in

        M_3 = mu(B_1 cap B_2 cap B_N)
              - p * sum_pairs mu(B_i cap B_j) + 2*p^3

    cancels every independent term.  What remains is exactly the even/odd
    expression returned below.  Since D and E are one-periodic, ``N*M_3`` is
    28-periodic before any finite computation is performed.
    """
    alpha = F(n, 28)
    d_alpha = integer_tooth_discrepancy(alpha)
    d_twice = integer_tooth_discrepancy(2 * alpha)
    if n % 2 == 0:
        return (5 * d_alpha - d_twice) / (7 * n)
    e_alpha = half_integer_tooth_discrepancy(alpha)
    return (6 * d_alpha - d_twice - e_alpha) / (7 * n)


def positive_coefficient(residue: int) -> int:
    """k(r) on 1 <= r <= 13, derived by substituting in D and E."""
    assert 1 <= residue <= 13
    if residue % 2 == 0:
        return 56 - 3 * residue
    if residue == 1:
        return 25
    if residue == 13:
        return 24
    return 70 - 3 * residue


def residue_coefficient(n: int) -> int:
    residue = n % 28
    if residue in (0, 14):
        return 0
    if residue < 14:
        return positive_coefficient(residue)
    return -positive_coefficient(28 - residue)


def closed_mass(n: int) -> F:
    assert n > 2
    return F(residue_coefficient(n), 686 * n)


def main() -> None:
    coefficients = tuple(residue_coefficient(residue) for residue in range(28))
    expected = (
        0,
        25,
        50,
        61,
        44,
        55,
        38,
        49,
        32,
        43,
        26,
        37,
        20,
        24,
        0,
        -24,
        -20,
        -37,
        -26,
        -43,
        -32,
        -49,
        -38,
        -55,
        -44,
        -61,
        -50,
        -25,
    )
    assert coefficients == expected

    # Formula versus the endpoint sweep on a long exact range.
    for n in range(3, 2001):
        sweep_mass = direct_mass(intersection_sweep, n)
        assert sweep_mass == discrepancy_mass(n) == closed_mass(n), n

    # Independent midpoint referee, including every face used by the recursion.
    for n in range(3, 301):
        triple = (1, 2, n)
        for size in (2, 3):
            for face in combinations(triple, size):
                face = tuple(face)
                assert intersection_sweep(face) == intersection_midpoints(face), face
        assert direct_mass(intersection_midpoints, n) == closed_mass(n), n

    # Exact quasiperiodicity and antipodal sign reversal.
    for n in range(3, 2001):
        assert n * closed_mass(n) == (n + 28) * closed_mass(n + 28)
    for residue in range(1, 14):
        assert coefficients[28 - residue] == -coefficients[residue]

    print("EXACT SUPPORT-TRIPLE SAWTOOTH LAW")
    print("M({1,2,N}) = k[N mod 28] / (686*N), N > 2")
    print(f"k[0..27] = {coefficients}")
    print("positive residues: 1..13 mod 28")
    print("zero residues:     0,14 mod 28")
    print("negative residues: 15..27 mod 28")
    print()
    print("COUNTEREXAMPLE FAMILY")
    print("M({1,2,15+28m}) = -12 / (343*(15+28m)) for every m >= 0")
    print(f"m=0: {closed_mass(15)}")
    print(f"m=1: {closed_mass(43)}")
    print(f"m=2: {closed_mass(71)}")
    print()
    print("INDEPENDENT EXACT REFEREE")
    print("endpoint sweep == tooth formula for N=3..2000")
    print("endpoint sweep == midpoint-cell count for every face, N=3..300")
    print("N*M(N) == (N+28)*M(N+28) for N=3..2000")
    print()
    print("TOURNAMENT / ALTERNATE-VERTEX AUDIT")
    print("faithful quotient: signed support hypergraph + 28-residue tooth-boundary walk")
    print("antipodal residue involution r -> 28-r negates the mass")
    print("runner-pair observable is symmetric, so no canonical tournament exists")
    print("score/cycle/SCC/flip/Hamiltonian fingerprints would be tie-break artifacts")
    print("PASS")


if __name__ == "__main__":
    main()
