#!/usr/bin/env python3
"""codex-S59: exact Fourier/residue law for the support mass M({a,2a,b}).

Put ``p=1/7``, ``f(x)=1_{||x||<1/14}``, and ``g=f-p``.  For a three-element
support the THM-935/948 subset-Moebius mass is the third centered moment

    M({a,2a,b}) = integral g(a*t) g(2*a*t) g(b*t) dt.

Haar invariance under ``t -> gcd(a,b)*t`` reduces to coprime ``A=a/d`` and
``B=b/d``.  If ``H(x)=g(x)g(2x)``, then, for nonzero n,

    ghat(n) = sin(pi*n/7)/(pi*n),
    Hhat(n) = [sin(pi*n/14) - (1/7)sin(pi*n/7)
               - (2/7)1_{2|n}sin(pi*n/14)]/(pi*n).

The Fourier supports of ``H(A*t)`` and ``g(B*t)`` meet only at ``A*B*k``.
The resulting series is absolutely convergent:

    M({A,2A,B}) = sum_{k != 0} Hhat(B*k) ghat(A*k).

Write

    J(x,y) = B2({(x-y)/2}) - B2({(x+y)/2}),
    B2(u) = u^2-u+1/6,

where braces mean fractional part.  The classical absolutely convergent
Bernoulli identity

    sum_{k != 0} sin(pi*k*x)sin(pi*k*y)/k^2 = pi^2 J(x,y)

then gives ``M=Q(A,B)/(A*B)``, where

    Q = (5/7)J(A/7,B/14) - (1/7)J(A/7,B/7)                 (B even),
    Q = J(A/7,B/14) - (1/14)J(2A/7,B/7)
          - (1/7)J(A/7,B/7)                               (B odd).

Thus ``C(r,s)=686*Q(r,s)`` is an integer depending only on
``r=A mod 14`` and ``s=B mod 28``, and

    M({a,2a,b}) = C(A mod 14, B mod 28)/(686*A*B).

The exact 14-by-28 table proves the strongest zero/sign statement:

    C(r,s)=0  iff  r in {0,7} or s in {0,14};

away from zero, its sign is the product of the signs of the centered residue
representatives of r modulo 14 and s modulo 28.  In particular, after gcd
reduction,

    M({a,2a,b})=0  iff  7 | A or 14 | B.

This script proves the finite residue claims with exact ``Fraction`` arithmetic,
and referees the formula against an exact endpoint sweep on a large coprime box
and an independent midpoint-cell engine on a smaller box.

Tournament analysis / challenged assumption
--------------------------------------------
The faithful combinatorial object is a signed 3-uniform support hyperedge
decorated by the two-frequency residue cell ``(A mod 14, B mod 28)``.  Equivalently,
it is the intersection of the Fourier lattices ``A*Z`` and ``B*Z`` with the
vanishing divisors of ``ghat`` and ``Hhat`` retained.  A runner tournament loses
the symmetric third cumulant, the Fourier indices, and both residue coordinates;
any orientation of its pair data is an arbitrary tie-break.  Hence score, cycle,
SCC, edge-flip, and Hamiltonian-path fingerprints cannot detect this zero law.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, groupby
from math import gcd


P = F(1, 7)


def lcm(values: tuple[int, ...]) -> int:
    result = 1
    for value in values:
        result = result * abs(value) // gcd(result, abs(value))
    return result


@lru_cache(maxsize=None)
def intersection_sweep(support: tuple[int, ...]) -> F:
    """Exact endpoint sweep for the measure of an intersection of bad arcs."""
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
    """Independent exact count on the common endpoint-cell decomposition."""
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


def direct_mass(intersection, a: int, b: int) -> F:
    """The exact support-three Moebius mass for the distinct set {a,2a,b}."""
    support = (a, 2 * a, b)
    assert len(set(support)) == 3
    pair_sum = sum(intersection(tuple(sorted(pair))) for pair in combinations(support, 2))
    return intersection(tuple(sorted(support))) - P * pair_sum + 2 * P**3


def fractional_part(x: F) -> F:
    return x - x.numerator // x.denominator


def periodic_bernoulli_two(x: F) -> F:
    x = fractional_part(x)
    return x * x - x + F(1, 6)


def sine_pair_sum(x: F, y: F) -> F:
    """The exact value of pi^-2 sum_{k != 0} sin(pi*k*x)sin(pi*k*y)/k^2."""
    return periodic_bernoulli_two((x - y) / 2) - periodic_bernoulli_two((x + y) / 2)


def normalized_fourier_mass(a_residue: int, b_residue: int) -> F:
    """Q with M=Q/(A*B); inputs may be any representatives of the residues."""
    x = F(a_residue, 7)
    y = F(b_residue, 14)
    if b_residue % 2 == 0:
        return F(5, 7) * sine_pair_sum(x, y) - F(1, 7) * sine_pair_sum(x, 2 * y)
    return (
        sine_pair_sum(x, y)
        - F(1, 14) * sine_pair_sum(2 * x, 2 * y)
        - F(1, 7) * sine_pair_sum(x, 2 * y)
    )


def residue_coefficient(a: int, b: int) -> int:
    """C(A mod 14,B mod 28)=686*Q, proved integral by the residue audit."""
    value = 686 * normalized_fourier_mass(a % 14, b % 28)
    assert value.denominator == 1
    return value.numerator


def closed_mass(a: int, b: int) -> F:
    """Closed formula, automatically reducing a common gcd."""
    divisor = gcd(a, b)
    A, B = a // divisor, b // divisor
    return F(residue_coefficient(A, B), 686 * A * B)


def centered_residue_sign(residue: int, modulus: int) -> int:
    residue %= modulus
    if residue in (0, modulus // 2):
        return 0
    return 1 if residue < modulus // 2 else -1


def residue_table() -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(residue_coefficient(a, b) for b in range(28)) for a in range(14))


def main() -> None:
    table = residue_table()

    # Exact finite proof of the zero locus and the complete sign law.
    for a_residue in range(14):
        for b_residue in range(28):
            coefficient = table[a_residue][b_residue]
            expected_zero = a_residue in (0, 7) or b_residue in (0, 14)
            assert (coefficient == 0) == expected_zero
            if not expected_zero:
                expected_sign = centered_residue_sign(a_residue, 14) * centered_residue_sign(
                    b_residue, 28
                )
                assert (coefficient > 0) - (coefficient < 0) == expected_sign

    # Formula versus the endpoint sweep on 5,000+ distinct coprime triples.
    sweep_checks = 0
    for A in range(1, 65):
        for B in range(1, 129):
            if gcd(A, B) != 1 or B in (A, 2 * A):
                continue
            assert direct_mass(intersection_sweep, A, B) == closed_mass(A, B)
            sweep_checks += 1

    # Independent midpoint referee for all faces in a smaller complete box.
    midpoint_checks = 0
    for A in range(1, 13):
        for B in range(1, 37):
            if gcd(A, B) != 1 or B in (A, 2 * A):
                continue
            support = (A, 2 * A, B)
            for size in (2, 3):
                for face in combinations(support, size):
                    face = tuple(sorted(face))
                    assert intersection_sweep(face) == intersection_midpoints(face)
            assert direct_mass(intersection_midpoints, A, B) == closed_mass(A, B)
            midpoint_checks += 1

    # Haar scaling invariance checked independently on representative dilates.
    dilation_checks = 0
    for A, B in ((1, 15), (3, 10), (5, 14), (7, 2), (8, 13), (11, 28)):
        assert gcd(A, B) == 1 and B not in (A, 2 * A)
        for dilation in range(1, 10):
            assert direct_mass(intersection_sweep, dilation * A, dilation * B) == closed_mass(A, B)
            dilation_checks += 1

    # The requested branch and its converse, sampled well beyond one residue period.
    for A in range(1, 200):
        for B in range(1, 300):
            if gcd(A, B) != 1 or B in (A, 2 * A):
                continue
            mass = closed_mass(A, B)
            assert (mass == 0) == (A % 7 == 0 or B % 14 == 0)

    print("EXACT SUPPORT MASS FOR {a,2a,b} (codex-S59)")
    print("d=gcd(a,b), A=a/d, B=b/d, gcd(A,B)=1")
    print("M({a,2a,b}) = C[A mod 14][B mod 28] / (686*A*B)")
    print("C=0 iff 7|A or 14|B")
    print("sign(C)=centered_sign_14(A)*centered_sign_28(B) away from zero")
    print()
    print("C RESIDUE TABLE (rows A mod 14; columns B mod 28)")
    for row_index, row in enumerate(table):
        print(f"{row_index:2}: {row}")
    print()
    print("EXACT REFEREES")
    print(f"Fourier/Bernoulli formula == endpoint sweep: {sweep_checks} coprime triples")
    print(f"endpoint sweep == midpoint cells on every face: {midpoint_checks} coprime triples")
    print(f"common-dilation invariance checks: {dilation_checks}")
    print("all 14*28 residue cells satisfy the exact zero and sign laws")
    print()
    print("TOURNAMENT / ALTERNATE-VERTEX AUDIT")
    print("faithful object: signed support hyperedge + (A mod 14,B mod 28) Fourier cell")
    print("Fourier-lattice zeros: 7|A kills ghat(Ak); 14|B kills Hhat(Bk)")
    print("runner-pair tournaments discard the third cumulant and both residue labels")
    print("score/cycle/SCC/flip/Hamiltonian fingerprints are therefore non-diagnostic")
    print("PASS")


if __name__ == "__main__":
    main()
