#!/usr/bin/env python3
"""Exact local-residue scout for Bagshaw's five 5-rank-five specializations.

This does not prove a sixth class.  It tests whether the five record points
carry a cheap coordinate lost by an unstratified height-box search: reduction
onto a bad fibre of the Schoof polynomial at the conductor primes 47 and 103.
"""

from fractions import Fraction
from math import comb, gcd

from flint import fmpz


H6_COEFFICIENTS = (47, 21, 598, 1561, 1198, 261, 47)
SEEDS = (
    (477, 1694, 23454009318604054148884180799),
    (1464, 2567, 3363903161767512730048342385039),
    (846, 1733, 95047361962872405516772909439),
    (3206, 5, 525774905878035876898636843167),
    (3419, 994, 3881041670601462517919585588159),
)
BAD_PRIMES = (47, 103)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def quadratic_factor(x: int, y: int) -> int:
    return x * x + x * y + y * y


def sextic_factor(x: int, y: int) -> int:
    return sum(
        coefficient * x ** (6 - index) * y ** index
        for index, coefficient in enumerate(H6_COEFFICIENTS)
    )


def homogeneous_value(x: int, y: int) -> int:
    return quadratic_factor(x, y) * sextic_factor(x, y)


def fundamental_discriminant_abs(numerator: int, denominator: int) -> tuple[int, int]:
    value = homogeneous_value(numerator, denominator)
    require(value > 0, "the five source specializations should be imaginary")
    squarefree_part = 1
    square_factor = 1
    for prime, exponent in fmpz(value).factor():
        squarefree_part *= int(prime) ** (exponent % 2)
        square_factor *= int(prime) ** (exponent // 2)
    negative_radicand = -squarefree_part
    discriminant = (
        negative_radicand
        if negative_radicand % 4 == 1
        else 4 * negative_radicand
    )
    return -discriminant, square_factor


def projective_point(numerator: int, denominator: int, prime: int) -> tuple[int, int]:
    if denominator % prime:
        return numerator * pow(denominator, -1, prime) % prime, 1
    require(numerator % prime != 0, "source fraction must remain primitive")
    return 1, 0


def projective_roots(prime: int, factor: str) -> tuple[str, ...]:
    points = [(value, 1) for value in range(prime)] + [(1, 0)]
    roots = []
    for x, y in points:
        value = (
            quadratic_factor(x, y)
            if factor == "quadratic"
            else sextic_factor(x, y)
        )
        if value % prime == 0:
            roots.append("infinity" if y == 0 else str(x))
    return tuple(roots)


def label_at(numerator: int, denominator: int, prime: int) -> str:
    x, y = projective_point(numerator, denominator, prime)
    labels = []
    if quadratic_factor(x, y) % prime == 0:
        labels.append("quadratic")
    if sextic_factor(x, y) % prime == 0:
        labels.append("sextic")
    return "+".join(labels) if labels else "good"


def main() -> None:
    print("SCHOOF 5-RANK-FIVE BAD-FIBRE SCOUT -- 2026-08-24")
    print("status=FINITE-EXACT_SCOUT;no_rank_six_claim")

    projective_hit_fractions = []
    for prime in BAD_PRIMES:
        quadratic_roots = projective_roots(prime, "quadratic")
        sextic_roots = projective_roots(prime, "sextic")
        require(set(quadratic_roots).isdisjoint(sextic_roots),
                f"factors share a root modulo {prime}")
        hit_count = len(quadratic_roots) + len(sextic_roots)
        hit_fraction = Fraction(hit_count, prime + 1)
        projective_hit_fractions.append(hit_fraction)
        print(
            f"prime={prime};quadratic_roots={quadratic_roots};"
            f"sextic_roots={sextic_roots};"
            f"projective_hit_fraction={hit_fraction}"
        )

    rows = []
    for numerator, denominator, expected_discriminant in SEEDS:
        require(gcd(numerator, denominator) == 1, "source t is not reduced")
        discriminant, square_factor = fundamental_discriminant_abs(
            numerator, denominator
        )
        require(discriminant == expected_discriminant,
                f"source discriminant mismatch at {numerator}/{denominator}")
        labels = tuple(
            (prime, label_at(numerator, denominator, prime))
            for prime in BAD_PRIMES
        )
        hit = any(label != "good" for _, label in labels)
        rows.append(hit)
        print(
            f"seed={numerator}/{denominator};D={discriminant};"
            f"square_factor={square_factor};labels={labels};bad_hit={hit};"
            f"factorization={fmpz(discriminant).factor()}"
        )

    require(sum(rows) == 4, "the source record should have four bad-fibre hits")
    first, second = projective_hit_fractions
    union_density = first + second - first * second
    tail = sum(
        Fraction(comb(5, count))
        * union_density ** count
        * (1 - union_density) ** (5 - count)
        for count in range(4, 6)
    )
    require(union_density == Fraction(5, 26), "bad-fibre density changed")
    print(f"crt_uniform_union_density={union_density}")
    print(f"formal_binomial_tail_at_least_4_of_5={tail}={float(tail):.12f}")
    print(
        "interpretation=post_selected_signal_only;"
        "stratify_future_search_by_47_103_bad_fibres_and_retest_out_of_sample"
    )
    print(
        "boundary=no_independence_gain_no_sixth_class_no_p_value;"
        "the_fifth_seed_is_a_good_fibre_hostile"
    )


if __name__ == "__main__":
    main()
