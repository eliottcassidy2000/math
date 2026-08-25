"""Independent finite audit for THM-4056.

This path uses a totient sieve, recursive divisor inversion, and direct clock
arrays rather than the primary script's per-integer arithmetic routines.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import gcd, lcm, prod


def require(flag: bool, payload: object) -> None:
    if not flag:
        raise AssertionError(payload)


def canonical_cf(p: int, q: int) -> tuple[int, ...]:
    out = []
    while q:
        a, r = divmod(p, q)
        out.append(a)
        p, q = q, r
    return tuple(out)


def projective_positive_word(word: tuple[int, ...]) -> tuple[int, ...]:
    return word[1:] if word[0] == 0 else word


def main() -> None:
    n_max = 512
    phi = list(range(n_max + 1))
    for p in range(2, n_max + 1):
        if phi[p] == p:
            for n in range(p, n_max + 1, p):
                phi[n] -= phi[n] // p

    # Direct quotient histogram at four unrelated cutoffs.
    fibre_checks = 0
    for cutoff in (17, 64, 191, 512):
        fibres: dict[tuple[int, int], int] = defaultdict(int)
        for a in range(1, cutoff + 1):
            for b in range(a + 1, cutoff + 1):
                g = gcd(a, b)
                fibres[a // g, b // g] += 1
        require(all(count == cutoff // q for (p, q), count in fibres.items()), cutoff)
        require(len(fibres) == sum(phi[2 : cutoff + 1]), cutoff)
        require(sum(fibres.values()) == cutoff * (cutoff - 1) // 2, cutoff)
        fibre_checks += len(fibres)

    # Build a finite divisor clock by adding each pulse train as an array.
    q_max = 15
    period = 1
    for q in range(1, q_max + 1):
        period = lcm(period, q)
    psi = {q: Fraction(q + 2, q * q + 1) for q in range(2, q_max + 1)}
    weights = {q: 2 * phi[q] * psi[q] for q in psi}
    lift_clock = [Fraction() for _ in range(period + 1)]
    phase_clock = [Fraction() for _ in range(period + 1)]
    for q, weight in weights.items():
        for n in range(q, period + 1, q):
            lift_clock[n] += weight
        for n in range(1, period + 1):
            if gcd(n, q) == 1:
                phase_clock[n] += 2 * psi[q]
    expected_mean = sum((weights[q] / q for q in weights), Fraction())
    require(sum(lift_clock[1:], Fraction()) / period == expected_mean, expected_mean)
    require(sum(phase_clock[1:], Fraction()) / period == expected_mean, expected_mean)
    first_difference = next(n for n in range(1, period + 1) if lift_clock[n] != phase_clock[n])

    phase_2 = [2 if gcd(n, 2) == 1 else 0 for n in range(1, 5)]
    phase_4 = [2 if gcd(n, 4) == 1 else 0 for n in range(1, 5)]
    lift_2 = [2 if n % 2 == 0 else 0 for n in range(1, 5)]
    lift_4 = [4 if n % 4 == 0 else 0 for n in range(1, 5)]
    require(phase_2 == phase_4 and lift_2 != lift_4, (phase_2, phase_4, lift_2, lift_4))

    # Recursive inversion: C(n)=sum_{q|n}b(q), so subtract proper divisors.
    recovered: dict[int, Fraction] = {1: Fraction()}
    for n in range(2, q_max + 1):
        value = lift_clock[n] - sum(
            (recovered[d] for d in range(1, n) if n % d == 0),
            Fraction(),
        )
        recovered[n] = value
        require(value == weights[n], (n, value, weights[n]))

    # Projective Euclidean coefficients are reciprocal-even, but the standard
    # finite Khinchin word always deletes a_0 and therefore is not.
    left = canonical_cf(3, 5)
    right = canonical_cf(5, 3)
    require(left == (0, 1, 1, 2) and right == (1, 1, 2), (left, right))
    require(projective_positive_word(left) == projective_positive_word(right), (left, right))
    require(left[1:] != right[1:], (left, right))
    half = projective_positive_word(canonical_cf(1, 2))
    four_fifths = projective_positive_word(canonical_cf(4, 5))
    require(prod(half) ** len(four_fifths) == prod(four_fifths) ** len(half), (half, four_fifths))

    # New edges at N split into phi(q) types for each q|N.
    incidence_checks = 0
    for n in range(2, n_max + 1):
        histogram: dict[int, int] = defaultdict(int)
        for a in range(1, n):
            histogram[n // gcd(a, n)] += 1
        require(all(n % q == 0 and count == phi[q] for q, count in histogram.items()), n)
        require(sum(histogram.values()) == n - 1, n)
        incidence_checks += len(histogram)

    print("status=INDEPENDENT FINITE-EXACT audit")
    print(f"fibre_types_checked={fibre_checks};incidence_shells_checked={incidence_checks}")
    print(f"clock_qmax={q_max};lcm_period={period};recursive_recoveries={q_max-1}")
    print(f"clock_mean={expected_mean}")
    print(f"phase_lift_pointwise_hostile=first_difference_at_t={first_difference}")
    print("phase_kernel_hostile=q2_and_q4_weights_have_same_phase_clock_but_different_lift_clocks")
    print("cf_normalization_hostile=3/5_vs_5/3;equal_projective_word_but_unequal_standard_Khinchin_word")
    print("PASS")


if __name__ == "__main__":
    main()
