#!/usr/bin/env python3
"""S653: Basel-type identities as a pi carrier family.

The Basel problem is the first member of two coupled infinite families:

1. power sums:
      zeta(2k) = rational * pi^(2k)

2. elementary/disjoint sums:
      zeta({2}^m) = sum_{n1>...>nm>=1} 1/(n1^2 ... nm^2)
                  = pi^(2m)/(2m+1)!

The sine product is the clean carrier.  Its coefficients are the elementary
symmetric sums of {1/n^2}; the log derivative turns them into the power sums
zeta(2k).  Newton identities are the exact algebraic bridge.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, factorial


def bernoulli_numbers(n: int) -> list[Fraction]:
    """Return B_0..B_n with the convention B_1=-1/2."""

    b = [Fraction(0) for _ in range(n + 1)]
    b[0] = Fraction(1)
    for m in range(1, n + 1):
        s = sum(Fraction(comb(m + 1, j)) * b[j] for j in range(m))
        b[m] = -s / Fraction(m + 1)
    return b


def primes_upto(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            for q in range(p * p, n + 1, p):
                sieve[q] = False
    return [i for i, ok in enumerate(sieve) if ok]


def von_staudt_primes(weight: int) -> list[int]:
    return [p for p in primes_upto(weight + 1) if weight % (p - 1) == 0]


def euler_zeta_coeff(k: int, bernoulli: list[Fraction]) -> Fraction:
    """Coefficient c_k in zeta(2k)=c_k*pi^(2k)."""

    sign = 1 if (k + 1) % 2 == 0 else -1
    return Fraction(sign) * bernoulli[2 * k] * Fraction(2 ** (2 * k - 1), factorial(2 * k))


def newton_power_coeffs(max_k: int) -> list[Fraction]:
    """Recover zeta(2k)/pi^(2k) from e_m=1/(2m+1)!."""

    e = [Fraction(0)] + [Fraction(1, factorial(2 * m + 1)) for m in range(1, max_k + 1)]
    p = [Fraction(0) for _ in range(max_k + 1)]
    for k in range(1, max_k + 1):
        total = Fraction(0)
        for i in range(1, k):
            total += Fraction((-1) ** (i - 1)) * e[i] * p[k - i]
        total += Fraction((-1) ** (k + 1)) * k * e[k]
        p[k] = total
    return p


def beats(scores: list[tuple[int, int, int]], i: int, j: int) -> bool:
    return scores[i] >= scores[j]


def tournament_fingerprints(labels: list[str], scores: list[tuple[int, int, int]]) -> dict[str, object]:
    n = len(labels)
    out = [0] * n
    edges = {}
    for i in range(n):
        for j in range(i + 1, n):
            if beats(scores, i, j):
                out[i] += 1
                edges[(i, j)] = i
            else:
                out[j] += 1
                edges[(i, j)] = j

    c3 = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                winners = {
                    (a, b): edges[(a, b)],
                    (a, c): edges[(a, c)],
                    (b, c): edges[(b, c)],
                }
                wins = Counter(winners.values())
                if sorted(wins.values()) == [1, 1, 1]:
                    c3 += 1

    # Hamiltonian path count by DP.
    def arc(i: int, j: int) -> bool:
        if i < j:
            return edges[(i, j)] == i
        return edges[(j, i)] == i

    dp: dict[tuple[int, int], int] = {}
    for last in range(n):
        dp[(1 << last, last)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if arc(last, nxt):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    h = sum(dp.get((full, last), 0) for last in range(n))
    return {
        "score_hist": dict(sorted(Counter(out).items())),
        "directed_3cycles": c3,
        "hamiltonian_paths": h,
        "top_order": [labels[i] for i in sorted(range(n), key=lambda i: out[i], reverse=True)],
    }


def main() -> None:
    max_k = 12
    bern = bernoulli_numbers(2 * max_k)
    newton = newton_power_coeffs(max_k)

    print("=" * 78)
    print("S653 Basel-type pi identities: power sums, disjoint sums, Bernoulli sieve")
    print("=" * 78)
    print()
    print("Core carrier identities")
    print("  sin(pi*x)/(pi*x) = product_{n>=1} (1 - x^2/n^2)")
    print("  zeta({2}^m)      = pi^(2m)/(2m+1)!        [elementary/disjoint face]")
    print("  zeta(2k)         = c_k * pi^(2k)           [power-sum/moment face]")
    print("  Newton identities recover c_k from e_m=1/(2m+1)!.")
    print()

    print("Table: even-zeta coefficients and von Staudt sieve")
    print(f"{'k':>2} {'zeta(2k)/pi^(2k)':>22} {'Newton':>22} {'B_2k':>16} {'VS primes':>24}")
    for k in range(1, max_k + 1):
        coeff = euler_zeta_coeff(k, bern)
        vs = von_staudt_primes(2 * k)
        marker = "OK" if coeff == newton[k] else "BAD"
        print(f"{k:>2} {str(coeff):>22} {str(newton[k]):>22} {str(bern[2*k]):>16} {str(vs):>24} {marker}")

    print()
    print("Elementary/disjoint Basel family")
    print(f"{'m':>2} {'zeta({2}^m)/pi^(2m)':>28} {'denominator':>20}")
    for m in range(1, max_k + 1):
        coeff = Fraction(1, factorial(2 * m + 1))
        print(f"{m:>2} {str(coeff):>28} {factorial(2*m+1):>20,}")

    print()
    print("Von Staudt denominator chain")
    x = 6
    seen = set()
    for step in range(8):
        vs = von_staudt_primes(x)
        nxt = 1
        for p in vs:
            nxt *= p
        fixed = " fixed" if nxt == x else ""
        print(f"  step {step}: weight={x:<6} primes={vs!s:<24} denom={nxt}{fixed}")
        if x in seen:
            break
        seen.add(x)
        x = nxt

    print()
    print("Carrier interpretation")
    print("  elementary face = disjoint packet counts, like OCF collections of disjoint odd cycles")
    print("  power-sum face  = moments/log-derivative, like scalar H or aggregate spectra")
    print("  Bernoulli face  = rational descent and p-adic/von-Staudt side-channel")
    print("  pi is the period that makes the elementary and power-sum faces agree.")

    labels = [
        "sine_product_elementary_packets",
        "newton_log_derivative_power_sums",
        "bernoulli_von_staudt_sieve",
        "fourier_parseval_boundary",
        "dirichlet_L_modular_periods",
        "p_adic_zeta_denominator_audit",
        "raw_numeric_pslq",
    ]
    scores = [
        (10, 10, 10),
        (10, 9, 10),
        (9, 10, 8),
        (8, 8, 7),
        (7, 9, 7),
        (7, 8, 6),
        (3, 3, 2),
    ]
    fp = tournament_fingerprints(labels, scores)
    print()
    print("Tournament Analysis over proof lenses")
    print(f"  top_order={fp['top_order']}")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")

    print()
    print("Hypothesis")
    print("  Treat Basel identities as a period-carrier duality:")
    print("  disjoint elementary packets <-> power-sum moments <-> Bernoulli/p-adic sieve.")
    print("  This is the analytic cousin of the repo's OCF identity H(T)=I(Omega(T),2).")


if __name__ == "__main__":
    main()
