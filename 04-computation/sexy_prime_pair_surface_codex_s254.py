#!/usr/bin/env python3
"""
sexy_prime_pair_surface_codex_s254.py

Finite scout for the repo question:

    how does the pair-cell / LRC finite-address work relate to proving the
    sexy prime conjecture?

Sexy primes are prime pairs (p, p+6).  In midpoint coordinates

    (p, p+6) = (m-3, m+3),

so they are the fixed half-gap ray h=3 in the prime-pair surface.  This script
records the local residue ledger, Hardy-Littlewood chord, observed counts, and
the exact proof obstruction: local admissibility and singular-series positivity
are easy; infinitude needs a lower-bound sieve/level-of-distribution mechanism
that defeats parity.
"""

from __future__ import annotations

from fractions import Fraction
from math import log

LIMIT = 1_000_000
SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
HALF_GAPS = [1, 2, 3, 5, 6, 10, 15, 21, 30]


def sieve(n: int) -> list[bool]:
    is_prime = [True] * (n + 1)
    if n >= 0:
        is_prime[0] = False
    if n >= 1:
        is_prime[1] = False
    p = 2
    while p * p <= n:
        if is_prime[p]:
            step = p
            start = p * p
            is_prime[start : n + 1 : step] = [False] * (((n - start) // step) + 1)
        p += 1
    return is_prime


def primes_from_bits(bits: list[bool]) -> list[int]:
    return [i for i, ok in enumerate(bits) if ok]


IS_PRIME = sieve(LIMIT + 2 * max(HALF_GAPS) + 10)
PRIMES = primes_from_bits(IS_PRIME)


def distinct_odd_prime_factors(n: int) -> list[int]:
    out: list[int] = []
    m = abs(n)
    while m % 2 == 0 and m:
        m //= 2
    p = 3
    while p * p <= m:
        if m % p == 0:
            out.append(p)
            while m % p == 0:
                m //= p
        p += 2
    if m > 1:
        out.append(m)
    return out


def midpoint_bad_residues(h: int, q: int) -> set[int]:
    if q == 2:
        # For odd prime endpoints m-h and m+h are both odd, so m and h have
        # opposite parity.  This is better treated as a support-spine label.
        return set()
    return {h % q, (-h) % q}


def finite_survival(h: int, moduli: list[int]) -> Fraction:
    prod = Fraction(1)
    for q in moduli:
        if q == 2:
            continue
        prod *= Fraction(q - len(midpoint_bad_residues(h, q)), q)
    return prod


def hardy_littlewood_chord(h: int) -> Fraction:
    chord = Fraction(1)
    for p in distinct_odd_prime_factors(h):
        chord *= Fraction(p - 1, p - 2)
    return chord


def pair_count(h: int, limit: int) -> int:
    gap = 2 * h
    total = 0
    for p in PRIMES:
        q = p + gap
        if q > limit:
            break
        if IS_PRIME[q]:
            total += 1
    return total


def twin_constant(bound: int) -> float:
    prod = 1.0
    for p in PRIMES:
        if p > bound:
            break
        if p > 2:
            prod *= 1.0 - 1.0 / ((p - 1) * (p - 1))
    return prod


def predicted_pairs(h: int, limit: int, c2: float) -> float:
    if limit < 10:
        return 0.0
    return float(hardy_littlewood_chord(h)) * 2.0 * c2 * limit / (log(limit) ** 2)


def print_gap_rows() -> None:
    print("Fixed-gap prime-pair rows through", LIMIT)
    print("Sexy primes are h=3, gap=6.")
    print()
    c2 = twin_constant(200_000)
    twin_count = pair_count(1, LIMIT)
    twin_pred = predicted_pairs(1, LIMIT, c2)
    print(f"C2 finite Euler product through 200000: {c2:.9f}")
    print(f"twin baseline count={twin_count}, predicted={twin_pred:.1f}")
    print()
    print(f"{'h':>4} {'gap':>5} {'odd_factors':>15} {'HL_chord':>10} "
          f"{'count':>8} {'count/twin':>11} {'pred/twin':>10} {'count/pred':>11}")
    for h in HALF_GAPS:
        count = pair_count(h, LIMIT)
        pred = predicted_pairs(h, LIMIT, c2)
        chord = hardy_littlewood_chord(h)
        print(f"{h:4d} {2*h:5d} {str(distinct_odd_prime_factors(h)):>15} "
              f"{str(chord):>10} {count:8d} {count/twin_count:11.3f} "
              f"{pred/twin_pred:10.3f} {count/pred:11.3f}")


def print_local_ledger() -> None:
    print()
    print("Local midpoint-sieve ledger")
    print()
    for h in [1, 3, 15]:
        label = "twin h=1" if h == 1 else ("sexy h=3" if h == 3 else "gap30 h=15")
        parts = []
        for q in SMALL_PRIMES:
            if q == 2:
                continue
            parts.append(f"q={q}:{sorted(midpoint_bad_residues(h, q))}")
        print(f"{label:<12} survival={finite_survival(h, SMALL_PRIMES)}")
        print("  " + "; ".join(parts))
    print()
    print("For h=3 the q=3 obstruction collapses from two bad residues to one:")
    print("  m = +/-3 mod 3 both mean m = 0 mod 3.")
    print("That gives the Hardy-Littlewood chord factor 2 relative to twins.")


def print_proof_ledger() -> None:
    print()
    print("Proof-interface ledger")
    print()
    rows = [
        ("local_admissibility", "closed", "gap 6 has no prime covering obstruction"),
        ("singular_series_chord", "closed", "HL factor relative to twins is 2"),
        ("finite_pair_counts", "evidence", "counts track the fixed-gap chord shape"),
        ("parity_barrier", "open", "upper-bound sieves cannot force both endpoints prime"),
        ("level_of_distribution", "open", "need enough uniformity for a lower-bound pair sieve"),
        ("finite_address_guardrail", "repo-transfer", "retain gap, midpoint residues, squarefree/prime-power sidecars"),
    ]
    for name, status, note in rows:
        print(f"{name:<25} {status:<14} {note}")


def main() -> None:
    print("Sexy prime pair surface scout (codex-S254)")
    print("=" * 72)
    print_gap_rows()
    print_local_ledger()
    print_proof_ledger()


if __name__ == "__main__":
    main()
