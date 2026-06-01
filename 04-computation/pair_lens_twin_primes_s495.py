#!/usr/bin/env python3
"""
pair_lens_twin_primes_s495.py

codex-2026-06-01 S495

See the additive-basis picture in terms of pairs.

Prime pairs are points (m,h) with entries m-h and m+h.  Goldbach fixes the
midpoint m=N/2 and varies h.  Twin primes fix h=1 and vary m.  Hardy-Littlewood
is the local product attached to these pair slices.  Zeckendorf carry debt is
a normal-form cost attached to the same pair edge.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import log


MAX_N = 200_000
POLY_LIMIT = 300


def sieve(n: int) -> list[int]:
    is_prime = [True] * (n + 1)
    if n >= 0:
        is_prime[0] = False
    if n >= 1:
        is_prime[1] = False
    p = 2
    while p * p <= n:
        if is_prime[p]:
            for q in range(p * p, n + 1, p):
                is_prime[q] = False
        p += 1
    return [i for i, ok in enumerate(is_prime) if ok]


PRIMES = sieve(MAX_N + 1000)
PRIME_SET = set(PRIMES)


def distinct_odd_prime_factors(n: int) -> list[int]:
    out = []
    m = n
    while m % 2 == 0 and m > 0:
        m //= 2
    d = 3
    while d * d <= m:
        if m % d == 0:
            out.append(d)
            while m % d == 0:
                m //= d
        d += 2
    if m > 1:
        out.append(m)
    return out


def twin_constant(bound: int = 100_000) -> float:
    prod = 1.0
    for p in sieve(bound):
        if p > 2:
            prod *= 1.0 - 1.0 / ((p - 1) * (p - 1))
    return prod


C2 = twin_constant()


def pair_chord(h: int) -> Fraction:
    """Hardy-Littlewood correction for gap 2h relative to twin gap 2."""
    corr = Fraction(1)
    for p in distinct_odd_prime_factors(h):
        corr *= Fraction(p - 1, p - 2)
    return corr


def prime_pairs_with_half_gap(h: int, limit: int) -> list[tuple[int, int]]:
    gap = 2 * h
    out = []
    for p in PRIMES:
        q = p + gap
        if q > limit:
            break
        if q in PRIME_SET:
            out.append((p, q))
    return out


def goldbach_pairs(n: int) -> list[tuple[int, int]]:
    out = []
    for p in PRIMES:
        if p > n // 2:
            break
        q = n - p
        if q in PRIME_SET:
            out.append((p, q))
    return out


def twin_prediction(limit: int) -> float:
    if limit < 3:
        return 0.0
    return 2.0 * C2 * limit / (log(limit) ** 2)


def pair_prediction(h: int, limit: int) -> float:
    return float(pair_chord(h)) * twin_prediction(limit)


def midpoint_forbidden_residue_count(h: int, q: int) -> int:
    residues = {h % q, (-h) % q}
    return len(residues)


def midpoint_survival_product(h: int, primes: list[int]) -> Fraction:
    """Finite local survival over odd primes in midpoint sieve."""
    prod = Fraction(1)
    for q in primes:
        if q == 2:
            continue
        prod *= Fraction(q - midpoint_forbidden_residue_count(h, q), q)
    return prod


def fibs_upto(n: int) -> list[int]:
    fibs = [1, 2]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs


FIBS = fibs_upto(2 * MAX_N + 10)


def zeckendorf_indices(n: int) -> tuple[int, ...]:
    out = []
    rem = n
    for i in range(len(FIBS) - 1, -1, -1):
        f = FIBS[i]
        if f <= rem:
            out.append(i + 1)
            rem -= f
        if rem == 0:
            break
    return tuple(sorted(out))


def pair_carry_score(a: int, b: int) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    raw = Counter(zeckendorf_indices(a))
    raw.update(zeckendorf_indices(b))
    target = Counter(zeckendorf_indices(a + b))
    all_idx = set(raw) | set(target)
    repeats = sum(max(0, c - 1) for c in raw.values())
    support = sorted(i for i, c in raw.items() if c > 0)
    adjacent = sum(1 for i in support if raw.get(i + 1, 0) > 0)
    l1 = sum(abs(raw.get(i, 0) - target.get(i, 0)) for i in all_idx)
    return l1 + repeats + adjacent, tuple(support), tuple(sorted(target))


def kgonal(k: int, n: int) -> int:
    return n * ((k - 2) * n - (k - 4)) // 2


def kgonal_values(k: int, limit: int) -> list[int]:
    vals = [0]
    n = 1
    while True:
        v = kgonal(k, n)
        if v > limit:
            break
        vals.append(v)
        n += 1
    return vals


def polygonal_pair_coverage(k: int, limit: int) -> tuple[int, int, list[int]]:
    vals = kgonal_values(k, limit)
    pair_sums = {a + b for a in vals for b in vals if a + b <= limit}
    misses = [n for n in range(1, limit + 1) if n not in pair_sums]
    return limit - len(misses), len(misses), misses[:10]


def header(title: str) -> None:
    print()
    print("=" * 112)
    print(title)
    print("=" * 112)


def section_pair_dictionary() -> None:
    header("PAIR LENS DICTIONARY")
    print("Every unordered pair (a,b) can be written as (m-h,m+h).")
    print()
    print(f"{'object':<20} {'pair slice':<28} {'meaning'}")
    rows = [
        ("Twin primes", "h=1, entries prime", "fixed-gap prime-pair ray"),
        ("Prime gap 2h", "h fixed", "same ray with local chord factor"),
        ("Goldbach N", "m=N/2 fixed", "fixed-sum prime-pair row"),
        ("Polygonal pair", "a,b k-gonal", "two-term bounded-cover layer"),
        ("Zeckendorf pair", "digit supports", "carry debt to canonical sum"),
        ("Tournament edge", "{i,j}", "the primitive pairwise relation"),
    ]
    for name, pair_slice, meaning in rows:
        print(f"{name:<20} {pair_slice:<28} {meaning}")


def section_twin_counts() -> None:
    header("TWIN PRIME COUNTS AS THE BASELINE PAIR RAY")
    print(f"C2 approximation from finite Euler product: {C2:.9f}")
    print(f"{'limit':>8} {'twins':>8} {'2*C2*x/log^2x':>18} {'ratio':>9} {'first/last'}")
    for limit in [100, 1_000, 10_000, 100_000, 200_000]:
        twins = prime_pairs_with_half_gap(1, limit)
        pred = twin_prediction(limit)
        ratio = len(twins) / pred if pred else 0.0
        edge = f"{twins[:2]} ... {twins[-2:]}" if len(twins) >= 4 else str(twins)
        print(f"{limit:8d} {len(twins):8d} {pred:18.2f} {ratio:9.3f} {edge}")
    print()
    print("Twin primes are the unboosted fixed-gap ray: h=1 has no odd prime factors.")


def section_gap_chords() -> None:
    header("PAIR CHORDS: GAP 2h RELATIVE TO TWINS")
    print(f"All prime pairs p,p+2h with p+2h <= {MAX_N}; not just consecutive prime gaps.")
    base = len(prime_pairs_with_half_gap(1, MAX_N))
    print(f"{'h':>4} {'gap':>5} {'odd p|h':>14} {'HL chord':>10} {'pairs':>8} {'ratio/twin':>12} {'predicted':>10}")
    for h in [1, 2, 3, 4, 5, 6, 7, 10, 12, 15, 21, 30]:
        pairs = prime_pairs_with_half_gap(h, MAX_N)
        chord = pair_chord(h)
        ratio = len(pairs) / base if base else 0.0
        print(
            f"{h:4d} {2*h:5d} {str(distinct_odd_prime_factors(h)):>14} "
            f"{str(chord):>10} {len(pairs):8d} {ratio:12.3f} {float(chord):10.3f}"
        )
    print()
    print("The pair product is local: for odd q, midpoint residues m=+/-h mod q are forbidden.")
    print("If q divides h, two forbidden residues collapse to one, creating the HL boost.")


def section_midpoint_sieve() -> None:
    header("MIDPOINT SIEVE: TWIN PRIMES ARE SURVIVORS OF +/-1")
    small_primes = [3, 5, 7, 11, 13]
    print(f"{'h':>4} {'forbidden residues m mod q':<62} {'finite survival'}")
    for h in [1, 2, 3, 5, 6, 15]:
        pieces = []
        for q in small_primes:
            residues = sorted({h % q, (-h) % q})
            pieces.append(f"q={q}:{residues}")
        surv = midpoint_survival_product(h, small_primes)
        print(f"{h:4d} {'; '.join(pieces):<62} {surv}")
    print()
    print("For twin primes beyond (3,5), h=1 forces m divisible by 2 and 3,")
    print("so the midpoint lives on the 6Z spine: (6k-1, 6k+1).")


def section_goldbach_orthogonal() -> None:
    header("GOLDBACH IS THE ORTHOGONAL PAIR SLICE")
    print(f"{'N':>7} {'m=N/2':>8} {'pairs':>7} {'h range':>18} {'best carry pair':>18} {'score'}")
    for n in [100, 200, 500, 1000, 2000, 5000, 10000]:
        pairs = goldbach_pairs(n)
        hs = [(b - a) // 2 for a, b in pairs]
        audits = [(pair_carry_score(a, b)[0], abs(b - a), (a, b)) for a, b in pairs]
        best = min(audits) if audits else (None, None, None)
        h_range = f"{min(hs)}..{max(hs)}" if hs else "-"
        print(f"{n:7d} {n//2:8d} {len(pairs):7d} {h_range:>18} {str(best[2]):>18} {best[0]}")
    print()
    print("Twin primes fix h and move in m. Goldbach fixes m and scans h.")
    print("They are perpendicular slices through the same prime-pair lattice.")


def section_twin_carry() -> None:
    header("ZECKENDORF CARRY DEBT ON TWIN PAIRS")
    twins = prime_pairs_with_half_gap(1, 50_000)
    scores = []
    zeroes = []
    for a, b in twins:
        score, support, target = pair_carry_score(a, b)
        scores.append(score)
        if score == 0 and len(zeroes) < 12:
            zeroes.append((a, b, a + b, target))
    hist = Counter(scores)
    print(f"Twin pairs scanned: {len(twins)} up to 50000")
    print("Carry-score histogram, first bins:")
    for score, count in sorted(hist.items())[:12]:
        print(f"  score {score:2d}: {count:4d}")
    print()
    print("First zero-carry twin pairs:")
    if zeroes:
        for a, b, s, target in zeroes:
            print(f"  ({a},{b}) sum={s} targetZ={target}")
    else:
        print("  none found in this scan")
    print()
    print("Zero carry means the two prime supports already union to the")
    print("canonical Zeckendorf support of their sum.  Twin-ness and")
    print("normal-form compatibility are independent pair filters.")


def section_polygonal_pairs() -> None:
    header("FERMAT POLYGONAL BEGINS WITH PAIRS BUT NEEDS HIGHER ARITY")
    print(f"Pair coverage by two k-gonal atoms for targets 1..{POLY_LIMIT}.")
    print(f"{'k':>3} {'pair-covered':>13} {'misses':>8} {'first misses'}")
    for k in range(3, 9):
        covered, misses, first = polygonal_pair_coverage(k, POLY_LIMIT)
        print(f"{k:3d} {covered:13d} {misses:8d} {first}")
    print()
    print("Fermat polygonal says the pair layer is only the first floor:")
    print("at most k summands give universal coverage.")


def section_synthesis() -> None:
    header("SYNTHESIS")
    print(
        """
Pair-first dictionary:
  - A tournament is a complete decision on every pair.
  - Goldbach is prime-pair coverage along fixed midpoint rows.
  - Twin primes are prime-pair persistence along the h=1 ray.
  - Hardy-Littlewood is the local residue product of a pair slice.
  - Zeckendorf carry is the normal-form debt of a pair edge.
  - Fermat polygonal says pairs are not always enough; higher arity repairs.

New slogan:
  pair = smallest relational molecule.

Twin primes are the cleanest test object because no summand choice remains:
the pair is forced to be (m-1,m+1).  The entire problem is whether infinitely
many midpoints survive all local pair sieves.

Repo target:
  Turn every additive or LRC repair search into a pair movie in (midpoint,
  half-gap, local product, carry debt), then use Tournament Analysis to orient
  competing pairs by which debt they reduce.
""".strip()
    )


def main() -> None:
    print("Pair lens and twin primes (codex-2026-06-01 S495)")
    section_pair_dictionary()
    section_twin_counts()
    section_gap_chords()
    section_midpoint_sieve()
    section_goldbach_orthogonal()
    section_twin_carry()
    section_polygonal_pairs()
    section_synthesis()


if __name__ == "__main__":
    main()
