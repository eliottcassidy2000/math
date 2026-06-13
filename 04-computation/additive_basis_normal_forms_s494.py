#!/usr/bin/env python3
"""
additive_basis_normal_forms_s494.py

codex-2026-06-01 S494

Connect Goldbach/Helfgott/Hardy-Littlewood, Fermat polygonal numbers, and
Zeckendorf as additive-basis normal-form problems.

The shared object is a representation hypergraph:

  atoms A, target n, hyperedges a_1+...+a_r=n.

Goldbach uses prime atoms and asks for positive two-edge mass on even targets.
Helfgott proves positive three-edge mass on odd targets. Hardy-Littlewood
predicts the two-prime edge count as archimedean volume times a p-adic singular
series. Fermat polygonal uses polynomial atoms and bounded arity. Zeckendorf
uses recurrence atoms plus a confluent local carry rule, giving a unique
independent-set normal form.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import log


MAX_N = 10_000
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


PRIMES = sieve(MAX_N)
PRIME_SET = set(PRIMES)


def twin_prime_constant(bound: int = 50_000) -> float:
    primes = sieve(bound)
    prod = 1.0
    for p in primes:
        if p > 2:
            prod *= 1.0 - 1.0 / ((p - 1) * (p - 1))
    return prod


C2_APPROX = twin_prime_constant()


def factor_distinct(n: int) -> list[int]:
    out = []
    m = n
    for p in PRIMES:
        if p * p > m:
            break
        if m % p == 0:
            out.append(p)
            while m % p == 0:
                m //= p
    if m > 1:
        out.append(m)
    return out


def goldbach_pairs(n: int) -> list[tuple[int, int]]:
    pairs = []
    for p in PRIMES:
        if p > n // 2:
            break
        q = n - p
        if q in PRIME_SET:
            pairs.append((p, q))
    return pairs


def hardy_littlewood_unordered(n: int) -> float:
    """Unordered Goldbach pair heuristic.

    The weighted singular series for even n is
      2*C2*prod_{p|n,p>2} (p-1)/(p-2).
    For unordered unweighted pairs, use roughly half of the weighted count
    and divide prime weights by log(n)^2.
    """
    if n <= 2 or n % 2:
        return 0.0
    correction = 1.0
    for p in factor_distinct(n):
        if p > 2:
            correction *= (p - 1) / (p - 2)
    return C2_APPROX * correction * n / (log(n) ** 2)


def ternary_prime_triples(n: int) -> list[tuple[int, int, int]]:
    triples = []
    for i, p in enumerate(PRIMES):
        if p > n // 3:
            break
        for q in PRIMES[i:]:
            if q > (n - p) // 2:
                break
            r = n - p - q
            if r >= q and r in PRIME_SET:
                triples.append((p, q, r))
    return triples


def kgonal(k: int, m: int) -> int:
    return m * ((k - 2) * m - (k - 4)) // 2


def polygonal_values(k: int, limit: int) -> list[int]:
    vals = []
    m = 1
    while True:
        v = kgonal(k, m)
        if v > limit:
            break
        vals.append(v)
        m += 1
    return vals


def polygonal_min_terms(k: int, limit: int) -> tuple[list[int], dict[int, list[int]]]:
    vals = polygonal_values(k, limit)
    inf = 10**9
    dp = [inf] * (limit + 1)
    rep: dict[int, list[int]] = {0: []}
    dp[0] = 0
    for n in range(1, limit + 1):
        best = inf
        best_rep = None
        for v in vals:
            if v > n:
                break
            if dp[n - v] + 1 < best:
                best = dp[n - v] + 1
                best_rep = rep[n - v] + [v]
        dp[n] = best
        if best_rep is not None:
            rep[n] = best_rep
    return dp, rep


def fibs_upto(n: int) -> list[int]:
    fibs = [1, 2]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs


FIBS = fibs_upto(MAX_N)
FIB_INDEX = {v: i + 1 for i, v in enumerate(FIBS)}


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


@dataclass(frozen=True)
class CarryAudit:
    pair: tuple[int, int]
    raw_terms: int
    target_terms: int
    repeats: int
    adjacencies: int
    l1: int
    score: int
    support: tuple[int, ...]
    target: tuple[int, ...]


def carry_audit_for_pair(n: int, pair: tuple[int, int]) -> CarryAudit:
    raw = Counter()
    for summand in pair:
        raw.update(zeckendorf_indices(summand))
    target = Counter(zeckendorf_indices(n))
    all_idx = set(raw) | set(target)
    repeats = sum(max(0, raw[i] - 1) for i in raw)
    positive = sorted(i for i, c in raw.items() if c > 0)
    adj = sum(1 for i in positive if raw.get(i + 1, 0) > 0)
    l1 = sum(abs(raw.get(i, 0) - target.get(i, 0)) for i in all_idx)
    score = l1 + repeats + adj
    return CarryAudit(
        pair=pair,
        raw_terms=sum(raw.values()),
        target_terms=sum(target.values()),
        repeats=repeats,
        adjacencies=adj,
        l1=l1,
        score=score,
        support=tuple(positive),
        target=tuple(sorted(target)),
    )


def best_goldbach_carry(n: int) -> tuple[CarryAudit, CarryAudit]:
    audits = [carry_audit_for_pair(n, pair) for pair in goldbach_pairs(n)]
    best_carry = min(audits, key=lambda a: (a.score, abs(a.pair[1] - a.pair[0]), a.pair))
    best_balance = min(audits, key=lambda a: (abs(a.pair[1] - a.pair[0]), a.score, a.pair))
    return best_carry, best_balance


def print_header(title: str) -> None:
    print()
    print("=" * 112)
    print(title)
    print("=" * 112)


def section_dictionary() -> None:
    print_header("ADDITIVE BASIS NORMAL-FORM DICTIONARY")
    rows = [
        (
            "Goldbach",
            "primes",
            "2 on even n",
            "open positivity",
            "HL singular series predicts count",
        ),
        (
            "Helfgott",
            "primes",
            "3 on odd n",
            "proved positivity",
            "third prime smooths minor arcs",
        ),
        (
            "Fermat polygonal",
            "k-gonal numbers",
            "<= k",
            "proved bounded cover",
            "polynomial atoms with finite arity",
        ),
        (
            "Zeckendorf",
            "Fibonacci numbers",
            "variable",
            "unique normal form",
            "path independence plus local carries",
        ),
    ]
    print(f"{'object':<18} {'atoms':<18} {'arity':<14} {'status':<22} control")
    for row in rows:
        print(f"{row[0]:<18} {row[1]:<18} {row[2]:<14} {row[3]:<22} {row[4]}")
    print()
    print("Shared object: target n is covered by hyperedges a_1+...+a_r=n.")
    print("Deeper axis: local residue/carry constraints decide which hyperedges are legal.")


def section_goldbach_hl() -> None:
    print_header("GOLDBACH PAIRS VS HARDY-LITTLEWOOD LOCAL PRODUCT")
    print(f"Twin-prime constant approximation C2={C2_APPROX:.9f}")
    print(f"{'n':>8} {'pairs':>8} {'HL_unordered':>14} {'ratio':>10} {'prime factors':>20}")
    for n in [50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        pairs = goldbach_pairs(n)
        hl = hardy_littlewood_unordered(n)
        ratio = len(pairs) / hl if hl else 0.0
        print(f"{n:8d} {len(pairs):8d} {hl:14.2f} {ratio:10.3f} {str(factor_distinct(n)):>20}")
    print()
    print("Reading: Hardy-Littlewood is not just 'many primes'.")
    print("It is archimedean area n/log^2(n) times a product of local residue bonuses.")


def section_helfgott() -> None:
    print_header("HELFGOTT TERNARY LIFT: ADD ONE PRIME, GAIN SMOOTHING")
    print(f"{'n':>8} {'triples':>10} first triples")
    for n in [7, 9, 15, 31, 99, 501, 999, 1999]:
        triples = ternary_prime_triples(n)
        print(f"{n:8d} {len(triples):10d} {str(triples[:4])}")
    scan_limit = 501
    counts = []
    for n in range(7, scan_limit + 1, 2):
        counts.append((len(ternary_prime_triples(n)), n))
    fewest = sorted(counts)[:8]
    print()
    print(f"Fewest ternary prime triples for odd 7..{scan_limit}:")
    for count, n in fewest:
        print(f"  n={n:3d}: triples={count}")
    print()
    print("Reading: the third variable turns the prime hypergraph from a thin edge")
    print("problem into a surface with enough averaging for Helfgott's theorem.")


def section_polygonal() -> None:
    print_header("FERMAT POLYGONAL: BOUNDED ADDITIVE COVER BY POLYNOMIAL ATOMS")
    print(f"DP check for n <= {POLY_LIMIT}; repetitions allowed; fewer than k terms allowed.")
    print(f"{'k':>3} {'max min-terms':>14} {'first hits max':<26} sample representation at first hit")
    for k in range(3, 9):
        dp, rep = polygonal_min_terms(k, POLY_LIMIT)
        mx = max(dp[1:])
        hits = [n for n in range(1, POLY_LIMIT + 1) if dp[n] == mx][:8]
        sample = rep[hits[0]] if hits else []
        print(f"{k:3d} {mx:14d} {str(hits):<26} {str(sample)}")
    print()
    print("Reading: Fermat polygonal is a uniform arity theorem.")
    print("The atom set is structured enough that a finite number of summands")
    print("beats all local obstructions.")


def section_zeckendorf_carry() -> None:
    print_header("ZECKENDORF CARRY DEBT ON GOLDBACH PAIRS")
    print("Fibonacci convention: 1,2,3,5,... with indices 1,2,3,4,...")
    print(
        f"{'n':>6} {'pairs':>6} {'target Z':<18} "
        f"{'best carry pair':<17} {'score':>5} {'raw support':<24} "
        f"{'balanced pair':<17} {'score':>5}"
    )
    for n in [50, 100, 200, 500, 1000, 2000, 5000]:
        pairs = goldbach_pairs(n)
        best, balanced = best_goldbach_carry(n)
        print(
            f"{n:6d} {len(pairs):6d} {str(best.target):<18} "
            f"{str(best.pair):<17} {best.score:5d} {str(best.support):<24} "
            f"{str(balanced.pair):<17} {balanced.score:5d}"
        )
    print()
    print("Carry score = L1(raw Fibonacci digits, target Zeckendorf digits)")
    print("              + repeated-index penalties + adjacent-index penalties.")
    print("Reading: Goldbach gives many prime hyperedges; Zeckendorf gives a")
    print("canonical normal form. Carry debt measures how violently a chosen")
    print("prime edge must normalize back to the target integer.")


def section_synthesis() -> None:
    print_header("SYNTHESIS: FOUR LEVELS OF THE SAME ADDITIVE-BASIS MACHINE")
    print(
        """
1. Coverage layer.
   Does the representation hypergraph have at least one edge over each target?
   Goldbach asks this for two prime atoms; Helfgott proves it for three.

2. Density layer.
   Hardy-Littlewood predicts not only existence but the number of prime edges:
   archimedean volume times a p-adic singular series.

3. Bounded-arity layer.
   Fermat polygonal says polynomial atoms become a universal basis once the
   allowed arity is k. This is Waring-like: enough summands erase obstruction.

4. Normal-form layer.
   Zeckendorf is the opposite extreme: it gives uniqueness, not abundance.
   The price is a local carry law, equivalently independence in a path graph.

Repo pattern:
   additive theorem = legal hyperedges + local obstruction product + carry
   normal form + entropy lower bound.

Possible new invariant:
   For any additive basis A, record (representation count, local singular
   product, normal-form carry debt). Goldbach-HL controls the first two;
   Zeckendorf supplies the third; Fermat polygonal is the bounded-cover
   sanity check.
""".strip()
    )


def main() -> None:
    print("Additive basis normal forms (codex-2026-06-01 S494)")
    print("Goldbach / Helfgott / Hardy-Littlewood / Fermat polygonal / Zeckendorf")
    section_dictionary()
    section_goldbach_hl()
    section_helfgott()
    section_polygonal()
    section_zeckendorf_carry()
    section_synthesis()


if __name__ == "__main__":
    main()
