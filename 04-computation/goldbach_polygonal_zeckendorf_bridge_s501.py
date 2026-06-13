#!/usr/bin/env python3
"""
goldbach_polygonal_zeckendorf_bridge_s501.py

codex-2026-06-01 S501

Connect four additive-number systems as representation hypergraphs:

* Goldbach: even N should have a prime-pair edge p+q=N.
* Helfgott / weak Goldbach: odd N>5 has a prime-triple hyperedge.
* Hardy-Littlewood: prime hyperedge counts should follow a local singular
  series times an archimedean density.
* Fermat polygonal theorem: s-gonal atoms are an order-s additive basis.
* Zeckendorf: Fibonacci atoms form a unique no-adjacent-carry normal form.

The point is not to prove Goldbach.  The point is to make the shared shape
visible: local admissibility, representation entropy, and whether a theorem
returns abundance, bounded coverage, or a canonical peel.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt, log, prod


LIMIT = 2_000
POLYGONAL_LIMIT = 600


def prime_sieve(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (((limit - start) // p) + 1)
    return [n for n in range(2, limit + 1) if sieve[n]]


PRIMES = prime_sieve(LIMIT)
PRIME_SET = set(PRIMES)


def factor_distinct(n: int) -> tuple[int, ...]:
    factors: list[int] = []
    remaining = n
    for p in PRIMES:
        if p * p > remaining:
            break
        if remaining % p == 0:
            factors.append(p)
            while remaining % p == 0:
                remaining //= p
    if remaining > 1:
        factors.append(remaining)
    return tuple(factors)


def twin_prime_constant_approx(limit: int = LIMIT) -> float:
    value = 1.0
    for p in PRIMES:
        if p > limit:
            break
        if p > 2:
            value *= 1.0 - 1.0 / ((p - 1) ** 2)
    return value


C2_APPROX = twin_prime_constant_approx()


def goldbach_pairs(n: int) -> list[tuple[int, int]]:
    return [(p, n - p) for p in PRIMES if p <= n - p and n - p in PRIME_SET]


def hardy_littlewood_goldbach(n: int) -> float:
    """Hardy-Littlewood binary Goldbach prediction in our unordered convention.

    The common heuristic is 2*C2*n/log(n)^2 times the singular correction over
    odd prime divisors of n.  This is close enough for a finite diagnostic; the
    script uses it only as a shape comparison, not as a theorem.
    """

    correction = 1.0
    for p in factor_distinct(n):
        if p > 2:
            correction *= (p - 1) / (p - 2)
    return 2.0 * C2_APPROX * n * correction / (log(n) ** 2)


def prime_triples(n: int) -> list[tuple[int, int, int]]:
    triples: list[tuple[int, int, int]] = []
    for i, p in enumerate(PRIMES):
        if p > n // 3:
            break
        for q in PRIMES[i:]:
            r = n - p - q
            if r < q:
                break
            if r in PRIME_SET:
                triples.append((p, q, r))
    return triples


def polygonal(sides: int, k: int) -> int:
    return ((sides - 2) * k * k - (sides - 4) * k) // 2


def polygonal_atoms(sides: int, limit: int) -> tuple[int, ...]:
    atoms: list[int] = []
    k = 0
    while True:
        value = polygonal(sides, k)
        if value > limit:
            break
        atoms.append(value)
        k += 1
    return tuple(atoms)


@dataclass(frozen=True)
class PolygonalDepth:
    sides: int
    max_depth: int
    first_witness: int
    missing: tuple[int, ...]


def polygonal_depth(sides: int, limit: int) -> PolygonalDepth:
    atoms = tuple(x for x in polygonal_atoms(sides, limit) if x > 0)
    inf = limit + 1
    dp = [0] + [inf] * limit
    for n in range(1, limit + 1):
        best = inf
        for atom in atoms:
            if atom > n:
                break
            best = min(best, 1 + dp[n - atom])
        dp[n] = best
    missing = tuple(n for n in range(1, limit + 1) if dp[n] > sides)
    max_depth = max(dp[1:])
    first = next(n for n in range(1, limit + 1) if dp[n] == max_depth)
    return PolygonalDepth(sides, max_depth, first, missing[:10])


def fibs_up_to(value: int) -> list[int]:
    fibs = [1, 2]
    while fibs[-1] <= value:
        fibs.append(fibs[-1] + fibs[-2])
    return fibs[:-1]


def zeck(value: int) -> tuple[tuple[int, int], ...]:
    if value <= 0:
        return tuple()
    remaining = value
    out: list[tuple[int, int]] = []
    fibs = fibs_up_to(value)
    for offset in range(len(fibs) - 1, -1, -1):
        fib = fibs[offset]
        if fib <= remaining:
            out.append((offset + 1, fib))
            remaining -= fib
        if remaining == 0:
            break
    return tuple(out)


def zeck_text(value: int) -> str:
    terms = zeck(value)
    if not terms:
        return "0"
    return "+".join(f"F{i}({v})" for i, v in terms)


def zeck_gaps(value: int) -> tuple[int, ...]:
    indices = [i for i, _v in zeck(value)]
    return tuple(indices[i] - indices[i + 1] for i in range(len(indices) - 1))


def print_goldbach_hl() -> None:
    print("Binary Goldbach as a prime-pair edge search")
    print(f"  C2 approximation from primes <= {LIMIT}: {C2_APPROX:.9f}")
    print("  N    pairs  HL-pred  ratio  odd-prime factors  Zeck(N)")
    for n in (100, 250, 500, 1000, 1500, 2000):
        count = len(goldbach_pairs(n))
        hl = hardy_littlewood_goldbach(n)
        ratio = count / hl if hl else 0.0
        factors = ",".join(str(p) for p in factor_distinct(n) if p > 2) or "-"
        print(f"  {n:4d} {count:6d} {hl:8.2f} {ratio:6.3f}  {factors:<17} {zeck_text(n)}")
    print()

    window = range(500, LIMIT + 1, 2)
    sparse = sorted((len(goldbach_pairs(n)), n) for n in window)[:8]
    print("  Sparsest even rows in 500..2000:")
    print("  " + " ".join(f"{n}:{count}" for count, n in sparse))
    print()


def print_ternary_helfgott() -> None:
    print("Ternary Goldbach as the Helfgott-covered hypergraph")
    print("  N    unordered prime triples  example")
    for n in (101, 301, 701, 1001, 1501, 1999):
        triples = prime_triples(n)
        example = "+".join(str(x) for x in triples[0]) if triples else "NONE"
        print(f"  {n:4d} {len(triples):24d}  {example}")
    print()

    for lo, hi in ((7, 101), (101, 501), (501, 1001), (1001, 2000)):
        rows = [(len(prime_triples(n)), n) for n in range(lo | 1, hi + 1, 2) if n > 5]
        count, witness = min(rows)
        print(f"  min triples for odd N in [{lo},{hi}]: {witness} has {count}")
    print()


def print_polygonal_depths() -> None:
    print("Fermat polygonal theorem as bounded additive basis depth")
    print("  sides  max depth <=600  first witness  first missing above side-bound")
    for sides in range(3, 9):
        depth = polygonal_depth(sides, POLYGONAL_LIMIT)
        miss = ",".join(str(n) for n in depth.missing) or "-"
        print(f"  {sides:5d} {depth.max_depth:16d} {depth.first_witness:14d}  {miss}")
    print()


def print_zeckendorf() -> None:
    print("Zeckendorf as a unique no-adjacent-carry cover")
    print("  N     summands  gaps       representation")
    for n in (100, 250, 500, 1000, 1500, 1729, 2000):
        terms = zeck(n)
        gaps = zeck_gaps(n)
        gap_text = ",".join(str(g) for g in gaps) or "-"
        print(f"  {n:4d} {len(terms):9d}  {gap_text:<9} {zeck_text(n)}")
    print()


def print_synthesis() -> None:
    print("Synthesis: one representation-hypergraph grammar")
    rows = [
        (
            "Goldbach",
            "prime atoms",
            "2-uniform edges on even N",
            "conjectural nonempty; HL predicts abundance",
        ),
        (
            "Helfgott",
            "prime atoms",
            "3-uniform hyperedges on odd N",
            "proved nonempty for N>5",
        ),
        (
            "Fermat polygonal",
            "s-gonal atoms",
            "at most s summands",
            "proved bounded basis; local quadratic residues absorbed",
        ),
        (
            "Zeckendorf",
            "Fibonacci atoms",
            "independent support in path",
            "proved unique normal form; entropy forced to 1",
        ),
    ]
    print("  system              atoms             rule                         output")
    for name, atoms, rule, output in rows:
        print(f"  {name:<19} {atoms:<17} {rule:<28} {output}")
    print()
    print("  Deep pattern:")
    print("    Hardy-Littlewood = abundant local-global coverage.")
    print("    Helfgott = enough analytic control to turn ternary abundance into proof.")
    print("    Fermat polygonal = bounded-depth coverage by a structured quadratic basis.")
    print("    Zeckendorf = coverage plus a local no-carry law that kills all redundancy.")
    print("    So the axis is not only additive vs multiplicative; it is")
    print("    redundancy -> bounded basis -> canonical peel.")


def main() -> None:
    print_goldbach_hl()
    print_ternary_helfgott()
    print_polygonal_depths()
    print_zeckendorf()
    print_synthesis()


if __name__ == "__main__":
    main()
