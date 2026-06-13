#!/usr/bin/env python3
"""Cycle formalization, extension, and investigation for operation digraphs.

Mode A, formalization input:
  Lean now proves that the additive one-shadow is strict order and the
  nonunit multiplicative one-shadow is proper divisibility.

Mode B, extension:
  Enumerate minimal product-sum witnesses through a larger arity range using
  the formal seed-defect identity: for nonunit seed F, D=prod(F)-sum(F) ones
  repair the additive fold.

Mode C, investigation:
  Compare binary/divisor-layer witnesses with multi-factor packings, and
  measure how the divisibility DAG sits inside the additive transitive
  completion.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isqrt


MAX_K = 120


@dataclass(frozen=True)
class SeedWitness:
    k: int
    product: int
    seed: tuple[int, ...]
    ones: int

    @property
    def length(self) -> int:
        return len(self.seed)

    def compact(self) -> str:
        seed_text = " ".join(map(str, self.seed))
        if self.ones == 0:
            return seed_text
        if self.ones == 1:
            return f"1 + {seed_text}"
        return f"1^{self.ones} + {seed_text}"


def divisors(n: int) -> list[int]:
    out = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def binary_witnesses(k: int) -> list[SeedWitness]:
    """Two nonunit factors at total arity k.

    With r=k-2 ones and factors a,b>1, the formalized theorem gives
    (a-1)(b-1)=k-1.
    """

    out = []
    target = k - 1
    for d in divisors(target):
        e = target // d
        a, b = d + 1, e + 1
        if a <= b:
            product = a * b
            out.append(SeedWitness(k=k, product=product, seed=(a, b), ones=k - 2))
    return sorted(out, key=lambda w: (w.product, w.seed))


def enumerate_witnesses(max_k: int) -> dict[int, list[SeedWitness]]:
    """Enumerate all seed-defect witnesses that can be minimal for k<=max_k.

    The trivial binary seed `(2,k)` proves m(k) <= 2k, so product <= 2*max_k
    is enough to find all minima up to max_k.
    """

    max_product = 2 * max_k
    by_k: dict[int, list[SeedWitness]] = defaultdict(list)

    def rec(start: int, product: int, total: int, seed: list[int]) -> None:
        if len(seed) >= 2 and product >= total:
            ones = product - total
            k = len(seed) + ones
            if 2 <= k <= max_k:
                by_k[k].append(
                    SeedWitness(k=k, product=product, seed=tuple(seed), ones=ones)
                )

        limit = max_product // product
        for f in range(start, limit + 1):
            rec(f, product * f, total + f, seed + [f])

    rec(2, 1, 0, [])

    for k, values in by_k.items():
        seen = set()
        dedup = []
        for w in sorted(values, key=lambda w: (w.product, w.length, w.seed)):
            key = (w.product, w.seed, w.ones)
            if key not in seen:
                seen.add(key)
                dedup.append(w)
        by_k[k] = dedup
    return dict(by_k)


def divisor_count_sieve(n: int) -> list[int]:
    tau = [0] * (n + 1)
    for d in range(1, n + 1):
        for m in range(d, n + 1, d):
            tau[m] += 1
    return tau


def prime_sieve(n: int) -> list[bool]:
    is_prime = [False, False] + [True] * max(0, n - 1)
    for p in range(2, isqrt(n) + 1):
        if is_prime[p]:
            for m in range(p * p, n + 1, p):
                is_prime[m] = False
    return is_prime


def operation_shadow_stats() -> None:
    print("MODE C FROM MODE A: SHADOW SPARSITY")
    print("=" * 78)
    print("Formal input: additive shadow is total order; product shadow is divisibility.")
    print("Finite check: product_edges_with_unit = sum_{z<=N} (tau(z)-1).")
    print()
    header = (
        "N   additive_edges   product_edges   density     no_unit_edges   "
        "prime_cover_edges"
    )
    print(header)
    print("-" * len(header))
    for n in [25, 50, 100, 250, 500, 1000]:
        tau = divisor_count_sieve(n)
        primes = prime_sieve(n)
        complete = n * (n - 1) // 2
        product_edges = sum(tau[z] - 1 for z in range(1, n + 1))
        no_unit_edges = product_edges - (n - 1)
        prime_covers = sum(1 for x in range(1, n + 1) for p in range(2, n // x + 1) if primes[p])
        print(
            f"{n:>4} {complete:>16} {product_edges:>15} "
            f"{product_edges / complete:>9.5f} {no_unit_edges:>17} "
            f"{prime_covers:>18}"
        )
    print()


def product_sum_atlas() -> dict[int, list[SeedWitness]]:
    print("MODE B FROM MODE A: PRODUCT-SUM ATLAS")
    print("=" * 78)
    print(f"Enumerating all potentially minimal seed-defect witnesses for k <= {MAX_K}.")
    by_k = enumerate_witnesses(MAX_K)
    minima = {k: values[0] for k, values in by_k.items()}

    print()
    header = "k  minP  seed_len  ones  binary_min  gain  compact_min_witness"
    print(header)
    print("-" * len(header))
    for k in range(2, 41):
        best = minima[k]
        bmin = binary_witnesses(k)[0].product
        gain = bmin - best.product
        print(
            f"{k:>2} {best.product:>5} {best.length:>9} {best.ones:>5} "
            f"{bmin:>11} {gain:>5}  {best.compact()}"
        )

    multi_wins = [
        (k, minima[k], binary_witnesses(k)[0])
        for k in range(2, MAX_K + 1)
        if minima[k].length > 2 and minima[k].product < binary_witnesses(k)[0].product
    ]
    print()
    print(f"Multi-factor strict wins over the binary divisor layer: {len(multi_wins)} / {MAX_K-1}")
    print("First 24 wins:")
    for k, best, bbest in multi_wins[:24]:
        print(
            f"  k={k:>3}: P={best.product:<4} seed={best.seed!s:<18} "
            f"beats binary P={bbest.product:<4} seed={bbest.seed}"
        )

    ties = [
        k
        for k in range(2, MAX_K + 1)
        if minima[k].product == binary_witnesses(k)[0].product
        and minima[k].length > 2
    ]
    print()
    print(f"Arity values where a multi-factor seed ties the binary minimum: {ties[:40]}")

    by_len = Counter(minima[k].length for k in range(2, MAX_K + 1))
    by_product = Counter(minima[k].product for k in range(2, MAX_K + 1))
    print()
    print("Minimal witness seed-length distribution:")
    for length, count in sorted(by_len.items()):
        print(f"  seed length {length}: {count}")

    print("Most reused minimal endpoints:")
    for product, count in by_product.most_common(12):
        ks = [k for k in range(2, MAX_K + 1) if minima[k].product == product]
        print(f"  P={product:>4}: {count:>2} arities {ks}")

    print()
    return by_k


def investigate_minima(by_k: dict[int, list[SeedWitness]]) -> None:
    print("MODE C FROM MODE B: INVESTIGATION LEDGER")
    print("=" * 78)
    minima = {k: values[0] for k, values in by_k.items()}

    record_endpoint = 0
    records = []
    for k in range(2, MAX_K + 1):
        best = minima[k]
        if best.product > record_endpoint:
            records.append((k, best))
            record_endpoint = best.product
    print("Record jumps of the minimal endpoint m(k):")
    for k, best in records[:32]:
        print(f"  k={k:>3}: m(k)={best.product:<4} seed={best.seed} ones={best.ones}")

    print()
    print("Local plateaus of m(k):")
    start = 2
    while start <= MAX_K:
        product = minima[start].product
        end = start
        while end + 1 <= MAX_K and minima[end + 1].product == product:
            end += 1
        if end > start:
            print(f"  k={start:>3}..{end:<3}: m(k)={product}")
        start = end + 1

    print()
    print("Critical-pair multiplicity at the minimal endpoint:")
    for k in range(2, 41):
        best_product = minima[k].product
        min_witnesses = [w for w in by_k[k] if w.product == best_product]
        if len(min_witnesses) > 1:
            seeds = ", ".join(str(w.seed) for w in min_witnesses[:8])
            print(f"  k={k:>2}, P={best_product:>3}: {len(min_witnesses)} minimal seeds: {seeds}")

    print()
    print("SYNTHESIS")
    print("-" * 78)
    print(
        "The formal theorems explain why ordinary edges are the wrong level: "
        "addition completes every comparison, while multiplication retains only "
        "divisibility.  The extension shows the product-sum sequence is a "
        "minimal-endpoint problem in that completion.  The investigation then "
        "points back to formalization: prove the binary divisor layer, the "
        "trivial bound m(k)<=2k, and a packing lemma explaining when several "
        "small nonunit factors beat a single binary factor."
    )


def main() -> None:
    operation_shadow_stats()
    by_k = product_sum_atlas()
    investigate_minima(by_k)


if __name__ == "__main__":
    main()
