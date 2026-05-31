#!/usr/bin/env python3
"""Natural-number operation digraphs and product-sum critical pairs.

The additive two-input gate has hyperedges {x,y} -> x+y.  Its simple
one-input shadow on [N] has an edge x -> z exactly when x < z, since
z = x + (z-x).  Thus the shadow is the transitive tournament.

The multiplicative two-input gate has hyperedges {x,y} -> xy.  Its simple
shadow has x -> z exactly when x divides z, so it is the divisibility DAG.

Product-sum equations are the critical pairs where an additive fold and a
multiplicative fold reach the same endpoint.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import prod


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def additive_edges(n: int) -> set[tuple[int, int]]:
    return {(x, z) for z in range(2, n + 1) for x in range(1, z)}


def product_edges(n: int, include_unit: bool = True) -> set[tuple[int, int]]:
    lo = 1 if include_unit else 2
    return {
        (x, z)
        for z in range(2, n + 1)
        for x in range(lo, z)
        if z % x == 0
    }


def product_cover_edges(n: int) -> set[tuple[int, int]]:
    return {
        (x, x * p)
        for x in range(1, n + 1)
        for p in range(2, n // x + 1)
        if is_prime(p)
    }


def additive_gate_count(n: int) -> int:
    # Unordered positive pairs x <= y with x+y <= n.
    return sum(1 for x in range(1, n + 1) for y in range(x, n + 1) if x + y <= n)


def product_gate_count(n: int, include_unit: bool = True) -> int:
    lo = 1 if include_unit else 2
    return sum(1 for x in range(lo, n + 1) for y in range(x, n + 1) if x * y <= n)


@dataclass(frozen=True)
class ProductSumSeed:
    k: int
    product: int
    seed: tuple[int, ...]
    ones: int

    @property
    def witness(self) -> tuple[int, ...]:
        return (1,) * self.ones + self.seed


def compact_witness(seed: ProductSumSeed) -> str:
    seed_text = " ".join(map(str, seed.seed))
    if seed.ones == 0:
        return seed_text
    if seed.ones == 1:
        return f"1 + {seed_text}"
    return f"1^{seed.ones} + {seed_text}"


def enumerate_product_sum_seeds(max_k: int, max_product: int) -> dict[int, list[ProductSumSeed]]:
    """Enumerate nondecreasing seeds >1.

    If seed F has product P and sum S, then adding P-S copies of 1 gives
    a product-sum solution with arity |F| + P - S, whenever P >= S.
    """

    by_k: dict[int, list[ProductSumSeed]] = defaultdict(list)

    def rec(start: int, factors: list[int], product_so_far: int, sum_so_far: int) -> None:
        if len(factors) >= 2:
            ones = product_so_far - sum_so_far
            if ones >= 0:
                k = len(factors) + ones
                if 2 <= k <= max_k:
                    by_k[k].append(
                        ProductSumSeed(
                            k=k,
                            product=product_so_far,
                            seed=tuple(factors),
                            ones=ones,
                        )
                    )

        max_next = max_product // product_so_far
        for f in range(start, max_next + 1):
            rec(f, factors + [f], product_so_far * f, sum_so_far + f)

    rec(2, [], 1, 0)

    for k, seeds in by_k.items():
        seeds.sort(key=lambda s: (s.product, len(s.seed), s.seed))
        # Remove duplicate witnesses caused by different recursion paths.  The
        # nondecreasing recursion should prevent them, but this keeps the table
        # robust if the generator is edited later.
        dedup: list[ProductSumSeed] = []
        seen: set[tuple[int, tuple[int, ...], int]] = set()
        for seed in seeds:
            key = (seed.product, seed.seed, seed.ones)
            if key not in seen:
                seen.add(key)
                dedup.append(seed)
        by_k[k] = dedup
    return dict(by_k)


def two_factor_solutions_for_arity(k: int) -> list[tuple[int, int, int]]:
    """All two-nonunit product-sum seeds at total arity k.

    With r = k-2 ones, r + a + b = ab is equivalent to
    (a-1)(b-1) = k-1.
    """

    out = []
    target = k - 1
    for d in divisors(target):
        e = target // d
        a, b = d + 1, e + 1
        if a <= b:
            out.append((a, b, a * b))
    return out


def print_shadow_stats() -> None:
    print("NATURAL OPERATION DIGRAPH SHADOWS")
    print("=" * 72)
    print()
    print("Simple one-input shadows on [N]:")
    print("  additive:       x -> z iff exists y>=1 with x+y=z")
    print("                  so x -> z iff x<z: the transitive tournament")
    print("  multiplicative: x -> z iff exists y>=1 with xy=z")
    print("                  so x -> z iff x|z: the divisibility DAG")
    print()
    header = (
        "N  complete  add_edges  prod_edges(+1)  prod_edges(no1)  "
        "prod_density(+1)  add_gates  prod_gates(+1)  prime_covers"
    )
    print(header)
    print("-" * len(header))
    for n in [8, 12, 20, 50, 100]:
        complete = n * (n - 1) // 2
        add = len(additive_edges(n))
        prod_with_1 = len(product_edges(n, include_unit=True))
        prod_no_1 = len(product_edges(n, include_unit=False))
        density = prod_with_1 / complete
        print(
            f"{n:>3} {complete:>9} {add:>10} {prod_with_1:>15} "
            f"{prod_no_1:>16} {density:>17.4f} "
            f"{additive_gate_count(n):>10} {product_gate_count(n):>15} "
            f"{len(product_cover_edges(n)):>13}"
        )
    print()
    for n in [12, 20]:
        missing = sorted(additive_edges(n) - product_edges(n, include_unit=True))
        sample = ", ".join(f"{a}->{b}" for a, b in missing[:12])
        print(f"N={n}: product shadow misses {len(missing)} additive comparisons; first misses: {sample}")
    print()


def print_gate_collision_stats() -> None:
    print("BINARY GATE COLLISIONS")
    print("=" * 72)
    collisions = []
    for x in range(1, 40):
        for y in range(x, 40):
            if x + y == x * y:
                collisions.append((x, y, x + y))
    print(f"Positive unordered x+y=xy collisions below 40: {collisions}")
    print("  Rearrangement: (x-1)(y-1)=1, so the only positive collision is (2,2)->4.")
    print()
    print("Two-factor product-sum layer by total arity k:")
    print("  r ones plus a,b>1 solve r+a+b=ab with k=r+2")
    print("  equivalently (a-1)(b-1)=k-1, so this layer is the divisor graph of k-1.")
    for k in range(2, 15):
        sols = two_factor_solutions_for_arity(k)
        text = ", ".join(f"{a},{b}->P={p}" for a, b, p in sols)
        print(f"  k={k:>2}: {text}")
    print()


def print_product_sum_table() -> None:
    print("PRODUCT-SUM CRITICAL PAIRS")
    print("=" * 72)
    max_k = 24
    by_k = enumerate_product_sum_seeds(max_k=max_k, max_product=500)
    print("For seed F of entries >1, let P=prod(F), S=sum(F), defect D=P-S.")
    print("D copies of 1 repair the additive fold: D + S = P.")
    print("Total arity is k=|F|+D.  Minimal P by k gives the product-sum number.")
    print()
    header = "k  min_P  compact witness        seed>1        ones  runner_up_Ps"
    print(header)
    print("-" * len(header))
    for k in range(2, max_k + 1):
        seeds = by_k.get(k, [])
        if not seeds:
            print(f"{k:>2}  --     --                      --            --    --")
            continue
        best = seeds[0]
        runner_products = []
        for seed in seeds[1:]:
            if seed.product not in runner_products:
                runner_products.append(seed.product)
            if len(runner_products) >= 4:
                break
        witness = compact_witness(best)
        seed_text = " ".join(map(str, best.seed))
        runners = ",".join(map(str, runner_products)) if runner_products else "-"
        print(
            f"{k:>2} {best.product:>6}  {witness:<23} "
            f"{seed_text:<13} {best.ones:>4}  {runners}"
        )
    print()
    print("First arity where a multi-factor seed beats every two-factor seed:")
    for k in range(2, max_k + 1):
        seeds = by_k.get(k, [])
        if not seeds:
            continue
        best = seeds[0]
        two_factor_min = min((p for _, _, p in two_factor_solutions_for_arity(k)), default=None)
        if two_factor_min is not None and len(best.seed) > 2 and best.product < two_factor_min:
            print(
                f"  k={k}: best seed {best.seed} gives P={best.product}, "
                f"beating two-factor minimum P={two_factor_min}."
            )
            break
    print()


def print_defect_histogram() -> None:
    print("DEFECT HISTOGRAM FOR SMALL MULTIPLICATIVE SEEDS")
    print("=" * 72)
    seeds = enumerate_product_sum_seeds(max_k=40, max_product=250)
    all_seeds = [seed for values in seeds.values() for seed in values]
    by_seed_length = Counter(len(seed.seed) for seed in all_seeds)
    by_ones = Counter(seed.ones for seed in all_seeds)
    print(f"Seeds with product <=250 and induced arity <=40: {len(all_seeds)}")
    print("By nonunit seed length:")
    for length, count in sorted(by_seed_length.items()):
        print(f"  length {length}: {count}")
    print("Most common one-defects D=P-S:")
    for ones, count in by_ones.most_common(12):
        print(f"  D={ones:>2}: {count}")
    print()


def main() -> None:
    print_shadow_stats()
    print_gate_collision_stats()
    print_product_sum_table()
    print_defect_histogram()
    print("SYNTHESIS")
    print("=" * 72)
    print(
        "The additive shadow by itself is too coarse: it is just the total order. "
        "The multiplicative shadow is the divisibility DAG sitting inside that "
        "transitive completion.  Product-sum equations are therefore best treated "
        "as critical pairs in the labeled gate hypergraph, not as ordinary edges."
    )


if __name__ == "__main__":
    main()
