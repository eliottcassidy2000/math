#!/usr/bin/env python3
"""
natural_operation_graphs_s365.py

codex-2026-05-31 S365

Explore two operation-generated oriented graphs on positive natural numbers:

  additive:       x -> z and y -> z when x + y = z
  multiplicative: x -> z and y -> z when x * y = z

If labels are forgotten, the additive graph on {1,...,N} is exactly the
transitive tournament/order: x -> z iff x < z.  The multiplicative graph,
after removing unit loops, is the proper-divisor DAG.  The interesting object
is therefore the labeled operation hypergraph, not the simple edge set.

The script also enumerates the collision equations

  x_1 + ... + x_k = x_1 * ... * x_k,

which are the arity-k places where additive and multiplicative operation
hyperedges land on the same natural-number node.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import prod


@dataclass(frozen=True)
class GraphStats:
    n: int
    additive_edges: int
    additive_expected: int
    additive_transitive_reduction: int
    multiplicative_edges_with_unit: int
    multiplicative_edges_without_unit: int
    multiplicative_hasse_edges: int
    additive_indegree_hist: tuple[tuple[int, int], ...]
    multiplicative_indegree_hist: tuple[tuple[int, int], ...]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def additive_simple_edges(n: int) -> set[tuple[int, int]]:
    """Unlabeled additive edges induced by positive complements."""

    return {(x, z) for z in range(1, n + 1) for x in range(1, z)}


def additive_labeled_fibers(n: int) -> dict[int, list[tuple[int, int]]]:
    fibers: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for x in range(1, n + 1):
        for y in range(1, n + 1):
            z = x + y
            if z <= n:
                fibers[z].append((x, y))
    return dict(fibers)


def multiplicative_simple_edges(
    n: int,
    *,
    include_unit_source: bool,
) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for z in range(1, n + 1):
        for x in range(1, z):
            if z % x == 0:
                if include_unit_source or x != 1:
                    edges.add((x, z))
    return edges


def multiplicative_labeled_fibers(n: int) -> dict[int, list[tuple[int, int]]]:
    fibers: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for x in range(1, n + 1):
        for y in range(1, n + 1):
            z = x * y
            if z <= n:
                fibers[z].append((x, y))
    return dict(fibers)


def multiplicative_hasse_edges(n: int) -> set[tuple[int, int]]:
    """Transitive reduction of divisibility: multiply by one prime."""

    edges: set[tuple[int, int]] = set()
    for x in range(1, n + 1):
        p = 2
        while x * p <= n:
            if is_prime(p):
                edges.add((x, x * p))
            p += 1
    return edges


def graph_stats(n: int) -> GraphStats:
    add_edges = additive_simple_edges(n)
    mul_unit = multiplicative_simple_edges(n, include_unit_source=True)
    mul_nonunit = multiplicative_simple_edges(n, include_unit_source=False)
    mul_hasse = multiplicative_hasse_edges(n)

    add_indegree = Counter()
    for _x, z in add_edges:
        add_indegree[z] += 1
    mul_indegree = Counter()
    for _x, z in mul_nonunit:
        mul_indegree[z] += 1

    return GraphStats(
        n=n,
        additive_edges=len(add_edges),
        additive_expected=n * (n - 1) // 2,
        additive_transitive_reduction=max(0, n - 1),
        multiplicative_edges_with_unit=len(mul_unit),
        multiplicative_edges_without_unit=len(mul_nonunit),
        multiplicative_hasse_edges=len(mul_hasse),
        additive_indegree_hist=tuple(sorted(Counter(add_indegree.values()).items())),
        multiplicative_indegree_hist=tuple(sorted(Counter(mul_indegree.values()).items())),
    )


def product_sum_solutions(k: int) -> list[tuple[int, ...]]:
    """Enumerate sorted positive k-tuples with sum(tuple)=product(tuple).

    The recursion is exact by the prefix-product gate.  If

        1 <= a_1 <= ... <= a_k
        a_1 + ... + a_k = a_1 * ... * a_k,

    then P = a_1 * ... * a_{k-1} <= k.  Indeed,

        P * a_k = sum_i a_i <= k * a_k.

    So the first k-1 entries have bounded product.  Once those entries have
    product P and sum S, the last entry is forced by

        x * (P - 1) = S.
    """

    out: list[tuple[int, ...]] = []

    def search_prefix(prefix: tuple[int, ...], start: int, p: int, s: int) -> None:
        if len(prefix) == k - 1:
            if p <= 1:
                return
            denom = p - 1
            if s % denom:
                return
            x = s // denom
            if x >= start and x > 0:
                candidate = prefix + (x,)
                if sum(candidate) == prod(candidate):
                    out.append(candidate)
            return

        for x in range(start, k // p + 1):
            search_prefix(prefix + (x,), x, p * x, s + x)

    search_prefix(tuple(), 1, 1, 0)
    return out


def reciprocal_denominators(solution: tuple[int, ...]) -> tuple[int, ...]:
    p = prod(solution)
    return tuple(sorted(p // x for x in solution))


def core_without_ones(solution: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x for x in solution if x != 1)


def core_defect(core: tuple[int, ...]) -> int:
    return prod(core) - sum(core)


def product_sum_table(max_k: int) -> dict[int, list[tuple[int, ...]]]:
    return {k: product_sum_solutions(k) for k in range(2, max_k + 1)}


def print_graph_stats() -> None:
    print("Operation graph collapse/sparsity")
    print("  Additive simple edge: x->z iff exists y>=1 with x+y=z.")
    print("  Multiplicative simple edge: x->z iff exists y>=2 with x*y=z.")
    print()
    print(
        "  "
        f"{'N':>4} {'A_edges':>8} {'A_expected':>10} {'A_Hasse':>8} "
        f"{'M_edges':>8} {'M+unit':>8} {'M_Hasse':>8} {'M/A':>8}"
    )
    for n in [8, 12, 20, 30, 50, 80]:
        stats = graph_stats(n)
        ratio = stats.multiplicative_edges_without_unit / stats.additive_edges
        print(
            "  "
            f"{n:4d} {stats.additive_edges:8d} {stats.additive_expected:10d} "
            f"{stats.additive_transitive_reduction:8d} "
            f"{stats.multiplicative_edges_without_unit:8d} "
            f"{stats.multiplicative_edges_with_unit:8d} "
            f"{stats.multiplicative_hasse_edges:8d} {ratio:8.4f}"
        )
    print()
    print("  Interpretation:")
    print("    - Additive edges are exactly the transitive tournament/order.")
    print("    - The additive transitive reduction is the successor chain.")
    print("    - Multiplicative edges are the divisor DAG; its Hasse edges are x->xp.")
    print("    - Multiplication keeps the prime lattice visible after forgetting labels.")
    print()


def print_fiber_examples(n: int = 18) -> None:
    add = additive_labeled_fibers(n)
    mul = multiplicative_labeled_fibers(n)
    print(f"Labeled operation fibers up to N={n}")
    for z in range(2, n + 1):
        add_f = add.get(z, [])
        mul_f = mul.get(z, [])
        if z <= 12 or mul_f:
            add_unordered = len({tuple(sorted(pair)) for pair in add_f})
            mul_unordered = len({tuple(sorted(pair)) for pair in mul_f})
            print(
                "  "
                f"z={z:2d}: add ordered={len(add_f):2d} unordered={add_unordered:2d}; "
                f"mul ordered={len(mul_f):2d} unordered={mul_unordered:2d}"
            )
    print()
    print("  The unlabeled additive edge x->z forgets the complement y=z-x.")
    print("  The labeled fiber remembers Goldbach/partition-style information.")
    print("  Multiplicative fibers remember divisors/factorizations.")
    print()


def print_collision_gate() -> None:
    print("Finite product-sum collision gate")
    print("  For a sorted collision a_1<=...<=a_k,")
    print("    a_1*...*a_{k-1} <= k.")
    print("  Therefore each fixed arity has finitely many positive solutions.")
    print("  After deleting q ones and calling the remaining core c,")
    print("    q = product(c) - sum(c).")
    print("  Thus 1 is a defect absorber: it changes the sum without changing product.")
    print()


def print_product_sum_solutions(max_k: int = 20) -> None:
    table = product_sum_table(max_k)
    print("Product-sum collision equations")
    print("  Sorted positive solutions of x_1+...+x_k = x_1*...*x_k.")
    print("  Reciprocally: sum_i 1/(product/x_i) = 1.")
    print()
    print(
        "  "
        f"{'k':>2} {'count':>5} {'ones_hist':>22} {'max_entry':>9} examples"
    )
    for k, rows in table.items():
        ones_hist = Counter(row.count(1) for row in rows)
        max_entry = max((max(row) for row in rows), default=0)
        examples = ", ".join(str(row) for row in rows[:4])
        if len(rows) > 4:
            examples += ", ..."
        print(
            "  "
            f"{k:2d} {len(rows):5d} {str(dict(sorted(ones_hist.items()))):>22} "
            f"{max_entry:9d} {examples}"
        )
    print()

    print("  Egyptian-fraction shadows of selected collisions")
    selected: list[tuple[int, ...]] = []
    for k in range(2, max_k + 1):
        selected.extend(table[k][:3])
    for row in selected[:18]:
        dens = reciprocal_denominators(row)
        print(f"    {row} -> denominators {dens}, sum 1/d = 1")
    print()


def print_collision_targets(max_k: int = 20) -> None:
    table = product_sum_table(max_k)
    print("Collision targets by arity")
    print("  Each target z is a natural-number mode hit by both k-fold sum and product.")
    print("  Counts in parentheses mean multiple sorted collision fibers hit the same z.")
    for k, rows in table.items():
        targets = Counter(sum(row) for row in rows)
        target_text = ", ".join(
            f"{target}({count})" if count > 1 else str(target)
            for target, count in sorted(targets.items())
        )
        print(f"    k={k:2d}: {target_text}")
    print()


def print_target_spectra(max_k: int = 20, top_n: int = 8) -> None:
    table = product_sum_table(max_k)
    spectrum: dict[int, list[tuple[int, tuple[int, ...]]]] = defaultdict(list)
    for k, rows in table.items():
        for row in rows:
            spectrum[sum(row)].append((k, core_without_ones(row)))

    ranked = sorted(spectrum.items(), key=lambda item: (-len(item[1]), item[0]))

    print("Target-mode collision spectra")
    print("  Ranked by number of collision fibers seen through the arity cutoff.")
    for z, records in ranked[:top_n]:
        record_text = "; ".join(
            f"k={k}, core={core}" for k, core in records
        )
        print(f"    z={z:2d}: {len(records):2d} fibers | {record_text}")
    print()


def print_core_transitions(max_k: int = 20) -> None:
    table = product_sum_table(max_k)
    core_hist: Counter[tuple[int, ...]] = Counter()
    for rows in table.values():
        for row in rows:
            core_hist[core_without_ones(row)] += 1

    print("Core patterns after deleting 1s")
    print("  A nonunit core c creates a solution exactly when defect(c)>=0.")
    print("  The solution is defect(c) copies of 1 plus c.")
    for core, count in core_hist.most_common(24):
        print(f"    core={core} defect={core_defect(core):2d} count={count}")
    print()
    print("  Two-nonunit cores have a shifted divisor law:")
    print("    core=(a,b) occurs at arity k iff (a-1)(b-1)=k-1.")
    for k, rows in table.items():
        pairs = [core_without_ones(row) for row in rows if len(core_without_ones(row)) == 2]
        if pairs:
            pair_text = ", ".join(
                f"{pair}->z={prod(pair)}" for pair in pairs
            )
            print(f"    k={k:2d}: {pair_text}")
    print()


def print_sum_product_bridge() -> None:
    print("Repo bridge")
    print("  Additive graph, unlabeled: complete transitive tournament.")
    print("  Multiplicative graph, unlabeled: incomplete oriented graph / divisor DAG.")
    print("  Additive graph, labeled: partitions of z.")
    print("  Multiplicative graph, labeled: factorizations of z.")
    print("  Product-sum tuples are operation-hyperedge collisions.")
    print()
    print("  This reframes earlier repo threads:")
    print("    - oriented graphs are incomplete tournaments; divisor DAGs are arithmetic examples;")
    print("    - Zeckendorf writes naturals as independence packets in a path graph;")
    print("    - Cayley rapidity linearizes multiplication but leaves addition external;")
    print("    - Egyptian fractions appear by dividing sum=product by the product;")
    print("    - sum-product tournament phenomena compare additive fibers to multiplicative fibers.")


def main() -> None:
    print("Natural operation graphs and product-sum collisions (codex-2026-05-31 S365)")
    print("=" * 78)
    print()
    print_graph_stats()
    print_fiber_examples()
    print_collision_gate()
    print_product_sum_solutions(max_k=20)
    print_collision_targets(max_k=20)
    print_target_spectra(max_k=20)
    print_core_transitions(max_k=20)
    print_sum_product_bridge()


if __name__ == "__main__":
    main()
