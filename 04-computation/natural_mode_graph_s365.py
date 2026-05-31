#!/usr/bin/env python3
"""
natural_mode_graph_s365.py

codex-2026-05-31 S365

Exploration of the "natural mode" digraphs from earlier repo history:

* additive mode graph: parents x,y point to z when x+y=z;
* multiplicative mode graph: parents x,y point to z when x*y=z;
* product-sum resonances: tuples whose additive and multiplicative outputs agree.

The earlier summand graph used distinct binary parents.  That convention is
important: it turns 4=2+2=2*2 into a hidden diagonal resonance and makes 6,
via 1+2+3=1*2*3, the first visible distinct resonance at higher arity.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations
from math import isqrt, prod


def additive_pairs(n: int, *, allow_equal: bool = False) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for a in range(1, n // 2 + 1):
        b = n - a
        if a > b:
            continue
        if a == b and not allow_equal:
            continue
        pairs.append((a, b))
    return pairs


def product_pairs(
    n: int,
    *,
    allow_equal: bool = False,
    include_one: bool = False,
) -> list[tuple[int, int]]:
    lo = 1 if include_one else 2
    pairs: list[tuple[int, int]] = []
    for a in range(lo, isqrt(n) + 1):
        if n % a:
            continue
        b = n // a
        if b < lo or a > b:
            continue
        if a == b and not allow_equal:
            continue
        pairs.append((a, b))
    return pairs


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


def product_atoms(limit: int) -> list[int]:
    return [n for n in range(2, limit + 1) if not product_pairs(n)]


def additive_closure(seeds: set[int], limit: int) -> set[int]:
    available = set(seeds)
    changed = True
    while changed:
        changed = False
        for n in range(1, limit + 1):
            if n in available:
                continue
            if any(a in available and b in available for a, b in additive_pairs(n)):
                available.add(n)
                changed = True
    return available


def product_closure(seeds: set[int], limit: int) -> set[int]:
    available = set(seeds)
    changed = True
    while changed:
        changed = False
        for n in range(2, limit + 1):
            if n in available:
                continue
            if any(a in available and b in available for a, b in product_pairs(n)):
                available.add(n)
                changed = True
    return available


def product_sum_solutions(arity: int) -> list[tuple[int, ...]]:
    """All nondecreasing positive tuples of fixed arity with sum=product."""

    solutions: set[tuple[int, ...]] = set()
    # In a product-sum tuple of arity k, every non-one factor is at most k.
    # For a two-factor core this is explicit:
    #   ab = a+b+k-2 -> b=(a+k-2)/(a-1) <= k.
    # Adding more factors only tightens the bound.
    max_factor = arity

    def rec(start: int, core: list[int], product: int, total: int) -> None:
        r = len(core)
        if r:
            defect = product - total
            if defect == arity - r:
                solutions.add(tuple([1] * defect + core))
            if r >= arity or defect > arity - r:
                return

        for x in range(start, max_factor + 1):
            next_product = product * x if core else x
            next_total = total + x
            rec(x, core + [x], next_product, next_total)

    rec(2, [], 1, 0)
    return sorted(solutions)


def distinct_product_sum_solutions(max_arity: int) -> list[tuple[int, ...]]:
    found: list[tuple[int, ...]] = []
    for arity in range(2, max_arity + 1):
        for tup in product_sum_solutions(arity):
            if len(set(tup)) == len(tup):
                found.append(tup)
    return found


def core_defect(core: tuple[int, ...]) -> int:
    return prod(core) - sum(core)


def cores_by_defect(max_defect: int, max_core_len: int = 4) -> dict[int, list[tuple[int, ...]]]:
    rows: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    max_factor = max_defect + 4

    def rec(start: int, core: list[int], product: int, total: int) -> None:
        if core:
            defect = product - total
            if len(core) >= 2 and 0 <= defect <= max_defect:
                rows[defect].append(tuple(core))
            if len(core) >= max_core_len or defect > max_defect:
                return
        for x in range(start, max_factor + 1):
            next_product = product * x if core else x
            rec(x, core + [x], next_product, total + x)

    rec(2, [], 1, 0)
    return {d: sorted(set(rows[d])) for d in sorted(rows)}


def print_binary_mode_table(limit: int) -> None:
    print("Binary mode graph comparison")
    print("  additive: distinct positive x<y, x+y=n")
    print("  product : distinct factors 2<=x<y, x*y=n")
    print()
    print(f"{'n':>3} {'add_in':>6} {'prod_in':>7} {'add_pairs':>18} {'prod_pairs':>18}")
    for n in range(1, 31):
        ap = additive_pairs(n)
        pp = product_pairs(n)
        print(f"{n:>3} {2*len(ap):>6} {2*len(pp):>7} {str(ap):>18} {str(pp):>18}")

    add_nodes = {n for n in range(1, limit + 1) if additive_pairs(n)}
    prod_nodes = {n for n in range(1, limit + 1) if product_pairs(n)}
    both = sorted(add_nodes & prod_nodes)
    print()
    print(f"Up to {limit}:")
    print(f"  additive_non_atoms={len(add_nodes)}")
    print(f"  product_non_atoms={len(prod_nodes)}")
    print(f"  both={len(both)} first={both[:24]}")
    print(f"  product_atoms={product_atoms(limit)[:40]}")
    print()


def print_product_sum_ledger() -> None:
    print("Product-sum resonances sum(tuple)=product(tuple)")
    print("  The tuple is written with all 1s first; the remaining entries are the core.")
    print()
    for arity in range(2, 11):
        rows = product_sum_solutions(arity)
        print(f"  arity={arity}: count={len(rows)}")
        for row in rows[:12]:
            ones = row.count(1)
            core = tuple(x for x in row if x != 1)
            print(
                "    "
                f"{row} core={core} defect={core_defect(core)} ones={ones}"
            )
        if len(rows) > 12:
            print(f"    ... {len(rows)-12} more")
    print()

    distinct = distinct_product_sum_solutions(max_arity=10)
    print("Distinct positive product-sum tuples from the exact arity ledgers through arity 10")
    print(f"  found={distinct}")
    print("  theorem_suggested=the only distinct positive solution is (1,2,3).")
    print()


def print_defect_table() -> None:
    print("Core-defect normal form")
    print("  For a core C with entries >=2, defect d=prod(C)-sum(C).")
    print("  Adding exactly d copies of 1 gives a product-sum tuple.")
    print()
    rows = cores_by_defect(12, max_core_len=4)
    for defect in range(0, 13):
        cores = rows.get(defect, [])
        shown = ", ".join(str(c) for c in cores[:10])
        suffix = "" if len(cores) <= 10 else f", ... {len(cores)-10} more"
        print(f"  defect={defect:>2}: {shown}{suffix}")
    print()
    print("  Two-core formula:")
    print("    core=(a,b), defect=ab-a-b=(a-1)(b-1)-1")
    print("    final arity = 2+defect = (a-1)(b-1)+1")
    print("    the spine core=(2,n) gives (1^(n-2),2,n), product=sum=2n.")
    print()


def print_hidden_resonances() -> None:
    print("Hidden low-dimensional resonances")
    binary_distinct = []
    binary_equal = []
    for a in range(1, 20):
        for b in range(a, 20):
            if a + b == a * b:
                target = a + b
                row = (a, b, target)
                if a == b:
                    binary_equal.append(row)
                else:
                    binary_distinct.append(row)
    print(f"  binary distinct x+y=x*y: {binary_distinct}")
    print(f"  binary with diagonal allowed: {binary_equal}")
    print("  interpretation: 4 is the diagonal binary resonance 2+2=2*2.")

    ternary_distinct = []
    for tup in combinations(range(1, 20), 3):
        if sum(tup) == prod(tup):
            ternary_distinct.append(tup)
    print(f"  ternary distinct x+y+z=xyz: {ternary_distinct}")
    print("  interpretation: 6 is the first visible distinct resonance, 1+2+3=1*2*3.")
    print()


def print_closure_comparison(limit: int) -> None:
    add_23 = additive_closure({2, 3}, limit)
    prod_23 = product_closure({2, 3}, limit)
    print("Closure comparison from the same seeds {2,3}")
    print(f"  additive missing up to {limit}: {[n for n in range(1, limit+1) if n not in add_23]}")
    print(f"  product reachable up to {limit}: {sorted(prod_23)}")
    print("  lesson=addition from {2,3} is cofinite; multiplication from {2,3} stays 2,3-smooth.")
    print()


def main() -> None:
    print("Natural mode graphs: additive, multiplicative, and product-sum (S365)")
    print("=" * 76)
    print()
    print_binary_mode_table(80)
    print_hidden_resonances()
    print_product_sum_ledger()
    print_defect_table()
    print_closure_comparison(80)
    print("Synthesis")
    print("  1. The additive graph is dense and cofinite from the seeds {2,3}.")
    print("  2. The product graph is sparse because every prime and prime square is an atom.")
    print("  3. The family sum=product is not random; it is exactly the defect-zero interface.")
    print("  4. Stripping 1s turns product-sum equations into the defect law prod(C)-sum(C)=#ones.")
    print("  5. The old additive module {1,4,6} is now legible:")
    print("       1 = identity/source, 4 = hidden binary diagonal, 6 = first distinct ternary resonance.")


if __name__ == "__main__":
    main()
