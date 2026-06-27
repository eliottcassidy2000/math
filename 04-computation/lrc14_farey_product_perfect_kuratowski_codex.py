#!/usr/bin/env python3
"""Farey product sets, perfect-number shadows, and Kuratowski ledgers.

For Farey order n, use the reduced terms a/b with 0 < a <= b <= n and attach
the product a*b.  Numerically this product is |E(K_{a,b})|, but graph-minor
obstructions are not edge-count averages.  This script keeps those roles apart.
"""

from __future__ import annotations

from collections import defaultdict
from math import gcd, isqrt


PERFECTS = [
    6,
    28,
    496,
    8128,
    33550336,
    8589869056,
    137438691328,
]


def divisors(n: int) -> list[int]:
    out: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def sigma(n: int) -> int:
    return sum(divisors(n))


def farey_product_witnesses(n: int) -> dict[int, list[tuple[int, int]]]:
    witnesses: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for b in range(1, n + 1):
        for a in range(1, b + 1):
            if gcd(a, b) == 1:
                witnesses[a * b].append((a, b))
    return dict(sorted(witnesses.items()))


def product_set(n: int) -> list[int]:
    return sorted(farey_product_witnesses(n))


def max_product(n: int) -> int:
    if n == 1:
        return 1
    return n * (n - 1)


def graph_tag(a: int, b: int) -> str:
    if min(a, b) >= 3:
        return "K33-minor"
    if min(a, b) == 2:
        return "planar-two-block"
    return "planar-star"


def first_product_value(value: int, limit: int = 200) -> tuple[int, list[tuple[int, int]]] | None:
    for n in range(1, limit + 1):
        witnesses = farey_product_witnesses(n)
        if value in witnesses:
            return n, witnesses[value]
    return None


def first_k33_minor(limit: int = 200) -> tuple[int, tuple[int, int]] | None:
    for n in range(1, limit + 1):
        for product, witnesses in farey_product_witnesses(n).items():
            for a, b in witnesses:
                if min(a, b) >= 3:
                    return n, (a, b)
    return None


def perfect_hits(limit: int = 1000) -> list[tuple[int, int, int, str]]:
    hits: list[tuple[int, int, int, str]] = []
    pset: set[int] = set()
    total = 0
    for n in range(1, limit + 1):
        # Incremental update by new denominator n.
        for a in range(1, n + 1):
            if gcd(a, n) == 1:
                product = a * n
                if product not in pset:
                    pset.add(product)
                    total += product
        if total in PERFECTS:
            hits.append((n, len(pset), total, "sum-perfect"))
        if total % 2 == 0 and total // 2 in PERFECTS:
            hits.append((n, len(pset), total, "half-sum-perfect"))
    return hits


def main() -> None:
    print("LRC14 FAREY PRODUCT / PERFECT / KURATOWSKI ATLAS")
    print("=" * 76)
    print()

    print("A. Product-set divisor closure")
    print("-" * 76)
    print("P_n = {a*b : 0 < a <= b <= n, gcd(a,b)=1}.")
    print("M_n = max(P_n) = n*(n-1), except M_1=1.")
    print("D(M_n) is always contained in P_n; equality is the special gate.")
    print()
    print(" n | |P_n| | M_n | sigma(M_n) | sum(P_n) | excess | closure | perfect tag")
    for n in range(1, 16):
        pset = product_set(n)
        m = max_product(n)
        divs = divisors(m)
        total = sum(pset)
        excess = total - sum(divs)
        closure = pset == divs
        tags: list[str] = []
        if m in PERFECTS:
            tags.append("M perfect")
        if total in PERFECTS:
            tags.append("sum perfect")
        if total % 2 == 0 and total // 2 in PERFECTS:
            tags.append("half-sum perfect")
        tag = ", ".join(tags) if tags else "-"
        print(
            f"{n:2d} | {len(pset):5d} | {m:5d} | {sigma(m):10d} | "
            f"{total:8d} | {excess:6d} | {str(closure):7s} | {tag}"
        )
    print()

    print("Exact closure proof sketch:")
    print("  Every divisor d of n*(n-1) splits as d1*d2 with d1|n, d2|(n-1),")
    print("  gcd(d1,d2)=1.  Then d = a*b for a=min(d1,d2), b=max(d1,d2),")
    print("  and b <= n, so d is in P_n.")
    print("  For n >= 5, (n-2)*(n-1) is in P_n via (n-2)/(n-1),")
    print("  but it does not divide n*(n-1), since n-2 does not divide n.")
    print("  Thus P_n = D(M_n) only for n <= 4.")
    print()

    print("B. The F3/F4 perfect gate")
    print("-" * 76)
    for n in [3, 4]:
        pset = product_set(n)
        m = max_product(n)
        print(f"F_{n}: P_{n}={pset}, M={m}, sigma(M)={sigma(m)}, sum(P)={sum(pset)}")
        witnesses = farey_product_witnesses(n)
        for product in pset:
            forms = " ".join(f"{a}/{b}:{graph_tag(a,b)}" for a, b in witnesses[product])
            print(f"  product {product:2d}: {forms}")
    print()
    print("Readout:")
    print("  F3 is planar and divisor-closed at M=6; sigma(6)=12=2*6.")
    print("  F4 is the last divisor-closed order and sum(P4)=28, the next perfect number.")
    print("  The same F4 step introduces 3/4 -> K_{3,4}, which contains K33.")
    print()

    print("C. Edge-count aliases are not obstruction cores")
    print("-" * 76)
    aliases = {
        "K33 edge count 9": 9,
        "K5 edge count 10": 10,
        "K5 disjoint K33 edge count 19": 19,
        "F4 max/divisor gate 12": 12,
        "perfect-number edge alias 28": 28,
    }
    for name, value in aliases.items():
        hit = first_product_value(value)
        if hit is None:
            print(f"{name:31s}: no Farey product <= 200")
            continue
        n, witnesses = hit
        forms = ", ".join(f"{a}/{b}->{graph_tag(a,b)}" for a, b in witnesses)
        print(f"{name:31s}: first at F_{n}: {forms}")
    k33 = first_k33_minor()
    if k33:
        n, (a, b) = k33
        print(f"first actual K33-minor term: F_{n}, {a}/{b}, product {a*b}")
    print()
    print("Kuratowski guardrail:")
    print("  K5 and K33 are a two-element obstruction set.")
    print("  K5 disjoint K33 is nonplanar because it already contains the cores.")
    print("  Edge-count aliases such as product 10 or 19 do not create new cores.")
    print()

    print("D. Mediant/minor noncommutation")
    print("-" * 76)
    left = (2, 3)
    right = (1, 1)
    med = (left[0] + right[0], left[1] + right[1])
    print(f"{left[0]}/{left[1]} -> {graph_tag(*left)}")
    print(f"{right[0]}/{right[1]} -> {graph_tag(*right)}")
    print(f"mediant {med[0]}/{med[1]} -> {graph_tag(*med)}")
    print("  Two planar Farey parents can mediant into the first K33-minor child.")
    print("  Minor/subdivision containment is transitive; mediant-taking is not that operation.")
    print()

    print("E. Perfect-sum finite scout")
    print("-" * 76)
    hits = perfect_hits(1000)
    for n, size, total, tag in hits:
        print(f"F_{n}: |P|={size}, sum(P)={total}, {tag}")
    print("No other perfect or half-perfect product sums occur for n <= 1000.")
    print()

    print("F. New hypotheses generated by this atlas")
    print("-" * 76)
    print("HYP-2944: F3/F4 form a perfect-product obstruction gate.")
    print("HYP-2945: Farey mediants and graph minors do not commute; K33 begins at 3/4.")
    print("HYP-2946: After F4, product-excess P_n \\ D(n(n-1)) is the divisor-leakage signal.")
    print()

    print("G. Tournament analysis")
    print("-" * 76)
    roles = [
        "exact LRC M/Farey branch",
        "AP/GW C27 shell labels",
        "Farey product divisor closure",
        "Kpq minor/subdivision obstruction",
        "post-F4 product-excess leakage",
        "perfect/aliquot shadow",
        "raw edge-count aliases",
        "mediant as numerical average",
    ]
    edges = len(roles) * (len(roles) - 1) // 2
    print("Relation: x -> y iff x retains more proof-critical structure than y.")
    print(f"vertices={len(roles)} edges={edges} c3=0 hp=1")
    for i, role in enumerate(roles):
        print(f"  score {len(roles)-1-i}: {role}")
    print()
    print("Verdict:")
    print("  F3 and F4 are not just small examples.  F3 is the planar perfect-max")
    print("  closure, and F4 is simultaneously the last divisor closure, the first")
    print("  K33-minor Farey product term, and the product-sum 28 perfect-number hit.")
    print("  From F5 onward the extra product layer is real leakage, so LRC14 should")
    print("  keep exact C27/Farey labels before using product/perfect analogies.")


if __name__ == "__main__":
    main()
