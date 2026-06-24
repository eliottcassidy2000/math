#!/usr/bin/env python3
"""Summand/multiplicand merge for the additive-basis Farey bridge.

This is a synthesis scout.  It does not prove Goldbach, Fermat polygonal,
Zeckendorf, or LRC statements.  It makes the shared packet schema executable:

* additive representations live on summand antidiagonals x+y=N;
* multiplicative representations live on hyperbolas xy=N;
* a Farey packet p/q has both shadows, p+q and p*q;
* the Fibonacci row 1,1,1+1,1+2,... is the d=2 Pascal/path-rank row.

Tournament Analysis uses proof currencies as vertices, not runners.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import comb, isqrt


def fibonacci(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(values[-1] + values[-2])
    return values


FIB = fibonacci(40)


def pascal_fib_row(n: int) -> list[int]:
    """Return the row C(n-1-k,k), whose sum is F_n for F_1=F_2=1."""
    return [comb(n - 1 - k, k) for k in range((n - 1) // 2 + 1)]


def prime_sieve(limit: int) -> list[bool]:
    prime = [False, False] + [True] * (limit - 1)
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            step_start = p * p
            prime[step_start : limit + 1 : p] = [False] * (((limit - step_start) // p) + 1)
    return prime


PRIME = prime_sieve(2000)


def goldbach_pairs(n: int) -> list[tuple[int, int]]:
    return [(p, n - p) for p in range(2, n // 2 + 1) if PRIME[p] and PRIME[n - p]]


def ternary_prime_triples(n: int) -> list[tuple[int, int, int]]:
    triples: list[tuple[int, int, int]] = []
    primes = [p for p in range(2, n + 1) if PRIME[p]]
    for i, p in enumerate(primes):
        for q in primes[i:]:
            r = n - p - q
            if r < q:
                break
            if r <= 2000 and PRIME[r]:
                triples.append((p, q, r))
    return triples


def polygonal(side: int, m: int) -> int:
    return ((side - 2) * m * m - (side - 4) * m) // 2


def polygonal_atoms(side: int, limit: int) -> list[int]:
    atoms = {0}
    m = 1
    while True:
        value = polygonal(side, m)
        if value > limit:
            break
        atoms.add(value)
        m += 1
    return sorted(atoms)


def max_min_polygonal_summands(side: int, limit: int) -> tuple[int, int]:
    atoms = polygonal_atoms(side, limit)
    inf = limit + 1
    dp = [0] + [inf] * limit
    for n in range(1, limit + 1):
        dp[n] = min(1 + dp[n - atom] for atom in atoms[1:] if atom <= n)
    max_need = max(dp[1:])
    first_at_max = next(n for n, value in enumerate(dp) if value == max_need)
    return max_need, first_at_max


def zeckendorf(n: int) -> list[int]:
    atoms = [f for f in FIB[2:] if f <= n]
    result: list[int] = []
    remaining = n
    for f in reversed(atoms):
        if f <= remaining:
            result.append(f)
            remaining -= f
    return list(reversed(result))


def zeckendorf_indices(values: list[int]) -> list[int]:
    index = {value: idx for idx, value in enumerate(FIB) if idx >= 2}
    return [index[value] for value in values]


def summand_pair_count(n: int) -> int:
    return (n - 1) // 2


def multiplicand_pair_count(n: int) -> int:
    count = 0
    for a in range(2, isqrt(n) + 1):
        if n % a == 0 and a < n // a:
            count += 1
    return count


def distinct_factor_pairs(n: int) -> list[tuple[int, int]]:
    return [(a, n // a) for a in range(2, isqrt(n) + 1) if n % a == 0 and a < n // a]


@dataclass(frozen=True)
class Carrier:
    name: str
    role: str
    preserves: tuple[str, ...]
    destroys: tuple[str, ...]
    vector: tuple[int, ...]


RETENTION_KEYS = (
    "exact_packet",
    "additive_fiber",
    "product_fiber",
    "carry_normal_form",
    "bounded_arity",
    "smoothing_entropy",
    "local_residue",
    "lrc_transfer",
    "anti_scalar_guard",
)


CARRIERS = (
    Carrier(
        "exact_farey_packet",
        "root address p/q, q, and excess before any mutation",
        ("p/q", "q", "Farey excess", "route labels"),
        ("nothing intentional; this is the base packet",),
        (5, 4, 4, 4, 3, 3, 4, 5, 5),
    ),
    Carrier(
        "summand_antidiagonal_fiber",
        "additive graph node p+q and all a+b=p+q witnesses",
        ("additive lane", "pinch denominators", "Pascal rows"),
        ("factor sparsity", "ordered magnitude"),
        (4, 5, 2, 3, 3, 4, 3, 4, 4),
    ),
    Carrier(
        "multiplicand_hyperbola_fiber",
        "product graph node p*q and factor/Kpq incidence",
        ("product lane", "Kpq incidence", "factor fibers"),
        ("additive branch multiplicity", "carry order"),
        (4, 2, 5, 2, 2, 2, 4, 4, 4),
    ),
    Carrier(
        "zeckendorf_path_carry",
        "no-adjacent Fibonacci support and confluent local carries",
        ("canonical support", "path conflict graph", "carry confluence"),
        ("representation abundance", "prime residue smoothing"),
        (3, 4, 2, 5, 3, 1, 3, 4, 5),
    ),
    Carrier(
        "fibonacci_pascal_row",
        "rank vector C(n-1-k,k), not just scalar F_n",
        ("row fiber", "path independent-set ranks", "carry width"),
        ("chosen normal form", "multiplicative incidence"),
        (3, 5, 1, 4, 3, 2, 2, 3, 4),
    ),
    Carrier(
        "fermat_polygonal_invoice",
        "bounded arity and residue absorption for s-gonal atoms",
        ("finite summand budget", "side-count invoice", "local residues"),
        ("uniqueness", "prime smoothing entropy"),
        (2, 4, 1, 2, 5, 2, 4, 3, 4),
    ),
    Carrier(
        "ternary_goldbach_smoothing",
        "prime triple hypergraph with enough smoothing dimension",
        ("many branches", "local residue rank", "circle-method shape"),
        ("canonical support", "bounded deterministic invoice"),
        (2, 4, 1, 1, 2, 5, 5, 2, 3),
    ),
    Carrier(
        "binary_goldbach_pair_sieve",
        "prime pair graph p+q=N with thin smoothing",
        ("pair graph", "singular series hints", "prime residue channels"),
        ("third-variable smoothing", "normal form"),
        (2, 4, 1, 1, 1, 4, 5, 2, 3),
    ),
    Carrier(
        "farey_power_stress",
        "q^p and p^q ordered magnitude stress lanes",
        ("magnitude flip detector",),
        ("proof scale", "fiber identity", "safe quotient data"),
        (1, 1, 1, 1, 1, 1, 1, 1, 5),
    ),
)


def winner(i: int, j: int) -> int:
    left = CARRIERS[i].vector
    right = CARRIERS[j].vector
    lw = sum(a > b for a, b in zip(left, right))
    rw = sum(b > a for a, b in zip(left, right))
    if lw > rw:
        return i
    if rw > lw:
        return j
    return min(i, j)


def build_edges() -> dict[int, set[int]]:
    edges = {i: set() for i in range(len(CARRIERS))}
    for i, j in combinations(range(len(CARRIERS)), 2):
        w = winner(i, j)
        edges[w].add(j if w == i else i)
    return edges


def count_3cycles(edges: dict[int, set[int]]) -> int:
    total = 0
    for i, j, k in combinations(range(len(CARRIERS)), 3):
        triples = ((i, j, k), (i, k, j), (j, i, k), (j, k, i), (k, i, j), (k, j, i))
        if any(b in edges[a] and c in edges[b] and a in edges[c] for a, b, c in triples):
            total += 1
    return total


def scc_sizes(edges: dict[int, set[int]]) -> list[int]:
    n = len(CARRIERS)
    reverse = {i: set() for i in range(n)}
    for src, dsts in edges.items():
        for dst in dsts:
            reverse[dst].add(src)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in edges[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: int, bucket: list[int]) -> None:
        seen.add(v)
        bucket.append(v)
        for nxt in reverse[v]:
            if nxt not in seen:
                rdfs(nxt, bucket)

    for v in reversed(order):
        if v not in seen:
            bucket: list[int] = []
            rdfs(v, bucket)
            sizes.append(len(bucket))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(edges: dict[int, set[int]]) -> int:
    n = len(CARRIERS)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if count == 0:
                continue
            for nxt in edges[last]:
                if not (mask & (1 << nxt)):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def print_header(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def main() -> None:
    print("S171 summand/multiplicand Farey additive-basis merge")
    print("=" * 62)

    print_header("[1] Fibonacci row arrangement")
    print("row formula: row(n,k)=C(n-1-k,k), sum_k row(n,k)=F_n")
    for n in range(1, 11):
        row = pascal_fib_row(n)
        print(f"n={n:2d} F_n={sum(row):3d} row={'+'.join(str(x) for x in row)}")
    print("readout: the scalar F_n is a shadow; the row vector is the packet.")

    print_header("[2] Additive-basis finite checks")
    for n in (20, 50, 100, 200):
        print(f"binary Goldbach pairs N={n:3d}: count={len(goldbach_pairs(n)):2d} sample={goldbach_pairs(n)[:5]}")
    for n in (101, 501, 999):
        triples = ternary_prime_triples(n)
        print(f"ternary Goldbach triples N={n:3d}: count={len(triples):4d} sample={triples[:4]}")
    for side in range(3, 9):
        need, first = max_min_polygonal_summands(side, 300)
        print(f"{side}-gonal bounded arity through 300: max_min_summands={need}, first_at_max={first}")
    for n in (21, 42, 89, 144, 233, 377):
        z = zeckendorf(n)
        idx = zeckendorf_indices(z)
        print(f"Zeckendorf N={n:3d}: atoms={z}, indices={idx}, digits={len(z)}")

    print_header("[3] Summand and multiplicand graph split")
    print("summand node N: antidiagonal x+y=N; pair count=floor((N-1)/2)")
    print("multiplicand node N: hyperbola xy=N; nonunit factor-pair count is sparse")
    for n in (14, 21, 27, 41, 54, 89, 123, 377):
        print(
            f"N={n:3d} summand_pairs={summand_pair_count(n):3d} "
            f"factor_pairs={multiplicand_pair_count(n):2d} factors={distinct_factor_pairs(n)[:5]}"
        )

    print_header("[4] Golden Farey spine")
    print("On p/q=F_i/F_{i+1}, p+q is the next Fibonacci number.")
    print("p*q is the K_{p,q} incidence area; powers are order-stress tests.")
    for i in range(2, 10):
        p, q = FIB[i], FIB[i + 1]
        product = p * q
        z = zeckendorf(p + q)
        power_winner = "q^p" if q**p > p**q else "p^q"
        print(
            f"F_{i}/F_{i+1}={p}/{q}: p+q={p+q}=F_{i+2}, "
            f"p*q={product}, sum_pairs={summand_pair_count(p+q)}, "
            f"factor_pairs(product)={multiplicand_pair_count(product)}, "
            f"Zeck(p+q)={z}, power_winner={power_winner}"
        )

    print_header("[5] LRC14 unit-excess Farey packet")
    print("For p/q=p/(14p-1), keep e=14p-q=1 before using sum/product/powers.")
    for p in range(1, 9):
        q = 14 * p - 1
        product = p * q
        branch = "q-parent/AP-GW" if p == 1 else "C27 petal/two-block" if p == 2 else "Kpq/K33 product wall"
        print(
            f"p={p:2d} q={q:3d} p+q={p+q:3d} p*q={product:4d} "
            f"sum_pairs={summand_pair_count(p+q):3d} "
            f"factor_pairs(product)={multiplicand_pair_count(product):2d} "
            f"route={branch}"
        )

    print_header("[6] Tournament Analysis")
    print("vertices: proof currencies / representation fibers, not runners")
    print("pairwise observable: coordinate-wise retention vector")
    print("switch/gauge: majority of better-retained proof coordinates")
    print("tie Hamiltonian path: " + " > ".join(carrier.name for carrier in CARRIERS))
    edges = build_edges()
    scores = {carrier.name: len(edges[i]) for i, carrier in enumerate(CARRIERS)}
    for name, score in sorted(scores.items(), key=lambda item: (-item[1], item[0])):
        carrier = next(c for c in CARRIERS if c.name == name)
        print(f"{name:34s} score={score} vector={carrier.vector}")
        print(f"  preserves: {', '.join(carrier.preserves)}")
        print(f"  destroys: {', '.join(carrier.destroys)}")
    hist = Counter(scores.values())
    print(f"score_hist={dict(sorted(hist.items()))}")
    print(f"directed_3cycles={count_3cycles(edges)}")
    print(f"scc_sizes={scc_sizes(edges)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(edges)}")

    print_header("[7] Merged packet rule")
    print("A representation/Farey quotient is theorem-safe only after declaring:")
    print("  root Farey packet: p/q, q, excess, route labels")
    print("  summand shadow: p+q as antidiagonal fiber and additive/pinch row")
    print("  multiplicand shadow: p*q as factor/Kpq incidence fiber")
    print("  carry law: Pascal row and Zeckendorf no-adjacent normal form when used")
    print("  proof economy: smoothing, bounded arity, normal form, product incidence, or stress test")
    print("If one of those coordinates is forgotten, it must be reconstructed,")
    print("dual-annihilated, exact/coboundary, descended, or routed to a named residual.")


if __name__ == "__main__":
    main()
