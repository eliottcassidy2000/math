#!/usr/bin/env python3
"""S642: Goldbach/Lemoine same-pair projection.

The user suggested reading even targets as unordered prime pairs and odd targets
as ordered prime-plus-doubled-prime pairs.  This script records the exact linear
carrier:

    E = p + q          (Goldbach/even projection, swap-blind)
    O = p + 2q         (Lemoine/odd projection, q is the doubled coordinate)

For odd primes p,q, E is even and O is odd.  The pair (E,O) reconstructs the
ordered pair:

    q = O - E,    p = 2E - O.

Thus "same even and odd use the same pair" is a simple bipartite graph between
even E and odd O.  The duplicate diagonal p=q maps p -> (2p,3p); in particular
p=7 gives the visible (14,21) shadow.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations_with_replacement


LIMIT = 300


def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        sieve[0] = 0
    if n >= 1:
        sieve[1] = 0
    p = 2
    while p * p <= n:
        if sieve[p]:
            step = p
            start = p * p
            sieve[start : n + 1 : step] = b"\x00" * (((n - start) // step) + 1)
        p += 1
    return [i for i, ok in enumerate(sieve) if ok]


PRIMES = primes_upto(2 * LIMIT)
PRIME_SET = set(PRIMES)
ODD_PRIMES = [p for p in PRIMES if p % 2 == 1]
ODD_PRIME_SET = set(ODD_PRIMES)


@dataclass(frozen=True)
class SamePairEdge:
    even: int
    odd: int
    p: int
    q: int
    diagonal: bool


@dataclass(frozen=True)
class Route:
    name: str
    preserves_pair: int
    exposes_symmetry: int
    handles_exceptions: int
    proof_transfer: int
    scalar_risk: int


ROUTES = [
    Route("invertible_EO_pair_projection", 5, 4, 4, 5, 1),
    Route("swap_reflection_O_to_3E_minus_O", 5, 5, 3, 4, 1),
    Route("diagonal_duplicate_fixed_locus", 4, 5, 3, 5, 1),
    Route("exceptional_prime_2_channel", 4, 3, 5, 4, 1),
    Route("degree_and_shadow_price_ledger", 4, 3, 4, 3, 2),
    Route("raw_goldbach_lemoine_counts", 2, 1, 2, 2, 4),
]


def goldbach_pairs(even: int, odd_only: bool = False) -> list[tuple[int, int]]:
    source = ODD_PRIMES if odd_only else PRIMES
    out = []
    for p in source:
        if p > even - p:
            break
        q = even - p
        if q in PRIME_SET and (not odd_only or q in ODD_PRIME_SET):
            out.append((p, q))
    return out


def lemoine_pairs(odd: int, odd_q_only: bool = False) -> list[tuple[int, int]]:
    out = []
    for q in PRIMES:
        if odd_q_only and q == 2:
            continue
        p = odd - 2 * q
        if p < 2:
            break
        if p in PRIME_SET:
            out.append((p, q))
    return out


def same_pair_edges(limit: int) -> list[SamePairEdge]:
    edges: list[SamePairEdge] = []
    for p, q in combinations_with_replacement(ODD_PRIMES, 2):
        even = p + q
        if even > limit:
            continue
        shadows = [(q + 2 * p, q, p), (p + 2 * q, p, q)]
        for odd, first, doubled in shadows:
            if odd > limit:
                continue
            if first <= 0 or doubled <= 0:
                continue
            edge = SamePairEdge(
                even=even,
                odd=odd,
                p=first,
                q=doubled,
                diagonal=first == doubled,
            )
            if edge not in edges:
                edges.append(edge)
    return sorted(edges, key=lambda x: (x.even, x.odd, x.p, x.q))


def same_pair_from_eo(even: int, odd: int) -> tuple[int, int] | None:
    q = odd - even
    p = 2 * even - odd
    if p in ODD_PRIME_SET and q in ODD_PRIME_SET:
        return (p, q)
    return None


def route_score(route: Route) -> int:
    return (
        3 * route.preserves_pair
        + 2 * route.exposes_symmetry
        + 2 * route.proof_transfer
        + route.handles_exceptions
        - 2 * route.scalar_risk
    )


def beats(a: Route, b: Route) -> bool:
    sa = route_score(a)
    sb = route_score(b)
    if sa != sb:
        return sa > sb
    return a.name < b.name


def hamiltonian_paths(routes: list[Route]) -> int:
    n = len(routes)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(routes[last], routes[nxt]):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, i), 0) for i in range(n))


def directed_triangles(routes: list[Route]) -> int:
    total = 0
    for a in range(len(routes)):
        for b in range(a + 1, len(routes)):
            for c in range(b + 1, len(routes)):
                wins = []
                for x, y in ((a, b), (a, c), (b, c)):
                    wins.append((x, y) if beats(routes[x], routes[y]) else (y, x))
                out = Counter(x for x, _ in wins)
                if sorted(out.values()) == [1, 1, 1]:
                    total += 1
    return total


def print_definitions() -> None:
    print("S642 Goldbach/Lemoine same-pair projection")
    print("==========================================")
    print()
    print("Definitions")
    print("-----------")
    print("For odd primes p,q:")
    print("  even projection E = p + q")
    print("  odd projection  O = p + 2q, where q is the doubled coordinate")
    print("The two projections invert the ordered pair:")
    print("  q = O - E")
    print("  p = 2E - O")
    print("So an even E and odd O share a pair exactly when O-E and 2E-O")
    print("are odd primes.  Necessarily E < O < 2E.")
    print()
    print("Swapping the pair fixes E and reflects odd shadows:")
    print("  (p,q) <-> (q,p) sends O <-> 3E - O.")
    print("The duplicate diagonal p=q is the fixed locus O = 3E/2.")
    print()


def print_duplicate_diagonal() -> None:
    print("Duplicate diagonal p -> (2p,3p)")
    print("--------------------------------")
    rows = []
    for p in ODD_PRIMES:
        if 3 * p > LIMIT:
            break
        mark = "  <== p=7 gives the 14/21 shadow" if p == 7 else ""
        rows.append(f"  p={p:<3} even=2p={2*p:<3} odd=3p={3*p:<3}{mark}")
    print("\n".join(rows[:18]))
    print()


def print_even_shadow_samples(edges: list[SamePairEdge]) -> None:
    print("Even targets and their odd same-pair shadows")
    print("--------------------------------------------")
    by_even: dict[int, list[SamePairEdge]] = defaultdict(list)
    for edge in edges:
        by_even[edge.even].append(edge)
    for even in range(6, 32, 2):
        shadows = by_even.get(even, [])
        if not shadows:
            continue
        pieces = []
        for edge in shadows:
            tag = "diag" if edge.diagonal else "off"
            pieces.append(f"O={edge.odd} via ({edge.p},{edge.q}) {tag}")
        print(f"  E={even:<2}: " + "; ".join(pieces))
    print()


def print_coverage(edges: list[SamePairEdge]) -> None:
    print("Coverage and exception channels through LIMIT=300")
    print("-------------------------------------------------")
    goldbach_missing = [n for n in range(4, LIMIT + 1, 2) if not goldbach_pairs(n)]
    lemoine_missing = [n for n in range(7, LIMIT + 1, 2) if not lemoine_pairs(n)]
    compatible_missing = [
        n for n in range(7, LIMIT + 1, 2) if not lemoine_pairs(n, odd_q_only=True)
    ]
    q2_only = [
        n
        for n in range(7, LIMIT + 1, 2)
        if lemoine_pairs(n) and not lemoine_pairs(n, odd_q_only=True)
    ]
    print(f"  Goldbach misses among even 4..{LIMIT}: {goldbach_missing}")
    print(f"  Lemoine misses among odd 7..{LIMIT}: {lemoine_missing}")
    print(f"  Odd-prime-compatible Lemoine misses: {compatible_missing[:20]}")
    print(f"  q=2-only odd exceptions: {q2_only[:20]}")
    print("The q=2 channel is the exceptional even-prime seam: it proves small")
    print("Lemoine rows but has no same odd-prime Goldbach pair E=p+q.")
    print()

    even_degree = Counter(edge.even for edge in edges)
    odd_degree = Counter(edge.odd for edge in edges)
    top_even = even_degree.most_common(8)
    top_odd = odd_degree.most_common(8)
    print("Same-pair bipartite degrees, using odd-prime pairs only:")
    print(f"  edges={len(edges)}")
    print(f"  top even degrees: {top_even}")
    print(f"  top odd degrees : {top_odd}")
    print()


def print_reflection_checks() -> None:
    print("Reflection checks")
    print("-----------------")
    for even, odd in [(14, 17), (14, 21), (14, 25), (22, 25), (22, 33)]:
        pair = same_pair_from_eo(even, odd)
        reflected = 3 * even - odd
        pair_ref = same_pair_from_eo(even, reflected)
        print(
            f"  E={even:<2}, O={odd:<2}: pair={pair}, "
            f"reflected O'={reflected:<2}, reflected_pair={pair_ref}"
        )
    print()


def print_tournament_analysis() -> None:
    print("Tournament Analysis over proof lenses")
    print("-------------------------------------")
    routes = list(ROUTES)
    scores = {r.name: sum(1 for other in routes if r is not other and beats(r, other)) for r in routes}
    hist = Counter(scores.values())
    ranking = sorted(routes, key=lambda r: (-scores[r.name], r.name))
    print(f"vertices={len(routes)}")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print(f"directed_3cycles={directed_triangles(routes)}")
    print(f"hamiltonian_paths={hamiltonian_paths(routes)}")
    print("tie Hamiltonian path:")
    for route in ranking:
        print(f"  score={scores[route.name]} {route.name}")
    print()
    print("Pairwise observable: which lens better preserves the same-pair predicate")
    print("between Goldbach and Lemoine without discarding swap, diagonal, or p=2")
    print("side channels.  The vertices are proof lenses, not numbers.")


def main() -> None:
    edges = same_pair_edges(LIMIT)
    print_definitions()
    print_duplicate_diagonal()
    print_even_shadow_samples(edges)
    print_coverage(edges)
    print_reflection_checks()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
