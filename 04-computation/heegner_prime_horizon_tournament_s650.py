"""S650: Heegner prime-polynomial horizons and tournament n-2 witnesses.

The user pointed at the prime-generating set {2,3,5,11,17,41}, the
Heegner/class-number-one set {1,2,3,7,11,19,43,67,163}, and the fact that
Euler's quadratic has a prime run through its first p-2 positive inputs.

This script keeps the roles typed:

* p is the constant term in Q_p(x)=x^2+x+p.
* d=4p-1 is the negative discriminant |disc(Q_p)|.
* x=1..p-2 are the interior inputs.
* x=p-1 is the forced boundary collapse Q_p(p-1)=p^2.
* a single reversed long edge in a p-vertex transitive tournament creates
  exactly p-2 cyclic triangles, by THM-410.
"""

from __future__ import annotations

from itertools import combinations
from math import isqrt


EULER_PRIMES = [2, 3, 5, 11, 17, 41]
HEEGNER = [1, 2, 3, 7, 11, 19, 43, 67, 163]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    r = isqrt(n)
    f = 3
    while f <= r:
        if n % f == 0:
            return False
        f += 2
    return True


def primes_upto(limit: int) -> list[int]:
    return [n for n in range(2, limit + 1) if is_prime(n)]


def q_value(p: int, x: int) -> int:
    return x * x + x + p


def first_composite_x(p: int) -> int:
    """First positive x with Q_p(x) composite."""

    x = 1
    while is_prime(q_value(p, x)):
        x += 1
    return x


def horizon(p: int) -> int:
    """Number of consecutive positive prime inputs x=1,2,..."""

    return first_composite_x(p) - 1


def interval_reversal_tournament(m: int) -> list[list[bool]]:
    """Transitive tournament on 0..m-1 with the long edge 0->m-1 reversed."""

    adj = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            adj[i][j] = True
    if m >= 2:
        adj[0][m - 1] = False
        adj[m - 1][0] = True
    return adj


def is_cyclic_triple(adj: list[list[bool]], triple: tuple[int, int, int]) -> bool:
    a, b, c = triple
    out = {a: 0, b: 0, c: 0}
    for u, v in combinations(triple, 2):
        if adj[u][v]:
            out[u] += 1
        else:
            out[v] += 1
    return sorted(out.values()) == [1, 1, 1]


def count_cyclic_triples(adj: list[list[bool]]) -> int:
    return sum(is_cyclic_triple(adj, t) for t in combinations(range(len(adj)), 3))


def hamiltonian_paths_dp(adj: list[list[bool]]) -> int:
    """Count Hamiltonian paths by Held-Karp DP. Used only for small controls."""

    n = len(adj)
    if n == 0:
        return 0
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            row = adj[last]
            for nxt in range(n):
                bit = 1 << nxt
                if mask & bit:
                    continue
                if row[nxt]:
                    dp[mask | bit][nxt] += count
    return sum(dp[-1])


def transitive_lens_tournament(labels: list[str]) -> dict[str, object]:
    n = len(labels)
    scores = {labels[i]: n - 1 - i for i in range(n)}
    score_hist = {s: 1 for s in range(n)}
    return {
        "vertices": labels,
        "score_hist": score_hist,
        "directed_3cycles": 0,
        "scc_count": n,
        "hamiltonian_paths": 1,
        "tie_path": " > ".join(labels),
        "scores": scores,
    }


def main() -> None:
    print("S650 Heegner prime-polynomial horizons and tournament n-2 witnesses")
    print("=" * 78)
    print("Refines HYP-2225 by adding the THM-410 long-edge witness split.")
    print()

    print("Euler/Rabinowitsch horizon check for Q_p(x)=x^2+x+p")
    print("-" * 78)
    hits = [p for p in primes_upto(500) if first_composite_x(p) == p - 1]
    print(f"Prime p <= 500 with first positive composite at x=p-1: {hits}")
    print(f"Expected Euler/Rabinowitsch set: {EULER_PRIMES}")
    print(f"Matches expected set? {hits == EULER_PRIMES}")
    print()

    print("Heegner discriminant map")
    print("-" * 78)
    mapped = [4 * p - 1 for p in EULER_PRIMES]
    print(f"p -> d=4p-1 gives: {mapped}")
    print(f"Heegner numbers:     {HEEGNER}")
    print(f"Heegner values not in the Euler-prime tail: {[d for d in HEEGNER if d not in mapped]}")
    print()

    print("Boundary/interior table")
    print("-" * 78)
    header = (
        "p  d=4p-1  heegner  positive_horizon  p-2  first_bad  "
        "Q(first_bad)  boundary_p^2  c3_long_edge  H_small"
    )
    print(header)
    for p in EULER_PRIMES:
        d = 4 * p - 1
        bad = first_composite_x(p)
        q_bad = q_value(p, bad)
        adj = interval_reversal_tournament(p)
        c3 = count_cyclic_triples(adj)
        h_small = hamiltonian_paths_dp(adj) if p <= 11 else "skip"
        print(
            f"{p:2d} {d:7d}  {str(d in HEEGNER):7s}"
            f" {horizon(p):16d} {p-2:4d} {bad:10d}"
            f" {q_bad:12d} {p*p:13d} {c3:13d} {h_small}"
        )
    print()

    print("Non-Heegner controls: the prime run breaks before the boundary")
    print("-" * 78)
    controls = [p for p in primes_upto(80) if p not in EULER_PRIMES][:12]
    print("p  d=4p-1  first_bad  Q(first_bad)  horizon  p-2")
    for p in controls:
        bad = first_composite_x(p)
        print(
            f"{p:2d} {4*p-1:7d} {bad:10d}"
            f" {q_value(p, bad):12d} {horizon(p):8d} {p-2:4d}"
        )
    print()

    print("THM-410 tournament witness interpretation")
    print("-" * 78)
    print("In a transitive p-vertex tournament, reverse only the long edge 0 -> p-1.")
    print("THM-410 says every interior vertex v with 0 < v < p-1 creates one")
    print("cyclic triangle 0 -> v -> p-1 -> 0, so #C3 = p-2.")
    print("For the reversed long edge p-1 -> 0, the pair-witness split is")
    print("sigma=0, lambda=p-2, delta=0: all p-2 witnesses are cyclic.")
    print("This is a structural horizon, not a Hamiltonian-path count.")
    print()

    print("Typed bridge")
    print("-" * 78)
    print("Q_p(0)=p is a boundary prime.")
    print("Q_p(1)..Q_p(p-2) are the p-2 interior primality witnesses.")
    print("Q_p(p-1)=p^2 is the forced boundary composite.")
    print("d=4p-1 is the class-number-one / unique-factorization side channel.")
    print("The tournament analogue is an interval with p-2 interior witnesses.")
    print()

    lenses = [
        "class_number_one_factorization",
        "euler_prime_horizon",
        "heegner_discriminant_map",
        "interval_reversal_long_edge",
        "pair_witness_n_minus_2",
        "first_composite_boundary",
        "raw_number_list",
    ]
    tournament = transitive_lens_tournament(lenses)
    print("Tournament Analysis")
    print("-" * 78)
    print("vertices =", ", ".join(tournament["vertices"]))
    print("pairwise observable = which lens preserves more typed information")
    print("switch/gauge = orient A -> B when A explains B without scalar collapse")
    print(f"tie Hamiltonian path = {tournament['tie_path']}")
    print(f"score_hist = {tournament['score_hist']}")
    print(f"directed_3cycles = {tournament['directed_3cycles']}")
    print(f"scc_count = {tournament['scc_count']}")
    print(f"hamiltonian_paths = {tournament['hamiltonian_paths']}")


if __name__ == "__main__":
    main()
