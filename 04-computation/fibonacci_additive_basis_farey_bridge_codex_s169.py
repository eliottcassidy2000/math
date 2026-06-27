#!/usr/bin/env python3
"""Fibonacci/additive-basis/Farey bridge scout.

This is a synthesis computation, not a proof.  It makes the user's Fibonacci
row arrangement exact, then places it beside three representation economies:
Goldbach smoothing, Fermat polygonal bounded arity, and Zeckendorf normal form.
It also records how the four Farey payloads from HYP-2984 behave on the
LRC14 unit-excess chain.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, permutations
from math import comb, isqrt, log10


def fibonacci_numbers(n: int) -> list[int]:
    """Return standard F_1..F_n with F_1=F_2=1."""
    if n <= 0:
        return []
    fib = [1]
    if n == 1:
        return fib
    fib.append(1)
    while len(fib) < n:
        fib.append(fib[-1] + fib[-2])
    return fib


def fibonacci_rank_row(n: int) -> list[int]:
    """Coefficients in F_n = sum_k binom(n-k-1,k)."""
    return [comb(n - k - 1, k) for k in range((n - 1) // 2 + 1)]


def zeckendorf(n: int) -> list[int]:
    """Greedy Zeckendorf values using 1,2,3,5,... atoms."""
    if n <= 0:
        return []
    fib = [1, 2]
    while fib[-1] + fib[-2] <= n:
        fib.append(fib[-1] + fib[-2])
    out: list[int] = []
    r = n
    for f in reversed(fib):
        if f <= r:
            out.append(f)
            r -= f
    return out


def polygonal(s: int, m: int) -> int:
    return ((s - 2) * m * m - (s - 4) * m) // 2


def polygonal_depths(limit: int = 250, sides=range(3, 9)) -> dict[int, tuple[int, list[int]]]:
    """Minimum number of s-gonal atoms needed up to limit."""
    ans: dict[int, tuple[int, list[int]]] = {}
    for s in sides:
        atoms = []
        m = 1
        while True:
            val = polygonal(s, m)
            if val > limit:
                break
            atoms.append(val)
            m += 1
        dp = [10**9] * (limit + 1)
        dp[0] = 0
        for n in range(1, limit + 1):
            dp[n] = min((dp[n - a] + 1 for a in atoms if a <= n), default=10**9)
        mx = max(dp[1:])
        witnesses = [i for i in range(1, limit + 1) if dp[i] == mx][:12]
        ans[s] = (mx, witnesses)
    return ans


def primes_upto(n: int) -> list[int]:
    sieve = [True] * (n + 1)
    sieve[:2] = [False, False]
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start : n + 1 : step] = [False] * (((n - start) // step) + 1)
    return [i for i, ok in enumerate(sieve) if ok]


def goldbach_counts(samples_even=(20, 50, 100, 200), samples_odd=(101, 501, 999)):
    maxn = max(max(samples_even), max(samples_odd))
    primes = primes_upto(maxn)
    ps = set(primes)
    binary = {}
    for n in samples_even:
        binary[n] = sum(1 for p in primes if p <= n - p and n - p in ps)
    ternary = {}
    for n in samples_odd:
        count = 0
        for i, p in enumerate(primes):
            if p > n:
                break
            for q in primes[i:]:
                r = n - p - q
                if r < q:
                    break
                if r in ps:
                    count += 1
        ternary[n] = count
    return binary, ternary


def digit_count(n: int) -> int:
    if n == 0:
        return 1
    return int(log10(n)) + 1


def farey_unit_excess_payloads(n: int = 14, pmax: int = 8) -> list[dict[str, object]]:
    rows = []
    for p in range(1, pmax + 1):
        q = n * p - 1
        rows.append(
            {
                "p": p,
                "q": q,
                "e": n * p - q,
                "sum": p + q,
                "product": p * q,
                "q_power_p_digits": digit_count(q**p),
                "p_power_q_digits": digit_count(p**q) if p > 1 else 1,
                "sum_zeck": zeckendorf(p + q),
                "product_zeck": zeckendorf(p * q),
                "route": "AP/GW parent" if p == 1 else "C27 petal/two-block" if p == 2 else "K33/state-lift/Fejer",
            }
        )
    return rows


def tournament_from_vectors(names: list[str], vectors: dict[str, tuple[int, ...]]):
    """Orient i->j by coordinate majority, with names order as Hamiltonian tie path."""
    n = len(names)
    tie_order = {name: idx for idx, name in enumerate(names)}
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            wins = sum(x > y for x, y in zip(vectors[a], vectors[b]))
            losses = sum(x < y for x, y in zip(vectors[a], vectors[b]))
            if wins > losses or (wins == losses and tie_order[a] < tie_order[b]):
                adj[i][j] = True
    return adj


def tournament_fingerprint(names: list[str], adj: list[list[bool]]) -> dict[str, object]:
    n = len(names)
    scores = [sum(row) for row in adj]
    score_hist = dict(sorted(Counter(scores).items()))
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1

    reach = [[adj[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    seen = [False] * n
    scc_sizes = []
    for i in range(n):
        if seen[i]:
            continue
        comp = [j for j in range(n) if reach[i][j] and reach[j][i]]
        for j in comp:
            seen[j] = True
        scc_sizes.append(len(comp))

    hp_count = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
            hp_count += 1

    leaders = sorted(zip(scores, names), reverse=True)
    return {
        "score_hist": score_hist,
        "directed_3cycles": c3,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_path_count": hp_count,
        "top_path_by_score": [name for _, name in leaders],
    }


def main() -> None:
    print("# Fibonacci / Additive-Basis / Farey Bridge Scout")
    print()
    print("## Fibonacci Rows")
    print("Identity: F_n = sum_k binom(n-k-1,k).")
    print("For n>=2 this is the rank vector of independent sets in path P_{n-2}.")
    print()
    fib = fibonacci_numbers(14)
    for n in range(1, 12):
        row = fibonacci_rank_row(n)
        expr = "+".join(str(x) for x in row)
        print(f"n={n:2d} F_n={fib[n-1]:3d} row={expr:<18} zeck(F_n)={zeckendorf(fib[n-1])}")
    print()

    print("## Fermat Polygonal Bounded-Arity Check")
    for s, (mx, witnesses) in polygonal_depths().items():
        print(f"{s}-gonal atoms: max minimum summands through 250 = {mx}; first max witnesses={witnesses}")
    print()

    print("## Goldbach / Ternary Goldbach Sample Entropy")
    binary, ternary = goldbach_counts()
    print(f"binary unordered prime-pair counts: {binary}")
    print(f"ternary unordered prime-triple counts: {ternary}")
    print()

    print("## LRC14 Unit-Excess Farey Payloads")
    print("Rows are p/q = p/(14p-1), so e=14p-q=1.")
    for row in farey_unit_excess_payloads():
        print(
            "p={p:2d} q={q:3d} p+q={sum:3d} pq={product:4d} "
            "digits(q^p)={q_power_p_digits:2d} digits(p^q)={p_power_q_digits:3d} "
            "route={route}; Zeck(p+q)={sum_zeck}; Zeck(pq)={product_zeck}".format(**row)
        )
    print()

    print("## Three Representation Currencies")
    print("Goldbach: high-entropy sieve branch count.")
    print("Fermat polygonal: bounded arity / residue absorber.")
    print("Zeckendorf/Fibonacci: path-rank normal form / no-adjacent carry.")
    print("Farey sum/product/power: affine check / incidence blow-up / magnitude stress.")
    print()

    print("## Tournament Analysis")
    print("Vertices are proof currencies, not runners.")
    print("Pairwise observable: retained currency vector")
    print("(quotient_safety, local_to_global, entropy_control, LRC_transfer, exact_normal_form).")
    print("Gauge: coordinate majority; tie Hamiltonian path is listed order below.")
    names = [
        "zeckendorf_normal_form",
        "fibonacci_path_rank",
        "farey_sum_affine_check",
        "fermat_polygonal_bounded_arity",
        "ternary_goldbach_smoothing",
        "farey_product_incidence",
        "binary_goldbach_sieve",
        "farey_power_stress_test",
    ]
    vectors = {
        "zeckendorf_normal_form": (5, 5, 1, 5, 5),
        "fibonacci_path_rank": (5, 4, 3, 5, 4),
        "farey_sum_affine_check": (5, 3, 2, 5, 4),
        "fermat_polygonal_bounded_arity": (4, 4, 3, 4, 3),
        "ternary_goldbach_smoothing": (3, 4, 5, 3, 2),
        "farey_product_incidence": (3, 3, 3, 4, 3),
        "binary_goldbach_sieve": (2, 3, 4, 2, 1),
        "farey_power_stress_test": (1, 1, 1, 3, 1),
    }
    adj = tournament_from_vectors(names, vectors)
    fp = tournament_fingerprint(names, adj)
    print(f"tie_path={names}")
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"score_path={fp['top_path_by_score']}")
    print()
    print("Main readout: Fibonacci rows are not merely another sequence.")
    print("They are the rank-polynomial layer between abundant hypergraphs and")
    print("Zeckendorf uniqueness.  In LRC/Farey language, p+q is the additive")
    print("safe quotient, pq is the incidence/product side channel, and powers")
    print("are useful mostly because they reveal magnitude-blind false quotients.")


if __name__ == "__main__":
    main()
