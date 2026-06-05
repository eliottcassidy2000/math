#!/usr/bin/env python3
"""S651: number-theory tournament carrier atlas.

This session follows S649/S650.  The goal is not to pretend that tournaments
solve hard number-theory problems by themselves.  The goal is to make precise
where they fit:

* Paley tournaments are literal finite-field/character objects.
* Sieve tournaments orient local primes by obstruction pressure.
* Horizon tournaments orient boundary/interior witness ledgers.
* Hard-problem transfer tournaments rank proof carriers by retained side
  channels, not by scalar numerology.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations
from math import isfinite, isqrt, log


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d = 5
    step = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += step
        step = 6 - step
    return True


def primes_upto(limit: int) -> list[int]:
    return [n for n in range(2, limit + 1) if is_prime(n)]


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                bit = 1 << nxt
                if mask & bit:
                    continue
                if adj[last][nxt]:
                    dp[mask | bit][nxt] += count
    return sum(dp[-1])


def tournament_fingerprints(adj: list[list[bool]]) -> dict[str, object]:
    n = len(adj)
    scores = [sum(row) for row in adj]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        out = [0, 0, 0]
        triple = [a, b, c]
        for i, j in combinations(range(3), 2):
            u, v = triple[i], triple[j]
            if adj[u][v]:
                out[i] += 1
            else:
                out[j] += 1
        if sorted(out) == [1, 1, 1]:
            c3 += 1

    # Strong components.
    def reach_from(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            u = q.popleft()
            for v in range(n):
                edge = adj[v][u] if reverse else adj[u][v]
                if edge and v not in seen:
                    seen.add(v)
                    q.append(v)
        return seen

    remaining = set(range(n))
    scc_sizes: list[int] = []
    while remaining:
        start = min(remaining)
        scc = reach_from(start) & reach_from(start, reverse=True)
        scc_sizes.append(len(scc))
        remaining -= scc

    h = hamiltonian_paths(adj) if n <= 11 else "skip"
    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_paths": h,
    }


def transitive_tournament_from_scores(scores: list[tuple[str, float]]) -> list[list[bool]]:
    """Orient higher score to lower score; deterministic name tiebreak."""
    ordered = sorted(scores, key=lambda item: (-item[1], item[0]))
    rank = {name: i for i, (name, _score) in enumerate(ordered)}
    n = len(scores)
    adj = [[False] * n for _ in range(n)]
    names = [name for name, _score in scores]
    for i, j in combinations(range(n), 2):
        if rank[names[i]] < rank[names[j]]:
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def paley_tournament(p: int) -> list[list[bool]]:
    residues = {x * x % p for x in range(1, p)}
    adj = [[False] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in residues:
                adj[i][j] = True
    return adj


@dataclass(frozen=True)
class LocalRow:
    prime: int
    forbidden: int
    survivors: int
    density: float
    weight: float
    note: str


def local_row(prime: int, forbidden_residues: set[int], note: str) -> LocalRow:
    forbidden = len(forbidden_residues)
    survivors = prime - forbidden
    density = survivors / prime
    weight = float("inf") if survivors == 0 else -log(density)
    return LocalRow(prime, forbidden, survivors, density, weight, note)


def twin_rows(limit: int) -> list[LocalRow]:
    rows = []
    for ell in primes_upto(limit):
        forbidden = {0 % ell, (-2) % ell}
        rows.append(local_row(ell, forbidden, "n and n+2 both prime"))
    return rows


def goldbach_rows(even_n: int, limit: int) -> list[LocalRow]:
    rows = []
    for ell in primes_upto(limit):
        forbidden = {0, even_n % ell}
        note = "one forbidden lane" if even_n % ell == 0 else "two forbidden lanes"
        rows.append(local_row(ell, forbidden, note))
    return rows


def polynomial_rows(p: int, limit: int) -> list[LocalRow]:
    rows = []
    for ell in primes_upto(limit):
        forbidden = {
            x for x in range(ell)
            if (x * x + x + p) % ell == 0
        }
        note = f"roots of x^2+x+{p} mod ell"
        rows.append(local_row(ell, forbidden, note))
    return rows


def obstruction_tournament(rows: list[LocalRow]) -> dict[str, object]:
    scores = [(str(row.prime), row.weight if isfinite(row.weight) else 999.0) for row in rows]
    adj = transitive_tournament_from_scores(scores)
    fp = tournament_fingerprints(adj)
    ranking = [
        row.prime for row in sorted(
            rows,
            key=lambda row: (-(row.weight if isfinite(row.weight) else 999.0), row.prime),
        )
    ]
    fp["ranking"] = ranking
    return fp


def edge_flip_count(a: list[LocalRow], b: list[LocalRow]) -> int:
    by_a = {row.prime: row.weight for row in a}
    by_b = {row.prime: row.weight for row in b}
    primes = sorted(set(by_a) & set(by_b))
    flips = 0
    for x, y in combinations(primes, 2):
        ax = by_a[x] > by_a[y] or (by_a[x] == by_a[y] and x < y)
        bx = by_b[x] > by_b[y] or (by_b[x] == by_b[y] and x < y)
        if ax != bx:
            flips += 1
    return flips


@dataclass(frozen=True)
class MethodLens:
    name: str
    local_witness: int
    side_channel: int
    finite_probe: int
    hard_problem_reach: int
    overclaim_risk: int


METHODS = (
    MethodLens("local_sieve_obstruction_tournaments", 5, 5, 5, 5, 1),
    MethodLens("character_phase_paley_tournaments", 5, 5, 4, 4, 1),
    MethodLens("endpoint_horizon_witness_tournaments", 5, 4, 5, 4, 1),
    MethodLens("valuation_drift_tournaments", 4, 4, 4, 4, 2),
    MethodLens("explicit_formula_balance_tournaments", 4, 5, 3, 5, 3),
    MethodLens("local_global_rank_tournaments", 3, 5, 3, 5, 3),
    MethodLens("raw_conjecture_numerology", 1, 1, 1, 1, 5),
)


def method_score(lens: MethodLens) -> int:
    return (
        3 * lens.local_witness
        + 3 * lens.side_channel
        + 2 * lens.finite_probe
        + 2 * lens.hard_problem_reach
        - 4 * lens.overclaim_risk
    )


def method_tournament() -> dict[str, object]:
    names = [lens.name for lens in METHODS]
    scores = [(lens.name, method_score(lens)) for lens in METHODS]
    adj = transitive_tournament_from_scores(scores)
    fp = tournament_fingerprints(adj)
    fp["ranking"] = [name for name, _ in sorted(scores, key=lambda item: (-item[1], item[0]))]
    return fp


def print_local_table(title: str, rows: list[LocalRow], top: int = 8) -> None:
    print(title)
    print("-" * 78)
    print("ell  forbidden  survivors  density   weight    note")
    for row in sorted(rows, key=lambda r: (-(r.weight if isfinite(r.weight) else 999.0), r.prime))[:top]:
        weight = "inf" if not isfinite(row.weight) else f"{row.weight:.4f}"
        print(
            f"{row.prime:3d} {row.forbidden:10d} {row.survivors:10d}"
            f" {row.density:8.4f} {weight:>8s}  {row.note}"
        )
    fp = obstruction_tournament(rows)
    print(f"ranking_by_obstruction={fp['ranking'][:top]}")
    print(f"score_hist={fp['score_hist']}, directed_3cycles={fp['directed_3cycles']}, "
          f"scc_sizes={fp['scc_sizes']}, H={fp['hamiltonian_paths']}")
    print()


def main() -> None:
    print("S651 number-theory tournament carrier atlas")
    print("=" * 78)
    print("Incoming signal: S650/HYP-2226 made the Heegner p-2 horizon exact via THM-410.")
    print("This session generalizes the move: tournamentize local witnesses, not raw integers.")
    print()

    print("Paley character tournaments: finite fields already are tournaments")
    print("-" * 78)
    print("p  residues        score_hist       c3   SCC       H")
    for p in (3, 7, 11, 19):
        adj = paley_tournament(p)
        fp = tournament_fingerprints(adj)
        residues = sorted({x * x % p for x in range(1, p)})
        print(
            f"{p:2d} {str(residues):14s} {str(fp['score_hist']):15s}"
            f" {fp['directed_3cycles']:4d} {str(fp['scc_sizes']):9s} {fp['hamiltonian_paths']}"
        )
    print()

    local_limit = 67
    twin = twin_rows(local_limit)
    goldbach_210 = goldbach_rows(210, local_limit)
    goldbach_2110 = goldbach_rows(2110, local_limit)
    poly_41 = polynomial_rows(41, local_limit)

    print_local_table("Twin-prime local obstruction tournament for gap 2", twin)
    print_local_table("Goldbach local obstruction tournament for N=210", goldbach_210)
    print_local_table("Goldbach local obstruction tournament for N=2110", goldbach_2110)
    print_local_table("Euler polynomial local obstruction tournament for p=41", poly_41)

    print("Cross-sieve edge flips among local-prime obstruction tournaments")
    print("-" * 78)
    print(f"twin gap 2 vs Goldbach N=210:  {edge_flip_count(twin, goldbach_210)} flips")
    print(f"twin gap 2 vs Goldbach N=2110: {edge_flip_count(twin, goldbach_2110)} flips")
    print(f"Goldbach 210 vs Goldbach 2110: {edge_flip_count(goldbach_210, goldbach_2110)} flips")
    print(f"twin gap 2 vs Euler p=41:      {edge_flip_count(twin, poly_41)} flips")
    print("Interpretation: the same local primes reorder when the retained side channel")
    print("changes.  That is the tournament signal, not noise.")
    print()

    print("Hard-problem morale ledger")
    print("-" * 78)
    print("problem        vertices to tournamentize          concrete finite morale move")
    print("Twin primes    local primes/residue lanes          rank bottleneck primes; compare gaps")
    print("Goldbach       local primes vs target N            track one-lane vs two-lane primes")
    print("RH             prime-power/zero phase terms        orient explicit-formula contributors")
    print("BSD            bad primes/root-number factors      orient local-global rank obligations")
    print("Collatz        residue classes/valuation jumps     orient drift witnesses by descent")
    print("Heegner        boundary/interior proof slots       use S650 endpoint horizon witnesses")
    print("LRC/UD         shell or spine/bulk side channels   keep quotient-loss labels attached")
    print()

    fp = method_tournament()
    print("Tournament Analysis")
    print("-" * 78)
    print("challenged vertices = integers, primes, residues, valuations, characters,")
    print("  Euler factors, zero terms, class groups, proof obligations, and side channels")
    print("chosen vertices = proof-carrier methods")
    print("pairwise observable = which method preserves more local witness information")
    print("switch/gauge = orient A -> B when A gives smaller finite probes with less scalar collapse")
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_paths={fp['hamiltonian_paths']}")
    print("tie Hamiltonian path:")
    for i, name in enumerate(fp["ranking"], start=1):
        print(f"  {i}. {name}")
    print()

    print("Hypothesis generated")
    print("  HYP-2227: number-theory tournaments are local-witness carrier quotients.")
    print("  They give morale progress by turning impossible global claims into finite")
    print("  side-channel ledgers: what is the vertex set, what edge observable is kept,")
    print("  and which quotient would destroy the obstruction?")


if __name__ == "__main__":
    main()
