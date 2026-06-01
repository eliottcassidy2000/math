#!/usr/bin/env python3
"""Whacky n=145 LRC reframes.

codex-2026-06-01 S518

The user asked for wild ways to reframe and maybe solve LRC at n=145.
This script treats n=145 as the total LRC denominator/runner count used in the
repo and builds a proof-atlas rather than attempting impossible enumeration.

The core arithmetic facts are:

    145 = 5 * 29
    moving runners = 144
    threshold = 1/145

For every q <= 145, THM-369 says that if no moving speed is divisible by q,
then a rational time with denominator q is lonely.  Thus a counterexample must
be a divisibility-covering object.  The weird n=145-specific idea is to read
the q=145 wall through CRT residues mod 5 and mod 29.  At t=a/145, every
nonzero residue mod 145 is safe; only residue 0 blocks.  If residue 0 appears,
then because there are only 144 moving runners, at least one antipodal boundary
residue pair is missing, creating a one-sided unit-wall aperture.  Whether the
zero-residue blockers can escape through that aperture is the remaining local
problem.

Tournament Analysis declaration:

Whacky route tournament
    vertices: proposed proof reframes for n=145.
    pairwise observable: exactness, n145 specificity, source-reachability fit,
        CRT leverage, pressure compatibility, novelty, and proof risk.
    switch/gauge: a route points to another when it wins the majority of these
        dimensions; proof risk is reversed because lower risk is better.
    tie Hamiltonian path: the listed route order.

Stored output:
    05-knowledge/results/lrc_n145_whacky_reframes_s518.out
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from functools import reduce
from itertools import combinations
from math import gcd
from typing import Callable


N = 145
MOVING = N - 1


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    value = n
    while d * d <= value:
        while value % d == 0:
            out[d] = out.get(d, 0) + 1
            value //= d
        d += 1
    if value > 1:
        out[value] = out.get(value, 0) + 1
    return out


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def phi(n: int) -> int:
    result = n
    for p in factor(n):
        result = result // p * (p - 1)
    return result


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def lcm_range(limit: int) -> int:
    return reduce(lcm, range(1, limit + 1), 1)


LCM_145 = lcm_range(N)


def int_summary(value: int) -> str:
    text = str(value)
    if len(text) <= 18:
        return text
    return f"{text[:8]}...{text[-8]} ({len(text)} digits)"


def initial() -> tuple[int, ...]:
    return tuple(range(1, N))


def unit_wall_with_zero() -> tuple[int, ...]:
    return tuple(range(1, N - 1)) + (N,)


def zero_heavy() -> tuple[int, ...]:
    return (1,) + tuple(N * k for k in range(1, MOVING))


def five_column() -> tuple[int, ...]:
    return (1,) + tuple(5 * k for k in range(1, MOVING))


def twentynine_column() -> tuple[int, ...]:
    return (1,) + tuple(29 * k for k in range(1, MOVING))


def lcm_spike_units() -> tuple[int, ...]:
    return (LCM_145,) + tuple(range(1, N - 1))


def crt_ladder_mix() -> tuple[int, ...]:
    nonunits = [r for r in range(1, N) if gcd(r, N) > 1]
    units = [r for r in range(1, N) if gcd(r, N) == 1]
    speeds: list[int] = []
    speeds.extend(nonunits)
    speeds.extend(units[:56])
    speeds.extend(u + N for u in units[:56])
    assert len(speeds) == MOVING
    assert len(set(speeds)) == MOVING
    return tuple(speeds)


FAMILIES: tuple[tuple[str, Callable[[], tuple[int, ...]]], ...] = (
    ("initial_1_to_144", initial),
    ("unit_wall_with_one_zero", unit_wall_with_zero),
    ("zero_heavy_one_unit", zero_heavy),
    ("five_column_one_unit", five_column),
    ("twentynine_column_one_unit", twentynine_column),
    ("lcm_spike_units", lcm_spike_units),
    ("crt_ladder_mix", crt_ladder_mix),
)


def residues(speeds: tuple[int, ...], modulus: int = N) -> Counter[int]:
    return Counter(speed % modulus for speed in speeds)


def divisibility_counts(speeds: tuple[int, ...], limit: int = N) -> dict[int, int]:
    return {q: sum(1 for speed in speeds if speed % q == 0) for q in range(2, limit + 1)}


def antipodal_profile(residue_counter: Counter[int]) -> dict[str, int]:
    full = half = empty = 0
    for r in range(1, (N + 1) // 2):
        a = residue_counter[r] > 0
        b = residue_counter[N - r] > 0
        if a and b:
            full += 1
        elif a or b:
            half += 1
        else:
            empty += 1
    return {"full": full, "half": half, "empty": empty}


def missing_residues(residue_counter: Counter[int]) -> list[int]:
    return [r for r in range(N) if residue_counter[r] == 0]


def rational_trap_summary(speeds: tuple[int, ...]) -> dict[str, object]:
    counts = divisibility_counts(speeds)
    min_count = min(counts.values())
    best_qs = [q for q, count in counts.items() if count == min_count]
    missing_qs = [q for q, count in counts.items() if count == 0]
    return {
        "min_count": min_count,
        "best_qs": best_qs[:12],
        "missing_qs": missing_qs[:12],
        "missing_total": len(missing_qs),
        "count_5": counts[5],
        "count_29": counts[29],
        "count_145": counts[145],
        "sieve_complete": not missing_qs,
    }


def family_profile(label: str, speeds: tuple[int, ...]) -> dict[str, object]:
    assert len(speeds) == MOVING
    primitive = reduce(gcd, speeds, 0) == 1
    res = residues(speeds)
    anti = antipodal_profile(res)
    missing = missing_residues(res)
    trap = rational_trap_summary(speeds)
    zero_count = res[0]
    unit_count = sum(count for residue, count in res.items() if gcd(residue, N) == 1)
    five_nonzero = sum(count for residue, count in res.items() if residue % 5 == 0 and residue != 0)
    twentynine_nonzero = sum(count for residue, count in res.items() if residue % 29 == 0 and residue != 0)
    shared_nonunit = MOVING - zero_count - unit_count - five_nonzero - twentynine_nonzero
    return {
        "label": label,
        "primitive": primitive,
        "max_speed": max(speeds),
        "max_speed_summary": int_summary(max(speeds)),
        "residue_support": sum(1 for count in res.values() if count),
        "zero_count": zero_count,
        "unit_count": unit_count,
        "five_nonzero": five_nonzero,
        "twentynine_nonzero": twentynine_nonzero,
        "shared_nonunit": shared_nonunit,
        "missing_residue_total": len(missing),
        "missing_residue_sample": missing[:12],
        "antipodal": anti,
        "trap": trap,
    }


def print_header() -> None:
    print("Whacky n=145 LRC reframes (S518)")
    print("=" * 88)
    print(f"n={N} = {factor(N)}; moving runners={MOVING}; threshold=1/{N}")
    print(f"divisors(n)={divisors(N)} phi(n)={phi(N)}")
    print(f"unit residues={phi(N)}; nonzero nonunits={N - 1 - phi(N)}")
    print(f"lcm(1..145)={int_summary(LCM_145)}")
    print()


def print_family_profiles() -> list[dict[str, object]]:
    print("Part A. Model families and exact n=145 trap data")
    print("-" * 88)
    profiles = [family_profile(label, make()) for label, make in FAMILIES]
    for profile in profiles:
        trap = profile["trap"]
        anti = profile["antipodal"]
        assert isinstance(trap, dict)
        assert isinstance(anti, dict)
        print(
            f"{profile['label']:<28} primitive={profile['primitive']} "
            f"max={profile['max_speed_summary']} support={profile['residue_support']} "
            f"zero={profile['zero_count']} units={profile['unit_count']} "
            f"5*={profile['five_nonzero']} 29*={profile['twentynine_nonzero']}"
        )
        print(
            f"  q-trap min blockers={trap['min_count']} best_qs={trap['best_qs']} "
            f"missing_qs={trap['missing_total']} sample={trap['missing_qs']}"
        )
        print(
            f"  q=5 blockers={trap['count_5']} q=29 blockers={trap['count_29']} "
            f"q=145 blockers={trap['count_145']} sieve_complete={trap['sieve_complete']}"
        )
        print(
            f"  antipodal pairs full/half/empty="
            f"{anti['full']}/{anti['half']}/{anti['empty']} "
            f"missing_residues={profile['missing_residue_total']} "
            f"sample={profile['missing_residue_sample']}"
        )
    print()
    return profiles


def print_unit_wall_lemma() -> None:
    print("Part B. The unit-wall aperture thought experiment")
    print("-" * 88)
    print("At t=a/145 with gcd(a,145)=1:")
    print("  * every nonzero residue mod 145 is at distance at least 1/145;")
    print("  * residue +a^{-1} and -a^{-1} are exactly on the two safe walls;")
    print("  * only speeds divisible by 145 block the observer.")
    print()
    print("Therefore:")
    print("  If no speed is 0 mod 145, THM-369/q=145 already gives a source wall.")
    print("  If a 0 mod 145 speed exists, the 144 moving runners cannot also")
    print("  occupy all 144 nonzero residues.  Hence some antipodal boundary pair")
    print("  {r,-r} has a missing side.  Choosing a with a^{-1}=r makes one wall")
    print("  unpinned: a one-sided aperture opens next to the unit source wall.")
    print()
    print("Whacky remaining local problem:")
    print("  Can the 0 mod 145 blockers be moved through this one-sided aperture")
    print("  before any pinned boundary runner enters the forbidden cap?")
    print("  This is a miniature observer-score/endpoint-pressure problem for the")
    print("  quotient speeds v/145 living above the unit wall.")
    print()


@dataclass(frozen=True)
class Route:
    name: str
    exactness: int
    n145: int
    source_fit: int
    crt: int
    pressure: int
    novelty: int
    risk: int


ROUTES = (
    Route("unit_wall_aperture", 5, 5, 5, 4, 3, 4, 2),
    Route("crt_5_29_two_moons", 4, 5, 4, 5, 3, 5, 3),
    Route("zero_residue_embryo", 4, 5, 5, 3, 5, 4, 3),
    Route("almost_source_tunnel", 5, 3, 5, 2, 4, 3, 2),
    Route("observer_score_descent", 5, 3, 5, 2, 5, 2, 2),
    Route("paley_failure_bipartite_qr", 2, 4, 2, 5, 1, 5, 5),
    Route("covering_code_dual", 3, 4, 3, 3, 4, 5, 4),
    Route("tropical_chip_firing", 2, 3, 3, 3, 5, 5, 5),
    Route("source_menu_percolation", 3, 3, 4, 2, 2, 4, 4),
)


DIMENSIONS: tuple[tuple[str, Callable[[Route], int], bool], ...] = (
    ("exactness", lambda route: route.exactness, True),
    ("n145", lambda route: route.n145, True),
    ("source_fit", lambda route: route.source_fit, True),
    ("crt", lambda route: route.crt, True),
    ("pressure", lambda route: route.pressure, True),
    ("novelty", lambda route: route.novelty, True),
    ("risk", lambda route: route.risk, False),
)


def route_tournament() -> list[list[int]]:
    adj = [[0] * len(ROUTES) for _ in ROUTES]
    for i, j in combinations(range(len(ROUTES)), 2):
        wins_i = wins_j = 0
        for _, getter, high_good in DIMENSIONS:
            vi = getter(ROUTES[i])
            vj = getter(ROUTES[j])
            if vi == vj:
                continue
            if (vi > vj) == high_good:
                wins_i += 1
            else:
                wins_j += 1
        winner = i if wins_i >= wins_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for start in range(n):
        dp[1 << start][start] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[full])


def triangle_count(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> tuple[int, ...]:
    n = len(adj)
    seen = [False] * n
    order: list[int] = []

    def dfs(vertex: int) -> None:
        seen[vertex] = True
        for nxt, bit in enumerate(adj[vertex]):
            if bit and not seen[nxt]:
                dfs(nxt)
        order.append(vertex)

    for vertex in range(n):
        if not seen[vertex]:
            dfs(vertex)
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    sizes: list[int] = []
    for start in reversed(order):
        if seen[start]:
            continue
        queue = deque([start])
        seen[start] = True
        size = 0
        while queue:
            vertex = queue.popleft()
            size += 1
            for nxt, bit in enumerate(rev[vertex]):
                if bit and not seen[nxt]:
                    seen[nxt] = True
                    queue.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def print_route_tournament() -> None:
    print("Part C. Whacky route Tournament Analysis")
    print("-" * 88)
    adj = route_tournament()
    scores = [sum(row) for row in adj]
    print(
        f"routes={len(ROUTES)} H={hamiltonian_paths(adj)} c3={triangle_count(adj)} "
        f"SCCs={scc_sizes(adj)} score_hist={dict(sorted(Counter(scores).items()))}"
    )
    ranked = sorted(((scores[idx], -idx, ROUTES[idx]) for idx in range(len(ROUTES))), reverse=True)
    for score, _, route in ranked:
        print(
            f"  {score:>2} {route.name:<28} exact={route.exactness} n145={route.n145} "
            f"source={route.source_fit} crt={route.crt} pressure={route.pressure} "
            f"novelty={route.novelty} risk={route.risk}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 88)
    print("1. The non-whacky base case is THM-369: any missing q<=145 gives a")
    print("   rational lonely time.  A hard n=145 system must be sieve-complete.")
    print("2. The whacky n=145-specific door is the unit wall t=a/145.  Nonzero")
    print("   residues are already safe there; only 145-divisible speeds block.")
    print("3. If a 145-divisible blocker exists, some antipodal nonzero residue pair")
    print("   is necessarily incomplete, because there are only 144 moving runners.")
    print("   That creates a one-sided aperture next to a unit source wall.")
    print("4. The attempted proof becomes: every zero-residue embryo can be pushed")
    print("   through some aperture, or else its labelled endpoint-pressure core")
    print("   violates the THM-380 peeling/cycle trilemma.")
    print("5. The wildest but most concrete reframing is therefore not full A000568")
    print("   enumeration at n=145.  It is a CRT aperture plus observer-score")
    print("   descent problem around denominator 145 walls.")


def main() -> None:
    print_header()
    print_family_profiles()
    print_unit_wall_lemma()
    print_route_tournament()
    print_synthesis()


if __name__ == "__main__":
    main()
