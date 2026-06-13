#!/usr/bin/env python3
"""LRC proof-route landscape and Tournament Analysis atlas.

This script records the 2025-2026 Lonely Runner proof surge as a finite
route-comparison object.  The tournament lift is deliberately meta-level:
vertices are proof routes, pairwise observables are route-feature comparisons,
and the gauges ask which route is more useful for either repo action or
mathematical credence.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Callable


@dataclass(frozen=True)
class Route:
    code: str
    year: int
    title: str
    status: str
    scope: str
    url: str
    total_runners: int
    trust: int
    direct_lrc: int
    generality: int
    finite_check: int
    sieve: int
    endpoint: int
    zonotope: int
    spectral: int
    mixed: int
    reproducible: int
    repo_transfer: int
    frontier: int
    note: str


ROUTES: list[Route] = [
    Route(
        "BS07",
        2008,
        "Barajas-Serra seven runners",
        "peer-reviewed proof",
        "fixed case: seven total runners",
        "https://arxiv.org/abs/0710.4495",
        7,
        5,
        5,
        2,
        2,
        1,
        2,
        0,
        0,
        0,
        3,
        4,
        1,
        "distance-graph/chromatic-number technology; important baseline",
    ),
    Route(
        "Tao18",
        2018,
        "Tao finite checking and Bohr-set gap",
        "peer-reviewed proof machinery",
        "general reduction plus asymptotic lower bound",
        "https://doi.org/10.55016/ojs/cdm.v13i2.62728",
        0,
        5,
        4,
        5,
        5,
        2,
        2,
        0,
        4,
        0,
        3,
        7,
        4,
        "finite checking at n^{O(n^2)} and log-improved loneliness gap",
    ),
    Route(
        "PS60",
        2024,
        "Perarnau-Serra survey",
        "survey",
        "panoramic map and open problems",
        "https://arxiv.org/abs/2409.20160",
        0,
        4,
        2,
        4,
        2,
        1,
        2,
        2,
        1,
        1,
        2,
        6,
        4,
        "best orientation source for terminology and open routes",
    ),
    Route(
        "MSS25",
        2025,
        "Malikiosis-Santos-Schymura finite checking",
        "peer-reviewed proof machinery",
        "general linearly-exponential finite reduction",
        "https://arxiv.org/abs/2411.06903",
        0,
        5,
        4,
        5,
        5,
        2,
        3,
        5,
        0,
        0,
        4,
        8,
        7,
        "zonotopal reduction to speeds up to binomial(n+1,2)^(n-1)",
    ),
    Route(
        "ACS25",
        2025,
        "Alcantara-Criado-Santos shifted five runners",
        "prepublished computational proof",
        "shifted variant: five total runners",
        "https://arxiv.org/abs/2506.13379",
        5,
        4,
        1,
        3,
        4,
        0,
        2,
        5,
        0,
        1,
        5,
        5,
        5,
        "dyadic fundamental-domain covering-radius certificates",
    ),
    Route(
        "R8",
        2025,
        "Rosenfeld eight runners",
        "arXiv computer-assisted proof",
        "fixed case: eight total runners",
        "https://arxiv.org/abs/2509.14111",
        8,
        3,
        5,
        3,
        5,
        5,
        4,
        1,
        0,
        0,
        4,
        9,
        6,
        "minimal counterexample product bound plus prime-divisibility forcing",
    ),
    Route(
        "Bedert25",
        2025,
        "Bedert Riesz-product wider gap",
        "arXiv proof machinery",
        "general lower bound, not full LRC",
        "https://arxiv.org/abs/2511.16636",
        0,
        3,
        3,
        5,
        1,
        0,
        3,
        0,
        5,
        0,
        3,
        8,
        7,
        "nonnegative Riesz-product tests improve 1/(2n) by n^{-5/3+o(1)}",
    ),
    Route(
        "T9T10",
        2025,
        "Trakulthongchai nine and ten runners",
        "arXiv; EJC to appear",
        "fixed cases: nine and ten total runners",
        "https://arxiv.org/abs/2511.22427",
        10,
        4,
        5,
        3,
        5,
        5,
        4,
        1,
        0,
        0,
        5,
        10,
        8,
        "fiber sieve makes Rosenfeld prime-divisibility checks tractable",
    ),
    Route(
        "R9",
        2025,
        "Rosenfeld nine runners",
        "arXiv computer-assisted proof",
        "fixed case: nine total runners",
        "https://arxiv.org/abs/2512.01912",
        9,
        3,
        5,
        3,
        5,
        5,
        4,
        1,
        0,
        0,
        4,
        8,
        6,
        "independent refinement with small prime powers and SAT/set-cover view",
    ),
    Route(
        "BCS26",
        2026,
        "Blanco-Criado-Santos shifted counterexamples",
        "arXiv warning theorem",
        "shifted LRC and lonely-vector property fail",
        "https://arxiv.org/abs/2603.24784",
        0,
        3,
        1,
        4,
        2,
        0,
        2,
        5,
        0,
        1,
        3,
        5,
        6,
        "shows shifted analogies cannot be imported into original LRC blindly",
    ),
    Route(
        "ST13",
        2026,
        "Sungkawichai-Trakulthongchai eleven to thirteen",
        "arXiv computer-assisted proof",
        "fixed cases: eleven to thirteen total runners",
        "https://arxiv.org/abs/2604.23906",
        13,
        3,
        5,
        3,
        5,
        5,
        4,
        1,
        0,
        0,
        5,
        10,
        10,
        "intermediate sieves plus polynomial method for the (1,...,k) residue",
    ),
    Route(
        "Jensen26",
        2026,
        "Jensen mixed thresholds",
        "arXiv proof machinery",
        "mixed-threshold variant",
        "https://arxiv.org/abs/2605.27941",
        0,
        2,
        2,
        4,
        1,
        1,
        4,
        1,
        4,
        5,
        2,
        7,
        8,
        "exact MLPS_2 polygon and unequal-threshold Fourier formulas",
    ),
    Route(
        "Repo26",
        2026,
        "Repo endpoint-pressure program",
        "internal proof route",
        "endpoint incidence, pressure DAGs, Tournament Analysis",
        "05-knowledge/hypotheses/HYP-1950-lrc-n14-n18-pressure-peelability.md",
        0,
        2,
        1,
        3,
        2,
        4,
        5,
        2,
        3,
        2,
        4,
        10,
        9,
        "turn external sieves into endpoint row-cover and pressure-peel certificates",
    ),
]


def derived(route: Route, key: str) -> int:
    if key == "scope_strength":
        return 2 * route.generality + min(route.total_runners, 13)
    if key == "boundary_overlap":
        return route.sieve + route.endpoint + route.zonotope + route.spectral + route.mixed
    if key == "proof_engine":
        return route.finite_check + route.sieve + route.spectral + route.zonotope
    if key == "frontier_signal":
        return route.frontier + (1 if route.year >= 2025 else 0)
    if key == "repo_signal":
        return route.repo_transfer + route.endpoint + route.sieve
    return getattr(route, key)


def compare_values(a: Route, b: Route, key: str) -> int:
    av = derived(a, key)
    bv = derived(b, key)
    return (av > bv) - (av < bv)


def majority_orientation(
    routes: list[Route],
    channels: list[str],
    tie_rank: dict[str, int],
) -> list[list[bool]]:
    n = len(routes)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a = routes[i]
        b = routes[j]
        wins = 0
        losses = 0
        for channel in channels:
            cmp_value = compare_values(a, b, channel)
            if cmp_value > 0:
                wins += 1
            elif cmp_value < 0:
                losses += 1
        if wins > losses:
            adj[i][j] = True
        elif losses > wins:
            adj[j][i] = True
        elif tie_rank[a.code] < tie_rank[b.code]:
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def score_sequence(adj: list[list[bool]]) -> list[int]:
    return [sum(row) for row in adj]


def directed_triangle_count(adj: list[list[bool]]) -> int:
    count = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            count += 1
    return count


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    index = 0
    stack: list[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    lows = [0] * n
    sizes: list[int] = []

    def strongconnect(v: int) -> None:
        nonlocal index
        indices[v] = index
        lows[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w, edge in enumerate(adj[v]):
            if not edge:
                continue
            if indices[w] == -1:
                strongconnect(w)
                lows[v] = min(lows[v], lows[w])
            elif on_stack[w]:
                lows[v] = min(lows[v], indices[w])
        if lows[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack[w] = False
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for v in range(n):
        if indices[v] == -1:
            strongconnect(v)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: list[list[bool]]) -> int:
    n = len(adj)

    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if prev_mask & (1 << prev) and adj[prev][last]:
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def topological_layers_if_dag(adj: list[list[bool]]) -> list[list[int]] | None:
    n = len(adj)
    indeg = [0] * n
    outgoing = [set() for _ in range(n)]
    for i in range(n):
        for j, edge in enumerate(adj[i]):
            if edge:
                indeg[j] += 1
                outgoing[i].add(j)
    remaining = set(range(n))
    layers: list[list[int]] = []
    while remaining:
        layer = sorted(v for v in remaining if indeg[v] == 0)
        if not layer:
            return None
        layers.append(layer)
        for v in layer:
            remaining.remove(v)
            for w in outgoing[v]:
                indeg[w] -= 1
    return layers


def fingerprint(name: str, routes: list[Route], adj: list[list[bool]]) -> None:
    scores = score_sequence(adj)
    hist = Counter(scores)
    ordered = sorted(zip(routes, scores), key=lambda item: (-item[1], item[0].year, item[0].code))
    print(f"\n{name} tournament")
    print("  score histogram:", dict(sorted(hist.items())))
    print("  SCC sizes:", scc_sizes(adj))
    print("  directed 3-cycles:", directed_triangle_count(adj))
    print("  Hamiltonian-path count:", hamiltonian_path_count(adj))
    layers = topological_layers_if_dag(adj)
    if layers is not None:
        print("  DAG layers:", [[routes[i].code for i in layer] for layer in layers])
    print("  top routes:")
    for route, score in ordered[:6]:
        print(f"    {route.code:9s} score={score:2d}  {route.title}")


def edge_flips(
    routes: list[Route], first: list[list[bool]], second: list[list[bool]]
) -> list[tuple[str, str]]:
    flips = []
    for i, j in combinations(range(len(routes)), 2):
        if first[i][j] != second[i][j]:
            winner_first = routes[i].code if first[i][j] else routes[j].code
            winner_second = routes[i].code if second[i][j] else routes[j].code
            flips.append((winner_first, winner_second))
    return flips


def print_route_table() -> None:
    print("LRC proof-route atlas (codex-2026-06-01-S505)")
    print("Convention: total_runners counts the stationary/reference runner when applicable.")
    print("\nRoutes:")
    for route in ROUTES:
        print(
            f"  {route.code:9s} {route.year} trust={route.trust} "
            f"direct={route.direct_lrc} repo={route.repo_transfer} "
            f"runners={route.total_runners or '-':>2}  {route.title}"
        )
        print(f"             {route.scope}; {route.note}")


def main() -> None:
    print_route_table()
    tie_order = sorted(ROUTES, key=lambda r: (r.year, r.code))
    tie_rank = {route.code: i for i, route in enumerate(tie_order)}
    print("\nTournament Analysis declaration")
    print("  Pairwise observable: feature-wise comparison between two proof routes.")
    print(
        "  Actionability gauge: majority over "
        "repo_signal, boundary_overlap, proof_engine, reproducible, frontier_signal."
    )
    print(
        "  Credence gauge: majority over "
        "trust, direct_lrc, scope_strength, reproducible, generality."
    )
    print(
        "  Tie Hamiltonian path:",
        " -> ".join(route.code for route in tie_order),
    )

    action_channels = [
        "repo_signal",
        "boundary_overlap",
        "proof_engine",
        "reproducible",
        "frontier_signal",
    ]
    credence_channels = [
        "trust",
        "direct_lrc",
        "scope_strength",
        "reproducible",
        "generality",
    ]
    action_adj = majority_orientation(ROUTES, action_channels, tie_rank)
    credence_adj = majority_orientation(ROUTES, credence_channels, tie_rank)

    fingerprint("Actionability", ROUTES, action_adj)
    fingerprint("Credence", ROUTES, credence_adj)

    flips = edge_flips(ROUTES, credence_adj, action_adj)
    print("\nGauge edge flips (credence winner -> actionability winner):")
    for before, after in flips[:24]:
        print(f"  {before:9s} -> {after:9s}")
    print(f"  total flips: {len(flips)} of {len(ROUTES) * (len(ROUTES) - 1) // 2}")

    buckets: dict[str, list[str]] = defaultdict(list)
    for route in ROUTES:
        if route.sieve >= 4:
            buckets["finite product-sieve proofs"].append(route.code)
        if route.spectral >= 4:
            buckets["Fourier/Riesz kernel tools"].append(route.code)
        if route.zonotope >= 4:
            buckets["zonotope/covering-radius tools"].append(route.code)
        if route.mixed >= 4:
            buckets["mixed-threshold tools"].append(route.code)
        if route.endpoint >= 4:
            buckets["endpoint-incidence compatible"].append(route.code)

    print("\nProof-currency buckets:")
    for bucket, codes in sorted(buckets.items()):
        print(f"  {bucket}: {', '.join(codes)}")

    print("\nSynthesis")
    print(
        "  The current accepted-looking fixed-case frontier is product-sieve plus "
        "finite-checking: R8 -> T9T10/R9 -> ST13."
    )
    print(
        "  The best general compression is MSS25; it is the bridge that turns "
        "divisibility-forcing into a contradiction."
    )
    print(
        "  Bedert25 and Jensen26 do not prove LRC, but they are the best matches "
        "for the repo's endpoint-kernel pressure hypothesis."
    )
    print(
        "  Shifted-LRC results are useful geometry but BCS26 warns that shifted "
        "analogies can be false even when the original LRC remains open."
    )
    print(
        "  Repo26 should next translate J(k,p)-emptiness into endpoint row-cover "
        "language and test whether pressure DAG leaves mirror improper-tuple sieves."
    )


if __name__ == "__main__":
    main()
