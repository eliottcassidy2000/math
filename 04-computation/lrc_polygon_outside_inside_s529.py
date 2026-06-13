#!/usr/bin/env python3
"""Outside/inside regular-polygon lens for LRC tournament geometry.

codex-2026-06-01 S529

The simplex view of a tournament has one binary relation per simplex edge.
The regular-polygon view starts with a cyclic outside necklace of gaps.  Every
inside chord orientation is then the half-turn threshold of a consecutive gap
sum.  The dihedral group acts on the outside necklace, while the hidden inside
arcs split into chord-length channels.

For LRC, the tight initial row at t=1/n is the regular n-gon with the observer
as a missing clasp vertex.  The source gaps adjacent to the clasp are exactly
1/n, so this is a compactified wall.  For even n, the clasp-deleted polygon
also contains hidden diameter chords among runners; their half-turn comparator
is a tie, so the boundary/tie path is not cosmetic.

Tournament Analysis declaration:
    vertices: clasp-deleted regular n-gons for n=4..18
    pairwise observable: hidden boundary burden profile
        (diameter tie count, channel count, hidden chord count, outside ratio)
    switch/gauge: lexicographic comparison of that profile
    tie Hamiltonian path: increasing n

Stored output:
    05-knowledge/results/lrc_polygon_outside_inside_s529.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def full_channel_count(n: int, d: int) -> int:
    """Number of unoriented chords of cyclic distance d in a full n-gon."""
    if n % 2 == 0 and d == n // 2:
        return n // 2
    return n


def clasp_deleted_channel_count(n: int, d: int) -> int:
    """Channel count after deleting vertex 0 as the LRC observer/clasp."""
    if n % 2 == 0 and d == n // 2:
        return n // 2 - 1
    return n - 2


def clasp_channel_inventory(n: int) -> Counter[int]:
    return Counter(
        {
            d: clasp_deleted_channel_count(n, d)
            for d in range(1, n // 2 + 1)
        }
    )


def hidden_diameter_count(n: int) -> int:
    if n % 2:
        return 0
    return n // 2 - 1


def occupied_vertices(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def half_turn_tournament_on_clasp_deleted_ngon(n: int) -> list[list[int]]:
    """Tie-resolved half-turn comparator on vertices 1..n-1 of the n-gon."""
    verts = occupied_vertices(n)
    idx = {v: i for i, v in enumerate(verts)}
    m = len(verts)
    adj = [[0] * m for _ in range(m)]
    for a, b in combinations(verts, 2):
        step = (b - a) % n
        if 2 * step < n:
            winner, loser = a, b
        elif 2 * step > n:
            winner, loser = b, a
        else:
            # Diameter wall: the outside polygon gives a hidden tie; use the
            # increasing Hamiltonian path as the compactification.
            winner, loser = min(a, b), max(a, b)
        adj[idx[winner]][idx[loser]] = 1
    return adj


def score_hist(adj: list[list[int]]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(sum(row) for row in adj).items()))


def count_c3(adj: list[list[int]]) -> int:
    n = len(adj)
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        edges = adj[a][b] + adj[b][c] + adj[c][a]
        if edges in (0, 3):
            c3 += 1
    return c3


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    amat = tuple(tuple(row) for row in adj)

    @lru_cache(None)
    def dp(mask: int, last: int) -> int:
        if mask == full:
            return 1
        total = 0
        for nxt in range(n):
            if mask & (1 << nxt):
                continue
            if amat[last][nxt]:
                total += dp(mask | (1 << nxt), nxt)
        return total

    return sum(dp(1 << start, start) for start in range(n))


@dataclass(frozen=True)
class PolygonProfile:
    n: int
    runner_vertices: int
    hidden_edges: int
    outside_gaps: int
    channels: int
    channel_counts: tuple[tuple[int, int], ...]
    diameter_ties: int
    outside_inside_ratio: Fraction
    lrc_clasp_gaps: tuple[Fraction, Fraction]
    score_hist: tuple[tuple[int, int], ...]
    c3: int
    hamiltonian_paths: int | None

    def burden(self) -> tuple:
        return (
            self.diameter_ties,
            self.channels,
            self.hidden_edges,
            self.outside_inside_ratio,
        )


def profile(n: int) -> PolygonProfile:
    m = n - 1
    inv = clasp_channel_inventory(n)
    adj = half_turn_tournament_on_clasp_deleted_ngon(n)
    h_paths = hamiltonian_paths(adj) if m <= 13 else None
    return PolygonProfile(
        n=n,
        runner_vertices=m,
        hidden_edges=m * (m - 1) // 2,
        outside_gaps=n,
        channels=len(inv),
        channel_counts=tuple(sorted(inv.items())),
        diameter_ties=hidden_diameter_count(n),
        outside_inside_ratio=Fraction(m * (m - 1), 2 * n),
        lrc_clasp_gaps=(Fraction(1, n), Fraction(1, n)),
        score_hist=score_hist(adj),
        c3=count_c3(adj),
        hamiltonian_paths=h_paths,
    )


def scc_sizes(adj: list[list[int]]) -> tuple[int, ...]:
    n = len(adj)
    graph = {i: [j for j in range(n) if adj[i][j]] for i in range(n)}
    rev = {i: [] for i in range(n)}
    for i, outs in graph.items():
        for j in outs:
            rev[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: int) -> int:
        seen.add(v)
        total = 1
        for w in rev[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def route_tournament(profiles: list[PolygonProfile]) -> list[list[int]]:
    n = len(profiles)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        left = profiles[i].burden()
        right = profiles[j].burden()
        if left > right:
            adj[i][j] = 1
        elif right > left:
            adj[j][i] = 1
        else:
            adj[i][j] = 1
    return adj


def print_profile(p: PolygonProfile) -> None:
    print(f"n={p.n} (observer + {p.runner_vertices} runners)")
    print(
        "  outside gaps="
        f"{p.outside_gaps}; hidden runner chords={p.hidden_edges}; "
        f"inside/outside={fmt_frac(p.outside_inside_ratio)}"
    )
    print(
        "  D_n chord channels after clasp deletion: "
        f"{p.channel_counts}; hidden diameter ties={p.diameter_ties}"
    )
    print(
        "  LRC AP witness clasp gaps="
        f"({fmt_frac(p.lrc_clasp_gaps[0])}, {fmt_frac(p.lrc_clasp_gaps[1])}) "
        "(closed source wall, not open)"
    )
    h = "-" if p.hamiltonian_paths is None else str(p.hamiltonian_paths)
    print(f"  tie-resolved inside tournament: score_hist={p.score_hist} c3={p.c3} H={h}")


def main() -> None:
    print("Polygon outside/inside arcs lens for LRC -- codex-2026-06-01 S529")
    print()
    print("MODEL")
    print("=" * 72)
    print("Simplex view: a tournament has C(m,2) independent-looking binary edges.")
    print("Polygon view: a cyclic outside gap necklace determines every hidden")
    print("inside chord by consecutive gap sums and the half-turn threshold.")
    print("D_n acts on the outside necklace; inside chords split into distance")
    print("channels.  LRC marks one outside vertex as the observer/clasp.")
    print()

    print("PART A: clasp-deleted regular n-gons")
    print("=" * 72)
    for n in range(4, 19):
        p = profile(n)
        diameter = f"{p.diameter_ties:2d}" if p.diameter_ties else " -"
        print(
            f"n={n:2d} m={p.runner_vertices:2d} hidden={p.hidden_edges:3d} "
            f"outside={p.outside_gaps:2d} ratio={fmt_frac(p.outside_inside_ratio):>5} "
            f"channels={p.channels:2d} diameter_ties={diameter}"
        )
    print()

    print("PART B: focused frontiers")
    print("=" * 72)
    for n in (6, 8, 14, 18):
        print_profile(profile(n))
        print()

    print("PART C: outside/inside lesson for n=14")
    print("=" * 72)
    p14 = profile(14)
    print("At the initial AP witness speeds {1..13}, t=1/14:")
    print("  - outside: the 14-gon is regular and the observer is the missing clasp;")
    print("  - the two clasp-adjacent gaps are exactly 1/14, so loneliness is a")
    print("    compactified boundary source, not an open source cell;")
    print(
        "  - inside: 78 runner-runner chords are forced from only 14 outside gaps,"
    )
    print(
        "    split into D_14 channels "
        f"{p14.channel_counts}, with {p14.diameter_ties} hidden diameter ties;"
    )
    print("  - those diameter ties are precisely boundary data: a tie Hamiltonian path")
    print("    is needed before the inside tournament is a genuine simplex relation.")
    print()
    print("So the hidden inside arcs are not extra free variables.  They are")
    print("shadow data of the outside gap necklace, except on the wall where")
    print("diameter/endpoint ties create compactified choices.  This is why the")
    print("LRC proof target should be a labelled boundary/endpoint-debt statement,")
    print("not a free simplex-edge argument.")
    print()

    print("PART D: Tournament Analysis over n=4..18 polygon profiles")
    print("=" * 72)
    profiles = [profile(n) for n in range(4, 19)]
    adj = route_tournament(profiles)
    print(f"vertices={[p.n for p in profiles]}")
    print(f"score_hist={score_hist(adj)}")
    print(f"c3={count_c3(adj)}")
    print(f"SCCs={scc_sizes(adj)}")
    print(f"Hamiltonian_paths={hamiltonian_paths(adj)}")
    ranking = sorted(profiles, key=lambda p: p.burden(), reverse=True)
    print("ranking by hidden boundary burden:")
    for p in ranking:
        print(
            f"  n={p.n}: burden=(diam={p.diameter_ties}, channels={p.channels}, "
            f"hidden={p.hidden_edges}, ratio={fmt_frac(p.outside_inside_ratio)})"
        )
    print()

    print("SYNTHESIS")
    print("=" * 72)
    print("1. A tournament is a simplex relation, but circular/LRC-realizable")
    print("   tournaments are not a simplex cube: they are images of an outside")
    print("   polygon gap necklace under interval-sum threshold maps.")
    print("2. The dihedral group sees the outside polygon; its hidden inside arcs")
    print("   are chord-length channel orbits.  These channels are the right")
    print("   coordinates for regular-polygon extremals.")
    print("3. LRC marks a clasp vertex.  Source means the clasp has an empty")
    print("   forbidden arc.  The AP/regular polygon has that empty arc only on")
    print("   the closed wall; for even n it also carries hidden diameter ties.")
    print("4. For n=14, the regular-polygon wall has 6 hidden diameter ties and")
    print("   78 inside chords determined by 14 outside gaps.  The next proof")
    print("   target is to show non-AP sweeps cannot keep all outside clasp gaps")
    print("   closed while cycling the hidden boundary debt forever.")


if __name__ == "__main__":
    main()
