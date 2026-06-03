#!/usr/bin/env python3
"""
unit_distance_n22_tournament_lrc_s614.py

S614 scout for the planar unit distance problem at n=22.

This is not a full SAT/SMT enumeration of all 61-edge candidates.  It records
the current published state, extracts graph-theoretic consequences from the
new exact n<=21 frontier, runs small lattice scouts, and applies the repo's
Tournament Analysis default to the route space.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from math import comb


U_EXACT = {
    0: 0,
    1: 0,
    2: 1,
    3: 3,
    4: 5,
    5: 7,
    6: 9,
    7: 12,
    8: 14,
    9: 18,
    10: 20,
    11: 23,
    12: 27,
    13: 30,
    14: 33,
    15: 37,
    16: 41,
    17: 43,
    18: 46,
    19: 50,
    20: 54,
    21: 57,
}

U22_LOWER = 60
U22_UPPER = 61

N21_GRAPH6 = [
    "Tsc@IGC@GD?R?S?Wd@A_CK@HG@VM??PRKOUZ",
    "TCKx?D?OI?OMCBSA_L?ApA_gEA\\EG?PBSCPV",
    "TCTWACAG@@CDKC?e?QgQA@OMOq]F??OUcCEj",
    "TCS`?H??XIZ?K_Co`CG@JO[?EOSCpGOTSCE\\",
    "T??_`OhSCSYA@I?c?OyWBEa@c?SIU?Aa[?el",
]


def graph6_decode(code: str) -> list[list[int]]:
    """Decode graph6 strings for n <= 62."""
    data = [ord(ch) - 63 for ch in code.strip()]
    if not data:
        raise ValueError("empty graph6 code")
    n = data[0]
    if n >= 63:
        raise ValueError("this scout only handles short graph6 codes")
    bits: list[int] = []
    for x in data[1:]:
        for shift in range(5, -1, -1):
            bits.append((x >> shift) & 1)
    adj = [[0] * n for _ in range(n)]
    p = 0
    for j in range(1, n):
        for i in range(j):
            if bits[p]:
                adj[i][j] = adj[j][i] = 1
            p += 1
    return adj


def edge_count(adj: list[list[int]]) -> int:
    return sum(adj[i][j] for i in range(len(adj)) for j in range(i + 1, len(adj)))


def degrees(adj: list[list[int]]) -> list[int]:
    return [sum(row) for row in adj]


def degree_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for d in degrees(adj):
        hist[d] = hist.get(d, 0) + 1
    return dict(sorted(hist.items()))


def count_triangles(adj: list[list[int]]) -> int:
    n = len(adj)
    return sum(
        adj[i][j] and adj[i][k] and adj[j][k]
        for i, j, k in combinations(range(n), 3)
    )


def count_4cycles(adj: list[list[int]]) -> int:
    n = len(adj)
    total = 0
    for a, b, c, d in combinations(range(n), 4):
        verts = [a, b, c, d]
        e = 0
        for i, j in combinations(verts, 2):
            e += adj[i][j]
        # A K4 has three 4-cycles, a 4-cycle with one diagonal has one.
        if e == 4:
            total += 1
        elif e == 5:
            total += 2
        elif e == 6:
            total += 3
    return total


TRI_DIRS = ((1, 0), (0, 1), (1, -1), (-1, 0), (0, -1), (-1, 1))


def add2(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (a[0] + b[0], a[1] + b[1])


def sub2(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (a[0] - b[0], a[1] - b[1])


def canon2(cluster: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    pts = sorted(cluster)
    base = pts[0]
    return tuple(sorted(sub2(p, base) for p in pts))


def tri_edges(cluster: tuple[tuple[int, int], ...]) -> int:
    s = set(cluster)
    return sum(1 for p in s for u in TRI_DIRS if add2(p, u) in s) // 2


def tri_frontier(cluster: tuple[tuple[int, int], ...]) -> set[tuple[int, int]]:
    s = set(cluster)
    out: set[tuple[int, int]] = set()
    for p in s:
        for u in TRI_DIRS:
            q = add2(p, u)
            if q not in s:
                out.add(q)
    return out


def compact_span2(cluster: tuple[tuple[int, int], ...]) -> int:
    qs = [p[0] for p in cluster]
    rs = [p[1] for p in cluster]
    ss = [-p[0] - p[1] for p in cluster]
    return (max(qs) - min(qs)) + (max(rs) - min(rs)) + (max(ss) - min(ss))


def triangular_beam(k: int = 22, width: int = 1500) -> tuple[int, tuple[tuple[int, int], ...]]:
    beam = {((0, 0),)}
    best: tuple[int, tuple[tuple[int, int], ...]] = (0, ((0, 0),))
    for size in range(2, k + 1):
        seen: dict[tuple[tuple[int, int], ...], int] = {}
        for cluster in beam:
            for q in tri_frontier(cluster):
                nxt = canon2(cluster + (q,))
                seen[nxt] = tri_edges(nxt)
        ranked = sorted(
            ((e, -compact_span2(c), c) for c, e in seen.items()), reverse=True
        )
        beam = {c for _, _, c in ranked[:width]}
        best = (ranked[0][0], ranked[0][2])
        if size in {7, 13, 19, 22}:
            print(
                f"triangular beam size={size}: states={len(seen)} "
                f"best_edges={ranked[0][0]} compact_span={-ranked[0][1]}"
            )
    return best


def moser_unit_vectors() -> tuple[tuple[int, int, int, int], ...]:
    units: list[tuple[int, int, int, int]] = []
    for a, b, c, d in product(range(-4, 5), repeat=4):
        if a == b == c == d == 0:
            continue
        if a * d != b * c:
            continue
        value = (
            6 * a * a
            + 6 * a * b
            + 10 * a * c
            + 5 * a * d
            + 6 * b * b
            + 5 * b * c
            + 10 * b * d
            + 6 * c * c
            + 6 * c * d
            + 6 * d * d
        )
        if value == 6:
            units.append((a, b, c, d))
    return tuple(sorted(units))


@dataclass(frozen=True)
class Route:
    name: str
    burden_bits: float
    side_info: int
    relation_rank: int
    n22_specificity: int
    coimage_loss: int
    geometry_contact: int


ROUTES = (
    Route("enumerate 61-edge F-free graphs", 14.0, 2, 2, 5, 1, 2),
    Route("extend 21-vertex max/near-max cores", 7.0, 4, 2, 5, 0, 3),
    Route("Moser-ring beam search to 61", 5.0, 3, 4, 4, 1, 4),
    Route("mine totally-unfaithful obstructions", 6.0, 5, 3, 4, 0, 2),
    Route("graph-only 62 upper-bound transfer", 3.0, 1, 1, 3, 5, 1),
    Route("triangular/square grid variants", 2.0, 0, 1, 1, 6, 4),
    Route("CM tower asymptotic transfer", 8.0, 5, 5, 1, 3, 5),
    Route("LRC coimage side-channel transfer", 6.0, 5, 5, 2, 1, 0),
)


def route_votes(a: Route, b: Route) -> tuple[int, int]:
    criteria = [
        (a.burden_bits < b.burden_bits, b.burden_bits < a.burden_bits),
        (a.side_info > b.side_info, b.side_info > a.side_info),
        (a.relation_rank > b.relation_rank, b.relation_rank > a.relation_rank),
        (a.n22_specificity > b.n22_specificity, b.n22_specificity > a.n22_specificity),
        (a.coimage_loss < b.coimage_loss, b.coimage_loss < a.coimage_loss),
        (a.geometry_contact > b.geometry_contact, b.geometry_contact > a.geometry_contact),
    ]
    av = sum(1 for x, y in criteria if x and not y)
    bv = sum(1 for x, y in criteria if y and not x)
    return av, bv


def route_tournament(routes: tuple[Route, ...]) -> list[list[int]]:
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            iv, jv = route_votes(routes[i], routes[j])
            if iv > jv or (iv == jv and i < j):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj


def tournament_score_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for row in adj:
        d = sum(row)
        hist[d] = hist.get(d, 0) + 1
    return dict(sorted(hist.items()))


def tournament_c3(adj: list[list[int]]) -> int:
    n = len(adj)
    total = 0
    for i, j, k in combinations(range(n), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            total += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(starts: list[int], forward: bool) -> set[int]:
        seen = set(starts)
        stack = list(starts)
        while stack:
            v = stack.pop()
            for w in range(n):
                edge = adj[v][w] if forward else adj[w][v]
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    left = set(range(n))
    sizes = []
    while left:
        v = next(iter(left))
        comp = reach([v], True) & reach([v], False)
        sizes.append(len(comp))
        left -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def print_published_state() -> None:
    print("Published state for planar unit distances at n=22")
    print(f"  Exact values through n=21 end at u(21)={U_EXACT[21]}.")
    print(f"  Current n=22 bounds: {U22_LOWER} <= u(22) <= {U22_UPPER}.")
    print("  The 62-edge F-free graph-only frontier has exactly two candidates;")
    print("  both are killed by totally unfaithful subgraphs, giving u(22)<=61.")
    print("  The 60-edge lower bound is realized by known Moser-ring database examples.")
    print("  Therefore the n=22 optimum is exactly one bit: either 60 or 61.")
    print()


def print_n21_core_data() -> None:
    print("n=21 exact-core deletion pressure")
    decoded = [graph6_decode(code) for code in N21_GRAPH6]
    for idx, adj in enumerate(decoded, 1):
        ds = degrees(adj)
        print(
            f"  core {idx}: edges={edge_count(adj)} degree_hist={degree_hist(adj)} "
            f"min={min(ds)} max={max(ds)} triangles={count_triangles(adj)} "
            f"fourcycle_count={count_4cycles(adj)}"
        )
    print("  Consequence for a hypothetical 61-edge UDG on 22 vertices:")
    print("    average degree = 122/22, so a min-degree vertex has degree <=5;")
    print("    deleting it leaves at least 56 edges on 21 vertices;")
    print("    u(21)=57 forces min degree >=4, so min degree is 4 or 5.")
    print("    Thus a 61 proof/hunt can focus on 21-vertex 56/57-edge cores")
    print("    plus one degree-5 or degree-4 embeddable extension vertex.")
    print("  Consequence for a 60-edge witness:")
    print("    min degree is at least 3, and deleting a min vertex leaves 55..57")
    print("    edges; known 60-edge examples should have dense 21-decks.")
    print()


def print_lattice_scouts() -> None:
    print("Small lattice scouts")
    best_tri, tri_cluster = triangular_beam()
    print(
        f"  Triangular-lattice connected beam best at n=22: {best_tri} edges "
        f"(cluster size {len(tri_cluster)})."
    )
    print(
        "  This is intentionally a grid sanity check, not a proof of the "
        "global UDG optimum; it shows why the 60-edge lower bound must use "
        "richer unit-vector structure than the plain triangular lattice."
    )
    units = moser_unit_vectors()
    print(f"  Moser lattice unit shell: {len(units)} directed unit vectors, {len(units)//2} antipodal pairs.")
    print(f"  Unit vectors in Z^4 basis: {units}")
    print(
        "  The Moser model is already a coimage: rank-4 integer coordinates "
        "plus an 18-vector unit shell project to planar geometry."
    )
    print()


def print_route_tournament() -> None:
    print("Tournament Analysis on n=22 proof routes")
    adj = route_tournament(ROUTES)
    scores = [sum(row) for row in adj]
    for route, score in sorted(zip(ROUTES, scores), key=lambda x: (-x[1], x[0].name)):
        print(
            f"  score={score} route={route.name}; "
            f"features=(burden={route.burden_bits}, side={route.side_info}, "
            f"rank={route.relation_rank}, spec={route.n22_specificity}, "
            f"loss={route.coimage_loss}, geom={route.geometry_contact})"
        )
    print(f"  score histogram: {tournament_score_hist(adj)}")
    print(f"  directed 3-cycles: {tournament_c3(adj)}")
    print(f"  SCC sizes: {scc_sizes(adj)}")
    print(f"  Hamiltonian path count: {hamiltonian_paths(adj)}")
    print(
        "  Reading: the best immediate route is not a raw graph-only quotient. "
        "It is an extension/obstruction search that retains side information "
        "about dense cores, Moser coordinates, and unfaithful subgraphs."
    )
    print()


def print_lrc_bridge() -> None:
    print("LRC/tournament bridge")
    print("  Shared pattern: exactness fails when the quotient forgets side channels.")
    print("  Unit-distance n=22: F-free graph coimage says 62, but geometry/unfaithfulness")
    print("    lowers it to <=61; deciding 61 needs retained embedding data.")
    print("  LRC n=14: Res_27/coarse tournament coimages shrink the state, but floor")
    print("    status needs owner labels, carry cocycles, and endpoint/pinch channels.")
    print("  Recent asymptotic disproof: the old grid intuition fixes a low-degree")
    print("    visible lattice; the counterexample fixes small primes and varies the")
    print("    CM field degree through a Golod-Shafarevich class-field tower, then")
    print("    projects one coordinate to the plane.")
    print("  Hypothesis HYP-2176: n=22 is the first small unit-distance coimage")
    print("    failure frontier, just as n=14 LRC is a small coimage-with-side-channel")
    print("    frontier. The right tournament vertices are proof obligations/cores,")
    print("    not only points or runners.")
    print()


def main() -> None:
    print("==== S614 unit-distance n=22 / tournament / LRC scout ====")
    print_published_state()
    print_n21_core_data()
    print_lattice_scouts()
    print_route_tournament()
    print_lrc_bridge()
    print("Sources summarized in the companion reflection:")
    print("  Alexeev-Tikhonov arXiv:2412.11914 exact n<=21 and n=22 bounds.")
    print("  Engel-Hammond-Lee-Su-Varga-Zsamboki arXiv:2406.15317 Moser-ring database.")
    print("  Alon-Bloom-Gowers-Litt-Sawin-Shankar-Tsimerman-Wang-Wood arXiv:2605.20695.")
    print("  Sawin arXiv:2605.20579 for explicit n^1.014 asymptotic lower bound.")


if __name__ == "__main__":
    main()
