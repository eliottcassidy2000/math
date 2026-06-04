#!/usr/bin/env python3
"""
unit_distance_hgap_carrier_s627.py

S627: Reframe the unit-distance tournament question through H-gap carriers.

The key distinction is that unit-distance edge counts are additive geometry
observables, while tournament H is a Hamiltonian-path/OCF evaluation with
strong-component multiplicative constraints.  So a unit-distance row can have
7 unit edges without producing a tournament with H=7.  The useful bridge is the
side-channel carrier split: spine, tile bulk, frontier shell, exact-vs-lattice
carrier, and H-gap guardrail.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import floor, sqrt


FORBIDDEN_H = {
    7: "tournament H=7 permanent gap",
    21: "tournament H=21 permanent gap",
}

# Exact planar values used by the S625 scout through n=14.
EXACT_U_SMALL = {
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
}

# The S599z unit-spine tournament H values.  These are not edge counts; they
# are directed Hamiltonian-path counts for the unit-flip tournament built from
# a unit Hamiltonian base path.
S599Z_UNIT_TOURNAMENT_H = {
    3: 1,
    4: 5,
    5: 15,
    6: 43,
    7: 141,
    8: 513,
    9: 1605,
    10: 4915,
}


def harborth_lattice_edges(n: int) -> int:
    return floor(3 * n - sqrt(12 * n - 3))


def echo_label(value: int) -> str:
    return FORBIDDEN_H.get(value, "")


def fmt(value: object) -> str:
    return str(value) if value is not None else "-"


@dataclass(frozen=True)
class EchoRow:
    n: int
    exact_u: int
    harborth_u: int
    spine: int
    exact_tiles: int
    harborth_tiles: int
    literal_h: int | None

    @property
    def exact_echo(self) -> str:
        return echo_label(self.exact_u)

    @property
    def lattice_echo(self) -> str:
        return echo_label(self.harborth_u)


def echo_rows() -> list[EchoRow]:
    rows = []
    for n in range(3, 15):
        exact = EXACT_U_SMALL[n]
        lattice = harborth_lattice_edges(n)
        spine = n - 1
        rows.append(
            EchoRow(
                n=n,
                exact_u=exact,
                harborth_u=lattice,
                spine=spine,
                exact_tiles=exact - spine,
                harborth_tiles=lattice - spine,
                literal_h=S599Z_UNIT_TOURNAMENT_H.get(n),
            )
        )
    return rows


@dataclass(frozen=True)
class Carrier:
    name: str
    traceability: int
    decomposition: int
    hgap: int
    scaling: int
    lrc_bridge: int
    nonlattice: int
    recursion: int
    risk: int
    note: str


CARRIERS = (
    Carrier(
        "unit-spine Hamiltonian path",
        5,
        4,
        3,
        3,
        2,
        2,
        5,
        1,
        "keeps the path predicate: the unit graph is traceable",
    ),
    Carrier(
        "tile/bulk unit-edge packet",
        3,
        5,
        4,
        4,
        2,
        3,
        4,
        1,
        "turns u(n) into spine plus non-spine unit tiles",
    ),
    Carrier(
        "boundary frontier shell",
        4,
        4,
        3,
        5,
        4,
        3,
        5,
        2,
        "records the +3 gain and perimeter deficit recursion",
    ),
    Carrier(
        "Moser/non-lattice exact carrier",
        3,
        3,
        4,
        5,
        2,
        5,
        4,
        2,
        "keeps the exact planar side beyond the triangular lattice",
    ),
    Carrier(
        "OCF H-gap guardrail",
        2,
        5,
        5,
        3,
        5,
        2,
        3,
        1,
        "blocks collapsing edge counts into forbidden tournament H values",
    ),
    Carrier(
        "round-LRC worry-set channel",
        2,
        4,
        5,
        3,
        5,
        1,
        4,
        1,
        "imports the known round-tournament H-gap propagation channel",
    ),
    Carrier(
        "literal unit-tournament H",
        4,
        3,
        4,
        3,
        3,
        2,
        2,
        1,
        "counts directed Hamiltonian paths; odd and avoids 7,21 in S599z",
    ),
    Carrier(
        "raw edge-count quotient",
        1,
        1,
        1,
        4,
        1,
        4,
        3,
        4,
        "forgets whether an edge count is spine, tile, or proof side-channel",
    ),
    Carrier(
        "equidecomposability ledger",
        3,
        5,
        5,
        4,
        4,
        4,
        4,
        2,
        "records the scissors split behind same-cardinality counts",
    ),
)


def geometry_gauge(c: Carrier) -> tuple[int, ...]:
    return (
        c.traceability,
        c.recursion,
        c.decomposition,
        c.scaling,
        c.hgap,
        -c.risk,
        c.name,
    )


def hgap_gauge(c: Carrier) -> tuple[int, ...]:
    return (
        c.hgap,
        c.decomposition,
        c.lrc_bridge,
        c.traceability,
        -c.risk,
        c.name,
    )


def scaling_gauge(c: Carrier) -> tuple[int, ...]:
    return (
        c.scaling,
        c.nonlattice,
        c.recursion,
        c.decomposition,
        c.traceability,
        -c.risk,
        c.name,
    )


GAUGES = {
    "geometry": geometry_gauge,
    "H-gap": hgap_gauge,
    "scaling": scaling_gauge,
}


def compare_by_gauge(i: int, j: int, gauge_name: str) -> int:
    gauge = GAUGES[gauge_name]
    a = gauge(CARRIERS[i])
    b = gauge(CARRIERS[j])
    if a > b:
        return i
    if b > a:
        return j
    return i if i < j else j


def tournament_from_gauge(gauge_name: str) -> list[list[int]]:
    n = len(CARRIERS)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        winner = compare_by_gauge(i, j, gauge_name)
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def majority_tournament() -> list[list[int]]:
    n = len(CARRIERS)
    adj = [[0] * n for _ in range(n)]
    names = tuple(GAUGES)
    for i, j in combinations(range(n), 2):
        wins_i = sum(1 for name in names if compare_by_gauge(i, j, name) == i)
        winner = i if wins_i >= 2 else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def score_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for row in adj:
        score = sum(row)
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def sccs(adj: list[list[int]]) -> list[list[int]]:
    n = len(adj)

    def reach(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w, edge in enumerate(graph[v]):
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    out: list[list[int]] = []
    while remaining:
        v = min(remaining)
        comp = reach(v, adj) & reach(v, rev)
        out.append(sorted(comp))
        remaining -= comp
    return out


def hamiltonian_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp[mask][last]
            if not ways:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += ways
    return sum(dp[-1])


def insertion_hp(adj: list[list[int]]) -> list[int]:
    path: list[int] = []
    for v in range(len(adj)):
        pos = 0
        while pos < len(path) and adj[path[pos]][v]:
            pos += 1
        path.insert(pos, v)
    return path


def edge_flips(a: list[list[int]], b: list[list[int]]) -> int:
    flips = 0
    for i, j in combinations(range(len(a)), 2):
        if a[i][j] != b[i][j]:
            flips += 1
    return flips


def print_echo_table() -> None:
    print("EDGE-COUNT ECHO TABLE")
    print("---------------------")
    print(
        "literal_H is the S599z unit-spine tournament count; it is a carrier "
        "diagnostic, not the exact planar edge count."
    )
    print("n exact_u harborth_u spine exact_tiles lattice_tiles literal_H echoes")
    for row in echo_rows():
        echoes = []
        if row.exact_echo:
            echoes.append(f"exact:{row.exact_echo}")
        if row.lattice_echo and row.lattice_echo != row.exact_echo:
            echoes.append(f"lattice:{row.lattice_echo}")
        elif row.lattice_echo:
            echoes.append(f"lattice:{row.lattice_echo}")
        literal = fmt(row.literal_h)
        if row.literal_h in FORBIDDEN_H:
            echoes.append("literal-H-FORBIDDEN")
        print(
            f"{row.n:2d} {row.exact_u:7d} {row.harborth_u:10d} "
            f"{row.spine:5d} {row.exact_tiles:11d} {row.harborth_tiles:13d} "
            f"{literal:>9s} {'; '.join(echoes) if echoes else '-'}"
        )
    print()


def print_forbidden_splits() -> None:
    print("FORBIDDEN-H CARRIER SPLITS")
    print("--------------------------")
    row5 = next(row for row in echo_rows() if row.n == 5)
    row11 = next(row for row in echo_rows() if row.n == 11)
    print(
        f"n=5 exact unit edges: {row5.exact_u} = spine {row5.spine} "
        f"+ tile/bulk {row5.exact_tiles}; literal unit-tournament H={row5.literal_h}."
    )
    print(
        "  Reading: the H=7 obstruction is not contradicted.  The scalar 7 is"
    )
    print(
        "  legal here because it is an additive edge count split across carriers,"
    )
    print(
        "  not a single tournament H evaluation."
    )
    print(
        f"n=11 triangular/Harborth edges: {row11.harborth_u} = spine {row11.spine} "
        f"+ lattice tile/bulk {row11.harborth_tiles}; exact planar u(11)={row11.exact_u}."
    )
    print(
        "  Reading: H=21 appears as a triangular-lattice lower-bound echo, while"
    )
    print(
        "  the exact small planar optimum has already jumped to 23.  Again the"
    )
    print(
        "  21 is an edge/carrier scalar, not a literal H=21 tournament."
    )
    print()


def print_tournament_analysis() -> None:
    print("TOURNAMENT ANALYSIS OVER CARRIER PACKETS")
    print("----------------------------------------")
    print("Vertices are carrier packets/proof obligations, not points.")
    print(
        "Pairwise observable: which carrier better preserves traceability, "
        "the H-gap guardrail, and the exact-vs-lattice scaling channel."
    )
    print("Switches: geometry, H-gap, and scaling gauges; ties use insertion HP.")
    print()
    gauge_adjs = {name: tournament_from_gauge(name) for name in GAUGES}
    for name, adj in gauge_adjs.items():
        scores = [sum(row) for row in adj]
        path = insertion_hp(adj)
        print(f"{name} gauge:")
        print(
            f"  score_hist={score_hist(adj)} directed_3cycles={directed_triangles(adj)} "
            f"SCCs={sccs(adj)} H={hamiltonian_path_count(adj)}"
        )
        print("  tie Hamiltonian path:")
        for idx in path:
            print(f"    score={scores[idx]} {CARRIERS[idx].name}: {CARRIERS[idx].note}")
        print()
    for a, b in combinations(gauge_adjs, 2):
        print(f"edge flips {a} vs {b}: {edge_flips(gauge_adjs[a], gauge_adjs[b])}/36")
    print()
    majority = majority_tournament()
    scores = [sum(row) for row in majority]
    print("majority carrier tournament:")
    print(
        f"  score_hist={score_hist(majority)} directed_3cycles={directed_triangles(majority)} "
        f"SCCs={sccs(majority)} H={hamiltonian_path_count(majority)}"
    )
    print("  insertion Hamiltonian path:")
    for idx in insertion_hp(majority):
        print(f"    score={scores[idx]} {CARRIERS[idx].name}: {CARRIERS[idx].note}")
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("--------------------")
    print(
        "Alternate vertices considered: points, unit pairs, nonunit pairs, "
        "distance classes, Hamiltonian paths, frontier additions, centered "
        "hex shell events, OCF cycle packets, LRC residues, and proof obligations."
    )
    print(
        "Chosen vertices: carrier packets/proof obligations.  This keeps the "
        "question close to LRC-style tournament analysis, where the vertex is "
        "the retained side channel rather than the raw object being counted."
    )
    print(
        "Preserved predicate: unit-spine traceability plus the distinction "
        "between additive edge counts and forbidden tournament H evaluations."
    )
    print(
        "Destroyed information: continuous embeddings, all optimal planar graph "
        "representatives, exact less-than/greater-than-one nonunit distances, "
        "and literal H values beyond the rows already computed."
    )
    print(
        "Challenged assumption: a scalar equality like u(5)=7 means the same "
        "thing as H(T)=7.  S627 treats it as an equidecomposability question: "
        "which pieces make up the same visible number?"
    )
    print()


def main() -> None:
    print("S627 UNIT-DISTANCE H-GAP CARRIER REFRAMING")
    print("===========================================")
    print()
    print(
        "Thesis: the meaningful H=7/H=21 relation is not a literal "
        "unit-distance tournament with 7 or 21 Hamiltonian paths.  It is a "
        "forbidden-collapse test: raw unit-edge counts may hit 7 or 21 only "
        "after their spine/tile/frontier carriers have been forgotten."
    )
    print()
    print_echo_table()
    print_forbidden_splits()
    print_tournament_analysis()
    print_assumption_challenge()
    print("NEXT PROBES")
    print("-----------")
    print("1. Compute literal unit-tournament H for exact Moser/planar rows n=9..14.")
    print("2. Test whether known non-lattice n>=22 candidates keep a unit-spine HP.")
    print("3. Build a shared carrier entropy variable for Sawin's n^1.014 side and")
    print("   the feasible OCF/proof-obligation side, avoiding raw scalar matching.")


if __name__ == "__main__":
    main()
