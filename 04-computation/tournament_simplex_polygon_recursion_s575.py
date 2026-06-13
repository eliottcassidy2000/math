#!/usr/bin/env python3
"""
tournament_simplex_polygon_recursion_s575.py

codex-2026-06-03 S575

Audit the user's "tournaments are both simplices and regular polygons" frame.

Dictionary:
  * Simplex/mesh view: a tournament orients every edge of the complete graph,
    the 1-skeleton of the (m-1)-simplex.  Boundary ties require interior/mesh
    choices to become a full tournament.
  * Polygon/outside view: a circular configuration gives a half-turn
    tournament from the outside order.  The regular m-gon is clean only for odd
    m; for even m it has antipodal tie diagonals.

LRC uses m=n-1 runner sub-tournaments.  Therefore clean regular-polygon
dihedral behavior occurs exactly at even LRC n, every other n.  The n -> n+2
recursion preserves this parity:
  * even n: m odd, outside determines R_m, with C_m automorphisms and m
    reflection anti-automorphisms (a dihedral automorphism+anti-automorphism
    object of size 2m);
  * odd n: m even, the regular outside is a wall with m/2 antipodal tie
    diagonals and 2^(m/2) labelled tie-resolution mesh choices.

Tournament Analysis:
  Vertices: LRC n-ladder entries.
  Pairwise observable: (clean polygon flag, extended dihedral size, tie pairs,
    round count, all-tournament count, regular 3-cycle mesh count).
  Switch/gauge: cleaner polygon quotient wins; within a parity ladder, larger
    extended dihedral size wins, then fewer tie choices, then smaller round/all
    fraction.  Tie Hamiltonian path: increasing n.

Assumption challenge:
  Vertices need not be tournaments only.  Consider simplex faces, polygon
  vertices, antipodal tie diagonals, boundary resolutions, dihedral orbits,
  SCC blocks, and LRC proof obligations.  This script preserves parity,
  outside-vs-mesh ambiguity, and dihedral extended symmetry; it destroys exact
  LRC runner ownership, wall endpoint provenance, and non-regular iso-class
  detail.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations
from math import factorial, gcd, log10


@dataclass(frozen=True)
class LadderRow:
    n_lrc: int
    m_runners: int
    simplex_dim: int
    simplex_edges: int
    polygon_state: str
    tie_pairs: int
    tie_resolution_choices: int
    cyclic_automorphisms: int
    reflection_anti_automorphisms: int
    extended_dihedral_size: int
    regular_c3: int | None
    round_classes: int
    all_tournament_classes: int


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def totient(n: int) -> int:
    return sum(1 for k in range(1, n + 1) if gcd(k, n) == 1)


def a000016_round(m: int) -> int:
    return sum(totient(d) * 2 ** (m // d) for d in divisors(m) if d % 2 == 1) // (2 * m)


def integer_partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if max_part is None or max_part > n:
        max_part = n
    if n == 0:
        return [()]
    out: list[tuple[int, ...]] = []
    for first in range(min(max_part, n), 0, -1):
        for rest in integer_partitions(n - first, first):
            out.append((first,) + rest)
    return out


def a000568_tournaments(n: int) -> int:
    # Unlabeled tournaments by Burnside over cycle types.  Only odd permutation
    # cycles contribute.
    total = 0
    for parts in integer_partitions(n):
        if any(p % 2 == 0 for p in parts):
            continue
        counts = Counter(parts)
        autom = 1
        for length, multiplicity in counts.items():
            autom *= (length ** multiplicity) * factorial(multiplicity)
        exponent = sum((p - 1) // 2 for p in parts)
        for i in range(len(parts)):
            for j in range(i + 1, len(parts)):
                exponent += gcd(parts[i], parts[j])
        total += factorial(n) // autom * 2**exponent
    return total // factorial(n)


def regular_round_adj(m: int) -> list[list[bool]]:
    h = (m - 1) // 2
    return [[i != j and ((j - i) % m) in range(1, h + 1) for j in range(m)] for i in range(m)]


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    c3 = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1
    return c3


def row_for_n(n_lrc: int) -> LadderRow:
    m = n_lrc - 1
    edges = m * (m - 1) // 2
    round_classes = a000016_round(m)
    all_classes = a000568_tournaments(m)
    if m % 2 == 1:
        regular_c3 = count_directed_3cycles(regular_round_adj(m))
        return LadderRow(
            n_lrc=n_lrc,
            m_runners=m,
            simplex_dim=m - 1,
            simplex_edges=edges,
            polygon_state="clean_odd_polygon",
            tie_pairs=0,
            tie_resolution_choices=1,
            cyclic_automorphisms=m,
            reflection_anti_automorphisms=m,
            extended_dihedral_size=2 * m,
            regular_c3=regular_c3,
            round_classes=round_classes,
            all_tournament_classes=all_classes,
        )
    tie_pairs = m // 2
    return LadderRow(
        n_lrc=n_lrc,
        m_runners=m,
        simplex_dim=m - 1,
        simplex_edges=edges,
        polygon_state="even_polygon_wall",
        tie_pairs=tie_pairs,
        tie_resolution_choices=2**tie_pairs,
        cyclic_automorphisms=m,
        reflection_anti_automorphisms=m,
        extended_dihedral_size=2 * m,
        regular_c3=None,
        round_classes=round_classes,
        all_tournament_classes=all_classes,
    )


def fmt_ratio(row: LadderRow) -> str:
    ratio = row.round_classes / row.all_tournament_classes
    if ratio >= 0.001:
        return f"{ratio:.6f}"
    return f"10^{log10(ratio):.2f}"


def tournament_fingerprint(rows: list[LadderRow]) -> dict[str, object]:
    def key(row: LadderRow) -> tuple[int, int, int, float, int]:
        clean = 1 if row.tie_pairs == 0 else 0
        return (
            clean,
            row.extended_dihedral_size,
            -row.tie_resolution_choices,
            -row.round_classes / row.all_tournament_classes,
            row.n_lrc,
        )

    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            if i == j:
                continue
            adj[i][j] = key(left) > key(right) or (key(left) == key(right) and i < j)

    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]):
            c3 += 1

    def reach(start: int) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v, edge in enumerate(adj[u]):
                if edge and v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reach(u)
        comp = {v for v in remaining if v in ru and u in reach(v)}
        sccs.append(len(comp))
        remaining -= comp

    return {
        "vertices": [f"n={row.n_lrc}" for row in rows],
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": 1 if c3 == 0 and len(set(scores)) == n else "not_counted",
    }


def main() -> None:
    rows = [row_for_n(n) for n in range(4, 19)]
    print("S575 tournament simplex/polygon n+2 recursion audit")
    print("=" * 78)
    print("LRC n uses runner tournament size m=n-1.")
    print("Odd m: regular polygon is a clean rotational tournament with dihedral")
    print("automorphism+anti-automorphism symmetry.  Even m: regular polygon is a")
    print("wall with antipodal tie diagonals, so the simplex mesh must resolve it.")
    print()

    print("n  m  simplex  edges  polygon_state        ties  tie_choices  D_ext  c3(R_m)  round/all")
    for row in rows:
        c3 = "-" if row.regular_c3 is None else str(row.regular_c3)
        print(
            f"{row.n_lrc:2d} {row.m_runners:2d} {row.simplex_dim:8d} {row.simplex_edges:6d} "
            f"{row.polygon_state:19s} {row.tie_pairs:4d} {row.tie_resolution_choices:12d} "
            f"{row.extended_dihedral_size:5d} {c3:7s} {fmt_ratio(row):>10s}"
        )

    print("\nEvery-other pattern")
    clean_ns = [row.n_lrc for row in rows if row.tie_pairs == 0]
    wall_ns = [row.n_lrc for row in rows if row.tie_pairs != 0]
    print(f"  clean dihedral polygon ladder: n={clean_ns}")
    print(f"  wall/mesh tie-resolution ladder: n={wall_ns}")
    print("  n -> n+2 preserves the ladder: m -> m+2, simplex dimension -> +2.")
    print("  On the clean ladder, D_ext grows by 4 each step.")
    print("  On the wall ladder, tie pairs grow by 1 and labelled mesh choices double.")

    print("\nTournament Analysis")
    print("  vertices: LRC n-ladder entries")
    print("  observable: clean polygon flag, D_ext, tie choices, round/all, regular c3")
    print("  switch: cleaner polygon quotient wins; within parity, larger D_ext wins")
    print(f"  fingerprints: {tournament_fingerprint(rows)}")

    print("\nSynthesis")
    print("  A tournament is the simplex mesh when every edge/face choice matters.")
    print("  It is the polygon outside when the circular order determines all arcs.")
    print("  The regular outside determines a tournament exactly for odd m=n-1.")
    print("  Hence LRC even n is the clean dihedral ladder, and n -> n+2 is the")
    print("  parity-preserving recursion for that ladder.  LRC odd n is the wall")
    print("  ladder where antipodal tie diagonals force extra mesh data.")
    print("  HYP-2086/2087 live on the clean ladder's reflection/self-converse side;")
    print("  HYP-2089 says the strong lens sees the same regular encircling block.")


if __name__ == "__main__":
    main()
