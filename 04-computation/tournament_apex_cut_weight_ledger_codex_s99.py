#!/usr/bin/env python3
"""
tournament_apex_cut_weight_ledger_codex_s99.py

Exact audit of the owner's n=8 apex-cut lead.

We work in the fixed-Hamiltonian-path cube.  Start with a rooted n=7
tournament on vertices 0..6 with fixed path 0->1->...->6.  Add the apex
vertex 7 at the end of the fixed path, so 6->7 is fixed and the six free
apex-row bits are the arcs between 7 and vertices 0..5.

For an old Hamiltonian path P, the apex contributes one unit for every legal
insertion slot:

    before P[0]      if 7 -> P[0],
    after P[-1]      if P[-1] -> 7,
    between a,b      if a -> 7 -> b.

Thus H(T+7) is the sum of a path-weight distribution over H(T) old paths.  The
raw cut cardinality w = #{i<=5 : 7 -> i} is useful, but the actual apex-tile
carrier is this insertion-weight profile.  The script audits all 2^15 rooted
n=7 bases and all 64 apex cuts, i.e. all 2^21 rooted n=8 rows, but reports the
two weights in the prompt most closely:

  * w=3, balanced apex cuts;
  * w=1, a single unbalanced apex win.

Tournament Analysis declaration.
  Vertices: proof carriers for the n=8 strong-value extension.
  Pairwise observable:
    (preserves H exactly, isolates apex tile, supplies n=8 strong values,
     separates balanced/unbalanced cuts, proof usefulness).
  Gauge: lexicographic comparison; ties follow the declared path.
  Fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian paths.

Assumption challenge.
  I considered runners, free arcs, raw cut cardinality, insertion slots,
  strong components, old H-atoms, apex tile orientation, even-graph holes,
  and proof obligations as vertices.  The chosen quotient is the apex
  insertion-weight profile.  It preserves H exactly under apex extension and
  destroys most old-row geometry; raw cut cardinality alone is too coarse.
"""

from __future__ import annotations

import importlib.util
from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
S96_PATH = HERE / "tournament_even_minor_obstruction_codex_s96.py"
SPEC = importlib.util.spec_from_file_location("s96_even_minor", S96_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {S96_PATH}")
s96 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(s96)

OLD_N = 7
APEX = 7
SELECTED_WEIGHTS = (1, 3)
SELECTED_VALUES = (49, 75)


def row_bits(adj: list[list[bool]]) -> list[int]:
    return [sum((1 << j) for j, bit in enumerate(row) if bit) for row in adj]


def is_strong_from_rows(rows: list[int]) -> bool:
    n = len(rows)
    reach = [rows[i] | (1 << i) for i in range(n)]
    for k in range(n):
        bit = 1 << k
        reach_k = reach[k]
        for i in range(n):
            if reach[i] & bit:
                reach[i] |= reach_k
    full = (1 << n) - 1
    return all(row == full for row in reach)


def hamiltonian_path_stats(adj: list[list[bool]]) -> tuple[list[int], list[int], list[list[int]], int]:
    n = len(adj)
    start = [0] * n
    end = [0] * n
    adjacent_pair = [[0] * n for _ in range(n)]
    used = [False] * n
    path: list[int] = []

    def dfs(last: int | None = None) -> None:
        if len(path) == n:
            start[path[0]] += 1
            end[path[-1]] += 1
            for a, b in zip(path, path[1:]):
                adjacent_pair[a][b] += 1
            return
        if last is None:
            for v in range(n):
                used[v] = True
                path.append(v)
                dfs(v)
                path.pop()
                used[v] = False
            return
        for v, edge in enumerate(adj[last]):
            if edge and not used[v]:
                used[v] = True
                path.append(v)
                dfs(v)
                path.pop()
                used[v] = False

    dfs()
    base_h = sum(end)
    return start, end, adjacent_pair, base_h


def apex_wins(cut: int) -> list[int]:
    return [(cut >> i) & 1 for i in range(OLD_N - 1)] + [0]


def h_from_stats(start: list[int], end: list[int], adjacent_pair: list[list[int]], cut: int) -> int:
    wins = apex_wins(cut)
    total = sum(start[v] * wins[v] for v in range(OLD_N))
    total += sum(end[v] * (1 - wins[v]) for v in range(OLD_N))
    for a in range(OLD_N):
        if wins[a]:
            continue
        for b in range(OLD_N):
            if wins[b]:
                total += adjacent_pair[a][b]
    return total


def insertion_distribution(adj: list[list[bool]], cut: int) -> Counter[int]:
    wins = apex_wins(cut)
    dist: Counter[int] = Counter()
    used = [False] * OLD_N
    path: list[int] = []

    def dfs(last: int | None = None) -> None:
        if len(path) == OLD_N:
            weight = int(bool(wins[path[0]])) + int(not wins[path[-1]])
            for a, b in zip(path, path[1:]):
                if not wins[a] and wins[b]:
                    weight += 1
            dist[weight] += 1
            return
        if last is None:
            for v in range(OLD_N):
                used[v] = True
                path.append(v)
                dfs(v)
                path.pop()
                used[v] = False
            return
        for v, edge in enumerate(adj[last]):
            if edge and not used[v]:
                used[v] = True
                path.append(v)
                dfs(v)
                path.pop()
                used[v] = False

    dfs()
    return dist


def extension_rows(base_rows: list[int], cut: int) -> list[int]:
    rows: list[int] = []
    for i, row in enumerate(base_rows):
        if i == OLD_N - 1 or not ((cut >> i) & 1):
            rows.append(row | (1 << APEX))
        else:
            rows.append(row)
    apex_row = 0
    for i in range(OLD_N - 1):
        if (cut >> i) & 1:
            apex_row |= 1 << i
    rows.append(apex_row)
    return rows


def free_reversals(n: int, mask: int) -> list[tuple[int, int]]:
    return [
        pair
        for bit, pair in zip(s96.bits_of_mask(n, mask), s96.free_pairs(n))
        if bit
    ]


def small_gaps(values: set[int], lo: int = 1, hi: int = 100) -> list[int]:
    return [v for v in range(lo, hi + 1) if v not in values]


def path_tournament(vertices: list[str], order: list[str]) -> dict[tuple[str, str], str]:
    rank = {name: i for i, name in enumerate(order)}
    winners: dict[tuple[str, str], str] = {}
    for a, b in combinations(vertices, 2):
        winners[(a, b)] = a if rank[a] < rank[b] else b
    return winners


def tournament_analysis() -> None:
    vertices = [
        "apex-insertion-profile",
        "labelled-strong-H-atom",
        "balanced-w3-cut",
        "single-unbalanced-w1-cut",
        "raw-cut-cardinality",
        "even-graph-hole-address",
        "runner/speed-subset",
    ]
    winners = path_tournament(vertices, vertices)
    scores = Counter(winners.values())
    for vertex in vertices:
        scores.setdefault(vertex, 0)

    def edge(a: str, b: str) -> bool:
        key = (a, b) if vertices.index(a) < vertices.index(b) else (b, a)
        return winners[key] == a

    def winner(a: str, b: str) -> str:
        key = (a, b) if vertices.index(a) < vertices.index(b) else (b, a)
        return winners[key]

    directed_triangles = 0
    for a, b, c in combinations(vertices, 3):
        out = Counter([winner(a, b), winner(a, c), winner(b, c)])
        if sorted(out.values()) == [1, 1, 1]:
            directed_triangles += 1
    hpaths = sum(
        1
        for path in permutations(vertices)
        if all(edge(path[i], path[i + 1]) for i in range(len(path) - 1))
    )
    print("TOURNAMENT ANALYSIS")
    print("  vertices: proof carriers for the apex-extension route, not runners")
    print(f"  Hamiltonian path: {' > '.join(vertices)}")
    print(f"  score histogram: {dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed 3-cycles: {directed_triangles}; SCCs: singleton; Hamiltonian-path count: {hpaths}")
    print("  challenged assumption: raw balanced/unbalanced cut size is the whole carrier")


def main() -> None:
    cuts_by_w = {
        w: [
            sum(1 << i for i in subset)
            for subset in combinations(range(OLD_N - 1), w)
        ]
        for w in range(OLD_N)
    }
    values_by_w: dict[int, set[int]] = defaultdict(set)
    strong_rows_by_w: Counter[int] = Counter()
    examples: dict[tuple[int, int], tuple[int, int, int, bool]] = {}
    base_h_sources: dict[tuple[int, int], Counter[int]] = defaultdict(Counter)
    base_strong_sources: dict[tuple[int, int], Counter[bool]] = defaultdict(Counter)

    base_count = 1 << comb(OLD_N - 1, 2)
    for mask in range(base_count):
        bits = s96.bits_of_mask(OLD_N, mask)
        adj = s96.tournament_from_bits(OLD_N, bits)
        base_rows = row_bits(adj)
        base_strong = is_strong_from_rows(base_rows)
        start, end, adjacent_pair, base_h = hamiltonian_path_stats(adj)
        for w, cuts in cuts_by_w.items():
            for cut in cuts:
                h = h_from_stats(start, end, adjacent_pair, cut)
                if is_strong_from_rows(extension_rows(base_rows, cut)):
                    values_by_w[w].add(h)
                    strong_rows_by_w[w] += 1
                    key = (w, h)
                    examples.setdefault(key, (mask, cut, base_h, base_strong))
                    base_h_sources[key][base_h] += 1
                    base_strong_sources[key][base_strong] += 1

    print("=" * 78)
    print("APEX CUT WEIGHT LEDGER n=7 -> n=8")
    print("=" * 78)
    print(f"rooted n=7 bases checked: {base_count}")
    print(f"apex cuts per base: 64")
    print(f"rooted n=8 rows represented: {base_count * 64}")
    print("fixed edge: 6 -> 7; cut weight w = #{i in 0..5 : 7 -> i}")
    print()
    print("Strong n=8 extensions by raw apex cut weight:")
    for w in range(OLD_N):
        values = values_by_w[w]
        print(
            f"  w={w}: cuts={len(cuts_by_w[w]):2d}, strong_rows={strong_rows_by_w[w]:7d}, "
            f"H_count={len(values):3d}, min={min(values) if values else None}, "
            f"max={max(values) if values else None}, "
            f"has49={49 in values}, has75={75 in values}"
        )
        if w in SELECTED_WEIGHTS:
            print(f"       gaps in [1,100]: {small_gaps(values)}")
    print()
    print("Selected value witnesses:")
    for w in SELECTED_WEIGHTS:
        for h in SELECTED_VALUES:
            mask, cut, base_h, base_strong = examples[(w, h)]
            adj = s96.tournament_from_bits(OLD_N, s96.bits_of_mask(OLD_N, mask))
            dist = insertion_distribution(adj, cut)
            beats = [i for i in range(OLD_N - 1) if (cut >> i) & 1]
            print(f"  w={w}, H8={h}: base_mask={mask}, cut={cut}, apex_beats={beats}")
            print(
                f"       base_H={base_h}, base_strong={base_strong}, "
                f"insertion_weight_distribution={dict(sorted(dist.items()))}"
            )
            print(f"       base free reversals={free_reversals(OLD_N, mask)}")
            print(f"       source base_H counts={dict(base_h_sources[(w, h)].most_common())}")
            print(f"       source base_strong counts={dict(base_strong_sources[(w, h)].most_common())}")
    print()
    print("Interpretation:")
    print("  Raw w=3 balanced cuts produce a wide block of strong n=8 H-values,")
    print("  but raw cut size does not isolate the 49/75 phenomenon: w=1 also supplies")
    print("  both, and rooted w=3 can supply them too.  HYP-2879's non-isomorphic")
    print("  strong-ear quotient is sharper: balanced w=3 misses exactly 49,75 there,")
    print("  and boundary w=1 fills them.  The rooted ledger explains why raw cut")
    print("  cardinality is only an address.  The sharper carrier is the insertion")
    print("  profile.  The single-unbalanced witnesses are nearly constant-weight:")
    print("    49 = 2*25 - 1 from distribution {1:1, 2:24};")
    print("    75 = 2*39 - 3 from distribution {1:3, 2:36}.")
    print("  That odd defect is the plausible apex-tile signal: one boundary row")
    print("  almost doubles the old H-atom, then subtracts a small source-sink defect.")
    print()
    tournament_analysis()
    print("\nDONE.")


if __name__ == "__main__":
    main()
