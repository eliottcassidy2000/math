#!/usr/bin/env python3
"""
tournament_even_minor_obstruction_codex_s96.py

Exact small audit for the Kuratowski/Wagner graph-minor analogy.

We work in the fixed-Hamiltonian-path tournament cube.  The path
0 -> 1 -> ... -> n-1 is fixed, and the free arcs are exactly the pairs
(i,j) with j > i+1.  There are C(n-1,2) of them, the same dimension as the
cycle space of K_n.  A concrete bijection sends the free arc (i,j) to the
triangle-basis coefficient (0, i+1, j) of an even graph on the same vertices.

This script checks what survives when the even-graph shadow is simplified by
degree-2 smoothing or by GF(2) edge contraction.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import comb


def free_pairs(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 2, n)]


def bits_of_mask(n: int, mask: int) -> tuple[int, ...]:
    return tuple((mask >> k) & 1 for k in range(len(free_pairs(n))))


def tournament_from_bits(n: int, bits: tuple[int, ...]) -> list[list[bool]]:
    adj = [[False] * n for _ in range(n)]

    for i in range(n - 1):
        adj[i][i + 1] = True
        adj[i + 1][i] = False

    for bit, (i, j) in zip(bits, free_pairs(n)):
        if bit == 0:
            adj[i][j] = True
            adj[j][i] = False
        else:
            adj[j][i] = True
            adj[i][j] = False
    return adj


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += value
    return sum(dp[full])


def strong_components(adj: list[list[bool]]) -> list[tuple[int, ...]]:
    n = len(adj)
    reach = [[adj[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen = [False] * n
    comps: list[tuple[int, ...]] = []
    for i in range(n):
        if seen[i]:
            continue
        comp = tuple(j for j in range(n) if reach[i][j] and reach[j][i])
        for j in comp:
            seen[j] = True
        comps.append(comp)
    return comps


def even_graph_from_bits(n: int, bits: tuple[int, ...]) -> frozenset[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()

    def toggle(a: int, b: int) -> None:
        if a > b:
            a, b = b, a
        edge = (a, b)
        if edge in edges:
            edges.remove(edge)
        else:
            edges.add(edge)

    for bit, (i, j) in zip(bits, free_pairs(n)):
        if bit == 0:
            continue
        a = i + 1
        b = j
        toggle(0, a)
        toggle(a, b)
        toggle(0, b)
    return frozenset(edges)


def bits_from_even_graph(n: int, edges: frozenset[tuple[int, int]]) -> tuple[int, ...]:
    edge_set = set(edges)
    bits = []
    for i, j in free_pairs(n):
        a = i + 1
        b = j
        if a > b:
            a, b = b, a
        bits.append(1 if (a, b) in edge_set else 0)
    return tuple(bits)


def degrees(n: int, edges: frozenset[tuple[int, int]]) -> list[int]:
    deg = [0] * n
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    return deg


def normalize_keep_vertices(n: int, edges: frozenset[tuple[int, int]]) -> tuple[int, tuple[tuple[int, int], ...]]:
    return n, tuple(sorted(edges))


def normalize_drop_isolates(n: int, edges: frozenset[tuple[int, int]]) -> tuple[int, tuple[tuple[int, int], ...]]:
    active = sorted({v for edge in edges for v in edge})
    if not active:
        return 0, ()
    relabel = {v: i for i, v in enumerate(active)}
    out = tuple(sorted((min(relabel[a], relabel[b]), max(relabel[a], relabel[b])) for a, b in edges))
    return len(active), out


def suppress_vertex_gf2(
    n: int, edges: frozenset[tuple[int, int]], v: int
) -> tuple[int, frozenset[tuple[int, int]]]:
    neigh = sorted((a if b == v else b) for a, b in edges if a == v or b == v)
    assert len(neigh) == 2
    a, b = neigh
    new_edges = {edge for edge in edges if v not in edge}
    chord = (min(a, b), max(a, b))
    if chord in new_edges:
        new_edges.remove(chord)
    else:
        new_edges.add(chord)

    def relabel(x: int) -> int:
        return x - 1 if x > v else x

    relabeled = frozenset(
        (min(relabel(x), relabel(y)), max(relabel(x), relabel(y)))
        for x, y in new_edges
    )
    return n - 1, relabeled


def greedy_suppression_core(
    state: tuple[int, tuple[tuple[int, int], ...]]
) -> frozenset[tuple[int, tuple[tuple[int, int], ...]]]:
    n, edge_tuple = state
    edges = frozenset(edge_tuple)
    while True:
        deg = degrees(n, edges)
        candidates = [v for v, d in enumerate(deg) if d == 2]
        if not candidates:
            return frozenset({normalize_drop_isolates(n, edges)})
        n, edges = suppress_vertex_gf2(n, edges, candidates[0])


def contract_edge_gf2(
    n: int, edges: frozenset[tuple[int, int]], edge: tuple[int, int]
) -> tuple[int, frozenset[tuple[int, int]]]:
    a, b = edge
    if a > b:
        a, b = b, a

    def merge_label(x: int) -> int:
        if x == b:
            x = a
        if x > b:
            x -= 1
        return x

    out: set[tuple[int, int]] = set()
    for x, y in edges:
        nx = merge_label(x)
        ny = merge_label(y)
        if nx == ny:
            continue
        new_edge = (min(nx, ny), max(nx, ny))
        if new_edge in out:
            out.remove(new_edge)
        else:
            out.add(new_edge)
    return n - 1, frozenset(out)


def format_core(core: tuple[int, tuple[tuple[int, int], ...]]) -> str:
    n, edges = core
    if not edges:
        return f"v{n}:empty"
    return f"v{n}:" + ",".join(f"{a}{b}" for a, b in edges)


def format_core_set(cores: frozenset[tuple[int, tuple[tuple[int, int], ...]]]) -> str:
    return "{" + "; ".join(format_core(c) for c in sorted(cores)) + "}"


def first_examples(records: list[dict]) -> tuple[dict | None, dict | None, dict | None]:
    by_core: dict[
        tuple[int, frozenset[tuple[int, tuple[tuple[int, int], ...]]]], list[dict]
    ] = defaultdict(list)
    for row in records:
        by_core[(row["n"], row["cores"])].append(row)

    same_core_diff_h = None
    for rows in by_core.values():
        h_values = sorted({row["H"] for row in rows})
        if len(h_values) >= 2:
            lo_h, hi_h = h_values[0], h_values[-1]
            lo = next(row for row in rows if row["H"] == lo_h)
            hi = next(row for row in rows if row["H"] == hi_h)
            candidate = {"lo": lo, "hi": hi, "h_values": h_values, "count": len(rows)}
            if same_core_diff_h is None:
                same_core_diff_h = candidate
            else:
                old_score = (len(same_core_diff_h["h_values"]), same_core_diff_h["hi"]["n"])
                new_score = (len(candidate["h_values"]), candidate["hi"]["n"])
                if new_score > old_score:
                    same_core_diff_h = candidate

    contraction_up = None
    contraction_down = None
    for row in records:
        n = row["n"]
        if n <= 3:
            continue
        for edge in row["edges"]:
            child_n, child_edges = contract_edge_gf2(n, row["edges"], edge)
            child_bits = bits_from_even_graph(child_n, child_edges)
            child_h = hamiltonian_paths(tournament_from_bits(child_n, child_bits))
            if child_h > row["H"] and contraction_up is None:
                contraction_up = {
                    "parent": row,
                    "edge": edge,
                    "child_n": child_n,
                    "child_edges": child_edges,
                    "child_H": child_h,
                    "child_bits": child_bits,
                }
            if child_h < row["H"] and contraction_down is None:
                contraction_down = {
                    "parent": row,
                    "edge": edge,
                    "child_n": child_n,
                    "child_edges": child_edges,
                    "child_H": child_h,
                    "child_bits": child_bits,
                }
            if contraction_up and contraction_down:
                return same_core_diff_h, contraction_up, contraction_down
    return same_core_diff_h, contraction_up, contraction_down


def main() -> None:
    all_records: list[dict] = []

    print("FIXED-PATH TOURNAMENT <-> EVEN-GRAPH MINOR AUDIT")
    print("free arc (i,j), j>i+1, maps to triangle-basis bit (0,i+1,j)")
    print()

    for n in range(3, 8):
        m = comb(n - 1, 2)
        total = 1 << m
        h_counter: Counter[int] = Counter()
        strong_h_counter: Counter[int] = Counter()
        core_to_h: dict[frozenset[tuple[int, tuple[tuple[int, int], ...]]], set[int]] = defaultdict(set)
        contraction_pairs = 0
        contraction_even_failures = 0

        for mask in range(total):
            bits = bits_of_mask(n, mask)
            adj = tournament_from_bits(n, bits)
            h = hamiltonian_paths(adj)
            edges = even_graph_from_bits(n, bits)
            cores = greedy_suppression_core(normalize_keep_vertices(n, edges))
            comps = strong_components(adj)

            h_counter[h] += 1
            if len(comps) == 1:
                strong_h_counter[h] += 1
            core_to_h[cores].add(h)

            for edge in edges:
                child_n, child_edges = contract_edge_gf2(n, edges, edge)
                if any(d % 2 for d in degrees(child_n, child_edges)):
                    contraction_even_failures += 1
                contraction_pairs += 1

            all_records.append(
                {
                    "n": n,
                    "mask": mask,
                    "bits": bits,
                    "H": h,
                    "edges": edges,
                    "cores": cores,
                    "scc_count": len(comps),
                }
            )

        multi_h_cores = sum(1 for hs in core_to_h.values() if len(hs) > 1)
        biggest_h_fiber = max((len(hs) for hs in core_to_h.values()), default=0)
        print(f"n={n}")
        print(f"  rooted tournaments = even graphs = 2^C(n-1,2) = {total}")
        print(f"  H spectrum min/max/count = {min(h_counter)}/{max(h_counter)}/{len(h_counter)}")
        print(f"  H=7 count {h_counter.get(7, 0)}; H=21 count {h_counter.get(21, 0)}")
        print(f"  strong H values <=25 = {[h for h in sorted(strong_h_counter) if h <= 25]}")
        print(
            "  greedy suppression-core buckets = "
            f"{len(core_to_h)}, buckets with >1 H = {multi_h_cores}, max distinct H per bucket = {biggest_h_fiber}"
        )
        print(
            "  GF(2) contractions checked = "
            f"{contraction_pairs}, even-degree failures = {contraction_even_failures}"
        )
        print()

    same_core, up, down = first_examples(all_records)

    print("EXAMPLE: degree-2 suppression is too coarse for H")
    if same_core is None:
        print("  No example found.")
    else:
        lo = same_core["lo"]
        hi = same_core["hi"]
        print(f"  deterministic terminal GF(2)-suppression core set: {format_core_set(lo['cores'])}")
        print(
            "  same core bucket contains "
            f"{same_core['count']} rooted tournaments and H-values {same_core['h_values']}"
        )
        print(f"  low-H witness:  n={lo['n']} mask={lo['mask']} bits={lo['bits']} H={lo['H']}")
        print(f"  high-H witness: n={hi['n']} mask={hi['mask']} bits={hi['bits']} H={hi['H']}")
        print()

    print("EXAMPLE: GF(2) edge contraction is not an H-monotone tournament operation")
    if up is None:
        print("  No increasing example found through n=7.")
    else:
        parent = up["parent"]
        print(
            "  increasing contraction: "
            f"n={parent['n']} mask={parent['mask']} H={parent['H']} "
            f"contract edge={up['edge']} -> n={up['child_n']} H={up['child_H']}"
        )
        print(f"    parent even edges={sorted(parent['edges'])}")
        print(f"    child even edges ={sorted(up['child_edges'])}")
    if down is None:
        print("  No decreasing example found through n=7.")
    else:
        parent = down["parent"]
        print(
            "  decreasing contraction: "
            f"n={parent['n']} mask={parent['mask']} H={parent['H']} "
            f"contract edge={down['edge']} -> n={down['child_n']} H={down['child_H']}"
        )
        print(f"    parent even edges={sorted(parent['edges'])}")
        print(f"    child even edges ={sorted(down['child_edges'])}")
    print()

    print("INTERPRETATION")
    print("  1. The fixed-path cube is the right labelled bridge: it sees every tournament")
    print("     up to relabelling because every tournament has a Hamiltonian path.")
    print("  2. The even-graph quotient is a GF(2) cycle-space object. Degree-2 smoothing")
    print("     cancels triangles and other cycles modulo 2, so it is not the same data")
    print("     as topological-minor suppression in Kuratowski's theorem.")
    print("  3. GF(2) contraction stays inside even graphs, but the induced tournament H")
    print("     is not monotone. Thus the forbidden values {7,21} are not behaving like")
    print("     forbidden minors of a minor-closed graph family.")
    print("  4. The surviving reduction is the known one: strong-component condensation")
    print("     gives a multiplicative semigroup for H. The {7,21} obstruction lives in")
    print("     irreducible conflict-graph/strong-tournament cores, not in graph minors.")


if __name__ == "__main__":
    main()
