#!/usr/bin/env python3
"""Independent-pair channel audit for the n=14 LRC multi-channel route.

codex-2026-06-01 S532

User insight: for the four-vertex tournament, two independent pair-arcs can
coordinate the isomorphism class when the remaining arcs keep a fixed scaffold.
This script verifies the K4 toy model and then measures the corresponding
independent-pair inventory in the clasp-deleted regular 14-gon channels.

Tournament Analysis declaration:
    vertices: chord channels of the clasp-deleted regular 14-gon
    pairwise observable: independent-pair inventory
        (maximum matching size, edge count, number of maximum matchings)
    switch/gauge: lexicographic comparison of that inventory
    tie Hamiltonian path: increasing chord distance
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, permutations, product


def canon(edges: set[tuple[int, int]], n: int) -> str:
    """Canonical labelled bitstring for an unlabeled tournament on n vertices."""
    best = None
    for perm in permutations(range(n)):
        bits = []
        for i, j in combinations(range(n), 2):
            a, b = perm[i], perm[j]
            bits.append("1" if (a, b) in edges else "0")
        word = "".join(bits)
        if best is None or word < best:
            best = word
    assert best is not None
    return best


def scores(edges: set[tuple[int, int]], n: int) -> tuple[int, ...]:
    return tuple(
        sorted(sum((i, j) in edges for j in range(n) if i != j) for i in range(n))
    )


def hamiltonian_paths(edges: set[tuple[int, int]], n: int) -> int:
    total = 0
    for perm in permutations(range(n)):
        if all((perm[i], perm[i + 1]) in edges for i in range(n - 1)):
            total += 1
    return total


def cyclic_triangles(edges: set[tuple[int, int]], n: int) -> int:
    total = 0
    for tri in combinations(range(n), 3):
        out = [
            sum((x, y) in edges for y in tri if y != x)
            for x in tri
        ]
        if sorted(out) == [1, 1, 1]:
            total += 1
    return total


def scc_sizes(edges: set[tuple[int, int]], n: int) -> tuple[int, ...]:
    reach = [[i == j for j in range(n)] for i in range(n)]
    for a, b in edges:
        reach[a][b] = True
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])

    seen: set[int] = set()
    sizes: list[int] = []
    for i in range(n):
        if i in seen:
            continue
        comp = [j for j in range(n) if reach[i][j] and reach[j][i]]
        seen.update(comp)
        sizes.append(len(comp))
    return tuple(sorted(sizes, reverse=True))


def k4_coordinate_rows() -> list[dict[str, object]]:
    """A K4 scaffold where two independent pair-arcs hit all four iso classes."""
    n = 4
    # Fixed scaffold: transitive spine 0->1->2->3 plus long outside edge 0->3.
    fixed = {(0, 1), (1, 2), (2, 3), (0, 3)}
    free = [(0, 2), (1, 3)]
    rows: list[dict[str, object]] = []
    for bits in product([0, 1], repeat=2):
        edges = set(fixed)
        for bit, (a, b) in zip(bits, free):
            edges.add((a, b) if bit else (b, a))
        rows.append(
            {
                "state": "".join(str(b) for b in bits),
                "edges": tuple(sorted(edges)),
                "canon": canon(edges, n),
                "scores": scores(edges, n),
                "H": hamiltonian_paths(edges, n),
                "c3": cyclic_triangles(edges, n),
                "scc": scc_sizes(edges, n),
            }
        )
    return rows


def k4_scaffold_census() -> Counter[int]:
    """Count how often two independent free arcs coordinate K4 iso classes."""
    n = 4
    all_pairs = list(combinations(range(n), 2))
    census: Counter[int] = Counter()
    for free in combinations(all_pairs, 2):
        if len(set(free[0] + free[1])) != 4:
            continue
        rest = [edge for edge in all_pairs if edge not in free]
        for fixed_bits in product([0, 1], repeat=4):
            fixed: set[tuple[int, int]] = set()
            for bit, (a, b) in zip(fixed_bits, rest):
                fixed.add((a, b) if bit else (b, a))
            classes = set()
            for bits in product([0, 1], repeat=2):
                edges = set(fixed)
                for bit, (a, b) in zip(bits, free):
                    edges.add((a, b) if bit else (b, a))
                classes.add(canon(edges, n))
            census[len(classes)] += 1
    return census


def clasp_deleted_edges(n: int, distance: int) -> list[tuple[int, int]]:
    """Unoriented chord channel after deleting observer 0 from an n-gon."""
    edges = []
    runners = set(range(1, n))
    for a, b in combinations(sorted(runners), 2):
        delta = (b - a) % n
        dist = min(delta, n - delta)
        if dist == distance:
            edges.append((a, b))
    return edges


def max_matchings(edge_list: list[tuple[int, int]]) -> tuple[int, int]:
    """Return (maximum matching size, number of maximum matchings)."""
    best_size = 0
    best_count = 0

    def rec(idx: int, used: set[int], size: int) -> None:
        nonlocal best_size, best_count
        if idx == len(edge_list):
            if size > best_size:
                best_size = size
                best_count = 1
            elif size == best_size:
                best_count += 1
            return
        # Skip current edge.
        rec(idx + 1, used, size)
        a, b = edge_list[idx]
        if a not in used and b not in used:
            used.add(a)
            used.add(b)
            rec(idx + 1, used, size + 1)
            used.remove(a)
            used.remove(b)

    rec(0, set(), 0)
    return best_size, best_count


def channel_stats(n: int = 14) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for distance in range(1, n // 2 + 1):
        edges = clasp_deleted_edges(n, distance)
        mm, count = max_matchings(edges)
        degree = Counter()
        for a, b in edges:
            degree[a] += 1
            degree[b] += 1
        rows.append(
            {
                "distance": distance,
                "edges": len(edges),
                "max_matching": mm,
                "max_matching_count": count,
                "degree_hist": tuple(sorted(Counter(degree.values()).items())),
                "all_edges_independent": mm == len(edges),
            }
        )
    return rows


def tournament_analysis(rows: list[dict[str, object]]) -> dict[str, object]:
    """Tournament fingerprint over channel inventories."""
    labels = [int(row["distance"]) for row in rows]
    profile = {
        int(row["distance"]): (
            int(row["max_matching"]),
            int(row["edges"]),
            int(row["max_matching_count"]),
        )
        for row in rows
    }
    outdegree = Counter({label: 0 for label in labels})
    edges = []
    for a, b in combinations(labels, 2):
        winner, loser = (a, b) if profile[a] > profile[b] else (b, a)
        outdegree[winner] += 1
        edges.append((winner, loser))
    return {
        "score_hist": tuple(sorted(Counter(outdegree.values()).items())),
        "edges": tuple(edges),
        "tie_hamiltonian_path": tuple(sorted(labels)),
    }


def main() -> None:
    print("K4 independent-pair coordinate")
    print("fixed scaffold: 0->1, 1->2, 2->3, 0->3")
    print("free independent pairs: (0,2), (1,3)")
    rows = k4_coordinate_rows()
    for row in rows:
        print(
            "state={state} canon={canon} scores={scores} H={H} c3={c3} "
            "scc={scc}".format(**row)
        )
    print("distinct K4 iso classes reached:", len({row["canon"] for row in rows}))

    census = k4_scaffold_census()
    print("K4 independent-pair scaffold census:")
    for class_count in sorted(census):
        print(f"  {class_count} classes reached: {census[class_count]} scaffolds")

    print()
    print("n=14 clasp-deleted chord channels")
    rows14 = channel_stats(14)
    for row in rows14:
        print(
            "d={distance} edges={edges} max_matching={max_matching} "
            "max_matching_count={max_matching_count} degree_hist={degree_hist} "
            "all_edges_independent={all_edges_independent}".format(**row)
        )
    ta = tournament_analysis(rows14)
    print("Tournament Analysis over channel inventories:")
    print("  score_hist:", ta["score_hist"])
    print("  tie_hamiltonian_path:", ta["tie_hamiltonian_path"])
    print("  edges:", ta["edges"])
    print()
    print("Interpretation:")
    print("  K4: two independent pair states are a complete iso-class coordinate.")
    print("  n=14: the diameter channel is six independent pairs plus one singleton;")
    print("        multi-channel parity should track matching state, not only width.")


if __name__ == "__main__":
    main()
