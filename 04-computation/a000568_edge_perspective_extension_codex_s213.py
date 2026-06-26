#!/usr/bin/env python3
"""S213: edge-perspective extension lift at the first A000568 failure.

S211 showed that the shifted equality

    rooted node perspectives on m vertices = U(m+1)

first fails at m=5: P(5)=48 while U(6)=56.  This script tests the next
controlled-forgetting step.  A rooted 5-tournament plus the incident word of a
new observer is exactly a 6-tournament with an ordered pair of marked vertices:
the old root and the new observer.  Forgetting old/new roles turns that into a
directed-edge perspective.

The useful new finite fact is that the ordered-pair sector deck at n=6 is a
compact sidecar.  Sector sizes and internal sector tournaments miss exactly
one converse pair, while adding cross-sector orientations separates all 56
six-vertex tournament classes.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations

import tournament_rigidity_cascade_s589 as tour


def extend_with_new_vertex(mask: int, n: int, word: int) -> int:
    """Add vertex n.  Bit i of word means the new vertex beats old vertex i."""
    out = 0
    for i, j in combinations(range(n), 2):
        out = tour.set_edge(out, n + 1, i, j, tour.edge(mask, n, i, j))
    for i in range(n):
        out = tour.set_edge(out, n + 1, n, i, bool((word >> i) & 1))
    return out


def generated_classes6() -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                tour.canonical(extend_with_new_vertex(rep, 5, word), 6)
                for rep in tour.classes(5)
                for word in range(1 << 5)
            }
        )
    )


def ordered_pair_canonical(mask: int, n: int, first: int, second: int) -> int:
    others = tuple(v for v in range(n) if v not in (first, second))
    return min(tour.relabel(mask, n, (first, second) + p) for p in permutations(others))


def unordered_pair_canonical(mask: int, n: int, a: int, b: int) -> int:
    return min(ordered_pair_canonical(mask, n, a, b), ordered_pair_canonical(mask, n, b, a))


def directed_edge_canonical(mask: int, n: int, tail: int, tip: int) -> int:
    return ordered_pair_canonical(mask, n, tail, tip)


def oriented_edges(mask: int, n: int) -> list[tuple[int, int]]:
    out = []
    for a, b in combinations(range(n), 2):
        out.append((a, b) if tour.edge(mask, n, a, b) else (b, a))
    return out


def rooted_parent_incident_state(parent: int, root: int, word: int) -> int:
    """Canonical 6-state with old root first and new observer second."""
    ext = extend_with_new_vertex(parent, 5, word)
    old_vertices = tuple(v for v in range(5) if v != root)
    return min(tour.relabel(ext, 6, (root, 5) + p) for p in permutations(old_vertices))


SECTORS = ((0, 0), (0, 1), (1, 0), (1, 1))


def induced_canonical(mask: int, n: int, vertices: list[int]) -> int:
    if len(vertices) <= 1:
        return 0
    return tour.canonical(tour.induced(mask, n, tuple(vertices)), len(vertices))


def sector_signature(mask: int, n: int, first: int, second: int, mode: str) -> tuple[object, ...]:
    """Four-sector signature around an ordered pair.

    The sector key is `(first beats x, second beats x)`.
    """
    groups: dict[tuple[int, int], list[int]] = defaultdict(list)
    for x in range(n):
        if x in (first, second):
            continue
        key = (int(tour.edge(mask, n, first, x)), int(tour.edge(mask, n, second, x)))
        groups[key].append(x)

    sizes = tuple((sector, len(groups[sector])) for sector in SECTORS)
    internal = tuple(
        (sector, induced_canonical(mask, n, groups[sector])) for sector in SECTORS
    )
    cross = []
    for i, a in enumerate(SECTORS):
        for b in SECTORS[i + 1 :]:
            wins = sum(1 for u in groups[a] for v in groups[b] if tour.edge(mask, n, u, v))
            cross.append((a, b, wins, len(groups[a]) * len(groups[b])))
    cross_tuple = tuple(cross)

    if mode == "size":
        return sizes
    if mode == "internal":
        return sizes, internal
    if mode == "cross":
        return sizes, cross_tuple
    if mode == "full":
        return sizes, internal, cross_tuple
    raise ValueError(mode)


def class_sector_deck(mask: int, n: int, mode: str) -> tuple[tuple[object, ...], ...]:
    return tuple(
        sorted(
            sector_signature(mask, n, a, b, mode)
            for a in range(n)
            for b in range(n)
            if a != b
        )
    )


def opposite(mask: int, n: int) -> int:
    out = 0
    for i, j in combinations(range(n), 2):
        out = tour.set_edge(out, n, i, j, not tour.edge(mask, n, i, j))
    return out


def cyclic_triples(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        total += int(
            (tour.edge(mask, n, a, b) and tour.edge(mask, n, b, c) and tour.edge(mask, n, c, a))
            or (
                tour.edge(mask, n, a, c)
                and tour.edge(mask, n, c, b)
                and tour.edge(mask, n, b, a)
            )
        )
    return total


def hamiltonian_paths(mask: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if seen & (1 << u):
                    continue
                if tour.edge(mask, n, v, u):
                    dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def directed_three_cycles_adj(adj: list[list[int]]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        total += int(adj[a][b] and adj[b][c] and adj[c][a])
        total += int(adj[a][c] and adj[c][b] and adj[b][a])
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u, bit in enumerate(graph[v]):
                if bit and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    reverse = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    sizes = []
    while remaining:
        v = min(remaining)
        comp = reach(v, adj) & reach(v, reverse)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths_adj(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if seen & (1 << u):
                    continue
                if adj[v][u]:
                    dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def carrier_tournament() -> tuple[list[tuple[str, dict[str, int]]], list[list[int]]]:
    """Tournament Analysis over forgetful carriers.

    Vertices are sidecar/proof carriers, not tournament runners.  The pairwise
    observable rewards retained observer role, incident word, edge-sector
    cross-coupling, chirality, class separation, and LRC owner compatibility;
    proof cost is the tie gauge.
    """
    carriers = [
        ("rooted_5_node_cache", dict(role=1, incident=0, cross=0, chirality=0, separates=0, owner=0, cost=1)),
        ("ordered_pair_extension", dict(role=2, incident=2, cross=1, chirality=0, separates=1, owner=0, cost=3)),
        ("directed_edge_perspective", dict(role=1, incident=2, cross=1, chirality=0, separates=1, owner=0, cost=2)),
        ("sector_size_deck", dict(role=1, incident=1, cross=0, chirality=0, separates=0, owner=0, cost=1)),
        ("sector_internal_deck", dict(role=1, incident=1, cross=0, chirality=0, separates=0, owner=0, cost=2)),
        ("sector_cross_deck", dict(role=1, incident=2, cross=3, chirality=1, separates=2, owner=0, cost=3)),
        ("sector_full_deck", dict(role=1, incident=2, cross=3, chirality=1, separates=2, owner=0, cost=4)),
        ("endpoint_owner_packet", dict(role=2, incident=2, cross=3, chirality=1, separates=2, owner=3, cost=5)),
    ]

    def score(data: dict[str, int]) -> tuple[int, int]:
        retained = (
            3 * data["role"]
            + 3 * data["incident"]
            + 3 * data["cross"]
            + 2 * data["chirality"]
            + 2 * data["separates"]
            + 4 * data["owner"]
        )
        return retained, -data["cost"]

    n = len(carriers)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        si = score(carriers[i][1])
        sj = score(carriers[j][1])
        if si > sj or (si == sj and carriers[i][0] < carriers[j][0]):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return carriers, adj


def main() -> None:
    classes5 = tour.classes(5)
    classes6 = generated_classes6()
    rooted5 = {
        tour.rooted_canonical(rep, 5, v)
        for rep in classes5
        for v in range(5)
    }

    ordered6 = {
        ordered_pair_canonical(rep, 6, a, b)
        for rep in classes6
        for a in range(6)
        for b in range(6)
        if a != b
    }
    unordered6 = {
        unordered_pair_canonical(rep, 6, a, b)
        for rep in classes6
        for a, b in combinations(range(6), 2)
    }
    directed6 = {
        directed_edge_canonical(rep, 6, tail, tip)
        for rep in classes6
        for tail, tip in oriented_edges(rep, 6)
    }

    extension_states = set()
    states_by_root: dict[int, set[int]] = defaultdict(set)
    targets_by_root: dict[int, set[int]] = defaultdict(set)
    roots_by_target: dict[int, set[int]] = defaultdict(set)
    for rep in classes5:
        for orbit in tour.vertex_orbits(rep, 5):
            root = orbit[0]
            root_key = tour.rooted_canonical(rep, 5, root)
            for word in range(1 << 5):
                state = rooted_parent_incident_state(rep, root, word)
                target = tour.canonical(extend_with_new_vertex(rep, 5, word), 6)
                extension_states.add(state)
                states_by_root[root_key].add(state)
                targets_by_root[root_key].add(target)
                roots_by_target[target].add(root_key)

    print("=" * 80)
    print("S213: A000568 edge-perspective extension lift")
    print("=" * 80)
    print()
    print("1. FIRST FAILURE AND EXACT EXTENSION LIFT")
    print(f"   U(5)={len(classes5)}  P(5)={len(rooted5)}  U(6)={len(classes6)}  gap={len(classes6)-len(rooted5)}")
    print(f"   U(6) generated from U(5)+all incident words: {len(classes6)}")
    print(f"   rooted 5-perspectives plus incident word states: {len(extension_states)}")
    print(f"   exact ordered-pair perspectives on U(6): {len(ordered6)}")
    print(f"   extension states equal ordered-pair perspectives: {extension_states == ordered6}")
    print(f"   exact directed-edge perspectives on U(6): {len(directed6)}")
    print(f"   exact unordered-pair perspectives on U(6): {len(unordered6)}")
    print(f"   directed-edge equals unordered-pair count: {len(directed6) == len(unordered6)}")
    print(
        "   per rooted 5-parent extension-state count hist: "
        f"{dict(sorted(Counter(len(v) for v in states_by_root.values()).items()))}"
    )
    print(
        "   target classes reached per rooted 5-parent hist: "
        f"{dict(sorted(Counter(len(v) for v in targets_by_root.values()).items()))}"
    )
    print(
        "   rooted 5-parents per target U(6) class hist: "
        f"{dict(sorted(Counter(len(v) for v in roots_by_target.values()).items()))}"
    )
    print()
    print("2. ORDERED-PAIR SECTOR DECKS")
    print("   Sector key for x is (old-root beats x, new-observer beats x).")
    print("   A deck is the multiset over all ordered pairs of a 6-tournament class.")
    for mode in ("size", "internal", "cross", "full"):
        all_sigs = {
            sector_signature(rep, 6, a, b, mode)
            for rep in classes6
            for a in range(6)
            for b in range(6)
            if a != b
        }
        decks = {class_sector_deck(rep, 6, mode) for rep in classes6}
        print(
            f"   {mode:8s}: individual_sigs={len(all_sigs):3d}  "
            f"unique_class_decks={len(decks):2d}/56"
        )

    internal_decks: dict[tuple[tuple[object, ...], ...], list[int]] = defaultdict(list)
    for rep in classes6:
        internal_decks[class_sector_deck(rep, 6, "internal")].append(rep)
    collisions = [reps for reps in internal_decks.values() if len(reps) > 1]
    print()
    print("3. THE ONLY SIZE/INTERNAL COLLISION")
    for reps in collisions:
        print(f"   collided masks={reps}; converse pair={tour.canonical(opposite(reps[0], 6), 6) == reps[1]}")
        for rep in reps:
            print(
                f"     rep={rep} score={tour.score_sequence(rep, 6)} "
                f"c3={cyclic_triples(rep, 6)} H={hamiltonian_paths(rep, 6)} "
                f"aut={len(tour.automorphisms(rep, 6))} "
                f"self_converse={tour.canonical(opposite(rep, 6), 6) == rep}"
            )
    print("   Cross-sector orientation separates this converse/chirality collision.")

    carriers, adj = carrier_tournament()
    out_scores = [sum(row) for row in adj]
    order = sorted(range(len(carriers)), key=lambda i: -out_scores[i])
    print()
    print("4. TOURNAMENT ANALYSIS OVER FORGETFUL CARRIERS")
    print("   Pairwise observable: retained observer role, incident word, cross-sector")
    print("   coupling, chirality, class separation, owner compatibility; cost breaks ties.")
    print(f"   vertices={[name for name, _ in carriers]}")
    print(f"   score_hist={dict(sorted(Counter(out_scores).items()))}")
    print(f"   directed_3_cycles={directed_three_cycles_adj(adj)}")
    print(f"   scc_sizes={scc_sizes(adj)}")
    print(f"   hamiltonian_paths={hamiltonian_paths_adj(adj)}")
    print("   one_hamiltonian_path=" + " -> ".join(carriers[i][0] for i in order))

    print()
    print("READING")
    print(
        "  A rooted node perspective is one controlled-forgetting rung.  Adding the "
        "observer's incident word lifts it exactly to an ordered-pair perspective "
        "on the next tournament size."
    )
    print(
        "  Directed-edge perspective is the old/new-role quotient of that ordered "
        "pair.  It is natural and dualistic, but it already forgets which endpoint "
        "was the observer unless a sidecar says so."
    )
    print(
        "  The four-sector deck is the first compact edge sidecar.  Sector sizes and "
        "internal sector types almost work at U(6), but they miss one converse pair; "
        "cross-sector orientation is exactly the chirality/coupling repair."
    )
    print(
        "  For LRC, this says the safe quotient is not raw A000568 -> rooted node.  "
        "It is rooted node -> incident word -> edge-sector cross deck -> endpoint-owner "
        "packet, with each forgotten coordinate explicitly accounted for."
    )
    print()
    print("ASSUMPTION CHALLENGE")
    print(
        "  Considered vertices: rooted nodes, ordered pairs, directed edges, sector "
        "cells, converse/chirality pairs, endpoint owners, packet sidecars, and proof "
        "obligations.  The selected Tournament Analysis vertices are carriers."
    )
    print(
        "  Preserved predicate: observer-source/incident-coupling data needed for a "
        "safe-box hit.  Destroyed information: old/new endpoint role, labelled runner "
        "identity, and full extension rows; those are acceptable only with a sidecar "
        "or named residual debt."
    )


if __name__ == "__main__":
    main()
