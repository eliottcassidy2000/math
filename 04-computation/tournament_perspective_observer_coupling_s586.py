#!/usr/bin/env python3
"""S586: tournament perspectives as an observer-coupling diagnostic.

This is a focused reread of the old perspective-counting curiosity:

    P(n) = sum_T (# vertex orbits of T)

where T ranges over unlabelled n-vertex tournaments.  P(n) is the number of
pointed/rooted tournament isomorphism classes.  The small coincidences
P(3)=T(4) and P(4)=T(5) are real counts, but they should not be read as a
natural one-vertex extension bijection.  The old n=5 -> 6 "gap 8" is better
viewed as the first visible symptom that a root alone does not carry the
incident-edge/threshold payload needed by the next level.

The LRC analogy: unmarked classes are observer-blind; rooted classes are the
first observer-coupled lift; LRC needs a further threshold/endpoints/owner
fiber above the root.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from itertools import combinations, permutations


KNOWN_UNROOTED = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
}


def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {(i, j): k for k, (i, j) in enumerate(combinations(range(n), 2))}


PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}
PERMS: dict[int, list[tuple[int, ...]]] = {}


def pairs(n: int) -> dict[tuple[int, int], int]:
    if n not in PAIR_INDEX:
        PAIR_INDEX[n] = pair_index(n)
    return PAIR_INDEX[n]


def perms(n: int) -> list[tuple[int, ...]]:
    if n not in PERMS:
        PERMS[n] = list(permutations(range(n)))
    return PERMS[n]


def edge(mask: int, n: int, i: int, j: int) -> bool:
    """Return True iff i beats j in mask."""
    if i == j:
        raise ValueError("no loops in tournaments")
    if i < j:
        return bool((mask >> pairs(n)[(i, j)]) & 1)
    return not edge(mask, n, j, i)


def with_edge(mask: int, n: int, i: int, j: int, i_beats_j: bool) -> int:
    """Set the i/j edge in mask and return the new mask."""
    if i > j:
        return with_edge(mask, n, j, i, not i_beats_j)
    bit = 1 << pairs(n)[(i, j)]
    if i_beats_j:
        return mask | bit
    return mask & ~bit


def relabel(mask: int, n: int, old_for_new: tuple[int, ...]) -> int:
    """Relabel a tournament by new vertex i = old vertex old_for_new[i]."""
    out = 0
    for i, j in combinations(range(n), 2):
        if edge(mask, n, old_for_new[i], old_for_new[j]):
            out = with_edge(out, n, i, j, True)
    return out


def canonical(mask: int, n: int) -> int:
    return min(relabel(mask, n, p) for p in perms(n))


def automorphisms(mask: int, n: int) -> list[tuple[int, ...]]:
    return [p for p in perms(n) if relabel(mask, n, p) == mask]


def vertex_orbits(mask: int, n: int) -> list[tuple[int, ...]]:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for p in automorphisms(mask, n):
        for v in range(n):
            union(v, p[v])

    out: dict[int, list[int]] = defaultdict(list)
    for v in range(n):
        out[find(v)].append(v)
    return tuple(tuple(vs) for vs in sorted(out.values(), key=lambda x: (len(x), x)))


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(sum(edge(mask, n, i, j) for j in range(n) if j != i) for i in range(n)))


def count_cyclic_triples(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        if edge(mask, n, a, b) and edge(mask, n, b, c) and edge(mask, n, c, a):
            total += 1
        elif edge(mask, n, a, c) and edge(mask, n, c, b) and edge(mask, n, b, a):
            total += 1
    return total


def ham_paths(mask: int, n: int) -> int:
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
                if edge(mask, n, v, u):
                    dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[(1 << n) - 1])


def classes(n: int) -> list[int]:
    seen = set()
    reps = []
    for mask in range(1 << (n * (n - 1) // 2)):
        c = canonical(mask, n)
        if c not in seen:
            seen.add(c)
            reps.append(c)
    return sorted(reps)


def extend_with_new_vertex(mask: int, n: int, word: int) -> int:
    """Add vertex n.  Bit i of word means the new vertex beats i."""
    out = 0
    for i, j in combinations(range(n), 2):
        out = with_edge(out, n + 1, i, j, edge(mask, n, i, j))
    for i in range(n):
        new_beats_i = bool((word >> i) & 1)
        out = with_edge(out, n + 1, n, i, new_beats_i)
    return out


def perspective_table(max_n: int = 6) -> dict[int, dict[str, object]]:
    table = {}
    for n in range(1, max_n + 1):
        reps = classes(n)
        orbit_counts = [len(vertex_orbits(rep, n)) for rep in reps]
        table[n] = {
            "classes": reps,
            "T": len(reps),
            "P": sum(orbit_counts),
            "orbit_distribution": dict(sorted(Counter(orbit_counts).items())),
            "orbit_counts": orbit_counts,
        }
    return table


def rooted_types(reps: list[int], n: int) -> list[dict[str, object]]:
    rooted = []
    for class_index, rep in enumerate(reps):
        for orbit in vertex_orbits(rep, n):
            root = min(orbit)
            rooted.append(
                {
                    "class_index": class_index,
                    "mask": rep,
                    "root": root,
                    "root_score": sum(edge(rep, n, root, j) for j in range(n) if j != root),
                    "score": score_sequence(rep, n),
                    "orbit": orbit,
                }
            )
    return rooted


def extension_stats(table: dict[int, dict[str, object]], n: int) -> dict[str, object]:
    reps = table[n]["classes"]
    target_reps = set(table[n + 1]["classes"])
    rooted = rooted_types(reps, n)
    image_by_root = []
    class_to_parent_roots: dict[int, set[int]] = defaultdict(set)

    for idx, rt in enumerate(rooted):
        image = set()
        for word in range(1 << n):
            target = canonical(extend_with_new_vertex(rt["mask"], n, word), n + 1)
            image.add(target)
            class_to_parent_roots[target].add(idx)
        image_by_root.append(frozenset(image))

    root_sensitive_classes = 0
    for class_index in range(len(reps)):
        images = {
            image_by_root[idx]
            for idx, rt in enumerate(rooted)
            if rt["class_index"] == class_index
        }
        if len(images) > 1:
            root_sensitive_classes += 1

    parent_counts = [len(class_to_parent_roots[c]) for c in target_reps]
    reach_counts = [len(image) for image in image_by_root]
    covered = sum(1 for c in target_reps if class_to_parent_roots[c])

    return {
        "n": n,
        "rooted_count": len(rooted),
        "target_count": len(target_reps),
        "gap_target_minus_rooted": len(target_reps) - len(rooted),
        "covered_target_classes": covered,
        "orphan_target_classes": len(target_reps) - covered,
        "root_sensitive_parent_classes": root_sensitive_classes,
        "parent_count_distribution": dict(sorted(Counter(parent_counts).items())),
        "reach_count_distribution": dict(sorted(Counter(reach_counts).items())),
    }


def directed_three_cycles(adj: list[list[int]]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach_from(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for u, bit in enumerate(graph[v]):
                if bit and u not in seen:
                    seen.add(u)
                    q.append(u)
        return seen

    remaining = set(range(n))
    sizes = []
    reverse = [[adj[j][i] for j in range(n)] for i in range(n)]
    while remaining:
        v = min(remaining)
        comp = reach_from(v, adj) & reach_from(v, reverse)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: list[list[int]]) -> int:
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


def route_tournament() -> dict[str, object]:
    """Tournament Analysis over possible quotient choices.

    Vertices are not runners.  They are proof-state models; this is the explicit
    assumption challenge required for LRC/tournament sessions.
    """
    routes = [
        {
            "name": "unmarked_A000568_class",
            "preserves_predicate": 1,
            "compression": 5,
            "computable": 5,
            "proof_power": 1,
            "extension_payload": 0,
            "risk": 5,
        },
        {
            "name": "observer_pointed_class",
            "preserves_predicate": 2,
            "compression": 4,
            "computable": 4,
            "proof_power": 2,
            "extension_payload": 1,
            "risk": 4,
        },
        {
            "name": "root_plus_incident_word",
            "preserves_predicate": 3,
            "compression": 3,
            "computable": 4,
            "proof_power": 3,
            "extension_payload": 4,
            "risk": 3,
        },
        {
            "name": "observer_source_threshold",
            "preserves_predicate": 5,
            "compression": 3,
            "computable": 4,
            "proof_power": 5,
            "extension_payload": 3,
            "risk": 1,
        },
        {
            "name": "endpoint_owner_sheaf",
            "preserves_predicate": 5,
            "compression": 2,
            "computable": 3,
            "proof_power": 5,
            "extension_payload": 5,
            "risk": 1,
        },
        {
            "name": "proof_obligation_automaton",
            "preserves_predicate": 4,
            "compression": 4,
            "computable": 5,
            "proof_power": 4,
            "extension_payload": 4,
            "risk": 2,
        },
    ]

    def strength(route: dict[str, object]) -> int:
        return (
            4 * int(route["preserves_predicate"])
            + 2 * int(route["proof_power"])
            + 2 * int(route["extension_payload"])
            + int(route["computable"])
            + int(route["compression"])
            - 2 * int(route["risk"])
        )

    scores = [strength(route) for route in routes]
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if scores[i] > scores[j]:
            adj[i][j] = 1
        elif scores[j] > scores[i]:
            adj[j][i] = 1
        else:
            # Tie Hamiltonian path: prefer the less compressed, more labelled model.
            if routes[i]["extension_payload"] >= routes[j]["extension_payload"]:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    out_scores = [sum(row) for row in adj]
    ranking = sorted(zip(scores, [r["name"] for r in routes]), reverse=True)
    return {
        "routes": routes,
        "strength_scores": dict((r["name"], s) for r, s in zip(routes, scores)),
        "score_histogram": dict(sorted(Counter(out_scores).items())),
        "directed_3_cycles": directed_three_cycles(adj),
        "scc_sizes": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_path_count(adj),
        "ranking": ranking,
    }


def print_occurrence_map() -> None:
    print("Repo occurrence map")
    print("  S25g draft/script: first literal wording, 3+1=4 and 2+2+4+4=12.")
    print("  T075 / HTML explorer: vertex-orbit child conjecture; works 3->4, 4->5, fails 5->6.")
    print("  S4/S92t scripts: exhaustive P(n) vs T(n+1), gap 8, extension-map attempt.")
    print("  perspective_class_map.out: all 56 six-classes have parents; there are no orphans.")
    print("  THM-260 rooted layer: P(n)=n*T(n)-D(n), D indexed by odd automorphism partitions.")
    print("  S369/S370: gap 8 tied to LRC alpha stencils, chirality, 12+44 and 14+42.")
    print("  S509/HYP-1977: unmarked and observer-pointed LRC classes can mix safe/unsafe cells.")
    print("  T600/T629: quotient stack needs observer/source, threshold, gap/apex, and endpoint labels.")


def main() -> None:
    print("S586 observer-coupled tournament perspectives")
    print("=" * 72)
    print_occurrence_map()

    table = perspective_table(6)
    print()
    print("Unrooted classes U(n), rooted perspectives P(n), and shifted gap")
    print("  n  U(n)  P(n)  U(n+1)  U(n+1)-P(n)  orbit-count distribution")
    for n in range(1, 7):
        next_u = KNOWN_UNROOTED.get(n + 1)
        gap = "?" if next_u is None else str(next_u - int(table[n]["P"]))
        print(
            f"  {n:1d}  {table[n]['T']:4d}  {table[n]['P']:4d}"
            f"  {str(next_u):>6}  {gap:>12}  {table[n]['orbit_distribution']}"
        )

    print()
    print("The user's small counts")
    print("  n=3 classes contribute orbit counts [3, 1], so P(3)=4=U(4).")
    print("  n=4 classes contribute orbit counts [4, 4, 2, 2], so P(4)=12=U(5).")
    print("  n=5 contributes distribution {1:1, 3:4, 5:7}, so P(5)=48, but U(6)=56.")
    print("  n=6 contributes distribution {2:5, 4:10, 6:41}, so P(6)=296, but U(7)=456.")

    print()
    print("Extension-map sanity check")
    print("  This uses the old operation: fix an n-class, add a new vertex with all 2^n incident words.")
    print("  Crucial diagnostic: the chosen root is not used by this operation.")
    for n in (3, 4, 5):
        stats = extension_stats(table, n)
        print(f"  n={n}->n+1={n+1}:")
        print(
            f"    rooted={stats['rooted_count']}, target={stats['target_count']}, "
            f"target-rooted gap={stats['gap_target_minus_rooted']}"
        )
        print(
            f"    covered targets={stats['covered_target_classes']}, "
            f"orphans={stats['orphan_target_classes']}, "
            f"root-sensitive parent classes={stats['root_sensitive_parent_classes']}"
        )
        print(f"    parents/class distribution={stats['parent_count_distribution']}")
        print(f"    reached-classes/root distribution={stats['reach_count_distribution']}")

    print()
    print("Interpretation")
    print("  1. The small equality P(n)=U(n+1) is a rank-shift count, not a natural extension bijection.")
    print("  2. The old 'orphan six-classes' reading is false for the full extension map: all 56 are covered.")
    print("  3. The real missing datum is the incident word/cocycle that couples the observer to the extension.")
    print("  4. In LRC, that cocycle is not merely orientation; it is threshold, endpoint, owner, and pressure data.")
    print("  5. Therefore observer-blind clocks are useful caches, but exact safe-box hitting must stay in a labelled fiber.")

    print()
    print("Safe-to-ignore versus must-keep clocks")
    print("  Ignore for exact safety: raw unmarked A000568 class, raw H/score/c3 alone, and root-unused extension maps.")
    print("  Keep as first lift: observer-pointed/rooted class, because the stationary runner is not symmetric noise.")
    print("  Keep for exact LRC: observer-source threshold arcs, endpoint/gap labels, owner residues, and proof automata.")
    print("  Remodel safe-box hitting as: quotient walk + labelled fiber + accepting observer-source state.")

    ta = route_tournament()
    print()
    print("Tournament Analysis over quotient choices")
    print(f"  route strengths={ta['strength_scores']}")
    print(f"  ranking={ta['ranking']}")
    print(f"  score_histogram={ta['score_histogram']}")
    print(f"  directed_3_cycles={ta['directed_3_cycles']}")
    print(f"  scc_sizes={ta['scc_sizes']}")
    print(f"  hamiltonian_paths={ta['hamiltonian_paths']}")
    print()
    print("Assumption challenge")
    print("  Vertices here were quotient models, not runners or tournament vertices.")
    print("  Preserved predicate: existence of an LRC safe-box hit after projection.")
    print("  Destroyed information: phase duration and exact arithmetic unless the endpoint/residue fiber is retained.")


if __name__ == "__main__":
    main()
