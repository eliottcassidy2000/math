#!/usr/bin/env python3
"""S589: local fixed-point rigidity and global cascade tests for tournaments.

This script is intentionally small and exact.  It treats rooted tournament
classes as a toy model for the LRC observer/basepoint problem:

* a source root is a rigid fixed point: delete it and the parent is canonical;
* an arbitrary root is only a perspective: deleting it or recording only the
  two side subtournaments forgets cross-coupling data;
* global deletion decks already collide through n <= 6, so unmarked global
  shadows are weaker than the LRC source/observer predicate.

Tournament Analysis is included with quotient lenses as vertices.  The pairwise
observable is "how much rooted information does this lens preserve?"  The
switch orients toward lower collision count, then lower max fiber, then higher
unique count; name order is only a deterministic tie breaker.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import combinations, permutations


PAIR_INDEX: dict[int, dict[tuple[int, int], int]] = {}
PERMS: dict[int, tuple[tuple[int, ...], ...]] = {}


def pairs(n: int) -> dict[tuple[int, int], int]:
    if n not in PAIR_INDEX:
        PAIR_INDEX[n] = {
            (i, j): k for k, (i, j) in enumerate(combinations(range(n), 2))
        }
    return PAIR_INDEX[n]


def perms(n: int) -> tuple[tuple[int, ...], ...]:
    if n not in PERMS:
        PERMS[n] = tuple(permutations(range(n)))
    return PERMS[n]


def edge(mask: int, n: int, i: int, j: int) -> bool:
    """Return True iff i beats j."""
    if i == j:
        raise ValueError("tournaments have no loops")
    if i < j:
        return bool((mask >> pairs(n)[(i, j)]) & 1)
    return not edge(mask, n, j, i)


def set_edge(mask: int, n: int, i: int, j: int, i_beats_j: bool) -> int:
    if i > j:
        return set_edge(mask, n, j, i, not i_beats_j)
    bit = 1 << pairs(n)[(i, j)]
    return (mask | bit) if i_beats_j else (mask & ~bit)


def relabel(mask: int, n: int, old_for_new: tuple[int, ...]) -> int:
    out = 0
    for i, j in combinations(range(n), 2):
        if edge(mask, n, old_for_new[i], old_for_new[j]):
            out = set_edge(out, n, i, j, True)
    return out


@lru_cache(maxsize=None)
def canonical(mask: int, n: int) -> int:
    return min(relabel(mask, n, p) for p in perms(n))


@lru_cache(maxsize=None)
def rooted_canonical(mask: int, n: int, root: int) -> int:
    others = tuple(v for v in range(n) if v != root)
    return min(relabel(mask, n, (root,) + p) for p in permutations(others))


@lru_cache(maxsize=None)
def automorphisms(mask: int, n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(p for p in perms(n) if relabel(mask, n, p) == mask)


def vertex_orbits(mask: int, n: int) -> tuple[tuple[int, ...], ...]:
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

    groups: dict[int, list[int]] = defaultdict(list)
    for v in range(n):
        groups[find(v)].append(v)
    return tuple(tuple(vs) for vs in sorted(groups.values(), key=lambda x: (len(x), x)))


def score_sequence(mask: int, n: int) -> tuple[int, ...]:
    return tuple(
        sorted(sum(edge(mask, n, i, j) for j in range(n) if j != i) for i in range(n))
    )


def induced(mask: int, n: int, verts: tuple[int, ...]) -> int:
    m = len(verts)
    out = 0
    for i, j in combinations(range(m), 2):
        if edge(mask, n, verts[i], verts[j]):
            out = set_edge(out, m, i, j, True)
    return out


def delete_vertex(mask: int, n: int, victim: int) -> int:
    verts = tuple(v for v in range(n) if v != victim)
    return induced(mask, n, verts)


@lru_cache(maxsize=None)
def classes(n: int) -> tuple[int, ...]:
    seen: set[int] = set()
    reps: list[int] = []
    for mask in range(1 << (n * (n - 1) // 2)):
        c = canonical(mask, n)
        if c not in seen:
            seen.add(c)
            reps.append(c)
    return tuple(sorted(reps))


def rooted_types(n: int) -> list[dict[str, object]]:
    out = []
    for rep in classes(n):
        for orbit in vertex_orbits(rep, n):
            root = min(orbit)
            out_neighbors = tuple(j for j in range(n) if j != root and edge(rep, n, root, j))
            in_neighbors = tuple(j for j in range(n) if j != root and edge(rep, n, j, root))
            out_class = canonical(induced(rep, n, out_neighbors), len(out_neighbors))
            in_class = canonical(induced(rep, n, in_neighbors), len(in_neighbors))
            parent = canonical(delete_vertex(rep, n, root), n - 1) if n > 1 else 0
            score = len(out_neighbors)
            out.append(
                {
                    "unrooted": rep,
                    "root": root,
                    "orbit": orbit,
                    "orbit_size": len(orbit),
                    "rooted": rooted_canonical(rep, n, root),
                    "score": score,
                    "source": score == n - 1,
                    "sink": score == 0,
                    "parent": parent,
                    "split_profile": (score, out_class, in_class),
                    "unrooted_score": (rep, score),
                    "score_sequence": score_sequence(rep, n),
                }
            )
    return out


def collision_stats(values: list[object]) -> dict[str, object]:
    fibers = Counter(values)
    ambiguous = [v for v, c in fibers.items() if c > 1]
    return {
        "unique": len(fibers),
        "collisions": len(values) - len(fibers),
        "ambiguous_values": len(ambiguous),
        "max_fiber": max(fibers.values()) if fibers else 0,
        "fiber_hist": dict(sorted(Counter(fibers.values()).items())),
    }


def deletion_deck(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(canonical(delete_vertex(mask, n, v), n - 1) for v in range(n)))


def deck_collision_stats(n: int) -> dict[str, object]:
    deck_to_classes: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for rep in classes(n):
        deck_to_classes[deletion_deck(rep, n)].append(rep)
    bad = [v for v in deck_to_classes.values() if len(v) > 1]
    return {
        "classes": len(classes(n)),
        "unique_decks": len(deck_to_classes),
        "collision_decks": len(bad),
        "collided_classes": sum(len(v) for v in bad),
        "max_fiber": max(len(v) for v in deck_to_classes.values()),
    }


def source_cascade_stats(n: int) -> dict[str, object]:
    rtypes = rooted_types(n)
    source_roots = [r for r in rtypes if r["source"]]
    parent_fibers = Counter(r["parent"] for r in source_roots)
    all_one = all(c == 1 for c in parent_fibers.values())
    source_fixed = all(r["orbit_size"] == 1 for r in source_roots)
    return {
        "n": n,
        "source_rooted": len(source_roots),
        "parent_classes": len(classes(n - 1)) if n > 1 else 1,
        "parent_values": len(parent_fibers),
        "delete_parent_collisions": len(source_roots) - len(parent_fibers),
        "all_parent_fibers_one": all_one,
        "all_sources_fixed_by_aut": source_fixed,
    }


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for u, yes in enumerate(graph[v]):
                if yes and u not in seen:
                    seen.add(u)
                    q.append(u)
        return seen

    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    out = []
    while remaining:
        v = next(iter(remaining))
        comp = reach(v, adj) & reach(v, radj)
        out.append(len(comp))
        remaining -= comp
    return sorted(out)


def directed_three_cycles(adj: list[list[int]]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def ham_paths_adj(adj: list[list[int]]) -> int:
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


def lens_tournament(n: int) -> tuple[list[tuple[str, dict[str, object]]], list[list[int]]]:
    rtypes = rooted_types(n)
    lens_values = {
        "full_rooted_class": [r["rooted"] for r in rtypes],
        "split_profile_no_cross": [r["split_profile"] for r in rtypes],
        "unrooted_plus_root_score": [r["unrooted_score"] for r in rtypes],
        "underlying_unrooted_class": [r["unrooted"] for r in rtypes],
        "delete_root_parent_only": [r["parent"] for r in rtypes],
        "score_sequence_only": [r["score_sequence"] for r in rtypes],
        "root_score_only": [r["score"] for r in rtypes],
    }
    stats = [(name, collision_stats(vals)) for name, vals in lens_values.items()]

    def rank(item: tuple[str, dict[str, object]]) -> tuple[int, int, int, str]:
        name, s = item
        return (
            int(s["collisions"]),
            int(s["max_fiber"]),
            -int(s["unique"]),
            name,
        )

    m = len(stats)
    adj = [[0] * m for _ in range(m)]
    for i, j in combinations(range(m), 2):
        if rank(stats[i]) < rank(stats[j]):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return stats, adj


def main() -> None:
    print("=" * 80)
    print("S589: TOURNAMENT RIGIDITY CASCADE")
    print("=" * 80)

    print("\n1. Source roots are exact local fixed points")
    for n in range(2, 7):
        s = source_cascade_stats(n)
        print(
            f"  n={n}: source_rooted={s['source_rooted']} "
            f"parent_classes=U({n-1})={s['parent_classes']} "
            f"parent_values={s['parent_values']} "
            f"delete_collisions={s['delete_parent_collisions']} "
            f"fibers_one={s['all_parent_fibers_one']} "
            f"sources_fixed={s['all_sources_fixed_by_aut']}"
        )

    print("\n2. Arbitrary rooted views have defect fibers")
    for n in range(3, 7):
        rtypes = rooted_types(n)
        full = len(rtypes)
        split = collision_stats([r["split_profile"] for r in rtypes])
        parent = collision_stats([r["parent"] for r in rtypes])
        unrooted = collision_stats([r["unrooted"] for r in rtypes])
        fixed = Counter(r["orbit_size"] == 1 for r in rtypes)
        print(
            f"  n={n}: rooted={full}, fixed_orbit_roots={fixed[True]}, "
            f"split_unique={split['unique']} split_collisions={split['collisions']} "
            f"split_max_fiber={split['max_fiber']}, "
            f"delete_parent_unique={parent['unique']} parent_max_fiber={parent['max_fiber']}, "
            f"unrooted_unique={unrooted['unique']} unrooted_max_fiber={unrooted['max_fiber']}"
        )

    print("\n3. Global deletion decks already have defect fibers")
    for n in range(3, 7):
        d = deck_collision_stats(n)
        print(
            f"  n={n}: classes={d['classes']} unique_decks={d['unique_decks']} "
            f"collision_decks={d['collision_decks']} collided_classes={d['collided_classes']} "
            f"max_fiber={d['max_fiber']}"
        )

    print("\n4. Lens collision table on rooted 6-tournaments")
    lens_stats, adj = lens_tournament(6)
    for name, s in sorted(
        lens_stats,
        key=lambda item: (
            int(item[1]["collisions"]),
            int(item[1]["max_fiber"]),
            -int(item[1]["unique"]),
            item[0],
        ),
    ):
        print(
            f"  {name:26s} unique={s['unique']:3d} "
            f"collisions={s['collisions']:3d} "
            f"ambiguous_values={s['ambiguous_values']:3d} "
            f"max_fiber={s['max_fiber']:3d} "
            f"fiber_hist={s['fiber_hist']}"
        )

    print("\n5. Tournament Analysis on quotient lenses")
    scores = [sum(row) for row in adj]
    print(f"  vertices={[name for name, _ in lens_stats]}")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3_cycles={directed_three_cycles(adj)}")
    print(f"  scc_sizes={scc_sizes(adj)}")
    print(f"  hamiltonian_paths={ham_paths_adj(adj)}")
    order = sorted(range(len(lens_stats)), key=lambda i: -scores[i])
    print("  one_hamiltonian_order=" + " -> ".join(lens_stats[i][0] for i in order))

    print("\nREADING")
    print(
        "  A source root is locally rigid: it is fixed by every automorphism, "
        "deletes to exactly one parent class, and every parent re-adds exactly "
        "one source-rooted child.  This is the tournament version of a safe "
        "observer/source certificate."
    )
    print(
        "  Arbitrary roots are not rigid in that sense.  Split profiles that "
        "remember only the two side subtournaments already collide by n=4 and "
        "badly by n=6; the forgotten data is precisely cross-coupling.  This is "
        "the same warning as the LRC observer-coupling defect."
    )
    print(
        "  Global deletion decks already collide in this small tournament range. "
        "This strengthens the warning: LRC needs a marked local predicate, so the "
        "usable rigidity is source/fold/endpoint rigidity that cascades without "
        "dropping observer labels, not unmarked deck reconstruction."
    )


if __name__ == "__main__":
    main()
