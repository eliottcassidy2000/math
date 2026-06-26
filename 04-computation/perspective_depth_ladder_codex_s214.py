#!/usr/bin/env python3
"""Perspective-depth scout for tournament rooted viewpoints.

Session: codex-2026-06-26-S214

The script compares the classic rooted vertex-perspective count

    P(m) = sum_T |V(T)/Aut(T)|

against A000568(m+1), then tests local k-depth node perspectives and several
non-node rooted carriers: directed edges, triples, directed 3-cycles, and
disjoint 3-cycle conflict pairs.

The k-depth node color is the directed 1-WL recursion:

    color_{k+1}(v) = (color_k(v),
                      multiset colors of out-neighbors,
                      multiset colors of in-neighbors).

This is deliberately a "controlled forgetting" view: each depth states exactly
which observer-rooted data was retained before quotienting.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations


A000568 = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
}


def edge_indices(n: int) -> dict[tuple[int, int], int]:
    idx = {}
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            idx[(i, j)] = k
            k += 1
    return idx


EDGE_INDEX_CACHE: dict[int, dict[tuple[int, int], int]] = {}
PERM_CACHE: dict[int, list[tuple[int, ...]]] = {}
PERM_TRANSFORM_CACHE: dict[int, list[tuple[tuple[int, ...], tuple[tuple[int, int, bool], ...]]]] = {}


def edge_index(n: int) -> dict[tuple[int, int], int]:
    if n not in EDGE_INDEX_CACHE:
        EDGE_INDEX_CACHE[n] = edge_indices(n)
    return EDGE_INDEX_CACHE[n]


def perms(n: int) -> list[tuple[int, ...]]:
    if n not in PERM_CACHE:
        PERM_CACHE[n] = list(permutations(range(n)))
    return PERM_CACHE[n]


def perm_transforms(n: int):
    """Return (perm, edge map) records for fast bit relabeling."""
    if n not in PERM_TRANSFORM_CACHE:
        idx = edge_index(n)
        records = []
        for p in perms(n):
            maps = []
            for i in range(n):
                pi = p[i]
                for j in range(i + 1, n):
                    pj = p[j]
                    if pi < pj:
                        maps.append((idx[(i, j)], idx[(pi, pj)], False))
                    else:
                        maps.append((idx[(i, j)], idx[(pj, pi)], True))
            records.append((p, tuple(maps)))
        PERM_TRANSFORM_CACHE[n] = records
    return PERM_TRANSFORM_CACHE[n]


def has_arc(bits: int, n: int, i: int, j: int) -> bool:
    if i == j:
        raise ValueError("no loops")
    idx = edge_index(n)
    if i < j:
        return bool((bits >> idx[(i, j)]) & 1)
    return not bool((bits >> idx[(j, i)]) & 1)


def permute_bits(bits: int, n: int, perm: tuple[int, ...]) -> int:
    """Relabel old tournament by new vertex i -> old vertex perm[i]."""
    for p, maps in perm_transforms(n):
        if p == perm:
            return transform_bits(bits, maps)
    raise ValueError("unknown permutation")


def transform_bits(bits: int, maps: tuple[tuple[int, int, bool], ...]) -> int:
    out = 0
    for new_idx, old_idx, invert in maps:
        bit = (bits >> old_idx) & 1
        if invert:
            bit ^= 1
        if bit:
            out |= 1 << new_idx
    return out


def canonical_bits(bits: int, n: int) -> int:
    return min(transform_bits(bits, maps) for _, maps in perm_transforms(n))


def reps_upto(nmax: int) -> dict[int, list[int]]:
    reps = {}
    for n in range(1, nmax + 1):
        m = n * (n - 1) // 2
        classes = {}
        for bits in range(1 << m):
            c = canonical_bits(bits, n)
            classes.setdefault(c, c)
        reps[n] = sorted(classes)
    return reps


def automorphisms(bits: int, n: int) -> list[tuple[int, ...]]:
    return [p for p, maps in perm_transforms(n) if transform_bits(bits, maps) == bits]


class DSU:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb

    def count(self) -> int:
        return len({self.find(x) for x in self.parent})


def orbit_count(items, auts, action):
    dsu = DSU(items)
    for p in auts:
        for item in items:
            dsu.union(item, action(item, p))
    return dsu.count()


def vertex_orbit_count(bits: int, n: int, auts) -> int:
    return orbit_count(tuple(range(n)), auts, lambda v, p: p.index(v))


def arc_items(bits: int, n: int):
    return tuple((i, j) for i in range(n) for j in range(n) if i != j and has_arc(bits, n, i, j))


def arc_orbit_count(bits: int, n: int, auts) -> int:
    items = arc_items(bits, n)
    return orbit_count(items, auts, lambda e, p: (p.index(e[0]), p.index(e[1])))


def triple_kind(bits: int, n: int, tri: tuple[int, int, int]) -> str:
    a, b, c = tri
    # A triple is cyclic iff the three orientations form a directed cycle.
    cyc = (
        has_arc(bits, n, a, b) and has_arc(bits, n, b, c) and has_arc(bits, n, c, a)
    ) or (
        has_arc(bits, n, a, c) and has_arc(bits, n, c, b) and has_arc(bits, n, b, a)
    )
    return "cyclic" if cyc else "transitive"


def triple_orbit_count(bits: int, n: int, auts, kind: str | None = None) -> int:
    items = tuple(
        tri
        for tri in combinations(range(n), 3)
        if kind is None or triple_kind(bits, n, tri) == kind
    )
    if not items:
        return 0
    return orbit_count(items, auts, lambda tri, p: tuple(sorted(p.index(v) for v in tri)))


def conflict_pair_orbit_count(bits: int, n: int, auts) -> int:
    cycles = [
        tri
        for tri in combinations(range(n), 3)
        if triple_kind(bits, n, tri) == "cyclic"
    ]
    pairs = []
    for i, a in enumerate(cycles):
        sa = set(a)
        for b in cycles[i + 1 :]:
            if sa.isdisjoint(b):
                pairs.append((a, b))
    items = tuple(pairs)
    if not items:
        return 0

    def act(pair, p):
        moved = tuple(tuple(sorted(p.index(v) for v in tri)) for tri in pair)
        return tuple(sorted(moved))

    return orbit_count(items, auts, act)


def relabel_colors(signatures):
    mapping = {}
    colors = []
    for sig in signatures:
        if sig not in mapping:
            mapping[sig] = len(mapping)
        colors.append(mapping[sig])
    return tuple(colors)


def node_depth_colors(bits: int, n: int, depth: int) -> tuple[int, ...]:
    colors = tuple(0 for _ in range(n))
    for _ in range(depth):
        sigs = []
        for v in range(n):
            out = sorted(colors[w] for w in range(n) if w != v and has_arc(bits, n, v, w))
            inn = sorted(colors[w] for w in range(n) if w != v and has_arc(bits, n, w, v))
            sigs.append((colors[v], tuple(out), tuple(inn)))
        colors = relabel_colors(sigs)
    return colors


def node_depth_count(bits: int, n: int, depth: int) -> int:
    return len(set(node_depth_colors(bits, n, depth)))


def edge_depth_count(bits: int, n: int, depth: int) -> int:
    node_colors = node_depth_colors(bits, n, max(depth - 1, 0))
    sigs = []
    for u, v in arc_items(bits, n):
        external = []
        for w in range(n):
            if w in (u, v):
                continue
            external.append(
                (
                    1 if has_arc(bits, n, u, w) else 0,
                    1 if has_arc(bits, n, v, w) else 0,
                    node_colors[w],
                )
            )
        sigs.append((node_colors[u], node_colors[v], tuple(sorted(external))))
    return len(set(sigs))


def triple_depth_signature(bits: int, n: int, tri: tuple[int, int, int], depth: int):
    node_colors = node_depth_colors(bits, n, max(depth - 1, 0))
    options = []
    # Canonicalize over all internal labelings of the rooted 3-subset.  This
    # makes the local signature invariant under tournament automorphisms; using
    # numeric tuple order would accidentally over-split transitive triples.
    for rooted in permutations(tri):
        internal = tuple(
            1 if has_arc(bits, n, rooted[i], rooted[j]) else 0
            for i in range(3)
            for j in range(i + 1, 3)
        )
        external = []
        for w in range(n):
            if w in tri:
                continue
            rel = tuple(1 if has_arc(bits, n, r, w) else 0 for r in rooted)
            external.append((rel, node_colors[w]))
        options.append((triple_kind(bits, n, tri), internal, tuple(sorted(external))))
    return min(options)


def triple_depth_count(bits: int, n: int, depth: int, kind: str | None = None) -> int:
    sigs = []
    for tri in combinations(range(n), 3):
        if kind is not None and triple_kind(bits, n, tri) != kind:
            continue
        sigs.append(triple_depth_signature(bits, n, tri, depth))
    return len(set(sigs))


def summarize(nmax: int = 5, depth_max: int = 4) -> str:
    reps = reps_upto(nmax)
    lines = []
    lines.append("Perspective-depth ladder scout")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Shift check: A000568(m+1) versus exact rooted node perspectives P(m)")
    lines.append("m classes(m) P(m) A000568(m+1) defect")

    exact_vertex_totals = {}
    aut_cache = {}
    for n in range(1, nmax + 1):
        total = 0
        for bits in reps[n]:
            auts = automorphisms(bits, n)
            aut_cache[(n, bits)] = auts
            total += vertex_orbit_count(bits, n, auts)
        exact_vertex_totals[n] = total
        target = A000568.get(n + 1)
        defect = None if target is None else target - total
        lines.append(f"{n} {len(reps[n])} {total} {target} {defect}")

    lines.append("")
    lines.append("Node k-depth totals (directed WL perspective colors, summed over iso classes)")
    header = "m exact " + " ".join(f"d{k}" for k in range(1, depth_max + 1))
    lines.append(header)
    for n in range(1, nmax + 1):
        row = [str(n), str(exact_vertex_totals[n])]
        for depth in range(1, depth_max + 1):
            row.append(str(sum(node_depth_count(bits, n, depth) for bits in reps[n])))
        lines.append(" ".join(row))

    lines.append("")
    lines.append("Exact non-node rooted carriers, summed over iso classes")
    lines.append("m arc_orbits triple_orbits transitive_triple_orbits cyclic_triple_orbits cycle_conflict_pair_orbits")
    for n in range(1, nmax + 1):
        arc = triple = trans = cyc = conflict = 0
        for bits in reps[n]:
            auts = aut_cache[(n, bits)]
            arc += arc_orbit_count(bits, n, auts)
            triple += triple_orbit_count(bits, n, auts)
            trans += triple_orbit_count(bits, n, auts, "transitive")
            cyc += triple_orbit_count(bits, n, auts, "cyclic")
            conflict += conflict_pair_orbit_count(bits, n, auts)
        lines.append(f"{n} {arc} {triple} {trans} {cyc} {conflict}")

    lines.append("")
    lines.append("Local depth totals for edge and triple carriers")
    lines.append("m edge_d1 edge_d2 edge_d3 all_triple_d1 all_triple_d2 cyclic_d1 cyclic_d2")
    for n in range(1, nmax + 1):
        vals = defaultdict(int)
        for bits in reps[n]:
            vals["edge_d1"] += edge_depth_count(bits, n, 1)
            vals["edge_d2"] += edge_depth_count(bits, n, 2)
            vals["edge_d3"] += edge_depth_count(bits, n, 3)
            vals["all_triple_d1"] += triple_depth_count(bits, n, 1)
            vals["all_triple_d2"] += triple_depth_count(bits, n, 2)
            vals["cyclic_d1"] += triple_depth_count(bits, n, 1, "cyclic")
            vals["cyclic_d2"] += triple_depth_count(bits, n, 2, "cyclic")
        lines.append(
            f"{n} {vals['edge_d1']} {vals['edge_d2']} {vals['edge_d3']} "
            f"{vals['all_triple_d1']} {vals['all_triple_d2']} "
            f"{vals['cyclic_d1']} {vals['cyclic_d2']}"
        )

    lines.append("")
    lines.append("Key readout")
    lines.append("- The classic equality P(m)=A000568(m+1) holds through m=4 and first fails at m=5.")
    lines.append("- In LRC indexing, this is the first next-level value n=6: A000568(6)=56 while P(5)=48.")
    lines.append("- Under this stricter directed-WL convention, node depth reaches exact rooted node orbits by depth 3 at m=5.")
    lines.append("- THM-260 gives the next exact rooted value P(6)=296 versus A000568(7)=456; extending the non-node carrier table to m=6 needs a faster class generator.")
    lines.append("- Edge, triple, cycle, and conflict perspectives grow as different marked carriers, not as replacements for A000568.")
    lines.append("- Controlled forgetting suggests using these carriers only after declaring which rooted payload is preserved.")
    return "\n".join(lines)


def main() -> None:
    print(summarize())


if __name__ == "__main__":
    main()
