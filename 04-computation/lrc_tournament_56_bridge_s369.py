#!/usr/bin/env python3
"""
lrc_tournament_56_bridge_s369.py

codex-2026-05-31 S369

Investigate the user's 56 numerology:

    S367 LRC missed cells = 56
    unlabeled tournaments on 6 vertices = 56
    56 - 12 = 44, where tournaments on 5 vertices = 12
    56 - 14 = 42, where fourteen runners are the active frontier

This script does two exact things:

1. Enumerate tournament isomorphism classes on 5 and 6 vertices by orbit
   generation, then measure the self-converse/chiral split at n=6.
2. Re-read the S367 56 missed cells and split their 8 alpha stencils into
   mirror pairs and the 14+42 outer/interior decomposition.

The goal is not to claim a bijection.  It is to isolate which decompositions
are structural enough to be worth a proof search.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
import importlib.util
import sys


@dataclass(frozen=True)
class TournamentClass:
    rep: int
    orbit_size: int
    aut_size: int
    self_converse: bool
    strongly_connected: bool
    vertex_orbit_count: int
    score_sequence: tuple[int, ...]


def edge_pairs(n: int) -> list[tuple[int, int]]:
    return list(combinations(range(n), 2))


def edge_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: idx for idx, pair in enumerate(edge_pairs(n))}


def permutation_maps(n: int) -> list[tuple[tuple[int, ...], tuple[tuple[int, bool], ...]]]:
    """Return relabel maps.

    A permutation p is interpreted as: new vertex i is old vertex p[i].
    For each new edge index, store (old_edge_index, invert_bit).
    """
    idx = edge_index(n)
    maps = []
    for p in permutations(range(n)):
        mapping = []
        for a, b in edge_pairs(n):
            old_a = p[a]
            old_b = p[b]
            if old_a < old_b:
                mapping.append((idx[(old_a, old_b)], False))
            else:
                mapping.append((idx[(old_b, old_a)], True))
        maps.append((p, tuple(mapping)))
    return maps


def transform_bits(bits: int, mapping: tuple[tuple[int, bool], ...]) -> int:
    out = 0
    for new_idx, (old_idx, invert) in enumerate(mapping):
        bit = (bits >> old_idx) & 1
        if invert:
            bit ^= 1
        if bit:
            out |= 1 << new_idx
    return out


def enumerate_classes(n: int) -> list[int]:
    maps = permutation_maps(n)
    total = 1 << (n * (n - 1) // 2)
    seen = bytearray(total)
    reps: list[int] = []
    for bits in range(total):
        if seen[bits]:
            continue
        orbit = {transform_bits(bits, mapping) for _, mapping in maps}
        rep = min(orbit)
        reps.append(rep)
        for item in orbit:
            seen[item] = 1
    return reps


def canonical(bits: int, n: int, maps=None) -> int:
    if maps is None:
        maps = permutation_maps(n)
    return min(transform_bits(bits, mapping) for _, mapping in maps)


def automorphisms(bits: int, n: int, maps=None) -> list[tuple[int, ...]]:
    if maps is None:
        maps = permutation_maps(n)
    return [p for p, mapping in maps if transform_bits(bits, mapping) == bits]


def vertex_orbit_count(auts: list[tuple[int, ...]], n: int) -> int:
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for p in auts:
        for i, image in enumerate(p):
            union(i, image)
    return len({find(i) for i in range(n)})


def converse(bits: int, n: int) -> int:
    all_mask = (1 << (n * (n - 1) // 2)) - 1
    return bits ^ all_mask


def adjacency(bits: int, n: int) -> list[list[int]]:
    adj = [[] for _ in range(n)]
    idx = edge_index(n)
    for i, j in edge_pairs(n):
        if (bits >> idx[(i, j)]) & 1:
            adj[i].append(j)
        else:
            adj[j].append(i)
    return adj


def strongly_connected(bits: int, n: int) -> bool:
    adj = adjacency(bits, n)
    radj = [[] for _ in range(n)]
    for i, outs in enumerate(adj):
        for j in outs:
            radj[j].append(i)

    def reach(graph: list[list[int]], start: int) -> int:
        seen = 1 << start
        stack = [start]
        while stack:
            v = stack.pop()
            for u in graph[v]:
                if not (seen >> u) & 1:
                    seen |= 1 << u
                    stack.append(u)
        return seen

    full = (1 << n) - 1
    return reach(adj, 0) == full and reach(radj, 0) == full


def score_sequence(bits: int, n: int) -> tuple[int, ...]:
    scores = [0] * n
    idx = edge_index(n)
    for i, j in edge_pairs(n):
        if (bits >> idx[(i, j)]) & 1:
            scores[i] += 1
        else:
            scores[j] += 1
    return tuple(sorted(scores))


def class_data(n: int) -> list[TournamentClass]:
    maps = permutation_maps(n)
    reps = enumerate_classes(n)
    out = []
    for rep in reps:
        orbit = {transform_bits(rep, mapping) for _, mapping in maps}
        auts = automorphisms(rep, n, maps)
        out.append(
            TournamentClass(
                rep=rep,
                orbit_size=len(orbit),
                aut_size=len(auts),
                self_converse=canonical(converse(rep, n), n, maps) == rep,
                strongly_connected=strongly_connected(rep, n),
                vertex_orbit_count=vertex_orbit_count(auts, n),
                score_sequence=score_sequence(rep, n),
            )
        )
    return out


def delete_vertex(bits: int, n: int, vertex: int) -> int:
    old_vertices = [v for v in range(n) if v != vertex]
    old_idx = edge_index(n)
    out = 0
    out_idx = edge_index(n - 1)
    for i, j in edge_pairs(n - 1):
        oi, oj = old_vertices[i], old_vertices[j]
        if oi < oj:
            bit = (bits >> old_idx[(oi, oj)]) & 1
        else:
            bit = ((bits >> old_idx[(oj, oi)]) & 1) ^ 1
        if bit:
            out |= 1 << out_idx[(i, j)]
    return out


def extend_by_vertex(bits: int, n: int, new_beats_mask: int) -> int:
    """Add vertex n.  Mask bit i=1 means new vertex beats old vertex i."""
    old_idx = edge_index(n)
    new_idx = edge_index(n + 1)
    out = 0
    for i, j in edge_pairs(n):
        if (bits >> old_idx[(i, j)]) & 1:
            out |= 1 << new_idx[(i, j)]
    for i in range(n):
        # Edge index is (i,n).  Bit 1 means i -> n; bit 0 means n -> i.
        if not ((new_beats_mask >> i) & 1):
            out |= 1 << new_idx[(i, n)]
    return out


def permute_subset(mask: int, p: tuple[int, ...]) -> int:
    out = 0
    for new_i, old_i in enumerate(p):
        if (mask >> old_i) & 1:
            out |= 1 << new_i
    return out


def subset_orbits(n: int, auts: list[tuple[int, ...]]) -> list[set[int]]:
    seen: set[int] = set()
    out: list[set[int]] = []
    for mask in range(1 << n):
        if mask in seen:
            continue
        orbit = {permute_subset(mask, p) for p in auts}
        seen |= orbit
        out.append(orbit)
    return out


def load_s367():
    path = "04-computation/lonely_runner_k13_scalar_gauge_s367.py"
    spec = importlib.util.spec_from_file_location("s367", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def lrc_missed_cell_ledger() -> None:
    s367 = load_s367()
    system = s367.build_pattern_system(14)
    vector = (0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0)
    missed = s367.missed_candidates(system, vector)
    by_pattern: dict[int, list[int]] = defaultdict(list)
    for shift, pattern_idx in missed:
        by_pattern[pattern_idx].append(shift)

    pattern_items = []
    for pattern_idx in sorted(by_pattern):
        pattern = system.patterns[pattern_idx]
        pattern_items.append((pattern_idx, pattern))

    mirror_pairs = []
    used = set()
    for p_idx, pattern in pattern_items:
        if p_idx in used:
            continue
        mirror_bins = tuple((system.n - 1 - value) % system.n for value in pattern.bins)
        mate = next(q for q, candidate in pattern_items if candidate.bins == mirror_bins)
        mirror_pairs.append((p_idx, mate, pattern.hi - pattern.lo))
        used.add(p_idx)
        used.add(mate)

    widths = Counter(pattern.hi - pattern.lo for _, pattern in pattern_items)
    outer_width = max(widths)
    outer_patterns = [p for p, pattern in pattern_items if pattern.hi - pattern.lo == outer_width]
    inner_patterns = [p for p, pattern in pattern_items if pattern.hi - pattern.lo != outer_width]

    print("S367 LRC missed-cell ledger")
    print(f"  missed_cells={len(missed)}")
    print(f"  odd_shifts={sorted({s for s, _ in missed})}")
    print(f"  alpha_stencils={len(pattern_items)}")
    print(f"  mirror_pairs={mirror_pairs}")
    print(f"  width_hist={[(s367.fmt_frac(width), count) for width, count in sorted(widths.items())]}")
    print(
        "  14+42 split="
        f"outer_pair {len(outer_patterns)} stencils * 7 shifts = {len(outer_patterns) * 7}; "
        f"inner {len(inner_patterns)} stencils * 7 shifts = {len(inner_patterns) * 7}"
    )
    print()


def main() -> None:
    print("LRC / tournament 56 bridge sprint (S369)")
    print()

    data5 = class_data(5)
    data6 = class_data(6)
    print("Tournament class ledger")
    print(f"  T(5)={len(data5)}")
    print(f"  T(6)={len(data6)}")

    self_converse6 = [item for item in data6 if item.self_converse]
    chiral6 = [item for item in data6 if not item.self_converse]
    print(
        "  T(6) self-converse/chiral split="
        f"{len(self_converse6)} + {len(chiral6)} = {len(data6)}"
    )
    print(f"  chiral converse pairs={len(chiral6)//2}")
    print(
        "  This realizes 56-12=44 structurally: "
        "the 12 is also the self-converse n=6 layer."
    )
    print(f"  strongly_connected_hist={Counter(item.strongly_connected for item in data6)}")
    print(
        "  converse_x_strong_hist="
        f"{Counter((item.self_converse, item.strongly_connected) for item in data6)}"
    )
    print(f"  automorphism_size_hist_n6={Counter(item.aut_size for item in data6)}")
    print(f"  vertex_orbit_hist_n5={Counter(item.vertex_orbit_count for item in data5)}")
    print(f"  sum_vertex_orbits_n5={sum(item.vertex_orbit_count for item in data5)}")
    print()

    maps5 = permutation_maps(5)
    maps6 = permutation_maps(6)
    rep5_to_idx = {item.rep: idx for idx, item in enumerate(data5)}
    extension_orbit_total = 0
    child_sets_by_parent = []
    for idx5, item in enumerate(data5):
        auts = automorphisms(item.rep, 5, maps5)
        orbits = subset_orbits(5, auts)
        extension_orbit_total += len(orbits)
        children = {
            canonical(extend_by_vertex(item.rep, 5, min(orbit)), 6, maps6)
            for orbit in orbits
        }
        child_sets_by_parent.append(children)

    child_support_hist = Counter()
    parent_support_hist = Counter()
    parent_edge_count = 0
    for item6 in data6:
        parents = set()
        for v in range(6):
            parent = canonical(delete_vertex(item6.rep, 6, v), 5, maps5)
            parents.add(rep5_to_idx[parent])
        parent_support_hist[len(parents)] += 1
        parent_edge_count += len(parents)
        for p in parents:
            child_support_hist[p] += 1

    print("5-to-6 extension/deletion ledger")
    print(f"  subset-extension orbits over all 5-classes={extension_orbit_total}")
    print(f"  distinct 6-children covered={len(set().union(*child_sets_by_parent))}")
    print(f"  parent_support_hist_for_6_classes={sorted(parent_support_hist.items())}")
    print(f"  class-level deletion incidence edges={parent_edge_count}")
    print(
        "  child counts by 5-parent="
        f"{[(idx, child_support_hist[idx]) for idx in range(len(data5))]}"
    )
    print(
        "  perspective gap="
        f"T(6) - sum_vertex_orbits(T5) = {len(data6)} - "
        f"{sum(item.vertex_orbit_count for item in data5)} = "
        f"{len(data6) - sum(item.vertex_orbit_count for item in data5)}"
    )
    print()

    lrc_missed_cell_ledger()

    print("Synthesis")
    print("  1. 56=12+44 is not only T(6)=T(5)+44.")
    print("     At n=6 it is exactly self-converse classes plus chiral classes.")
    print("  2. The LRC missed stencils also come in mirror pairs, so chirality is")
    print("     the first tournament-side structure worth comparing, not raw classes.")
    print("  3. 56=14+42 is visible inside the missed cells: one mirror stencil pair")
    print("     accounts for 14 cells, and the six interior stencils account for 42.")
    print("     On the tournament side, 42 is also twice the 21 non-strong n=6 classes.")
    print("  4. The old perspective conjecture fails by 8 at n=5->6; this matches the")
    print("     eight LRC alpha stencils and deserves a direct stencil/perspective test.")


if __name__ == "__main__":
    main()
