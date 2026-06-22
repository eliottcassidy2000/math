#!/usr/bin/env python3
"""
e7_c5_h7_obstruction_map_codex_s100.py

Audit whether the C5 holes in the E7 even-graph metagraph literally realize the
H=7 / Omega=K3 obstruction under the fixed-path tournament <-> even-graph
cycle-space bijection.

Gauge used here:
  - Fix the Hamiltonian path 0 -> 1 -> ... -> 6.
  - A free arc (i,j), j>i+1, is a bit.
  - bit=1 means the arc is reversed, j -> i; this creates the directed cycle
    i -> i+1 -> ... -> j -> i.
  - The corresponding even-graph basis vector is the undirected fundamental
    cycle on the same path interval plus chord (i,j).

This is the path-fundamental-cycle form of the even-graph/tournament-cycle
bijection, matching the E7 metagraph construction.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations, product


N = 7


def chords(n: int = N) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 2, n)]


CHORDS = chords(N)
CHORD_INDEX = {edge: idx for idx, edge in enumerate(CHORDS)}


def fundamental_cycle_edges(chord: tuple[int, int]) -> frozenset[tuple[int, int]]:
    i, j = chord
    edges = {(i, j)} | {(k, k + 1) for k in range(i, j)}
    return frozenset((min(a, b), max(a, b)) for a, b in edges)


FC = [fundamental_cycle_edges(c) for c in CHORDS]


K7_EDGES = [(i, j) for i in range(N) for j in range(i + 1, N)]
K7_EDGE_INDEX = {edge: idx for idx, edge in enumerate(K7_EDGES)}


def edges_to_bitset(edges: set[tuple[int, int]] | frozenset[tuple[int, int]]) -> int:
    out = 0
    for a, b in edges:
        if a > b:
            a, b = b, a
        out |= 1 << K7_EDGE_INDEX[(a, b)]
    return out


FC_BITS = [edges_to_bitset(cyc) for cyc in FC]


def evengraph_bits(mask: int) -> int:
    bits = 0
    for b, cyc_bits in enumerate(FC_BITS):
        if (mask >> b) & 1:
            bits ^= cyc_bits
    return bits


def bitset_to_edges(bits: int) -> frozenset[tuple[int, int]]:
    return frozenset(edge for idx, edge in enumerate(K7_EDGES) if (bits >> idx) & 1)


def evengraph(mask: int) -> frozenset[tuple[int, int]]:
    return bitset_to_edges(evengraph_bits(mask))


def graph_degrees(bits: int) -> tuple[int, ...]:
    deg = [0] * N
    for idx, (a, b) in enumerate(K7_EDGES):
        if (bits >> idx) & 1:
            deg[a] += 1
            deg[b] += 1
    return tuple(deg)


def permute_graph_bits(bits: int, mapping: tuple[int, ...]) -> int:
    out = 0
    work = bits
    while work:
        lsb = work & -work
        idx = lsb.bit_length() - 1
        a, b = K7_EDGES[idx]
        na, nb = mapping[a], mapping[b]
        if na > nb:
            na, nb = nb, na
        out |= 1 << K7_EDGE_INDEX[(na, nb)]
        work ^= lsb
    return out


def color_signature(bits: int) -> tuple[int, ...]:
    """1-WL color refinement, enough to make n=7 canonicalization cheap."""
    colors = graph_degrees(bits)
    while True:
        descriptors = []
        for v in range(N):
            neigh = []
            for u in range(N):
                if u == v:
                    continue
                a, b = (u, v) if u < v else (v, u)
                if (bits >> K7_EDGE_INDEX[(a, b)]) & 1:
                    neigh.append(colors[u])
            descriptors.append((colors[v], tuple(sorted(neigh))))
        uniq = {d: i for i, d in enumerate(sorted(set(descriptors)))}
        new_colors = tuple(uniq[d] for d in descriptors)
        if new_colors == colors:
            return colors
        colors = new_colors


def canonical_graph_bits(bits: int) -> int:
    colors = color_signature(bits)
    groups: list[list[int]] = []
    for color in sorted(set(colors)):
        groups.append([v for v, c in enumerate(colors) if c == color])

    target_blocks = []
    start = 0
    for group in groups:
        block = list(range(start, start + len(group)))
        target_blocks.append(block)
        start += len(group)

    best = None
    for assignments in product(
        *(permutations(block) for block in target_blocks)
    ):
        mapping = [None] * N
        for group, images in zip(groups, assignments):
            for old, new in zip(group, images):
                mapping[old] = new
        image = permute_graph_bits(bits, tuple(mapping))  # type: ignore[arg-type]
        if best is None or image < best:
            best = image
    assert best is not None
    return best


def evengraph_edges_from_bits(mask: int) -> frozenset[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    for b, cyc in enumerate(FC):
        if (mask >> b) & 1:
            edges ^= set(cyc)
    return frozenset(edges)


def mask_from_even_graph(edges: frozenset[tuple[int, int]]) -> int:
    """Inverse in the path-fundamental basis: chord coefficients are chord edges."""
    out = 0
    edge_set = set(edges)
    for idx, edge in enumerate(CHORDS):
        if edge in edge_set:
            out |= 1 << idx
    assert evengraph(out) == edges
    return out


def build_e7_classes():
    cls_of = [0] * (1 << len(CHORDS))
    canon_to_id: dict[int, int] = {}
    rep_bits: list[int] = []

    for mask in range(1 << len(CHORDS)):
        canon = canonical_graph_bits(evengraph_bits(mask))
        cid = canon_to_id.get(canon)
        if cid is None:
            cid = len(rep_bits)
            canon_to_id[canon] = cid
            rep_bits.append(canon)
        cls_of[mask] = cid

    meta_adj = [set() for _ in rep_bits]
    for mask in range(1 << len(CHORDS)):
        a = cls_of[mask]
        for b in range(len(CHORDS)):
            nb = cls_of[mask ^ (1 << b)]
            if nb != a:
                meta_adj[a].add(nb)
                meta_adj[nb].add(a)

    return cls_of, rep_bits, meta_adj


def chordless_c5s(meta_adj: list[set[int]]) -> list[tuple[int, ...]]:
    holes = []
    for verts in combinations(range(len(meta_adj)), 5):
        deg = {v: 0 for v in verts}
        edge_count = 0
        for a, b in combinations(verts, 2):
            if b in meta_adj[a]:
                deg[a] += 1
                deg[b] += 1
                edge_count += 1
        if edge_count == 5 and all(d == 2 for d in deg.values()):
            holes.append(verts)
    return holes


def tournament_adj(mask: int) -> list[list[int]]:
    adj = [[0] * N for _ in range(N)]
    for i in range(N - 1):
        adj[i][i + 1] = 1
    for idx, (i, j) in enumerate(CHORDS):
        if (mask >> idx) & 1:
            adj[j][i] = 1
        else:
            adj[i][j] = 1
    return adj


TRIANGLES = list(combinations(range(N), 3))
FIVE_SUBSETS = list(combinations(range(N), 5))
SEVEN_SUBSETS = [tuple(range(N))]
FIVE_CYCLES = [
    (subset[0],) + perm
    for subset in FIVE_SUBSETS
    for perm in permutations(subset[1:])
]
SEVEN_CYCLES = [
    (0,) + perm
    for perm in permutations(range(1, N))
]


def directed_cycle(adj: list[list[int]], cyc: tuple[int, ...]) -> bool:
    return all(adj[cyc[i]][cyc[(i + 1) % len(cyc)]] for i in range(len(cyc)))


def triangle_cycles(adj: list[list[int]]) -> list[tuple[int, int, int]]:
    out = []
    for a, b, c in TRIANGLES:
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            out.append((a, b, c))
    return out


def cycle_stats(mask: int) -> dict[str, object]:
    adj = tournament_adj(mask)
    tris = triangle_cycles(adj)
    c5 = sum(1 for cyc in FIVE_CYCLES if directed_cycle(adj, cyc))
    c7 = sum(1 for cyc in SEVEN_CYCLES if directed_cycle(adj, cyc))
    tri_disjoint = 0
    for a, b in combinations(tris, 2):
        if set(a).isdisjoint(b):
            tri_disjoint += 1
    alpha1 = len(tris) + c5 + c7
    alpha2 = tri_disjoint
    return {
        "c3": len(tris),
        "c5": c5,
        "c7": c7,
        "alpha1": alpha1,
        "alpha2": alpha2,
        "near_h7_p3": alpha1 == 3 and alpha2 == 1,
        "forbidden_h7_point": alpha1 == 3 and alpha2 == 0,
        # Three directed triangles, no disjoint triangle pair, and a forced
        # directed pentagon: the local THM-200/K3-forces-C5 shadow.
        "k3_forces_pentagon": len(tris) == 3 and tri_disjoint == 0 and c5 > 0,
    }


def class_of_cycle_support(cls_of: list[int], cycle_vertices: tuple[int, ...]) -> int:
    edges = set()
    for i in range(len(cycle_vertices)):
        a = cycle_vertices[i]
        b = cycle_vertices[(i + 1) % len(cycle_vertices)]
        edges.add((min(a, b), max(a, b)))
    mask = mask_from_even_graph(frozenset(edges))
    return cls_of[mask]


def pct(num: int, den: int) -> str:
    if den == 0:
        return "n/a"
    return f"{100.0 * num / den:.2f}%"


def main() -> None:
    cls_of, rep_bits, meta_adj = build_e7_classes()
    assert len(rep_bits) == 54

    holes5 = chordless_c5s(meta_adj)
    assert len(holes5) == 1496, len(holes5)

    class_masks: dict[int, list[int]] = defaultdict(list)
    for mask, cid in enumerate(cls_of):
        class_masks[cid].append(mask)

    summary = {
        cid: {
            "masks": len(masks),
            "any_near_h7_p3": False,
            "any_forbidden_h7_point": False,
            "any_k3_forces_pentagon": False,
            "min_alpha1": 10**9,
            "max_alpha1": -1,
            "min_c3": 10**9,
            "max_c3": -1,
            "min_c5": 10**9,
            "max_c5": -1,
        }
        for cid, masks in class_masks.items()
    }

    global_counts = Counter()
    for mask, cid in enumerate(cls_of):
        st = cycle_stats(mask)
        for key in ["near_h7_p3", "forbidden_h7_point", "k3_forces_pentagon"]:
            if st[key]:
                global_counts[key] += 1
                summary[cid][f"any_{key}"] = True
        summary[cid]["min_alpha1"] = min(summary[cid]["min_alpha1"], st["alpha1"])
        summary[cid]["max_alpha1"] = max(summary[cid]["max_alpha1"], st["alpha1"])
        summary[cid]["min_c3"] = min(summary[cid]["min_c3"], st["c3"])
        summary[cid]["max_c3"] = max(summary[cid]["max_c3"], st["c3"])
        summary[cid]["min_c5"] = min(summary[cid]["min_c5"], st["c5"])
        summary[cid]["max_c5"] = max(summary[cid]["max_c5"], st["c5"])

    hole_vertices = Counter(v for h in holes5 for v in h)
    hole_vertex_set = set(hole_vertices)
    near_classes = {cid for cid, s in summary.items() if s["any_near_h7_p3"]}
    forbidden_classes = {cid for cid, s in summary.items() if s["any_forbidden_h7_point"]}
    k3_forced_classes = {cid for cid, s in summary.items() if s["any_k3_forces_pentagon"]}

    def hole_counts(flagged: set[int]) -> tuple[int, int, int]:
        full = sum(1 for h in holes5 if all(v in flagged for v in h))
        anyhit = sum(1 for h in holes5 if any(v in flagged for v in h))
        incid = sum(1 for h in holes5 for v in h if v in flagged)
        return full, anyhit, incid

    # The literal directed pentagon support is a single E7 vertex class, not a
    # metagraph C5 hole.  Compute its class and hole incidence.
    pentagon_support_classes = {
        class_of_cycle_support(cls_of, cyc)
        for cyc in combinations(range(N), 5)
    }

    print("=== E7 C5 holes vs H=7/K3 obstruction audit (codex S100) ===")
    print(f"E7 classes: {len(rep_bits)}")
    print(f"E7 metagraph edges: {sum(len(n) for n in meta_adj) // 2}")
    print(f"Chordless C5 holes: {len(holes5)}")
    print(f"C5 holes touch classes: {len(hole_vertex_set)}/54")
    print()

    print("Fixed-path tournament flags over all 2^15 masks:")
    for key in ["forbidden_h7_point", "near_h7_p3", "k3_forces_pentagon"]:
        print(f"  {key}: {global_counts[key]} masks")
    print(f"  classes with forbidden_h7_point: {len(forbidden_classes)}")
    print(f"  classes with near_h7_p3: {len(near_classes)}")
    print(f"  classes with k3_forces_pentagon: {len(k3_forced_classes)}")
    print()

    print("C5-hole alignment with obstruction flags:")
    for name, flagged in [
        ("forbidden_h7_point", forbidden_classes),
        ("near_h7_p3", near_classes),
        ("k3_forces_pentagon", k3_forced_classes),
    ]:
        full, anyhit, incid = hole_counts(flagged)
        print(f"  {name}:")
        print(f"    flagged classes touched by C5 holes: {len(flagged & hole_vertex_set)}/{len(flagged) or 1}")
        print(f"    C5 holes fully inside flagged classes: {full}/{len(holes5)} ({pct(full, len(holes5))})")
        print(f"    C5 holes with at least one flagged class: {anyhit}/{len(holes5)} ({pct(anyhit, len(holes5))})")
        print(f"    flagged vertex incidences in C5 holes: {incid}/{5 * len(holes5)} ({pct(incid, 5 * len(holes5))})")
    print()

    print("Directed pentagon support under the cycle-space map:")
    print(f"  undirected 5-cycle support classes: {sorted(pentagon_support_classes)}")
    for cid in sorted(pentagon_support_classes):
        print(
            "  class"
            f" {cid}: edge_count={rep_bits[cid].bit_count()},"
            f" C5-hole incidence={hole_vertices.get(cid, 0)},"
            f" any_k3_forces_pentagon={summary[cid]['any_k3_forces_pentagon']},"
            f" any_near_h7_p3={summary[cid]['any_near_h7_p3']}"
        )
    print()

    top_hole_classes = hole_vertices.most_common(10)
    print("Top C5-hole class incidences:")
    for cid, count in top_hole_classes:
        s = summary[cid]
        print(
            f"  class {cid:2d}: incidence={count:4d},"
            f" edges={rep_bits[cid].bit_count():2d},"
            f" alpha1=[{s['min_alpha1']},{s['max_alpha1']}],"
            f" c3=[{s['min_c3']},{s['max_c3']}],"
            f" c5=[{s['min_c5']},{s['max_c5']}],"
            f" k3_forces={s['any_k3_forces_pentagon']}"
        )
    print()

    print("Verdict:")
    print("  The literal H=7 point alpha=(3,0) has zero masks and zero classes, as expected.")
    print("  E7 C5 metagraph holes are not the same object as the directed pentagon support:")
    print("  the directed C5 support maps to a single E7 vertex class, while a metagraph C5")
    print("  is a five-class quotient cycle.  The useful bridge is therefore not an exact")
    print("  identification of E7 C5 holes with Omega=K3.  It is an obstruction-layer")
    print("  analogy: both are first odd cycles created by cycle-space closure after scalar")
    print("  data forgets incidence.")
    print()

    print("Tournament Analysis:")
    print("  vertices = {E7_C5_hole_classes, directed_C5_support_class, near_h7_P3,")
    print("              k3_forces_pentagon, strong_component_atom, euler_totient_packet}")
    print("  challenged assumption = the E7 C5 hole and the H=7 pentagon are the same")
    print("  object under the bijection; audit says they live at different quotient levels.")


if __name__ == "__main__":
    main()
