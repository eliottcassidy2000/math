#!/usr/bin/env python3
"""
codex-S105: Clebsch, Bruhat/truncated-octahedral, and unital-style carriers for
the LRC(14) signed residual.

The goal is not to prove LRC14 by importing graph facts.  It is to test whether
the user's three prompts supply usable quotient objects for the current
HYP-2889/HYP-2890 proof state:

  * Clebsch graph: folded 5-cube on 16 vertices.  Use it as a quotient of the
    64 residual masks after choosing a tangent/anchor sector.
  * Truncated octahedral graph: Bruhat graph B_4 = Cayley graph of S_4 by
    adjacent transpositions.  Use it as the local compression/swap carrier.
  * Unital: design language.  Use tangent/secant discipline as the warning that
    a useful quotient should have controlled pair incidence, not only a scalar
    energy.

Tournament Analysis discipline:
  vertices are proof carriers (Clebsch quotient, Bruhat compression, unital
  incidence, scalar energy), not runners.  The observable is
  (preserves signed labels, finite exactness, compression relevance,
  scalar-loss warning).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, permutations

import networkx as nx


def bits(x: int, n: int) -> tuple[int, ...]:
    return tuple(i for i in range(n) if (x >> i) & 1)


def popcount(x: int) -> int:
    return x.bit_count()


def canonical_fold5(x: int) -> int:
    comp = x ^ 31
    return min(x, comp)


def clebsch_graph() -> nx.Graph:
    verts = sorted({canonical_fold5(x) for x in range(32)})
    g = nx.Graph()
    g.add_nodes_from(verts)
    for a, b in combinations(verts, 2):
        d = popcount(a ^ b)
        # Folded 5-cube adjacency: distance 1 after antipodal quotient.
        if min(d, 5 - d) == 1:
            g.add_edge(a, b)
    return g


def strongly_regular_params(g: nx.Graph) -> tuple[bool, dict[str, object]]:
    n = g.number_of_nodes()
    degs = sorted(dict(g.degree()).values())
    k = degs[0] if len(set(degs)) == 1 else None
    adj_common = Counter()
    non_common = Counter()
    for a, b in combinations(g.nodes, 2):
        c = len(set(g.neighbors(a)) & set(g.neighbors(b)))
        if g.has_edge(a, b):
            adj_common[c] += 1
        else:
            non_common[c] += 1
    ok = len(set(degs)) == 1 and adj_common == Counter({0: g.number_of_edges()})
    ok = ok and len(non_common) == 1 and next(iter(non_common)) == 2
    return ok, {
        "n": n,
        "m": g.number_of_edges(),
        "degree_hist": Counter(degs),
        "adjacent_common_neighbors": dict(adj_common),
        "nonadjacent_common_neighbors": dict(non_common),
    }


def closed_neighborhood_design(g: nx.Graph) -> dict[str, object]:
    blocks = {v: frozenset([v, *g.neighbors(v)]) for v in g.nodes}
    block_sizes = Counter(len(b) for b in blocks.values())
    point_degrees = Counter()
    pair_degrees = Counter()
    for block in blocks.values():
        for v in block:
            point_degrees[v] += 1
        for a, b in combinations(sorted(block), 2):
            pair_degrees[(a, b)] += 1
    return {
        "num_blocks": len(blocks),
        "block_size_hist": dict(block_sizes),
        "point_rep_hist": dict(Counter(point_degrees.values())),
        "pair_rep_hist": dict(Counter(pair_degrees.values())),
    }


def residual_mask_to_clebsch(mask: int, tangent: int = 5) -> tuple[int, int, int]:
    """Project a 6-sector residual mask to a tangent sector plus folded 5-cube class.

    Returns (tangent_bit, folded_5_class, residual_depth).
    """
    tangent_bit = (mask >> tangent) & 1
    x = 0
    out = 0
    for i in range(6):
        if i == tangent:
            continue
        if (mask >> i) & 1:
            x |= 1 << out
        out += 1
    return tangent_bit, canonical_fold5(x), popcount(mask)


def residual_quotient_stats(tangent: int = 5) -> dict[str, object]:
    by_class: dict[int, list[int]] = defaultdict(list)
    by_tangent_class: dict[tuple[int, int], list[int]] = defaultdict(list)
    for mask in range(64):
        tb, cls, depth = residual_mask_to_clebsch(mask, tangent)
        by_class[cls].append(depth)
        by_tangent_class[(tb, cls)].append(depth)
    class_depths = {cls: Counter(ds) for cls, ds in by_class.items()}
    tangent_depths = {key: Counter(ds) for key, ds in by_tangent_class.items()}
    mixed_classes = sum(1 for c in class_depths.values() if len(c) > 1)
    return {
        "classes": len(by_class),
        "preimage_size_hist": dict(Counter(len(v) for v in by_class.values())),
        "mixed_depth_classes": mixed_classes,
        "sample_depth_profiles": dict(list(class_depths.items())[:6]),
        "tangent_refined_classes": len(by_tangent_class),
        "tangent_refined_preimage_hist": dict(Counter(len(v) for v in by_tangent_class.values())),
        "sample_tangent_depth_profiles": dict(list(tangent_depths.items())[:6]),
    }


def adjacent_swap(p: tuple[int, ...], i: int) -> tuple[int, ...]:
    q = list(p)
    q[i], q[i + 1] = q[i + 1], q[i]
    return tuple(q)


def bruhat_s4_graph() -> nx.Graph:
    verts = list(permutations(range(4)))
    g = nx.Graph()
    g.add_nodes_from(verts)
    for p in verts:
        for i in range(3):
            q = adjacent_swap(p, i)
            if p < q:
                g.add_edge(p, q, label=i + 1)
    return g


def inversion_count(p: tuple[int, ...]) -> int:
    return sum(1 for i, j in combinations(range(len(p)), 2) if p[i] > p[j])


def bruhat_faces(g: nx.Graph) -> dict[str, object]:
    # Commutation squares: cosets of <s1,s3>.
    squares = set()
    for p in g.nodes:
        cyc = frozenset(
            {
                p,
                adjacent_swap(p, 0),
                adjacent_swap(p, 2),
                adjacent_swap(adjacent_swap(p, 0), 2),
            }
        )
        if len(cyc) == 4:
            squares.add(cyc)

    # Braid hexagons: cosets of adjacent S3 subgroups <s1,s2> and <s2,s3>.
    hexagons = set()
    for p in g.nodes:
        for pair in ((0, 1), (1, 2)):
            seen = {p}
            frontier = [p]
            while frontier:
                cur = frontier.pop()
                for i in pair:
                    nxt = adjacent_swap(cur, i)
                    if nxt not in seen:
                        seen.add(nxt)
                        frontier.append(nxt)
            if len(seen) == 6:
                hexagons.add(frozenset(seen))
    return {"commutation_squares": len(squares), "braid_hexagons": len(hexagons)}


def edge_orbits(g: nx.Graph) -> dict[str, object]:
    matcher = nx.algorithms.isomorphism.GraphMatcher(g, g)
    autos = list(matcher.isomorphisms_iter())
    edge_set = {frozenset(e) for e in g.edges}
    unseen = set(edge_set)
    orbits = []
    while unseen:
        e = unseen.pop()
        orb = set()
        for a in autos:
            image = frozenset(a[v] for v in e)
            orb.add(image)
        orbits.append(orb)
        unseen -= orb
    label_counts = []
    for orb in orbits:
        c = Counter()
        for e in orb:
            u, v = tuple(e)
            c[g.edges[u, v]["label"]] += 1
        label_counts.append(dict(c))
    return {
        "automorphisms": len(autos),
        "edge_orbit_sizes": sorted(len(o) for o in orbits),
        "edge_orbit_label_counts": label_counts,
    }


def bruhat_stats() -> dict[str, object]:
    g = bruhat_s4_graph()
    degs = Counter(dict(g.degree()).values())
    inv_hist = Counter(inversion_count(p) for p in g.nodes)
    edge_labels = Counter(nx.get_edge_attributes(g, "label").values())
    return {
        "vertices": g.number_of_nodes(),
        "edges": g.number_of_edges(),
        "degree_hist": dict(degs),
        "bipartition_by_inversion_parity": dict(Counter(inversion_count(p) % 2 for p in g.nodes)),
        "inversion_rank_hist": dict(inv_hist),
        "edge_label_counts": dict(edge_labels),
        **bruhat_faces(g),
        **edge_orbits(g),
    }


def tournament_lens() -> list[tuple[str, tuple[int, int, int, int]]]:
    """Hand-scored proof-carrier tournament; higher tuple is better."""
    lenses = [
        ("AP_Jensen_labelled_functional", (4, 4, 4, 3)),
        ("Clebsch_tangent_residual_design", (3, 4, 3, 3)),
        ("Bruhat_S4_compression_carrier", (3, 3, 4, 2)),
        ("unital_tangent_secant_incidence", (2, 4, 2, 3)),
        ("scalar_additive_energy", (1, 1, 2, 1)),
    ]
    return sorted(lenses, key=lambda x: x[1], reverse=True)


def main() -> None:
    print("=" * 78)
    print("S105: Clebsch / Bruhat / unital-style carriers for LRC14 residuals")
    print("=" * 78)

    cleb = clebsch_graph()
    ok, params = strongly_regular_params(cleb)
    print("\n[1] Clebsch folded-5-cube carrier")
    print(f"  strongly_regular_(16,5,0,2) verified: {ok}")
    print(f"  params: {params}")
    print(f"  closed-neighborhood design: {closed_neighborhood_design(cleb)}")
    print(f"  residual-mask quotient by tangent sector 5: {residual_quotient_stats(5)}")
    print("  readout: folded residual classes preserve one-flip adjacency but destroy missed-depth.")
    print("           Therefore Clebsch is a carrier for signed covariance, not a scalar q_t proof.")

    print("\n[2] Truncated-octahedral / Bruhat B4 compression carrier")
    bst = bruhat_stats()
    print(f"  stats: {bst}")
    print("  readout: adjacent swaps form the S4 permutohedron with square and braid faces.")
    print("           Local compression failures in HYP-2889 should be sorted by face type;")
    print("           vertex-transitive but non-edge-transitive means edge labels cannot be erased.")

    print("\n[3] Unital-style tangent/secant lesson")
    print("  MathWorld's geometric-combinatorics unital is a 2-(q^3+1,q+1,1) design.")
    print("  LRC transfer: choose a tangent sector/anchor, then demand controlled pair incidence")
    print("  before scalarizing.  Clebsch closed neighborhoods give a nearby exact 2-(16,6,2)")
    print("  design; this is the finite covariance object, while scalar additive energy forgets")
    print("  tangent/secant labels.")

    print("\n[4] Tournament-analysis proof-lens ranking")
    for name, score in tournament_lens():
        print(f"  {name:35s} score={score}")
    print("  Hamiltonian path:")
    print("    AP_Jensen_labelled_functional > Clebsch_tangent_residual_design >")
    print("    Bruhat_S4_compression_carrier > unital_tangent_secant_incidence > scalar_additive_energy")

    print("\nProof target emitted:")
    print("  Decompose the S104 residual leak over tangent Clebsch classes and Bruhat S4")
    print("  compression faces.  Prove square/commuting faces are nonpositive by design")
    print("  balance, then isolate braid hexagons as the finite AP/Freiman atlas.")


if __name__ == "__main__":
    main()
