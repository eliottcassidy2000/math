#!/usr/bin/env python3
"""cgp_planar_witness_kps_S128c109.py -- kind-pasteur-2026-07-20-S128c109

THE CGP QUESTION, IN THE READING THE OWNER'S SENTENCE FORCES.

  "a point of K_n is a vertex, its clique is a STAR, stars are INCIDENCE ROWS --
   whence duality."

Take that literally at the TILE level.  Let the faces be the TILES (non-base-path
arcs of K_n) and let two tiles be adjacent when they share a vertex of K_n.  Then
the clique at the point v is exactly the star S_v, and the natural bipartite witness
is forced:

    B  =  the TILE-VERTEX INCIDENCE graph  of  K_n minus the base path,
    tiles on one side, points (vertices of K_n) on the other.

Chen-Grigni-Papadimitriou: G is a map graph iff G is the half-square of a PLANAR
bipartite graph.  The half-square of B on the tile side is exactly the tile-adjacency
graph, so the CGP question becomes a planarity question about ONE explicit graph --
and the incidence graph of a graph H is planar precisely when H is planar (it is a
subdivision of H).  Hence

    tile-adjacency graph is a map graph VIA THIS WITNESS
       <=>  K_n minus a Hamiltonian path is planar.

That is decidable outright, so this reading of the question closes.  Note what it
does NOT settle: CGP asks whether SOME planar bipartite witness exists, so a
non-planar B rules out only this witness, not every witness.  Reported honestly.

Also computed: the alternating-vertex clique K_ceil(n/2) inside K_n minus the path
(vertices 0,2,4,... are pairwise non-consecutive so all their edges survive), which
forces non-planarity for all large n independently of any planarity routine.
"""
import sys
from itertools import combinations

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10

try:
    import networkx as nx
    HAVE_NX = True
except Exception:
    HAVE_NX = False


def kn_minus_path(n):
    """K_n minus the base Hamiltonian path n-1 -> ... -> 0."""
    E = set()
    for a, b in combinations(range(n), 2):
        if b - a >= 2:
            E.add((a, b))
    return E


print("=" * 78)
print("IS  K_n MINUS A HAMILTONIAN PATH  PLANAR?")
print("  (equivalently: is the forced tile-vertex incidence witness B planar?)")
print("=" * 78)
print("  networkx available for exact planarity: %s" % HAVE_NX)
print()
print("  %-4s %-8s %-8s %-10s %-14s %s"
      % ("n", "|V|", "|E|", "3n-6", "alt-clique", "planar?"))
first_bad = None
for n in range(3, NMAX + 1):
    E = kn_minus_path(n)
    e = len(E)
    bound = 3 * n - 6 if n >= 3 else 3
    altc = (n + 1) // 2          # vertices 0,2,4,... are pairwise non-adjacent on the path
    if HAVE_NX:
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(E)
        planar = nx.check_planarity(G)[0]
    else:
        planar = None
    verdict = ("YES" if planar else "NO") if planar is not None else "?"
    if planar is False and first_bad is None:
        first_bad = n
    flag = ""
    if e > bound and n >= 3:
        flag = "  (exceeds 3n-6, non-planar by Euler)"
    if altc >= 5:
        flag += "  (contains K_%d, so K_5 minor)" % altc
    print("  %-4d %-8d %-8d %-10d K_%-12d %s%s" % (n, n, e, bound, altc, verdict, flag))

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
if first_bad:
    print("  B is planar for n < %d and NON-planar from n = %d onward." % (first_bad, first_bad))
else:
    print("  planarity did not fail in the tested range (or networkx unavailable).")
print()
print("  The alternating-vertex argument is independent of any planarity routine:")
print("  vertices 0,2,4,... are pairwise NON-consecutive, so every edge among them")
print("  survives the removal of the base path, giving a clique of size ceil(n/2).")
print("  That is a K_5 as soon as ceil(n/2) >= 5, i.e. n >= 9, so B is non-planar")
print("  for ALL n >= 9 with no computation at all.")
print()
print("  CONSEQUENCE for the CGP question, stated at its true strength:")
print("   * In the TILE reading the witness is FORCED, and it is planar only in a")
print("     small initial range -- so the tile-adjacency graph is NOT a map graph")
print("     via its natural witness beyond that range.")
print("   * This does NOT answer CGP for G^(<=k) on ISO CLASSES, which is a different")
print("     graph, nor does it rule out some OTHER planar witness: CGP is an")
print("     existential over all bipartite witnesses.  A non-planar natural witness")
print("     is evidence, not a proof.")
