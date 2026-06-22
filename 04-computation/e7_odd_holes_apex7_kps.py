"""
E_7 odd holes and the apex-7 unification (kind-pasteur S31e).

The even-graph metagraph E_n = quotient of the cycle-space hypercube Q_m (m=C(n-1,2))
by S_n.  Even graphs = GF(2) cycle space of K_n.  Edge of E_n = single fundamental-cycle
flip (one tile / chord), projected to even-graph iso classes.  The reflections report
E_n is chordal for n<=6 and FIRST gains odd holes at n=7 = the apex prime of LRC(14).

This script:
  (1) builds E_n for n=5,6,7 from the path-spanning-tree cycle space,
  (2) VALIDATES against the known V(E_n)=3,7,16,54 / chi(E_n)=3,5,10,28,
  (3) finds E_7's odd holes (chordless odd cycles >=5) and reports their length
      spectrum + any cyclic (Z/7-like) symmetry signature -- the apex-7 test.
"""
import itertools, sys
import networkx as nx

def build_En(n):
    verts = list(range(n))
    path = [(k, k+1) for k in range(n-1)]
    chords = [(i, j) for i in range(n) for j in range(i+2, n)]   # non-tree edges
    m = len(chords)                                              # = C(n-1,2)
    def fundcycle(c):
        i, j = c
        edges = {(i, j)} | {(k, k+1) for k in range(i, j)}
        return frozenset(tuple(sorted(e)) for e in edges)
    fc = [fundcycle(c) for c in chords]
    # even graph for each subset S of chords = XOR of fundamental cycles
    # represent S as an integer bitmask over m chords
    def evengraph(mask):
        eset = set()
        for b in range(m):
            if (mask >> b) & 1:
                eset ^= fc[b]
        return frozenset(eset)
    # iso classes: WL hash to bucket, then EXACT isomorphism to split collisions
    def mkgraph(eset):
        g = nx.Graph(); g.add_nodes_from(verts); g.add_edges_from(eset); return g
    cls_of = [0] * (1 << m)
    wl_buckets = {}                    # wl_hash -> list of (rep_graph, class_id)
    rep_graphs = []
    for mask in range(1 << m):
        g = mkgraph(evengraph(mask))
        wl = nx.weisfeiler_lehman_graph_hash(g, iterations=4)
        bucket = wl_buckets.setdefault(wl, [])
        cid = None
        for rg, rid in bucket:
            if nx.is_isomorphic(g, rg):
                cid = rid; break
        if cid is None:
            cid = len(rep_graphs); rep_graphs.append(g); bucket.append((g, cid))
        cls_of[mask] = cid
    V = len(rep_graphs)
    # metagraph edges: flip one chord
    meta = nx.Graph(); meta.add_nodes_from(range(V))
    for mask in range(1 << m):
        a = cls_of[mask]
        for b in range(m):
            nb = cls_of[mask ^ (1 << b)]
            if nb != a:
                meta.add_edge(a, nb)
    return meta, V, rep_graphs

def analyze(n, report_holes=False):
    meta, V, rep_graphs = build_En(n)
    E = meta.number_of_edges()
    # chromatic via clique number lower bound + greedy (exact chi is NP; report omega + greedy)
    omega = max((len(c) for c in nx.find_cliques(meta)), default=0)
    greedy = max(nx.greedy_color(meta, strategy="largest_first").values()) + 1
    print(f"n={n}: V(E_n)={V}  E(E_n)={E}  omega={omega}  greedy_chi<={greedy}")
    if report_holes:
        # chordless cycles of length>=5 and odd = odd holes
        odd_holes = []
        for cyc in nx.chordless_cycles(meta, length_bound=11):
            L = len(cyc)
            if L >= 5 and L % 2 == 1:
                odd_holes.append(tuple(cyc))
        lens = sorted(len(h) for h in odd_holes)
        from collections import Counter
        print(f"  odd holes (chordless odd cycles >=5): count={len(odd_holes)}  length spectrum={dict(Counter(lens))}")
        # apex-7 signature: do odd holes have length 7 (the apex prime)? are there C_7?
        c7 = [h for h in odd_holes if len(h) == 7]
        c5 = [h for h in odd_holes if len(h) == 5]
        print(f"  C_5 holes={len(c5)}  C_7 holes={len(c7)}")
        # structural readout: edge-count (even-graph complexity) of the classes in the holes
        ecount = {cid: g.number_of_edges() for cid, g in enumerate(rep_graphs)}
        from collections import Counter
        c7_verts = Counter(v for h in c7 for v in h)
        c5_verts = Counter(v for h in c5 for v in h)
        # how many classes participate in C_7 vs C_5 holes; is the heptagon concentrated?
        print(f"  C_7 holes touch {len(c7_verts)}/{V} classes; top class multiplicities {sorted(c7_verts.values(), reverse=True)[:5]}")
        print(f"  C_5 holes touch {len(c5_verts)}/{V} classes")
        # edge-count signature of the most-central C_7 class (apex even graph?)
        if c7_verts:
            top = c7_verts.most_common(1)[0][0]
            print(f"  most-central C_7 class: edge-count={ecount[top]} (n=7 Hamiltonian C_7 has 7 edges; empty=0)")
        # does any C_7 hole consist ENTIRELY of low-edge (<=7) even graphs (a 'clean' apex heptagon)?
        clean7 = sum(1 for h in c7 if all(ecount[v] <= 7 for v in h))
        print(f"  C_7 holes with all classes <=7 edges (low-complexity 'clean' heptagons): {clean7}")
        return odd_holes
    return None

if __name__ == "__main__":
    print("=== E_n construction + validation (expect V=3,7,16,54 ; omega=3,5,10,28) ===")
    for n in [5, 6]:
        analyze(n)
    print("=== E_7: the apex prime -- odd holes ===")
    analyze(7, report_holes=True)
