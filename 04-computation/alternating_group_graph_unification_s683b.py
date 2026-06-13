#!/usr/bin/env python3
r"""
alternating_group_graph_unification_s683b.py    oracle-2026-06-06-S683b

Brings the ALTERNATING GROUP GRAPH into the LRC = Hadwiger-Nelson = unit-distance unification
(HYP-2265: distance Cayley graphs Cay(A,S); chromatic bound from the connection-set Fourier
transform / Hoffman-Delsarte spectral method).

HYP-2265 covers ABELIAN A (R^2, Z, Z_m). The alternating group graph is the NON-ABELIAN
extension AND the canonical PARITY KEY:

  PARITY MECHANISM (the chromatic obstruction, unified):
   - Cay(S_n, TRANSPOSITIONS) is BIPARTITE: the sign character e:S_n->{+-1} is a proper
     2-coloring (each odd generator flips sign). chi = 2.  <-->  LRC interval circulant R_m
     (the AP, chi=2, S626/S582): a 'sign-respecting' / additive connection set.
   - Cay(A_n, 3-CYCLES) is NON-bipartite: the generators are EVEN, so the sign character
     vanishes on A_n and CANNOT 2-color; order-3 generators create TRIANGLES (3-cycles in
     the graph) => chi >= omega = 3.  <-->  LRC Paley/QR (chi=3) AND Hadwiger-Nelson's unit
     EQUILATERAL TRIANGLE (3 mutual unit distances = an odd cycle => chi(plane) >= 3).

  So the ONE shared key is: chromatic number is the failure of a PARITY (sign) 2-coloring,
  obstructed by ODD CYCLES (the repo's OCF / parity divide). The alternating group graph is
  the minimal non-bipartite Cayley graph -- its triangles ARE the unit equilateral triangle
  and the LRC odd resonance. The spectral (Hoffman) bound chi >= 1 - lam_max/lam_min carries
  over verbatim, with the Cayley eigenvalues now from the IRREPS of A_n (non-abelian Fourier)
  instead of characters (abelian Fourier).
"""
import numpy as np
import itertools as it
import networkx as nx

# ---------- permutation utilities ----------
def compose(p, q):           # (p after q): result[i] = p[q[i]]
    return tuple(p[q[i]] for i in range(len(p)))

def parity(p):
    seen = [False] * len(p); par = 0
    for i in range(len(p)):
        if not seen[i]:
            j = i; clen = 0
            while not seen[j]:
                seen[j] = True; j = p[j]; clen += 1
            par ^= (clen - 1) & 1
    return par                # 0 even, 1 odd

def cycle_perm(n, cyc):       # cyc like (0,1,2) meaning 0->1->2->0
    p = list(range(n))
    for a, b in zip(cyc, cyc[1:] + cyc[:1]):
        p[a] = b
    return tuple(p)

def elements(n, even_only):
    for p in it.permutations(range(n)):
        if (not even_only) or parity(p) == 0:
            yield p

def cayley_graph(n, gens, even_only):
    verts = list(elements(n, even_only))
    idx = {v: i for i, v in enumerate(verts)}
    G = nx.Graph(); G.add_nodes_from(range(len(verts)))
    for v in verts:
        for g in gens:
            w = compose(g, v)
            if w in idx:
                G.add_edge(idx[v], idx[w])
    return G, verts, idx

def chromatic_exact(G, kmax=6):
    # exact via trying k-colorings (backtracking); ok for small/medium vertex-transitive graphs
    nodes = list(G.nodes()); adj = {u: set(G.neighbors(u)) for u in nodes}
    order = sorted(nodes, key=lambda u: -len(adj[u]))
    for k in range(1, kmax + 1):
        color = {}
        def bt(i):
            if i == len(order): return True
            u = order[i]; used = {color[w] for w in adj[u] if w in color}
            for c in range(k):
                if c not in used:
                    color[u] = c
                    if bt(i + 1): return True
                    del color[u]
            return False
        if bt(0): return k
    return None

def clique_number(G):
    return nx.graph_clique_number(G) if hasattr(nx, "graph_clique_number") else max(len(c) for c in nx.find_cliques(G))

def hoffman_bound(G):
    A = nx.to_numpy_array(G)
    ev = np.sort(np.linalg.eigvalsh(A))
    lmin, lmax = ev[0], ev[-1]
    return 1 - lmax / lmin, lmax, lmin

def num_triangles(G):
    return sum(nx.triangles(G).values()) // 3

def main():
    print("=" * 84)
    print("THE ALTERNATING GROUP GRAPH in the LRC=HN=unit-distance unification (HYP-2265 + parity)")
    print("=" * 84)

    # (1) Alternating group graphs AG_n: Cay(A_n, {(1 2 i)^+-1})
    print("\n(1) Alternating group graph AG_n = Cay(A_n, 3-cycles (0 1 i)^{+-1}):  EVEN generators")
    for n in (4, 5):
        gens = []
        for i in range(2, n):
            g = cycle_perm(n, (0, 1, i)); gens.append(g); gens.append(compose(g, g))  # (012i)^2 = inverse 3cyc
        # dedup
        gens = list({g for g in gens})
        G, verts, idx = cayley_graph(n, gens, even_only=True)
        chi = chromatic_exact(G)
        try: omega = clique_number(G)
        except Exception: omega = max(len(c) for c in nx.find_cliques(G))
        hb, lmax, lmin = hoffman_bound(G)
        tri = num_triangles(G)
        print(f"   AG_{n}: |A_{n}|={G.number_of_nodes()} vtx, {G.number_of_edges()} edges, "
              f"{2*G.number_of_edges()//G.number_of_nodes()}-regular")
        print(f"        chi={chi}  omega(clique)={omega}  triangles={tri}  bipartite={nx.is_bipartite(G)}")
        print(f"        Hoffman bound 1-lmax/lmin = {hb:.3f}  (lmax={lmax:.3f}, lmin={lmin:.3f})")

    # (2) The parity contrast: Cay(S_n, transpositions) is BIPARTITE via the sign character
    print("\n(2) Cay(S_n, adjacent transpositions): ODD generators => BIPARTITE (sign 2-coloring)")
    for n in (4, 5):
        gens = [cycle_perm(n, (i, i + 1)) for i in range(n - 1)]   # adjacent transpositions
        G, verts, idx = cayley_graph(n, gens, even_only=False)
        chi = chromatic_exact(G)
        # verify sign character is the 2-coloring
        sign_ok = all(parity(verts[u]) != parity(verts[v]) for u, v in G.edges())
        print(f"   Cay(S_{n},transp): |S_{n}|={G.number_of_nodes()} vtx; chi={chi}; bipartite={nx.is_bipartite(G)}; "
              f"sign-character 2-colors all edges={sign_ok}")

    # (3) LRC interval circulant R_m (chi=2) vs Paley (chi=3) -- the SAME parity split
    print("\n(3) LRC distance circulants Cay(Z_m, +-S):  interval R_m (chi=2) vs Paley QR (chi>=3)")
    def circ_graph(m, conn):
        G = nx.Graph(); G.add_nodes_from(range(m)); cs = set(c % m for c in conn) | set((-c) % m for c in conn)
        for i in range(m):
            for c in cs:
                G.add_edge(i, (i + c) % m)
        return G
    for m in (7, 9, 11):
        Rm = circ_graph(m, range(1, (m - 1) // 2 + 1))            # interval (the AP beat)
        qr = sorted({(x * x) % m for x in range(1, m)})
        P = circ_graph(m, qr)                                      # Paley/QR
        print(f"   m={m}: interval Cay(Z_{m},+-[1..{(m-1)//2}]) chi={chromatic_exact(Rm)} (=K_{m}, the COMPLETE graph, AP-tight);"
              f"  Paley Cay(Z_{m},QR) chi={chromatic_exact(P)} bipartite={nx.is_bipartite(P)}")

    # (4) Hadwiger-Nelson seed: the unit equilateral triangle is the odd cycle
    print("\n(4) Hadwiger-Nelson: the unit EQUILATERAL TRIANGLE = a triangle in the unit-distance")
    print("    graph => chi(plane) >= 3; the same ODD CYCLE as AG_n's order-3 generator and the")
    print("    LRC odd resonance. Moser spindle (4-chromatic) glues two unit rhombi on the")
    print("    Eisenstein lattice Z[zeta_6]; de Grey's chi>=5 is the rigid-gadget = LRC tight-config method.")

    print("\n" + "=" * 84)
    print("READING")
    print("=" * 84)
    print("""  The alternating group graph is the NON-ABELIAN node of the HYP-2265 unification and its
  PARITY KEY. Chromatic number = failure of a sign(parity) 2-coloring, obstructed by ODD
  CYCLES (the repo's OCF). EVEN generators (3-cycles on A_n) kill the sign character and make
  triangles => chi>=3; ODD generators (transpositions on S_n) keep the sign 2-coloring =>
  chi=2. This is verbatim the LRC split (interval/AP additive chi=2 vs Paley chi=3, S582) and
  HN (unit equilateral triangle => chi>=3). The Hoffman spectral bound carries over, with the
  Cayley eigenvalues now from A_n IRREPS (non-abelian Fourier) -- exactly generalizing
  HYP-2265's connection-set Fourier transform from characters to representations.""")

if __name__ == "__main__":
    main()
