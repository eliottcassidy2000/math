#!/usr/bin/env python3
"""collatz_ck_fixed_points_20260930_dihedral.py -- identify the connected non-circulant C_k-fixed Cayley graphs of the census
(thirteenth note, addendum): Cay(D_6, S) with k = 4 (12 vertices, degree 4), Cay(D_8, S) with k = 4 (16 vertices, degree 4),
Cay(D_8, S) with k = 5 (16 vertices, degree 6).  Compared with the fixed circulants of the same order and degree, with the
Shrikhande graph (the 4 x 4 triangular torus) and the 4 x 4 rook graph, and described by their induced k-cycles.
"""
import itertools
from collatz_ck_fixed_points_20260930 import small_groups, symmetric_subsets, cayley, local_count, cycle_operator, is_isomorphic, circulant, abelian, induced_cycles, is_locally_ck, connected


def fixed_examples(name, degree, k, want=3):
    els, mul, inv, ident = small_groups()[name]
    out = []
    for S in symmetric_subsets(els, mul, inv, ident, degree):
        if len(S) != degree:
            continue
        G = cayley(els, mul, inv, S)
        if not connected(G) or local_count(G, k, k) != k:
            continue
        Ck, cyc = cycle_operator(G, k, limit=len(els))
        if Ck is None or len(cyc) != len(els):
            continue
        if is_isomorphic(Ck, G):
            if not any(is_isomorphic(G, H) for _, H in out):
                out.append((S, G))
            if len(out) >= want:
                break
    return out


def describe(G, k):
    n = len(G)
    cyc = induced_cycles(G, k)
    tri = sum(1 for u in G[0] for w in G[0] if u < w and w in G[u])
    # common neighbours profile (lambda for adjacent, mu for non-adjacent pairs)
    lam = sorted(set(len(G[0] & G[v]) for v in G[0])); mu = sorted(set(len(G[0] & G[v]) for v in range(1, n) if v not in G[0]))
    return "n=%d degree=%d induced %d-cycles=%d triangles/vertex=%d lambda=%s mu=%s locally-C_k=%s" % (n, len(G[0]), k, len(cyc), tri, lam, mu, is_locally_ck(G, k))


shrikhande = cayley(*abelian(4, 4)[:3], [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)])
rook = cayley(*abelian(4, 4)[:3], [(1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3)])
print("Shrikhande (4x4 triangular torus):", describe(shrikhande, 5), "; C_5-fixed:", (lambda C: C is not None and len(C) == 16 and is_isomorphic(C, shrikhande))(cycle_operator(shrikhande, 5, limit=16)[0]))

for name, degree, k, cands in (("D_6", 4, 4, [("C_12(1,4)", circulant(12, [1, 4, 8, 11])), ("C_12(1,5)", circulant(12, [1, 5, 7, 11]))]),
                              ("D_8", 4, 4, [("C_16(1,6)", circulant(16, [1, 6, 10, 15])), ("C_16(2,3)", circulant(16, [2, 3, 13, 14])), ("C_16(1,4)", circulant(16, [1, 4, 12, 15]))]),
                              ("D_8", 6, 5, [("C_16(1,2,3)", circulant(16, [1, 2, 3, 13, 14, 15])), ("C_16(1,2,5)", circulant(16, [1, 2, 5, 11, 14, 15])), ("C_16(1,3,4)", circulant(16, [1, 3, 4, 12, 13, 15])), ("Shrikhande", shrikhande), ("rook 4x4", rook)])):
    ex = fixed_examples(name, degree, k)
    print("== Cay(%s, S), degree %d, k = %d: %d non-isomorphic connected fixed graphs found" % (name, degree, k, len(ex)))
    for S, G in ex:
        print("  ", describe(G, k))
        matches = [nm for nm, H in cands if len(H) == len(G) and is_isomorphic(G, H)]
        print("   isomorphic to:", matches if matches else "none of the candidates")
        if k == 5 and not matches:
            # is it a circulant at all? test against all 6-regular circulants on 16 vertices
            hits = []
            for a, b, c in itertools.combinations(range(1, 8), 3):
                H = circulant(16, [a, b, c, 16 - a, 16 - b, 16 - c])
                if is_isomorphic(G, H):
                    hits.append((a, b, c))
            for a, b in itertools.combinations(range(1, 8), 2):
                H = circulant(16, [a, b, 8, 16 - a, 16 - b])
                if len(H[0]) == 5:
                    continue
            print("   6-regular circulants on 16 vertices isomorphic to it:", hits if hits else "none (a genuinely non-circulant C_5 fixed point)")
