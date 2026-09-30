#!/usr/bin/env python3
"""collatz_pentagon_operator_20260930.py -- the pentagon graph operator C_5 (Gervacio-Maehara-Ramos, arXiv:2604.18984) on the
repo's graphs, and the carry cocycle with its hexagon identities (session collatz-posets-zeta5-20260927, opus, 2026-09-30,
twelfth note).

 (1) C_5(G): vertices = induced 5-cycles of G, adjacent iff they share an edge.  Trajectories (vertex counts) for C_5, K_5,
     K_(3,3), K_6, Petersen, dodecahedron, icosahedron, the Paley graphs P_13, P_17, the Clebsch graph, the Seidel doubles
     D(P_5), D^2(P_5) of the eighth note, the Heawood and Desargues graphs, and the functional graph of the 5x+1 trivial cycle.
     Identifications: C_5(D) = I = C_5(I) (the paper), C_5(Petersen) = ?, and whether C_5 commutes with the antipodal folds
     D -> Petersen, I -> K_6.
 (2) The carry cocycle S_(uv) = 3^(p_v) S_u + 2^(A_u) S_v on the free monoid of words, the clock cocycle
     D_(uv) = 2^(A_u) D_v + 3^(p_v) D_u, the commutation defect beta(u,v) = S_(uv) - S_(vu) = D_u D_v (x_v - x_u), its
     antisymmetry, and the two hexagon identities beta(u, vw) = 3^(p_w) beta(u,v) + 2^(A_v) beta(u,w) and
     beta(uv, w) = 3^(p_v) beta(u,w) + 2^(A_u) beta(v,w); the coboundary reading of cycles (S_w = x_w D_w) and
     H^1(<w>; Z) = Z/(2^A - 3^p); checked on random words and on the splits of the -17 cycle word.
 (3) Catalan: critical spine blocks at the formal slope 2 (Dyck words) against the vertices of the associahedra (1, 2, 5, 14, 42).
Usage: python3 collatz_pentagon_operator_20260930.py
"""
import itertools, random
from fractions import Fraction


# ---------- graphs ----------
def graph_from_edges(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    return adj


def complete(n):
    return graph_from_edges(n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def complete_bipartite(a, b):
    return graph_from_edges(a + b, [(i, a + j) for i in range(a) for j in range(b)])


def cycle(n):
    return graph_from_edges(n, [(i, (i + 1) % n) for i in range(n)])


def petersen():
    e = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    return graph_from_edges(10, e)


def dodecahedron():
    # standard construction: outer pentagon 0-4, middle 5-9 (spokes), inner ring 10-14 and 15-19 pentagon
    e = [(i, (i + 1) % 5) for i in range(5)] + [(i, 5 + i) for i in range(5)]
    e += [(5 + i, 10 + i) for i in range(5)] + [(5 + i, 10 + (i - 1) % 5) for i in range(5)]
    e += [(10 + i, 15 + i) for i in range(5)] + [(15 + i, 15 + (i + 1) % 5) for i in range(5)]
    return graph_from_edges(20, e)


def icosahedron():
    # 0 top, 1-5 upper ring, 6-10 lower ring, 11 bottom
    e = [(0, i) for i in range(1, 6)] + [(i, 1 + i % 5) for i in range(1, 6)]
    e += [(11, i) for i in range(6, 11)] + [(6 + i, 6 + (i + 1) % 5) for i in range(5)]
    e += [(1 + i, 6 + i) for i in range(5)] + [(1 + i, 6 + (i + 1) % 5) for i in range(5)]
    return graph_from_edges(12, e)


def paley(q):
    sq = {(x * x) % q for x in range(1, q)}
    return graph_from_edges(q, [(i, j) for i in range(q) for j in range(i + 1, q) if (j - i) % q in sq])


def clebsch():
    # vertices: 5-bit vectors of even weight (16); adjacent iff they differ in exactly 4 bits ... use the standard
    # folded 5-cube: vertices 0..15 as 4-bit vectors, adjacent iff differ in 1 bit or are complements (differ in 4 bits)
    e = [(i, j) for i in range(16) for j in range(i + 1, 16) if bin(i ^ j).count("1") in (1, 4)]
    return graph_from_edges(16, e)


def heawood():
    e = []
    for i in range(14):
        e.append((i, (i + 1) % 14))
        if i % 2 == 0:
            e.append((i, (i + 5) % 14))
    return graph_from_edges(14, e)


def desargues():
    # generalized Petersen GP(10, 3)
    e = [(i, (i + 1) % 10) for i in range(10)] + [(i, 10 + i) for i in range(10)] + [(10 + i, 10 + (i + 3) % 10) for i in range(10)]
    return graph_from_edges(20, e)


def seidel_double_graph(adj):
    """Seidel matrix S (0 diag, -1 adjacent, +1 non-adjacent); S' = [[S, S+I],[S+I, -S]]; back to a graph."""
    n = len(adj)
    S = [[0 if i == j else (-1 if j in adj[i] else 1) for j in range(n)] for i in range(n)]
    N = 2 * n
    S2 = [[0] * N for _ in range(N)]
    for i in range(n):
        for j in range(n):
            S2[i][j] = S[i][j]
            S2[i][j + n] = S[i][j] + (1 if i == j else 0)
            S2[i + n][j] = S[i][j] + (1 if i == j else 0)
            S2[i + n][j + n] = -S[i][j]
    return graph_from_edges(N, [(i, j) for i in range(N) for j in range(i + 1, N) if S2[i][j] == -1])


class TooMany(Exception):
    pass


def induced_pentagons(adj, limit=60000):
    n = len(adj); out = set()
    for a in range(n):
        if len(out) > limit:
            raise TooMany(len(out))
        for b in adj[a]:
            if b < a:
                continue
            for c in adj[b]:
                if c <= a or c == b or c in adj[a]:
                    continue
                for d in adj[c]:
                    if d <= a or d in (b, c) or d in adj[a] or d in adj[b]:
                        continue
                    for e in adj[d]:
                        if e <= a or e in (b, c, d) or e not in adj[a] or e in adj[b] or e in adj[c]:
                            continue
                        # induced: a-b-c-d-e-a with no chords: checked ac, ad, bd, be, ce
                        cyc = (a, b, c, d, e)
                        # canonical: start at a (min), the two directions: (a,b,c,d,e) and (a,e,d,c,b)
                        key = min(cyc, (a, e, d, c, b))
                        out.add(key)
    return sorted(out)


def pentagon_graph(adj):
    pents = induced_pentagons(adj)
    edges_of = []
    for p in pents:
        es = set()
        for i in range(5):
            x, y = p[i], p[(i + 1) % 5]
            es.add((min(x, y), max(x, y)))
        edges_of.append(es)
    m = len(pents)
    e = [(i, j) for i in range(m) for j in range(i + 1, m) if edges_of[i] & edges_of[j]]
    return graph_from_edges(m, e), pents


def degree_sequence(adj):
    return sorted(len(a) for a in adj)


def is_iso_small(adj1, adj2):
    """brute-force isomorphism for very small graphs (n <= 12) via degree-refined backtracking."""
    n = len(adj1)
    if n != len(adj2) or degree_sequence(adj1) != degree_sequence(adj2):
        return False
    d1 = [len(a) for a in adj1]; d2 = [len(a) for a in adj2]
    order = sorted(range(n), key=lambda v: -d1[v])
    mapping = {}
    used = set()

    def bt(k):
        if k == n:
            return True
        v = order[k]
        for w in range(n):
            if w in used or d2[w] != d1[v]:
                continue
            ok = all(((u in adj1[v]) == (mapping[u] in adj2[w])) for u in mapping)
            if ok:
                mapping[v] = w; used.add(w)
                if bt(k + 1):
                    return True
                del mapping[v]; used.discard(w)
        return False
    return bt(0)


def trajectory(name, adj, steps=5, cap=60):
    """iterate C_5 while the current graph has at most `cap` vertices and at most 60000 induced pentagons (the pentagon
    graphs are pentagon-rich and dense; C_5(P_13) on 39 vertices already has 4563 induced pentagons)."""
    sizes = [len(adj)]; cur = adj; graphs = [adj]
    for _ in range(steps):
        if len(cur) == 0:
            break
        if len(cur) > cap:
            sizes.append("(>%d vertices: not iterated)" % cap); break
        try:
            nxt, _ = pentagon_graph(cur)
        except TooMany as e:
            sizes.append("(>%s induced pentagons: not iterated)" % e); break
        sizes.append(len(nxt)); graphs.append(nxt); cur = nxt
    return sizes, graphs


def part1():
    print("== (1) the pentagon graph operator on the repo's graphs ==")
    I = icosahedron(); D = dodecahedron(); P = petersen()
    tests = [("C_5 (= Paley P_5)", cycle(5)), ("K_5", complete(5)), ("K_(3,3)", complete_bipartite(3, 3)), ("K_6", complete(6)),
             ("Petersen (= D/+-1)", P), ("dodecahedron", D), ("icosahedron", I), ("Paley P_13", paley(13)), ("Paley P_17", paley(17)),
             ("Clebsch", clebsch()), ("Seidel double D(P_5)", seidel_double_graph(cycle(5))), ("Seidel double D^2(P_5)", seidel_double_graph(seidel_double_graph(cycle(5)))),
             ("Heawood", heawood()), ("Desargues GP(10,3)", desargues()), ("5x+1 trivial cycle 1-3-8-4-2", cycle(5))]
    for name, g in tests:
        sizes, graphs = trajectory(name, g)
        edges = [sum(len(a) for a in h) // 2 for h in graphs]
        print(" %-32s |V(C_5^k)| = %s   edges %s" % (name, sizes, edges), flush=True)
    C5D, _ = pentagon_graph(D); C5I, _ = pentagon_graph(I); C5P, pentsP = pentagon_graph(P)
    print(" C_5(D) iso icosahedron: %s; C_5(I) iso icosahedron: %s" % (is_iso_small(C5D, I), is_iso_small(C5I, I)))
    print(" Petersen: %d induced pentagons; C_5(Petersen) has %d vertices, degree sequence %s" % (len(pentsP), len(C5P), degree_sequence(C5P)))
    # structure of C_5(Petersen): complement's components
    m = len(C5P)
    comp_adj = [set(j for j in range(m) if j != i and j not in C5P[i]) for i in range(m)]
    print("  complement of C_5(Petersen): degree sequence %s" % degree_sequence(comp_adj))
    print("  C_5(Petersen) iso icosahedron: %s; iso K_12 minus a perfect matching: %s" % (is_iso_small(C5P, I), all(len(a) == 1 for a in comp_adj)))
    # which pentagons of Petersen are faces of the hemi-dodecahedron (lift to faces of D) vs Petrie: count shared-edge multiplicities
    shared = {}
    edges_of = [set((min(p[i], p[(i + 1) % 5]), max(p[i], p[(i + 1) % 5])) for i in range(5)) for p in pentsP]
    for i in range(m):
        for j in range(i + 1, m):
            k = len(edges_of[i] & edges_of[j]); shared[k] = shared.get(k, 0) + 1
    print("  pairs of Petersen pentagons by number of shared edges: %s" % dict(sorted(shared.items())))
    C5C5P, _ = pentagon_graph(C5P)
    print("  C_5^2(Petersen) has %d vertices -> Petersen is pentagon-vanishing in %d steps" % (len(C5C5P), 2 if len(C5C5P) == 0 else -1))
    # the fold check: C_5 of the folds vs folds of C_5
    print(" folds: D -> Petersen (C_5: I vs %d-vertex graph), I -> K_6 (C_5: I vs empty): the operator sees the double cover, not the fold" % len(C5P))


# ---------- the carry cocycle ----------
def carry(w):
    S = 0; d = 0
    for v in w:
        S = 3 * S + 2 ** d; d += v
    return S


def clock(w):
    return 2 ** sum(w) - 3 ** len(w)


def fixed_point(w):
    return Fraction(carry(w), clock(w))


def part2():
    print("== (2) the carry cocycle and its hexagon identities ==")
    random.seed(7)
    ok_cocycle = ok_clock = ok_beta = ok_anti = ok_hex1 = ok_hex2 = True
    for _ in range(400):
        u = tuple(random.randint(1, 4) for _ in range(random.randint(1, 5)))
        v = tuple(random.randint(1, 4) for _ in range(random.randint(1, 5)))
        w = tuple(random.randint(1, 4) for _ in range(random.randint(1, 5)))
        pu, Au = len(u), sum(u); pv, Av = len(v), sum(v); pw, Aw = len(w), sum(w)
        ok_cocycle &= carry(u + v) == 3 ** pv * carry(u) + 2 ** Au * carry(v)
        ok_clock &= clock(u + v) == 2 ** Au * clock(v) + 3 ** pv * clock(u)
        beta = lambda a, b: carry(a + b) - carry(b + a)
        ok_beta &= beta(u, v) == clock(u) * clock(v) * (fixed_point(v) - fixed_point(u))
        ok_anti &= beta(u, v) == -beta(v, u)
        ok_hex1 &= beta(u, v + w) == 3 ** pw * beta(u, v) + 2 ** Av * beta(u, w)
        ok_hex2 &= beta(u + v, w) == 3 ** pv * beta(u, w) + 2 ** Au * beta(v, w)
    print(" cocycle S_(uv) = 3^(p_v) S_u + 2^(A_u) S_v: %s; clock cocycle D_(uv) = 2^(A_u) D_v + 3^(p_v) D_u: %s" % (ok_cocycle, ok_clock))
    print(" commutation defect beta(u,v) = D_u D_v (x_v - x_u): %s; antisymmetric: %s" % (ok_beta, ok_anti))
    print(" hexagon 1  beta(u, vw) = 3^(p_w) beta(u,v) + 2^(A_v) beta(u,w): %s" % ok_hex1)
    print(" hexagon 2  beta(uv, w) = 3^(p_v) beta(u,w) + 2^(A_u) beta(v,w): %s" % ok_hex2)
    w17 = (1, 1, 1, 2, 1, 1, 4)
    print(" the -17 word: fixed point %s, carry %d = x_w * D_w (coboundary over Z): %s" % (fixed_point(w17), carry(w17), carry(w17) == fixed_point(w17) * clock(w17)))
    rows = []
    for k in range(1, 7):
        u, v = w17[:k], w17[k:]
        b = carry(u + v) - carry(v + u)
        rows.append((k, str(fixed_point(u)), str(fixed_point(v)), b, b % 139 == 0))
    print(" splits u|v of the -17 word: (k, x_u, x_v, beta(u,v), 139 | beta):", rows)
    # H^1(<w>; M) = M/(3^p - 2^A) M for M = Z, Z_2, Z_3: the class of the carry
    for w in ((1, 2), (1, 1, 1, 2, 1, 1, 4), (2, 2, 1), (1, 3)):
        D = abs(clock(w)); print("  word %s: clock %d, H^1(<w>; Z) = Z/%d, class of the carry = %d mod %d -> %s" % (w, clock(w), D, carry(w) % D if D > 1 else 0, D, "coboundary: integer cycle" if carry(w) % D == 0 else "not a coboundary: no integer cycle (2-adically and 3-adically the class vanishes: the clock is a unit there)"))


def part3():
    print("== (3) Catalan: critical spine blocks at the formal slope 2 (Dyck words) and the associahedra ==")
    def dyck_blocks(k):
        # words (v_1..v_k) in {1,2,3}: steps v - 2 in {-1,0,+1}? use the seventh note's convention: a block is a rise (v=1 -> step +1)
        # followed by a Dyck path in steps v-2 ... here simply count Dyck words of semilength k (Catalan) by DP
        c = [1] + [0] * k
        for n in range(1, k + 1):
            c[n] = sum(c[i] * c[n - 1 - i] for i in range(n))
        return c[k]
    print(" Catalan numbers C_k = %s: vertices of the associahedra K_(k+2) (K_4 is the pentagon of Mac Lane's axiom, 5 bracketings of 4 letters; K_5 has 14)" % [dyck_blocks(k) for k in range(1, 7)])
    print(" the seventh note counted the critical spine blocks at q = 4 by the same numbers: the pentagon axiom's five bracketings are in bijection with the five critical blocks of length 3 (Dyck paths of semilength 3); a bijection, not a mechanism")


def induced_cycles_k(adj, k):
    """induced k-cycles (k = 3 or 4) as canonical tuples."""
    n = len(adj); out = set()
    if k == 3:
        for a in range(n):
            for b in adj[a]:
                for c in adj[b]:
                    if a < b < c and c in adj[a]:
                        out.add((a, b, c))
    elif k == 4:
        for a in range(n):
            for b in adj[a]:
                if b < a:
                    continue
                for c in adj[b]:
                    if c <= a or c == b or c in adj[a]:
                        continue
                    for d in adj[c]:
                        if d <= a or d in (b, c) or d not in adj[a] or d in adj[b]:
                            continue
                        out.add(min((a, b, c, d), (a, d, c, b)))
    return sorted(out)


def cycle_graph_k(adj, k):
    cyc = induced_cycles_k(adj, k)
    edges_of = [set((min(c[i], c[(i + 1) % k]), max(c[i], c[(i + 1) % k])) for i in range(k)) for c in cyc]
    m = len(cyc)
    return graph_from_edges(m, [(i, j) for i in range(m) for j in range(i + 1, m) if edges_of[i] & edges_of[j]])


def octahedron():
    return graph_from_edges(6, [(i, j) for i in range(6) for j in range(i + 1, 6) if j != i + 3 or i >= 3])


def cube():
    return graph_from_edges(8, [(i, j) for i in range(8) for j in range(i + 1, 8) if bin(i ^ j).count("1") == 1])


def part4():
    print("== (4) the Platonic solids under the induced-cycle operators C_3, C_4, C_5 ==")
    solids = {"tetrahedron K_4": complete(4), "octahedron": octahedron(), "cube": cube(), "icosahedron": icosahedron(), "dodecahedron": dodecahedron()}
    names = list(solids)
    def ident(g):
        for nm, h in solids.items():
            if len(h) == len(g) and is_iso_small(g, h):
                return nm
        if len(g) == 0:
            return "empty"
        if all(len(a) == 0 for a in g):
            return "%d K_1" % len(g)
        return "%d vertices, degrees %s" % (len(g), sorted(set(len(a) for a in g)))
    for nm in names:
        g = solids[nm]
        row = []
        for k in (3, 4, 5):
            h = cycle_graph_k(g, k) if k < 5 else pentagon_graph(g)[0]
            row.append("C_%d -> %s" % (k, ident(h)))
        print(" %-16s %s" % (nm, "; ".join(row)))
    print(" fixed points among the Platonic graphs: C_3(K_4) = K_4 (the self-dual tetrahedron), C_5(I) = I (the vertex links of the icosahedron are its twelve pentagons); the dual pairs appear as operator steps C_k(P) = P* when the faces are induced k-gons; the octahedron/cube pair vanishes under C_3 and C_4")
    print(" clocks mod 6: 2^A - 3^p = (-1)^A mod 6 for all A, p >= 1: %s (every clock is a Syracuse-core residue +-1 mod 6, the sign being the parity of the halving count)" % all((2 ** A - 3 ** p - (-1) ** A) % 6 == 0 for A in range(1, 40) for p in range(1, 40)))


if __name__ == "__main__":
    part1(); part2(); part3(); part4()
