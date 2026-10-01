#!/usr/bin/env python3
"""procgen_petersen_20261001_orchestrator_check.py -- the orchestrator's independent audit of the petersen lane (2026-10-01).

Written from the statements in the lane note; the lane's code was not read. Exact arc-HP counts come from the orchestrator's
own C engine (procgen_petersen_20261001_orchestrator_check.c), compiled into a temp dir.

Checks:
  A. QR_127 restricted to mu_14 = +-<2>: a tournament, H = 24540117, all 91 arc counts odd, |Aut| = 7 (multipliers <2>).
  B. Anti-circulant tournaments T_s on Z/2m (m = 3, 5, 7): every sign pattern gives a tournament with x -> x+1 an
     anti-automorphism; Theorem 6.1 (antipodal arcs odd) on every pattern; the number of all-odd isomorphism classes is 1
     for N = 6, 10, 14; QR_q - 0 (q = 7, 11) and QR_127[mu_14] are anti-circulant.
  C. QR_p[mu_6] all-odd iff p = 7 (mod 8); QR_p[mu_10] all-odd iff p = 3 (mod 8) (primes p < 3000 with 2m | p-1, p = 3 mod 4).
  D. Theorem 5.1: random Cayley digraphs (any connection set) of odd-order abelian groups have every arc on an even
     number of Hamiltonian paths.
  E. Prop. 3.2: a non-complete graph has an orientation with H = 0 (construction); complete graphs: H odd (Redei; spot).
  F. Petersen family: the Delta-Y / Y-Delta closure of K6 has 7 graphs, all with 15 edges, Aut orders
     {720, 36, 72, 72, 8, 12, 120}, Delta-Y Hasse diagram a tree with 6 arrows; the Heawood family (closure of K7) has 20
     graphs, 14 of them Delta-Y descendants of K7.
  G. Paley coordinates: Delta-Y on the 7 translates of {1,2,4} in K7 = Heawood graph (Fano incidence graph); Delta-Y on the
     4 translates avoiding 0 in K7 - 0 = Petersen graph, the 3 translates through 0 leave the matching {u, 3u};
     Delta-Y on {1,2,4} and {3,5,6} in K7 - 0 = K_{4,4} - e; Aut(P7 - 0) = <x -> 2x> (order 3).
  H. Theorem L (linear K6): for random integer configurations of 6 points in R^3 (general position), the number of linked
     complementary triangle pairs is 1 or 3; the Gale tournament T* (exact rational Gale vectors) is locally transitive;
     #linked = 3 iff T* = C3[TT2,TT2,TT2]; the piercing rule (number of B-edges piercing A = number of uniform-opposite
     b in B) holds for every split.
"""
import itertools
import math
import os
import random
import re
import subprocess
import sys
import tempfile
import time
from fractions import Fraction

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
OKS = []


def ok(cond, msg):
    OKS.append(bool(cond))
    print(('[OK] ' if cond else '[FAIL] ') + msg, flush=True)


ENGINE = None


def arcs(out):
    """exact arc-HP counts of a digraph given as list of out-neighbour sets"""
    N = len(out)
    masks = ['%x' % sum(1 << j for j in out[i]) for i in range(N)]
    r = subprocess.run([ENGINE, str(N)] + masks, capture_output=True, text=True, check=True).stdout.split('\n')
    H = int(r[0].split()[1])
    c = {}
    for line in r[1:]:
        if line.strip():
            u, v, x = map(int, line.split())
            c[(u, v)] = x
    return H, c


def is_tournament(out):
    N = len(out)
    return all((j in out[i]) != (i in out[j]) for i in range(N) for j in range(N) if i != j)


def qr_set(p):
    return {(x * x) % p for x in range(1, p)}


def restricted_paley(p, verts):
    Q = qr_set(p)
    idx = {v: i for i, v in enumerate(verts)}
    return [[idx[w] for w in verts if w != v and (w - v) % p in Q] for v in verts]


def tournament_canon(out):
    """canonical form of a tournament by brute force (small N only)"""
    N = len(out)
    adj = [[j in out[i] for j in range(N)] for i in range(N)]
    best = None
    for perm in itertools.permutations(range(N)):
        key = tuple(adj[perm[i]][perm[j]] for i in range(N) for j in range(N))
        if best is None or key < best:
            best = key
    return best


def aut_count_tournament(out):
    G = nx.DiGraph()
    G.add_nodes_from(range(len(out)))
    for i, s in enumerate(out):
        for j in s:
            G.add_edge(i, j)
    return sum(1 for _ in nx.algorithms.isomorphism.DiGraphMatcher(G, G).isomorphisms_iter())


def check_A():
    p = 127
    mu = sorted({pow(2, k, p) for k in range(7)} | {(-pow(2, k, p)) % p for k in range(7)})
    out = restricted_paley(p, mu)
    ok(len(mu) == 14 and is_tournament(out), 'A: QR_127 on mu_14 = +-<2> is a tournament on 14 vertices')
    H, c = arcs(out)
    vals = sorted(set(c.values()))
    ok(H == 24540117, f'A: H = {H} (lane: 24540117)')
    ok(len(c) == 91 and all(x % 2 == 1 for x in c.values()), f'A: all 91 arc counts odd; distinct values {vals}')
    ok(vals == [3085307, 3125731, 3361361, 3510057, 3604781, 4056367, 4087295], 'A: the 7 values match the lane')
    a = aut_count_tournament(out)
    ok(a == 7, f'A: |Aut| = {a}')


def anti_tournament(m, s):
    """T_s on Z/2m: x -> y iff s[(y-x) mod 2m] == (-1)^x"""
    N = 2 * m
    return [[y for y in range(N) if y != x and s[(y - x) % N] == (1 if x % 2 == 0 else -1)] for x in range(N)]


def all_patterns(m):
    N = 2 * m
    reps = [d for d in range(1, N) if d < N - d] + [m]
    reps = sorted(set(reps))
    for bits in itertools.product([1, -1], repeat=len(reps)):
        s = [0] * N
        for d, b in zip(reps, bits):
            s[d] = b
            if d != m:
                s[N - d] = -((-1) ** d) * b
        # check constraint s(-d) = -(-1)^d s(d) for all d
        assert all(s[(N - d) % N] == -((-1) ** d) * s[d] for d in range(1, N))
        yield s


def check_B():
    for m in (3, 5, 7):
        N = 2 * m
        pats = list(all_patterns(m))
        good_t = good_anti = thm61 = 0
        odd_classes = {}
        for s in pats:
            out = anti_tournament(m, s)
            if is_tournament(out):
                good_t += 1
            # x -> x+1 anti-automorphism: x->y iff (y+1) -> (x+1)
            if all((y in out[x]) == (((x + 1) % N) in out[(y + 1) % N]) for x in range(N) for y in range(N) if x != y):
                good_anti += 1
            H, c = arcs(out)
            if all(c[(x, (x + m) % N)] % 2 == 1 for x in range(N) if (x, (x + m) % N) in c):
                thm61 += 1
            if all(v % 2 == 1 for v in c.values()):
                key = (H, tuple(sorted(c.values())))
                odd_classes.setdefault(key, s)
        ok(good_t == len(pats) == 2 ** m, f'B: N={N}: all {len(pats)} sign patterns give tournaments')
        ok(good_anti == len(pats), f'B: N={N}: x -> x+1 is an anti-automorphism of every T_s')
        ok(thm61 == len(pats), f'B: N={N}: Theorem 6.1 (antipodal arcs odd) on every pattern')
        # isomorphism classes of all-odd: group by invariant, then confirm isomorphism inside each group via networkx
        classes = []
        for key, s in odd_classes.items():
            out = anti_tournament(m, s)
            G = nx.DiGraph([(i, j) for i in range(N) for j in out[i]])
            if not any(nx.is_isomorphic(G, G2) for G2 in classes):
                classes.append(G)
        ok(len(classes) == 1, f'B: N={N}: all-odd anti-circulant isomorphism classes = {len(classes)} (lane: 1)')
    # QR_q - 0 and QR_127[mu14] are anti-circulant: multiplication by a non-residue generator is a cyclic anti-automorphism
    for q in (7, 11, 19):
        g = next(x for x in range(2, q) if len({pow(x, k, q) for k in range(q - 1)}) == q - 1)
        verts = [pow(g, k, q) for k in range(q - 1)]
        out = restricted_paley(q, verts)
        N = q - 1
        anti = all((y in out[x]) == (((x + 1) % N) in out[(y + 1) % N]) for x in range(N) for y in range(N) if x != y)
        ok(anti, f'B: QR_{q} - 0 in generator order (g = {g}) is anti-circulant (x -> x+1 reverses every arc)')
    p = 127
    zeta = next(z for z in range(2, p) if pow(z, 14, p) == 1 and all(pow(z, k, p) != 1 for k in (1, 2, 7)))
    verts = [pow(zeta, k, p) for k in range(14)]
    out = restricted_paley(p, verts)
    anti = all((y in out[x]) == (((x + 1) % 14) in out[(y + 1) % 14]) for x in range(14) for y in range(14) if x != y)
    ok(sorted(verts) == sorted({pow(2, k, p) for k in range(7)} | {(-pow(2, k, p)) % p for k in range(7)}) and anti,
       f'B: QR_127[mu_14] in the order of a primitive 14th root (zeta = {zeta}) is anti-circulant')


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def check_C():
    for m, residue in ((3, 7), (5, 3)):
        N = 2 * m
        tested = mism = 0
        for p in range(3, 3000):
            if not is_prime(p) or p % 4 != 3 or (p - 1) % N != 0:
                continue
            zeta = next(z for z in range(2, p) if pow(z, N, p) == 1 and all(pow(z, N // r, p) != 1 for r in (2, m) if N % r == 0)
                        and all(pow(z, d, p) != 1 for d in range(1, N)))
            verts = [pow(zeta, k, p) for k in range(N)]
            out = restricted_paley(p, verts)
            assert is_tournament(out)
            H, c = arcs(out)
            allodd = all(v % 2 == 1 for v in c.values())
            tested += 1
            if allodd != (p % 8 == residue):
                mism += 1
        ok(tested > 20 and mism == 0, f'C: QR_p[mu_{N}] all-odd iff p = {residue} (mod 8): {tested} primes < 3000, {mism} mismatches')


def check_D():
    rng = random.Random(20261001)
    groups = [('Z9', [(i,) for i in range(9)], lambda a, b: ((a[0] - b[0]) % 9,)),
              ('Z15', [(i,) for i in range(15)], lambda a, b: ((a[0] - b[0]) % 15,)),
              ('Z3xZ3', [(i, j) for i in range(3) for j in range(3)], lambda a, b: ((a[0] - b[0]) % 3, (a[1] - b[1]) % 3)),
              ('Z5x Z3', [(i, j) for i in range(5) for j in range(3)], lambda a, b: ((a[0] - b[0]) % 5, (a[1] - b[1]) % 3)),
              ('Z7', [(i,) for i in range(7)], lambda a, b: ((a[0] - b[0]) % 7,)),
              ('Z11', [(i,) for i in range(11)], lambda a, b: ((a[0] - b[0]) % 11,))]
    total = bad = nontour = 0
    for name, G, diff in groups:
        zero = tuple(0 for _ in G[0])
        nonzero = [g for g in G if g != zero]
        for trial in range(6):
            C = {g for g in nonzero if rng.random() < 0.5}
            if not C:
                continue
            out = [[j for j, y in enumerate(G) if j != i and diff(y, x) in C] for i, x in enumerate(G)]
            if not is_tournament(out):
                nontour += 1
            H, c = arcs(out)
            total += 1
            if any(v % 2 for v in c.values()):
                bad += 1
    ok(total >= 30 and bad == 0 and nontour > 20,
       f'D: Theorem 5.1 on {total} random Cayley digraphs ({nontour} non-tournaments): every arc even ({bad} failures)')


def check_E():
    # Prop 3.2: non-complete graph -> orientation with H = 0 (u, v non-adjacent, all edges at u, v oriented away)
    rng = random.Random(7)
    fails = 0
    for trial in range(40):
        n = rng.randint(4, 9)
        edges = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.8]
        E = set(edges)
        non = [(i, j) for i in range(n) for j in range(i + 1, n) if (i, j) not in E]
        if not non:
            continue
        u, v = non[0]
        out = [set() for _ in range(n)]
        for (a, b) in edges:
            if a in (u, v):
                out[a].add(b)
            elif b in (u, v):
                out[b].add(a)
            elif rng.random() < 0.5:
                out[a].add(b)
            else:
                out[b].add(a)
        H, c = arcs([sorted(s) for s in out])
        if H != 0:
            fails += 1
    ok(fails == 0, 'E: Prop 3.2 construction gives H = 0 on 40 random non-complete graphs')


def dy(G, tri, newnode):
    a, b, c = tri
    H = G.copy()
    H.remove_edges_from([(a, b), (b, c), (a, c)])
    H.add_node(newnode)
    H.add_edges_from([(newnode, a), (newnode, b), (newnode, c)])
    return H


def yd(G, y):
    nb = list(G.neighbors(y))
    H = G.copy()
    H.remove_node(y)
    H.add_edges_from([(nb[0], nb[1]), (nb[1], nb[2]), (nb[0], nb[2])])
    return H


def closure(G0):
    """Delta-Y / Y-Delta closure (simple graphs; moves that would create multi-edges are skipped)"""
    fam = [G0]
    frontier = [G0]
    arrows = set()
    def index(G):
        for i, F in enumerate(fam):
            if F.number_of_nodes() == G.number_of_nodes() and nx.is_isomorphic(F, G):
                return i
        return None
    while frontier:
        new = []
        for G in frontier:
            gi = index(G)
            for tri in [t for t in itertools.combinations(G.nodes(), 3) if G.has_edge(t[0], t[1]) and G.has_edge(t[1], t[2]) and G.has_edge(t[0], t[2])]:
                H = dy(G, tri, max(G.nodes()) + 1)
                hi = index(H)
                if hi is None:
                    fam.append(H); new.append(H); hi = len(fam) - 1
                arrows.add((gi, hi))
            for y in [v for v in G.nodes() if G.degree(v) == 3]:
                nb = list(G.neighbors(y))
                if any(G.has_edge(nb[i], nb[j]) for i in range(3) for j in range(i + 1, 3)):
                    continue  # would create a multi-edge
                H = nx.convert_node_labels_to_integers(yd(G, y))
                hi = index(H)
                if hi is None:
                    fam.append(H); new.append(H); hi = len(fam) - 1
                arrows.add((hi, gi))
        frontier = new
    return fam, arrows


def aut_order(G):
    return sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())


def check_F():
    fam, arrows = closure(nx.complete_graph(6))
    ok(len(fam) == 7 and all(G.number_of_edges() == 15 for G in fam), f'F: closure of K6 has {len(fam)} graphs, all with 15 edges')
    auts = sorted(aut_order(G) for G in fam)
    ok(auts == sorted([720, 36, 72, 72, 8, 12, 120]), f'F: Aut orders {auts}')
    U = nx.Graph()
    U.add_nodes_from(range(len(fam)))
    U.add_edges_from(arrows)
    ok(len(arrows) == 6 and nx.is_tree(U), f'F: Delta-Y Hasse diagram has {len(arrows)} arrows and is a tree')
    pet = nx.petersen_graph()
    ok(any(nx.is_isomorphic(G, pet) for G in fam), 'F: the Petersen graph is in the family')
    k331 = nx.complete_multipartite_graph(3, 3, 1)
    # Delta-Y descendants of K6
    desc = {0}
    changed = True
    while changed:
        changed = False
        for (a, b) in arrows:
            if a in desc and b not in desc:
                desc.add(b); changed = True
    idx331 = next(i for i, G in enumerate(fam) if nx.is_isomorphic(G, k331))
    ok(len(desc) == 6 and idx331 not in desc, 'F: every member except K3,3,1 is a Delta-Y descendant of K6')
    famH, arrowsH = closure(nx.complete_graph(7))
    descH = {0}
    changed = True
    while changed:
        changed = False
        for (a, b) in arrowsH:
            if a in descH and b not in descH:
                descH.add(b); changed = True
    ok(len(famH) == 20 and len(descH) == 14, f'F: Heawood family (closure of K7) has {len(famH)} graphs, {len(descH)} Delta-Y descendants of K7')


def check_G():
    lines = [tuple(sorted(((x + t) % 7) for x in (1, 2, 4))) for t in range(7)]
    G = nx.complete_graph(7)
    nxt = 7
    for L in lines:
        G = dy(G, L, nxt); nxt += 1
    heawood = nx.heawood_graph()
    ok(nx.is_isomorphic(G, heawood) and nx.is_bipartite(G), 'G: Delta-Y on the 7 translates of {1,2,4} in K7 is the Heawood graph')
    K6 = nx.complete_graph(range(1, 7))
    avoid = [L for L in lines if 0 not in L]
    through = [tuple(x for x in L if x != 0) for L in lines if 0 in L]
    ok(len(avoid) == 4 and sorted(tuple(sorted(e)) for e in through) == sorted([(1, 3), (2, 6), (4, 5)]) and
       all((3 * e[0]) % 7 == e[1] or (3 * e[1]) % 7 == e[0] for e in through),
       'G: 4 translates avoid 0; the 3 through 0 leave the matching {1,3},{2,6},{4,5} = {u, 3u}')
    P = K6.copy()
    nxt = 7
    for L in avoid:
        P = dy(P, L, nxt); nxt += 1
    ok(nx.is_isomorphic(P, nx.petersen_graph()), 'G: Delta-Y on the 4 translates avoiding 0 in K6 = K7 - 0 is the Petersen graph')
    K = K6.copy()
    K = dy(K, (1, 2, 4), 7)
    K = dy(K, (3, 5, 6), 8)
    k44e = nx.complete_bipartite_graph(4, 4)
    k44e.remove_edge(0, 4)
    ok(nx.is_isomorphic(K, k44e), 'G: Delta-Y on {1,2,4} and {3,5,6} in K6 is K_{4,4} - e')
    verts = [1, 2, 3, 4, 5, 6]
    out = restricted_paley(7, verts)
    a = aut_count_tournament(out)
    two = all(((verts.index((2 * verts[w]) % 7)) in out[verts.index((2 * verts[v]) % 7)]) == (w in out[v])
              for v in range(6) for w in range(6) if v != w)
    ok(a == 3 and two, f'G: |Aut(P7 - 0)| = {a}, and x -> 2x is an automorphism')


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))


def orient(p, q, r, s):
    return det3([q[i] - p[i] for i in range(3)], [r[i] - p[i] for i in range(3)], [s[i] - p[i] for i in range(3)])


def seg_pierces_tri(a, b, t0, t1, t2):
    """does the open segment ab cross the triangle t0 t1 t2 (general position)"""
    s1 = orient(t0, t1, t2, a)
    s2 = orient(t0, t1, t2, b)
    if s1 * s2 >= 0:
        return False
    o0 = orient(a, b, t0, t1)
    o1 = orient(a, b, t1, t2)
    o2 = orient(a, b, t2, t0)
    return (o0 > 0 and o1 > 0 and o2 > 0) or (o0 < 0 and o1 < 0 and o2 < 0)


def gale_vectors(P):
    """basis of {lam : sum lam_i = 0, sum lam_i p_i = 0} as a 6x2 rational matrix"""
    rows = [[Fraction(1)] * 6] + [[Fraction(P[i][k]) for i in range(6)] for k in range(3)]
    # reduce rows to RREF
    M = [r[:] for r in rows]
    piv = []
    r = 0
    for c in range(6):
        pr = next((i for i in range(r, 4) if M[i][c] != 0), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(4):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(6) if c not in piv]
    basis = []
    for fc in free:
        v = [Fraction(0)] * 6
        v[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            v[pc] = -M[i][fc]
        basis.append(v)
    assert len(basis) == 2
    return [(basis[0][i], basis[1][i]) for i in range(6)]


def c3_tt2_canon():
    # blocks {0,1},{2,3},{4,5}: 0->1, 2->3, 4->5; block i beats block i+1
    out = [set() for _ in range(6)]
    for b in range(3):
        a0, a1 = 2 * b, 2 * b + 1
        out[a0].add(a1)
        nb = (b + 1) % 3
        for x in (a0, a1):
            for y in (2 * nb, 2 * nb + 1):
                out[x].add(y)
    return tournament_canon([sorted(s) for s in out])


def locally_transitive(out):
    N = len(out)
    adj = [[j in out[i] for j in range(N)] for i in range(N)]
    def transitive(S):
        S = list(S)
        return all(not (adj[a][b] and adj[b][c] and adj[c][a]) for a in S for b in S for c in S if len({a, b, c}) == 3)
    return all(transitive(out[v]) and transitive([u for u in range(N) if v in out[u]]) for v in range(N))


def check_H():
    rng = random.Random(4524)
    C3T = c3_tt2_canon()
    trials = gp = 0
    counts = {}
    bad_count = bad_iff = bad_pierce = bad_lt = 0
    while gp < 1500:
        P = [tuple(rng.randint(-60, 60) for _ in range(3)) for _ in range(6)]
        trials += 1
        if any(orient(*[P[i] for i in q]) == 0 for q in itertools.combinations(range(6), 4)):
            continue
        gp += 1
        g = gale_vectors(P)
        out = [[j for j in range(6) if j != i and g[i][0] * g[j][1] - g[i][1] * g[j][0] > 0] for i in range(6)]
        if not locally_transitive(out):
            bad_lt += 1
        linked = 0
        for A in itertools.combinations(range(6), 3):
            if 0 not in A:
                continue
            B = tuple(x for x in range(6) if x not in A)
            pierce = sum(1 for (x, y) in itertools.combinations(B, 2) if seg_pierces_tri(P[x], P[y], P[A[0]], P[A[1]], P[A[2]]))
            if pierce == 1:
                linked += 1
            # piercing rule: number of b in B uniform-opposite in T*
            uo = 0
            for b in B:
                rest = [x for x in B if x != b]
                beats_A = all(a in out[b] for a in A)
                A_beats = all(b in out[a] for a in A)
                rest_to_b = all(b in out[x] for x in rest)
                b_to_rest = all(x in out[b] for x in rest)
                if (beats_A and rest_to_b) or (A_beats and b_to_rest):
                    uo += 1
            if uo != pierce:
                bad_pierce += 1
        counts[linked] = counts.get(linked, 0) + 1
        if linked not in (1, 3):
            bad_count += 1
        iso = tournament_canon(out) == C3T
        if (linked == 3) != iso:
            bad_iff += 1
    ok(bad_count == 0, f'H: {gp} random general-position configurations: linked-pair counts {dict(sorted(counts.items()))} (always 1 or 3)')
    ok(bad_lt == 0, 'H: the Gale tournament T* is locally transitive (circular) in every configuration')
    ok(bad_iff == 0 and counts.get(3, 0) > 0, 'H: #linked = 3 iff T* is isomorphic to C3[TT2,TT2,TT2]')
    ok(bad_pierce == 0, 'H: piercing rule (# B-edges piercing A = # uniform-opposite b in B) for all 10 splits of every configuration')


def main():
    global ENGINE
    t0 = time.time()
    with tempfile.TemporaryDirectory() as d:
        ENGINE = os.path.join(d, 'arcs')
        subprocess.run(['cc', '-O2', '-o', ENGINE, os.path.join(HERE, 'procgen_petersen_20261001_orchestrator_check.c')], check=True)
        for name, fn in [('A', check_A), ('B', check_B), ('C', check_C), ('D', check_D), ('E', check_E), ('F', check_F),
                         ('G', check_G), ('H', check_H)]:
            print(f'==== {name} ====', flush=True)
            fn()
    print(f'elapsed {time.time() - t0:.0f} s')
    print('ALL CHECKS PASSED' if all(OKS) else 'SOME CHECK FAILED')


if __name__ == '__main__':
    main()
