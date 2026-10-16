#!/usr/bin/env python3
"""S38 audit (independent referee): exhaustive check that no 19-vertex graph can be a vertex-minimal failed
six-colour exterior in the ABSTRACT model of THM-4661, using none of the theorem's list-colouring lemmas.

Abstract model (weaker than the geometry, so a negative answer here implies the geometric statement):
  * H is a finite simple graph, triangle-free (I2) and any two distinct vertices have at most two common
    neighbours (I3);
  * every vertex b has a list L(b) = [6] or [6] minus one colour (I1);
  * H is not L-colourable and every proper induced subgraph is (vertex-minimality), hence deg(b) >= |L(b)| >= 5.

Steps (each by SAT with my own encodings; no use of Borodin / Erdos-Rubin-Taylor, of the Gallai count (1), or of
the parity step (v) of THM-4661):
  A1a  no vertex of degree >= 7 on 19 vertices (triangle-free, codegree <= 2, min degree 5);
  A1b  rooted at a degree-6 vertex z (one exists: 19 * 5 is odd), the structure is tight: the six neighbours have
       degree 5, every other vertex has exactly two neighbours in N(z);
  A1c  the neighbour pairs of the twelve second neighbours are distinct and form a 4-regular graph on six indices,
       hence K_{2,2,2} (complement of a perfect matching); fixed WLOG;
  A2   enumerate every completion (edges inside the second layer), reduce modulo Aut(K_{2,2,2}) (order 48);
  A3   for every completion: three-colourability, and a CEGAR search for a list assignment (one forbidden colour
       per vertex, the hardest case) that makes it non-L-colourable.
Also: the same pipeline at 16 vertices (5-regular, tight: Clebsch type) and an existence test at 18 vertices.
Output: s38_audit_six_colour_20261016/audit_nineteen_vertex_exhaustive.out (stdout) and .json.
"""
import itertools, json, os, sys, time
from pysat.formula import IDPool, CNF
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_nineteen_vertex_exhaustive.json')
RES = {}


def log(*a):
    print(*a, flush=True)


class GraphSAT:
    """Edge variables e(i,j) on n vertices with triangle-free and codegree <= 2 constraints."""

    def __init__(self, n, solver='cadical195'):
        self.n = n
        self.pool = IDPool()
        self.cls = []
        self.solver_name = solver
        for i, j, k in itertools.combinations(range(n), 3):
            self.cls.append([-self.e(i, j), -self.e(j, k), -self.e(i, k)])
        for i, j in itertools.combinations(range(n), 2):
            cs = []
            for k in range(n):
                if k in (i, j):
                    continue
                c = self.pool.id(('c', i, j, k))
                self.cls.append([-self.e(i, k), -self.e(j, k), c])
                cs.append(c)
            self.card(cs, 2, 'atmost')

    def e(self, i, j):
        if i > j:
            i, j = j, i
        assert i != j
        return self.pool.id(('e', i, j))

    def card(self, lits, k, kind, cond=None):
        if kind == 'atmost':
            enc = CardEnc.atmost(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter)
        elif kind == 'atleast':
            enc = CardEnc.atleast(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter)
        else:
            enc = CardEnc.equals(lits=lits, bound=k, vpool=self.pool, encoding=EncType.seqcounter)
        for cl in enc.clauses:
            self.cls.append(cl + ([-cond] if cond is not None else []))

    def nbrs(self, v):
        return [self.e(v, u) for u in range(self.n) if u != v]

    def solve(self, extra=(), assumptions=()):
        with Solver(name=self.solver_name, bootstrap_with=self.cls + list(extra)) as S:
            ok = S.solve(assumptions=list(assumptions))
            return ok, (S.get_model() if ok else None)


def graph_from_model(gs, model):
    pos = set(l for l in model if l > 0)
    return [(i, j) for i, j in itertools.combinations(range(gs.n), 2) if gs.e(i, j) in pos]


# ------------------------------------------------------------------------------------------------ A1
def step_A1c(n=19):
    """Tight structure: two second neighbours on the same pair of N(z) would give that pair three common neighbours."""
    t0 = time.time()
    gs = GraphSAT(n)
    for v in range(n):
        gs.card(gs.nbrs(v), 5, 'atleast')
        gs.card(gs.nbrs(v), 6, 'atmost')
    extra = [[gs.e(0, u)] for u in range(1, 7)] + [[-gs.e(0, w)] for w in range(7, n)]
    for u in range(1, 7):
        gs.card(gs.nbrs(u), 5, 'equals')
    for w in range(7, n):
        gs.card([gs.e(w, u) for u in range(1, 7)], 2, 'equals')
    sel = []
    for w1, w2 in itertools.combinations(range(7, n), 2):
        for u1, u2 in itertools.combinations(range(1, 7), 2):
            s = gs.pool.id(('same', w1, w2, u1, u2))
            extra += [[-s, gs.e(w1, u1)], [-s, gs.e(w1, u2)], [-s, gs.e(w2, u1)], [-s, gs.e(w2, u2)]]
            sel.append(s)
    extra.append(sel)
    ok, _ = gs.solve(extra=extra)
    log(f'A1c n={n}: tight structure with two second neighbours on the same pair: SAT={ok} ({time.time()-t0:.1f}s)')
    assert not ok
    return dict(A1c_repeated_pair_sat=ok)


def step_A1(n=19):
    out = {}
    # A1a: vertex 0 of degree >= 7 (min degree 5 everywhere)
    t0 = time.time()
    gs = GraphSAT(n)
    for v in range(n):
        gs.card(gs.nbrs(v), 5, 'atleast')
    gs.card(gs.nbrs(0), 7, 'atleast')
    ok, _ = gs.solve()
    out['A1a_degree_ge7_sat'] = ok
    log(f'A1a n={n}: a vertex of degree >= 7 with min degree 5, triangle-free, codegree<=2: SAT={ok} '
        f'({time.time()-t0:.1f}s)')
    # A1b: root z=0 adjacent exactly to 1..6; a violation of tightness
    t0 = time.time()
    gs = GraphSAT(n)
    for v in range(n):
        gs.card(gs.nbrs(v), 5, 'atleast')
        gs.card(gs.nbrs(v), 6, 'atmost')          # justified by A1a (by symmetry, any vertex)
    extra = [[gs.e(0, u)] for u in range(1, 7)] + [[-gs.e(0, w)] for w in range(7, n)]
    sel = []
    for u in range(1, 7):                          # u has degree >= 6
        a = gs.pool.id(('selU', u))
        gs.card(gs.nbrs(u), 6, 'atleast', cond=a)
        sel.append(a)
    for w in range(7, n):                          # w has <= 1 or >= 3 neighbours in N(z)
        b = gs.pool.id(('selLo', w)); c = gs.pool.id(('selHi', w))
        gs.card([gs.e(w, u) for u in range(1, 7)], 1, 'atmost', cond=b)
        gs.card([gs.e(w, u) for u in range(1, 7)], 3, 'atleast', cond=c)
        sel += [b, c]
    extra.append(sel)
    ok, model = gs.solve(extra=extra)
    out['A1b_tightness_violation_sat'] = ok
    log(f'A1b n={n}: rooted at a degree-6 vertex, some violation of tightness: SAT={ok} ({time.time()-t0:.1f}s)')
    # A1c: tight structure; are two second neighbours ever attached to the same pair? is the pair graph 4-regular?
    t0 = time.time()
    gs = GraphSAT(n)
    for v in range(n):
        gs.card(gs.nbrs(v), 5, 'atleast')
        gs.card(gs.nbrs(v), 6, 'atmost')
    extra = [[gs.e(0, u)] for u in range(1, 7)] + [[-gs.e(0, w)] for w in range(7, n)]
    for u in range(1, 7):
        gs.card(gs.nbrs(u), 5, 'equals')
    for w in range(7, n):
        gs.card([gs.e(w, u) for u in range(1, 7)], 2, 'equals')
    sel = []
    for w1, w2 in itertools.combinations(range(7, n), 2):
        for u1, u2 in itertools.combinations(range(1, 7), 2):
            s = gs.pool.id(('same', w1, w2, u1, u2))
            extra += [[-s, gs.e(w1, u1)], [-s, gs.e(w1, u2)], [-s, gs.e(w2, u1)], [-s, gs.e(w2, u2)]]
            sel.append(s)
    extra.append(sel)
    ok, _ = gs.solve(extra=extra)
    out['A1c_repeated_pair_sat'] = ok
    log(f'A1c n={n}: tight structure with two second neighbours on the same pair: SAT={ok} ({time.time()-t0:.1f}s)')
    return out


# ------------------------------------------------------------------------------------------------ A2
def octahedron_pairs():
    match = [(1, 2), (3, 4), (5, 6)]
    pairs = [p for p in itertools.combinations(range(1, 7), 2) if p not in match]
    assert len(pairs) == 12 and all(sum(1 for p in pairs if u in p) == 4 for u in range(1, 7))
    return pairs, match


def aut_k222(match):
    grp = []
    for perm in itertools.permutations(range(1, 7)):
        pi = dict(zip(range(1, 7), perm))
        if all(tuple(sorted((pi[a], pi[b]))) in match for a, b in match):
            grp.append(pi)
    assert len(grp) == 48
    return grp


def step_A2(n=19):
    pairs, match = octahedron_pairs()
    W = list(range(7, 19))
    pair_of = dict(zip(W, pairs))
    gs = GraphSAT(n)
    fixed = []
    for i, j in itertools.combinations(range(n), 2):
        if i == 0:
            fixed.append([gs.e(i, j)] if 1 <= j <= 6 else [-gs.e(i, j)])
        elif 1 <= i <= 6 and 1 <= j <= 6:
            fixed.append([-gs.e(i, j)])
        elif 1 <= i <= 6 and j >= 7:
            fixed.append([gs.e(i, j)] if i in pair_of[j] else [-gs.e(i, j)])
    for w in W:
        gs.card(gs.nbrs(w), 5, 'atleast')
        gs.card(gs.nbrs(w), 6, 'atmost')
    WW = [(a, b) for a, b in itertools.combinations(W, 2)]
    sols = []
    t0 = time.time()
    with Solver(name='cadical195', bootstrap_with=gs.cls + fixed) as S:
        while S.solve():
            m = set(l for l in S.get_model() if l > 0)
            es = frozenset((a, b) for a, b in WW if gs.e(a, b) in m)
            sols.append(es)
            S.add_clause([(-gs.e(a, b) if (a, b) in es else gs.e(a, b)) for a, b in WW])
    log(f'A2: {len(sols)} labelled completions of the tight octahedral structure ({time.time()-t0:.1f}s)')
    # orbit reduction under Aut(K222)
    grp = aut_k222(match)
    idx = {p: w for w, p in pair_of.items()}

    def act(pi, es):
        out = []
        for a, b in es:
            pa = tuple(sorted((pi[pair_of[a][0]], pi[pair_of[a][1]])))
            pb = tuple(sorted((pi[pair_of[b][0]], pi[pair_of[b][1]])))
            x, y = idx[pa], idx[pb]
            out.append((min(x, y), max(x, y)))
        return tuple(sorted(out))

    canon = {}
    for es in sols:
        c = min(act(pi, es) for pi in grp)
        canon.setdefault(c, 0)
        canon[c] += 1
    reps = sorted(canon)
    log(f'A2: {len(reps)} orbits under Aut(K222) (order 48); orbit sizes {sorted(canon.values())}')
    graphs = []
    for c in reps:
        E = [(0, u) for u in range(1, 7)] + [(u, w) for w in W for u in pair_of[w]] + list(c)
        graphs.append(E)
    # sanity: every solution graph satisfies the constraints (independent recheck with networkx)
    for E in graphs:
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
        assert nx.number_of_selfloops(G) == 0
        assert sum(nx.triangles(G).values()) == 0
        assert all(5 <= d <= 6 for _, d in G.degree())
        for a, b in itertools.combinations(range(n), 2):
            assert len(set(G[a]) & set(G[b])) <= 2
    # isomorphism classes of the full graphs
    classes = []
    for E in graphs:
        G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
        if not any(nx.is_isomorphic(G, H) for H in classes):
            classes.append(G)
    log(f'A2: {len(classes)} isomorphism classes of 19-vertex graphs (triangle-free, codegree <= 2, degrees 5/6)')
    return sols, graphs, classes


# ------------------------------------------------------------------------------------------------ A3
def colourable(G, k, forbid=None, solver='glucose4'):
    nodes = list(G.nodes())
    var = {(v, c): i + 1 for i, (v, c) in enumerate(itertools.product(nodes, range(k)))}
    cls = [[var[(v, c)] for c in range(k)] for v in nodes]
    for a, b in G.edges():
        for c in range(k):
            cls.append([-var[(a, c)], -var[(b, c)]])
    if forbid:
        for v, c in forbid.items():
            cls.append([-var[(v, c)]])
    with Solver(name=solver, bootstrap_with=cls) as S:
        if not S.solve():
            return None
        m = set(l for l in S.get_model() if l > 0)
    col = {v: next(c for c in range(k) if var[(v, c)] in m) for v in nodes}
    assert all(col[a] != col[b] for a, b in G.edges())
    if forbid:
        assert all(col[v] != c for v, c in forbid.items())
    return col


def blocking_list_assignment(G, q=6, max_iter=200000):
    """CEGAR: is there f: V -> [q] such that no proper q-colouring avoids f everywhere?  Returns (f or None, iters).
    Each found colouring c is added together with its q! colour permutations (colourings are permutation-closed)."""
    nodes = list(G.nodes())
    pool = IDPool()
    fv = lambda v, c: pool.id(('f', v, c))
    S = Solver(name='cadical195')
    for v in nodes:
        enc = CardEnc.equals(lits=[fv(v, c) for c in range(q)], bound=1, vpool=pool, encoding=EncType.pairwise)
        for cl in enc.clauses:
            S.add_clause(cl)
    perms = list(itertools.permutations(range(q)))
    it = 0
    while it < max_iter:
        it += 1
        if not S.solve():
            S.delete()
            return None, it
        m = set(l for l in S.get_model() if l > 0)
        f = {v: next(c for c in range(q) if fv(v, c) in m) for v in nodes}
        col = colourable(G, q, forbid=f, solver='cadical195')
        if col is None:
            S.delete()
            return f, it
        for p in perms:
            S.add_clause([fv(v, p[col[v]]) for v in nodes])
    S.delete()
    raise RuntimeError('CEGAR did not converge')


def step_A3(classes, tag):
    out = []
    for idx, G in enumerate(classes):
        t0 = time.time()
        degs = sorted(d for _, d in G.degree())
        c3 = colourable(G, 3)
        f, iters = blocking_list_assignment(G)
        rec = dict(index=idx, edges=G.number_of_edges(), degree_sequence_counts={d: degs.count(d) for d in set(degs)},
                   three_colourable=c3 is not None, blocking_list_assignment_exists=f is not None,
                   cegar_iterations=iters, seconds=round(time.time() - t0, 1),
                   edge_list=sorted(tuple(sorted(e)) for e in G.edges()))
        out.append(rec)
        log(f'A3 {tag} #{idx}: edges {rec["edges"]}, degrees {rec["degree_sequence_counts"]}, 3-colourable '
            f'{rec["three_colourable"]}, blocking list assignment exists: {rec["blocking_list_assignment_exists"]} '
            f'(CEGAR {iters} rounds, {rec["seconds"]}s)')
    return out


# ------------------------------------------------------------------------------------------------ 16 and 18
def step_16():
    """5-regular on 16 vertices; a degree-6 vertex would need 1 + 6 + 12 = 19 vertices.  Tight at a root:
    N(z) = 1..5, ten second neighbours on the ten pairs of 1..5 (distinct pairs by codegree <= 2)."""
    n = 16
    pairs = list(itertools.combinations(range(1, 6), 2))
    W = list(range(6, 16))
    pair_of = dict(zip(W, pairs))
    gs = GraphSAT(n)
    fixed = []
    for i, j in itertools.combinations(range(n), 2):
        if i == 0:
            fixed.append([gs.e(i, j)] if 1 <= j <= 5 else [-gs.e(i, j)])
        elif 1 <= i <= 5 and 1 <= j <= 5:
            fixed.append([-gs.e(i, j)])
        elif 1 <= i <= 5 and j >= 6:
            fixed.append([gs.e(i, j)] if i in pair_of[j] else [-gs.e(i, j)])
    for w in W:
        gs.card(gs.nbrs(w), 5, 'equals')
    WW = list(itertools.combinations(W, 2))
    graphs = []
    with Solver(name='cadical195', bootstrap_with=gs.cls + fixed) as S:
        while S.solve():
            m = set(l for l in S.get_model() if l > 0)
            es = [(a, b) for a, b in WW if gs.e(a, b) in m]
            E = [(0, u) for u in range(1, 6)] + [(u, w) for w in W for u in pair_of[w]] + es
            G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(E)
            graphs.append(G)
            S.add_clause([(-gs.e(a, b) if (a, b) in es else gs.e(a, b)) for a, b in WW])
    classes = []
    for G in graphs:
        if not any(nx.is_isomorphic(G, H) for H in classes):
            classes.append(G)
    clebsch = nx.Graph()
    import numpy as np
    verts = list(itertools.product([0, 1], repeat=4))
    for a, b in itertools.combinations(range(16), 2):
        d = sum(x != y for x, y in zip(verts[a], verts[b]))
        if d == 1 or d == 4:
            clebsch.add_edge(a, b)
    log(f'16: {len(graphs)} labelled tight completions, {len(classes)} isomorphism class(es); '
        f'Clebsch: {[nx.is_isomorphic(G, clebsch) for G in classes]}')
    return classes


def step_18_existence():
    n = 18
    t0 = time.time()
    gs = GraphSAT(n)
    for v in range(n):
        gs.card(gs.nbrs(v), 5, 'equals')
    extra = [[gs.e(0, u)] for u in range(1, 6)] + [[-gs.e(0, w)] for w in range(6, n)]
    ok, model = gs.solve(extra=extra)
    log(f'18: a 5-regular triangle-free codegree<=2 graph on 18 vertices exists: {ok} ({time.time()-t0:.1f}s)')
    return ok, (graph_from_model(gs, model) if ok else None)


if __name__ == '__main__':
    T = time.time()
    # A1a/A1b (no degree >= 7; tightness at a degree-6 root) are pure double counting of pigeonhole type, which is
    # exponentially hard for resolution-based SAT (a first run of step_A1 with the full 19-vertex encoding did not
    # finish in several minutes and was stopped).  They are recorded as arithmetic; the substantive steps A2/A3
    # are exhaustive SAT enumeration.  step_A1 is kept for reference (run with --a1 to try it).
    n = 19
    RES['A1'] = dict(
        degree7_lower_bound=1 + 7 + (7 * 4 + 1) // 2,
        degree6_lower_bound=1 + 6 + 6 * 4 // 2,
        five_regular_parity_impossible=(19 * 5) % 2 == 1,
        tight_incidences_rows=6 * 4, tight_incidences_capacity=12 * 2)
    assert RES['A1']['degree7_lower_bound'] == 22 > n and RES['A1']['degree6_lower_bound'] == 19
    log('A1 (arithmetic): degree 7 forces >= 22 vertices; degree 6 forces >= 19 with equality only if tight '
        '(6*4 = 24 incidences into 12 second neighbours of capacity 2); 19*5 is odd, so a degree-6 vertex exists')
    if '--a1' in sys.argv:
        RES['A1_sat'] = step_A1(19)
    RES['A1c'] = step_A1c(19)
    sols, graphs, classes = step_A2(19)
    RES['A2'] = dict(labelled_completions=len(sols), orbits=len(graphs), isomorphism_classes=len(classes))
    RES['A3_19'] = step_A3(classes, 'n=19')
    c16 = step_16()
    RES['n16'] = step_A3(c16, 'n=16')
    # 18 vertices: moved to audit_eighteen_vertex_five_cores.py (the unstructured SAT existence query
    # step_18_existence did not finish within 15 minutes and was stopped).
    RES['seconds'] = round(time.time() - T, 1)
    json.dump(RES, open(OUTJ, 'w'), indent=1, default=str)
    log('wrote', OUTJ)
