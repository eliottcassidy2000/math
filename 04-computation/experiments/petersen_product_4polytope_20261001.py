"""Is the Cartesian product of two Petersen graphs the graph of a 4-polytope?  (Ziegler's question, dimension 4.)

Companion to 05-knowledge/results/petersen_product_4polytope_attempt_20261001.md.

Everything here is a NECESSARY condition for a graph G to be the graph of a 4-polytope Pi, encoded in SAT, with
lazily added exact cuts.  "UNSAT" therefore means: no 4-polytope with graph G whose 2-faces are among the listed
candidate cycles.

Model (explicit candidate 2-faces = induced cycles of G):
  F_Q        2-face Q is chosen;
  A(x;a,b)   the angle a-x-b lies in a chosen 2-face            (A <-> OR of faces through the angle);
  M(x,y;z,w) the 3-path z-x-y-w lies in a chosen 2-face         (M <-> OR of faces through the 3-path);
  C(x,y;z,z') z,z' are consecutive around y in Gamma_x            (Tutte-face DNF of the A's at x).
Static constraints:
  (I) two chosen faces never share a NON-ADJACENT pair of vertices (exact characterisation of improper meeting of two
      induced cycles): at-most-one per non-adjacent pair;
  (V) every vertex figure Gamma_x (the graph of covered angles at x) is polyhedral = 3-connected planar
      (CNF verified against brute force for 4, 5, 6 neighbours: 1, 25, 1227 labelled graphs);
  (R) edge figures are polygons seen identically from both ends: two faces through the edge xy are consecutive around
      y in Gamma_x iff they are consecutive around x in Gamma_y.
Lazy constraints (exact cuts, see the note for the soundness arguments):
  (H) the chosen faces span the Z/2 cycle space of G (H_1(S^3; Z/2) = 0): cut = "some face with odd w-weight";
  (F) facets: germs = faces of Gamma_x; germs glue across an edge when they share two consecutive 2-faces; every glued
      class must be a 3-polytope: one germ per vertex (F1), induced (F4), Euler characteristic 2 (F3), 3-connected
      (F5), and two facets meet in a common face (F6).  Cuts use germ paths (short witnesses).
Optional:
  (CGS) Conway-Gordon-Sachs: in each Petersen fibre some complementary pentagon pair has neither member capped by
      2-faces avoiding the other (otherwise the linking sum would be even).

Usage:  python petersen_product_4polytope_20261001.py [verify|sanity|structure|k33c3|pc3 L|pp-packing|
                                                       pp-prismgerm|pp-nobent-norot|pp-packing-norot|orientation|
                                                       pp-levels k [m]|pp-levels-orient k|all|all-long]
        'all' takes about 10 minutes (it includes the second-solver re-solve for P x C3 with faces <= 10); 'all-long'
        adds the 9-fibre relaxation, the 8-fibre relaxation with orientability, and P x C3 with all 7681 induced cycles
        (about 45 minutes; pass --crosscheck for the Glucose re-solve of the final formula, another 20 minutes).
"""
import itertools
import sys
import time
from collections import Counter, deque

import networkx as nx
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

SOLVER = 'cadical195'


# ============================================================================ graphs
def petersen():
    G = nx.Graph()
    for i in range(5):
        G.add_edge(i, (i + 1) % 5); G.add_edge(i, i + 5); G.add_edge(5 + i, 5 + (i + 2) % 5)
    return G


def product(F1, F2, n=10):
    """Cartesian product with vertex (u, v) encoded as n*u + v; 'H' edges change u, 'V' edges change v."""
    G = nx.Graph()
    for u in F1.nodes():
        for v in F2.nodes():
            G.add_node(n * u + v)
            for w in F1.neighbors(u):
                G.add_edge(n * u + v, n * w + v)
            for w in F2.neighbors(v):
                G.add_edge(n * u + v, n * u + w)
    return G


def induced_cycles(G, L):
    """All induced cycles with at most L vertices (each once), as vertex tuples."""
    nodes = sorted(G.nodes())
    order = {x: i for i, x in enumerate(nodes)}
    adj = {x: set(G.neighbors(x)) for x in nodes}
    out = []

    def dfs(path, onp):
        last, s = path[-1], path[0]
        for w in adj[last]:
            if w in onp or order[w] < order[s]:
                continue
            if any(w in adj[q] for q in path[1:-1]):
                continue
            if len(path) >= 2 and s in adj[w]:
                if len(path) + 1 >= 3 and order[path[1]] < order[w]:
                    out.append(tuple(path + [w]))
                continue
            if len(path) + 1 >= L:
                continue
            onp.add(w); path.append(w)
            dfs(path, onp)
            path.pop(); onp.discard(w)

    for s in nodes:
        dfs([s], {s})
    return out


# ============================================================================ polyhedral vertex figures
def polyhedral_cnf(d):
    """CNF over the C(d,2) pair variables: graph is 3-connected and planar (d <= 6)."""
    assert d <= 6, 'the planarity clauses are complete only for at most 6 vertices'
    pairs = list(itertools.combinations(range(d), 2))
    pid = {p: i for i, p in enumerate(pairs)}
    P = lambda a, b: pid[(min(a, b), max(a, b))]
    cls = []
    for r in range(0, 3):
        for R in itertools.combinations(range(d), r):
            W = [w for w in range(d) if w not in R]
            for k in range(1, len(W)):
                for S in itertools.combinations(W, k):
                    if W[0] not in S:
                        continue
                    T = [w for w in W if w not in S]
                    cls.append([+(P(s, t) + 1) for s in S for t in T])
    for Q in itertools.combinations(range(d), 5):
        cls.append([-(P(a, b) + 1) for a, b in itertools.combinations(Q, 2)])
    if d == 6:
        for c in range(6):
            rest = [w for w in range(6) if w != c]
            for a, b in itertools.combinations(rest, 2):
                E = [(s, t) for s, t in itertools.combinations(rest, 2) if {s, t} != {a, b}] + [(c, a), (c, b)]
                cls.append([-(P(s, t) + 1) for s, t in E])
        for S in itertools.combinations(range(6), 3):
            if 0 not in S:
                continue
            T = [w for w in range(6) if w not in S]
            cls.append([-(P(s, t) + 1) for s in S for t in T])
    return pairs, cls


def polyhedral_labelled(d):
    pairs = list(itertools.combinations(range(d), 2))
    out = []
    for mask in range(1 << len(pairs)):
        E = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
        H = nx.Graph(); H.add_nodes_from(range(d)); H.add_edges_from(E)
        if min(dd for _, dd in H.degree()) < 3 or not nx.is_connected(H) or nx.node_connectivity(H) < 3:
            continue
        ok, emb = nx.check_planarity(H)
        if not ok:
            continue
        consec = set()
        for y in range(d):
            nb = list(emb.neighbors_cw_order(y))
            for i in range(len(nb)):
                z1, z2 = nb[i], nb[(i + 1) % len(nb)]
                consec.add((y, min(z1, z2), max(z1, z2)))
        out.append((mask, set(E), consec))
    return out


def consec_terms(d, y, z1, z2):
    """DNF: z1, z2 consecutive around y in a polyhedral graph on d <= 6 vertices <=> the path z1-y-z2 lies on a face;
    faces = induced non-separating cycles (Tutte).  Terms are lists of ((a,b), polarity)."""
    assert d <= 6, 'the DNF is exact only for at most 6 vertices'
    W = [w for w in range(d) if w not in (y, z1, z2)]
    E = lambda a, b: (min(a, b), max(a, b))
    base = [(E(y, z1), True), (E(y, z2), True)]
    terms = []
    if d == 6:
        for a, b, c in [(W[0], W[1], W[2]), (W[1], W[2], W[0]), (W[2], W[0], W[1])]:
            terms.append(base + [(E(z1, z2), True), (E(a, b), True), (E(b, c), True)])
    elif d == 5:
        terms.append(base + [(E(z1, z2), True), (E(W[0], W[1]), True)])
    elif d == 4:
        terms.append(base + [(E(z1, z2), True)])
    for w in W:
        rest = [r for r in W if r != w]
        t = base + [(E(z1, z2), False), (E(y, w), False), (E(z1, w), True), (E(z2, w), True)]
        if len(rest) == 2:
            t = t + [(E(rest[0], rest[1]), True)]
        terms.append(t)
    for w1, w2 in itertools.permutations(W, 2):
        if len(W) - 2 > 1:
            continue
        terms.append(base + [(E(z1, w1), True), (E(w1, w2), True), (E(w2, z2), True), (E(z1, z2), False),
                             (E(y, w1), False), (E(y, w2), False), (E(z1, w2), False), (E(z2, w1), False)])
    return terms


# ============================================================================ explicit-face model
class Model:
    def __init__(self, G, faces, rotation=True, X=None, verbose=True):
        self.G = G
        self.nodes = sorted(G.nodes())
        self.N = {x: sorted(G.neighbors(x)) for x in self.nodes}
        self.adj = {x: set(self.N[x]) for x in self.nodes}
        self.X = set(self.nodes) if X is None else set(X)
        self.faces = [tuple(f) for f in faces]
        self.fidx = {f: i for i, f in enumerate(self.faces)}
        self.pool = IDPool()
        self.cls = []
        self.pending = []
        self.extra_static = []
        t0 = time.time()
        self.build(rotation)
        if verbose:
            print(f'  model: {len(self.faces)} candidate faces, {self.pool.top} vars, {len(self.cls)} clauses '
                  f'({time.time()-t0:.1f}s)', flush=True)

    def F(self, i): return self.pool.id(('F', i))

    def A(self, x, a, b):
        if a > b: a, b = b, a
        return self.pool.id(('A', x, a, b))

    def M(self, x, y, z, w):
        return self.pool.id(('M', x, y, z, w)) if x < y else self.pool.id(('M', y, x, w, z))

    def C(self, x, y, z1, z2):
        if z1 > z2: z1, z2 = z2, z1
        return self.pool.id(('C', x, y, z1, z2))

    def build(self, rotation):
        cls, N, adj = self.cls, self.N, self.adj
        angle_f, path_f, pair_f = {}, {}, {}
        for i, f in enumerate(self.faces):
            n = len(f)
            for t in range(n):
                a, x, b = f[t - 1], f[t], f[(t + 1) % n]
                angle_f.setdefault((x, min(a, b), max(a, b)), []).append(i)
                z, x2, y, w = f[t - 1], f[t], f[(t + 1) % n], f[(t + 2) % n]
                path_f.setdefault((x2, y, z, w) if x2 < y else (y, x2, w, z), []).append(i)
            for s, t in itertools.combinations(range(n), 2):
                if f[t] not in adj[f[s]]:
                    pair_f.setdefault((min(f[s], f[t]), max(f[s], f[t])), []).append(i)
        self.angle_index = angle_f
        for fl in pair_f.values():                                   # (I)
            if len(fl) > 1:
                cls.extend(CardEnc.atmost([self.F(i) for i in fl], bound=1, vpool=self.pool,
                                          encoding=EncType.ladder).clauses)
        for x in self.nodes:                                         # A
            for a, b in itertools.combinations(N[x], 2):
                fl = angle_f.get((x, a, b), [])
                v = self.A(x, a, b)
                cls.append([-v] + [self.F(i) for i in fl])
                for i in fl: cls.append([-self.F(i), v])
        for x in self.nodes:                                         # M
            for y in N[x]:
                if not x < y: continue
                for z in N[x]:
                    if z == y: continue
                    for w in N[y]:
                        if w == x: continue
                        fl = path_f.get((x, y, z, w), [])
                        v = self.M(x, y, z, w)
                        cls.append([-v] + [self.F(i) for i in fl])
                        for i in fl: cls.append([-self.F(i), v])
        cache = {}
        for x in self.nodes:                                         # (V)
            if x not in self.X: continue
            d = len(N[x])
            if d not in cache: cache[d] = polyhedral_cnf(d)
            pairs, pc = cache[d]
            var = [self.A(x, N[x][a], N[x][b]) for a, b in pairs]
            for c in pc:
                cls.append([var[abs(l) - 1] if l > 0 else -var[abs(l) - 1] for l in c])
        if not rotation:
            return
        for x in self.nodes:                                         # C via Tutte DNF
            if x not in self.X: continue
            d = len(N[x]); loc = N[x]
            for yi in range(d):
                for zi, zj in itertools.combinations([q for q in range(d) if q != yi], 2):
                    c = self.C(x, loc[yi], loc[zi], loc[zj])
                    auxs = []
                    for ti, term in enumerate(consec_terms(d, yi, zi, zj)):
                        t = self.pool.id(('T', x, yi, zi, zj, ti)); auxs.append(t)
                        lits = [self.A(x, loc[a], loc[b]) if pol else -self.A(x, loc[a], loc[b]) for (a, b), pol in term]
                        for l in lits: cls.append([-t, l])
                        cls.append([t] + [-l for l in lits])
                        cls.append([-t, c])
                    cls.append([-c] + auxs)
        for x in self.nodes:                                         # (R)
            for y in N[x]:
                if not x < y or x not in self.X or y not in self.X: continue
                Zs = [z for z in N[x] if z != y]; Ws = [w for w in N[y] if w != x]
                for z1, z2 in itertools.combinations(Zs, 2):
                    for w1 in Ws:
                        for w2 in Ws:
                            if w1 == w2: continue
                            m1, m2 = self.M(x, y, z1, w1), self.M(x, y, z2, w2)
                            cls.append([-m1, -m2, -self.C(x, y, z1, z2), self.C(y, x, w1, w2)])
                            cls.append([-m1, -m2, self.C(x, y, z1, z2), -self.C(y, x, w1, w2)])

    # ---------------------------------------------------------------- (H) homology cut
    def homology_cut(self, chosen):
        edges = sorted(set((min(a, b), max(a, b)) for a, b in self.G.edges()))
        eid = {e: k for k, e in enumerate(edges)}
        m = len(edges)

        def vec(f):
            v = 0
            for t in range(len(f)):
                a, b = f[t], f[(t + 1) % len(f)]
                v ^= 1 << eid[(min(a, b), max(a, b))]
            return v
        piv = {}
        for i in chosen:
            v = vec(self.faces[i])
            while v:
                h = v.bit_length() - 1
                if h in piv: v ^= piv[h]
                else: piv[h] = v; break
        if len(piv) == m - len(self.nodes) + nx.number_connected_components(self.G):
            return None
        keys = sorted(piv)
        for h in keys:
            for h2 in keys:
                if h2 != h and (piv[h2] >> h) & 1:
                    piv[h2] ^= piv[h]
        cob = {}
        for x in self.nodes:
            v = 0
            for y in self.N[x]: v ^= 1 << eid[(min(x, y), max(x, y))]
            while v:
                h = v.bit_length() - 1
                if h in cob: v ^= cob[h]
                else: cob[h] = v; break

        def in_cob(w):
            while w:
                h = w.bit_length() - 1
                if h not in cob: return False
                w ^= cob[h]
            return True
        for fk in range(m):
            if fk in piv: continue
            w = 1 << fk
            for h, r in piv.items():
                if (r >> fk) & 1: w |= 1 << h
            if not in_cob(w):
                return [self.F(i) for i, f in enumerate(self.faces) if bin(vec(f) & w).count('1') % 2 == 1]
        raise RuntimeError('no homology cut found')


# ============================================================================ (F) facets
def faces_of_planar(H):
    ok, emb = nx.check_planarity(H)
    assert ok
    seen, out = set(), []
    for u, v in emb.edges():
        if (u, v) in seen: continue
        out.append(emb.traverse_face(u, v, mark_half_edges=seen))
    return out


def extract(G, chosen):
    angle = {}
    for qi, f in enumerate(chosen):
        n = len(f)
        for t in range(n):
            a, x, b = f[t - 1], f[t], f[(t + 1) % n]
            angle[(x, min(a, b), max(a, b))] = qi
    byx = {}
    for (x, a, b), q in angle.items():
        byx.setdefault(x, {})[(a, b)] = q
    germs = []
    for x in sorted(G.nodes()):
        fa = byx.get(x, {})
        H = nx.Graph(); H.add_edges_from(fa.keys())
        for cyc in faces_of_planar(H):
            k = len(cyc)
            qs = [fa[(min(cyc[i], cyc[(i + 1) % k]), max(cyc[i], cyc[(i + 1) % k]))] for i in range(k)]
            germs.append((x, tuple(cyc), qs))
    par = list(range(len(germs)))

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    pairidx = {}
    for gi, (x, cyc, qs) in enumerate(germs):
        for i in range(len(cyc)):
            q1, q2 = qs[i - 1], qs[i]
            pairidx[(x, cyc[i], min(q1, q2), max(q1, q2))] = gi
    problems, glue = [], []
    for (x, y, q1, q2), gi in pairidx.items():
        gj = pairidx.get((y, x, q1, q2))
        if gj is None:
            problems.append(('rotation-mismatch', x, y)); continue
        par[find(gi)] = find(gj)
        if x < y: glue.append((gi, gj, q1, q2))
    classes = {}
    for gi in range(len(germs)):
        classes.setdefault(find(gi), []).append(gi)
    facets = []
    for gl in classes.values():
        fq = set()
        for g in gl: fq.update(germs[g][2])
        facets.append({'germs': gl, 'verts': [germs[g][0] for g in gl], 'faces': sorted(fq)})
    return germs, facets, problems, glue


def check_facets(G, chosen, germs, facets):
    adj = {x: set(G.neighbors(x)) for x in G.nodes()}
    bad, vsets = [], []
    for fi, Fc in enumerate(facets):
        vs = Fc['verts']; reasons = []
        if len(set(vs)) != len(vs): reasons.append('F1')
        Vset = set(vs); E = set()
        for q in Fc['faces']:
            f = chosen[q]
            for t in range(len(f)):
                a, b = f[t], f[(t + 1) % len(f)]
                E.add((min(a, b), max(a, b)))
        if len(Vset) - len(E) + len(Fc['faces']) != 2: reasons.append('F3')
        if set((min(a, b), max(a, b)) for a in Vset for b in adj[a] if b in Vset) != E: reasons.append('F4')
        H = nx.Graph(); H.add_edges_from(E)
        if H.number_of_nodes() >= 4 and nx.node_connectivity(H) < 3: reasons.append('F5')
        if reasons: bad.append((fi, reasons))
        vsets.append((Vset, E, set(Fc['faces'])))
    inter, byv, done = [], {}, set()
    for fi, (Vs, E, Fs) in enumerate(vsets):
        for x in Vs: byv.setdefault(x, []).append(fi)
    for x, fl in byv.items():
        for i, j in itertools.combinations(fl, 2):
            if (i, j) in done: continue
            done.add((i, j))
            com = vsets[i][0] & vsets[j][0]
            if len(com) == 1: continue
            if len(com) == 2:
                a, b = tuple(com); e = (min(a, b), max(a, b))
                if e in vsets[i][1] and e in vsets[j][1]: continue
            if any(set(chosen[q]) == com for q in vsets[i][2] & vsets[j][2]): continue
            inter.append((i, j))
    return bad, inter


def canon_cycle(c):
    k = len(c)
    return min(tuple(c[(s + d * i) % k] for i in range(k)) for s in range(k) for d in (1, -1))


def germ_var(model, x, cyc):
    """Variable <-> 'cyc is a face of Gamma_x' (induced and non-separating, Tutte); defined on demand."""
    key = ('G', x, canon_cycle(list(cyc)))
    if key in model.pool.obj2id:
        return model.pool.id(key)
    g = model.pool.id(key)
    k = len(cyc)
    lits = [model.A(x, cyc[i], cyc[(i + 1) % k]) for i in range(k)]
    for i, j in itertools.combinations(range(k), 2):
        if (j - i) % k not in (1, k - 1):
            lits.append(-model.A(x, cyc[i], cyc[j]))
    rest = [r for r in model.N[x] if r not in cyc]
    assert len(rest) <= 3, 'germ_var is exact only for vertex degree <= 6'
    cl = model.pending
    for l in lits: cl.append([-g, l])
    if len(rest) == 3:
        rp = [model.A(x, a, b) for a, b in itertools.combinations(rest, 2)]
        for p, q in itertools.combinations(rp, 2):
            cl.append([-g, p, q]); cl.append([g] + [-l for l in lits] + [-p, -q])
    elif len(rest) == 2:
        r = model.A(x, rest[0], rest[1])
        cl.append([-g, r]); cl.append([g] + [-l for l in lits] + [-r])
    else:
        cl.append([g] + [-l for l in lits])
    return g


def facet_literals(model, chosen, germs, Fc):
    return [model.F(model.fidx[chosen[q]]) for q in Fc['faces']] + \
           [germ_var(model, germs[g][0], germs[g][1]) for g in Fc['germs']]


def germ_path_literals(model, chosen, germs, adjg, src, dst):
    prev = {src: None}; dq = deque([src])
    while dq:
        g = dq.popleft()
        if g == dst: break
        for h, q1, q2 in adjg.get(g, []):
            if h not in prev:
                prev[h] = (g, q1, q2); dq.append(h)
    lits, g = [], dst
    while True:
        lits.append(germ_var(model, germs[g][0], germs[g][1]))
        if prev[g] is None: break
        h, q1, q2 = prev[g]
        lits += [model.F(model.fidx[chosen[q1]]), model.F(model.fidx[chosen[q2]])]
        g = h
    return lits


def solve_full(model, solver, maxit=100000, report=50, facets_on=True, homology_on=True, crosscheck=None):
    """Lazy loop.  Returns (chosen faces, facets) if everything checks, None if UNSAT.
    crosscheck: name of a second SAT solver that re-solves the final formula (base + every added clause) from scratch."""
    t0 = time.time()
    st = Counter()
    added = []
    _add = solver.add_clause

    def add_clause(c):
        added.append(list(c)); _add(c)
    solver.add_clause = add_clause
    for it in range(1, maxit + 1):
        if not solver.solve():
            print(f'  UNSAT after {it} iterations ({time.time()-t0:.1f}s), cuts {dict(st)}', flush=True)
            if crosscheck:
                t1 = time.time()
                s2 = Solver(name=crosscheck, bootstrap_with=model.cls + model.extra_static + added)
                ok2 = s2.solve()
                print(f'  cross-check of the final formula ({len(model.cls) + len(model.extra_static) + len(added)} '
                      f'clauses) with {crosscheck}: {"SAT (DISAGREEMENT)" if ok2 else "UNSAT (agrees)"} '
                      f'({time.time()-t1:.1f}s)', flush=True)
                assert not ok2
            return None
        vs = set(l for l in solver.get_model() if l > 0)
        ch = [i for i in range(len(model.faces)) if model.F(i) in vs]
        if homology_on:
            cut = model.homology_cut(ch)
            if cut is not None:
                solver.add_clause(cut); st['H'] += 1
                continue
        if not facets_on:
            print(f'  SAT (2-face level) after {it} iterations ({time.time()-t0:.1f}s)', flush=True)
            return [model.faces[i] for i in ch], None
        chosen = [model.faces[i] for i in ch]
        germs, facets, problems, glue = extract(model.G, chosen)
        assert not problems, problems[:3]
        adjg = {}
        for gi, gj, q1, q2 in glue:
            adjg.setdefault(gi, []).append((gj, q1, q2)); adjg.setdefault(gj, []).append((gi, q1, q2))
        bad, inter = check_facets(model.G, chosen, germs, facets)
        if not bad and not inter:
            print(f'  SAT (all checks) after {it} iterations ({time.time()-t0:.1f}s), cuts {dict(st)}', flush=True)
            return chosen, facets
        new = []
        for fi, reasons in bad:
            Fc = facets[fi]; byx = {}
            for g in Fc['germs']: byx.setdefault(germs[g][0], []).append(g)
            rep = [x for x, l in byx.items() if len(l) > 1]
            if rep:
                g1, g2 = byx[rep[0]][:2]
                new.append([-l for l in germ_path_literals(model, chosen, germs, adjg, g1, g2)]); st['F1'] += 1
                continue
            if 'F4' in reasons:
                found = None
                for g in Fc['germs']:
                    a, cyc, _ = germs[g]
                    for b in model.adj[a]:
                        if b in byx and b not in cyc:
                            found = (g, byx[b][0]); break
                    if found: break
                new.append([-l for l in germ_path_literals(model, chosen, germs, adjg, *found)]); st['F4'] += 1
                continue
            new.append([-l for l in facet_literals(model, chosen, germs, Fc)]); st['Fglobal'] += 1
        for i, j in inter:
            Fi, Fj = facets[i], facets[j]
            bi = {germs[g][0]: g for g in Fi['germs']}; bj = {germs[g][0]: g for g in Fj['germs']}
            com = [x for x in bi if x in bj]
            shared = [set(chosen[q]) for q in set(Fi['faces']) & set(Fj['faces'])]
            pair = next(((a, b) for a, b in itertools.combinations(com, 2)
                         if b not in model.adj[a] and not any(a in sq and b in sq for sq in shared)), None)
            if pair is None:
                new.append([-l for l in facet_literals(model, chosen, germs, Fi)] +
                           [-l for l in facet_literals(model, chosen, germs, Fj)]); st['F6global'] += 1
                continue
            a, b = pair
            lits = germ_path_literals(model, chosen, germs, adjg, bi[a], bi[b]) + \
                germ_path_literals(model, chosen, germs, adjg, bj[a], bj[b])
            ca, cb = germs[bi[a]][1], germs[bj[a]][1]
            ea = set(frozenset((ca[t], ca[(t + 1) % len(ca)])) for t in range(len(ca)))
            eb = set(frozenset((cb[t], cb[(t + 1) % len(cb)])) for t in range(len(cb)))
            esc = []
            for sh in ea & eb:
                z1, z2 = sorted(sh)
                esc += [model.F(qi) for qi in model.angle_index.get((a, z1, z2), []) if b in model.faces[qi]]
            new.append([-l for l in lits] + esc); st['F6'] += 1
        for c in model.pending: solver.add_clause(c)
        model.pending = []
        for c in new: solver.add_clause(c)
        if it % report == 0:
            print(f'    it {it}: cuts {dict(st)}  ({time.time()-t0:.0f}s)', flush=True)
    return 'maxit'


# ============================================================================ (O) orientability (optional layer)
def add_orientation(model, X=None):
    """Static clauses: the vertex links are coherently oriented, i.e. the 3-manifold is orientable (S^3 is).
    Su(x,y,z1,z2): in the oriented rotation system of Gamma_x, z2 is the successor of z1 around y.
      (1) successors are present and consecutive (C); (2) around each y the successor relation is one cycle on the
      present neighbours (exactly one successor and one predecessor, no 2-cycles); (3) every face phi of Gamma_x is
      traced coherently: Su(n_i; n_(i-1) -> n_(i+1)) for all i, or the reverse for all i (face tracing a->b->sigma_b(a));
      (4) across an edge xy, the cyclic order of the 2-faces around xy seen from x is the reverse of the order seen
      from y: for faces Q1 = (z1 at x, w1 at y), Q2 = (z2, w2):  Su(x,y,z1,z2) <-> Su(y,x,w2,w1)."""
    N, pool, cls = model.N, model.pool, []
    X = set(model.nodes) if X is None else set(X)
    Su = lambda x, y, z1, z2: pool.id(('Su', x, y, z1, z2))
    for x in model.nodes:
        if x not in X: continue
        loc = N[x]
        assert len(loc) <= 6
        for y in loc:
            others = [z for z in loc if z != y]
            for z1, z2 in itertools.permutations(others, 2):
                v = Su(x, y, z1, z2)
                cls += [[-v, model.A(x, y, z1)], [-v, model.A(x, y, z2)], [-v, model.C(x, y, z1, z2)]]
            for z in others:
                succ = [Su(x, y, z, w) for w in others if w != z]
                pred = [Su(x, y, w, z) for w in others if w != z]
                cls.append([-model.A(x, y, z)] + succ); cls.append([-model.A(x, y, z)] + pred)
                for i, j in itertools.combinations(range(len(succ)), 2):
                    cls += [[-succ[i], -succ[j]], [-pred[i], -pred[j]]]
            for z1, z2 in itertools.combinations(others, 2):
                cls.append([-Su(x, y, z1, z2), -Su(x, y, z2, z1)])
        seen = set()
        for k in (3, 4, 5):
            for S in itertools.combinations(loc, k):
                for perm in itertools.permutations(S[1:]):
                    c = canon_cycle((S[0],) + perm)
                    if c in seen: continue
                    seen.add(c)
                    model.pending = []
                    g = germ_var(model, x, list(c))
                    cls += model.pending; model.pending = []
                    dr = pool.id(('Dr', x, c))
                    for i in range(k):
                        a, b, cc = c[i - 1], c[i], c[(i + 1) % k]
                        cls += [[-g, -dr, Su(x, b, a, cc)], [-g, dr, Su(x, b, cc, a)]]
    for x in model.nodes:
        for y in N[x]:
            if not x < y or x not in X or y not in X: continue
            Zs = [z for z in N[x] if z != y]; Ws = [w for w in N[y] if w != x]
            for z1, z2 in itertools.permutations(Zs, 2):
                for w1 in Ws:
                    for w2 in Ws:
                        if w1 == w2: continue
                        m1, m2 = model.M(x, y, z1, w1), model.M(x, y, z2, w2)
                        cls += [[-m1, -m2, -Su(x, y, z1, z2), Su(y, x, w2, w1)],
                                [-m1, -m2, Su(x, y, z1, z2), -Su(y, x, w2, w1)]]
    return cls


def orientation_conflicts(G, faces, X=None):
    """Independent check on a found 2-face system: number of BFS conflicts when trying to orient the vertex figures
    coherently (0 = orientable).  Uses networkx embeddings, not the SAT encoding."""
    X = set(G.nodes()) if X is None else set(X)
    angle = {}
    for qi, f in enumerate(faces):
        n = len(f)
        for t in range(n):
            a, x, b = f[t - 1], f[t], f[(t + 1) % n]
            angle[(x, min(a, b), max(a, b))] = qi
    rot = {}
    for x in X:
        H = nx.Graph(); H.add_edges_from((a, b) for (y, a, b) in angle if y == x)
        ok, emb = nx.check_planarity(H)
        rot[x] = {y: list(emb.neighbors_cw_order(y)) for y in H.nodes()}
    sgn = {}
    for x in X:
        for y in G.neighbors(x):
            if y not in X or not x < y: continue
            fx = [angle[(x, min(y, z), max(y, z))] for z in rot[x][y]]
            fy = [angle[(y, min(x, w), max(x, w))] for w in rot[y][x]]
            k = len(fx); i = fy.index(fx[0])
            sgn[(x, y)] = 1 if all(fy[(i - j) % k] == fx[j] for j in range(k)) else -1
    adj = {}
    for (x, y), v in sgn.items():
        adj.setdefault(x, []).append((y, v)); adj.setdefault(y, []).append((x, v))
    o, conflicts = {}, 0
    for x0 in X:
        if x0 in o: continue
        o[x0] = 1; dq = deque([x0])
        while dq:
            x = dq.popleft()
            for y, v in adj.get(x, []):
                if y not in o: o[y] = o[x] * v; dq.append(y)
                elif o[y] != o[x] * v: conflicts += 1
    return conflicts // 2


# ============================================================================ helpers for P x P
def typ(a, b, n=10):
    return 'H' if a % n == b % n else 'V'


def kind(f, n=10):
    k = len(f)
    t = [typ(f[i], f[(i + 1) % k], n) for i in range(k)]
    turns = sum(1 for i in range(k) if t[i] != t[i - 1])
    return 'fibre' if turns == 0 else ('square' if k == 4 and turns == 4 else 'bent')


def pp_graph():
    return product(petersen(), petersen())


# ============================================================================ experiments
def exp_verify():
    print('== verify: polyhedral CNF and rotation DNF against brute force')
    for d in (4, 5, 6):
        pairs, cls = polyhedral_cnf(d)
        lab = polyhedral_labelled(d)
        masks = set(m for m, _, _ in lab)
        m = len(pairs)
        agree = all((all(any((l > 0) == bool(mask >> (abs(l) - 1) & 1) for l in c) for c in cls)) == (mask in masks)
                    for mask in range(1 << m))
        rot_ok = True
        for mask, E, consec in lab:
            for y in range(d):
                for z1, z2 in itertools.combinations([q for q in range(d) if q != y], 2):
                    val = any(all(((p in E) == pol) for p, pol in t) for t in consec_terms(d, y, z1, z2))
                    rot_ok &= (val == ((y, z1, z2) in consec))
        print(f'  d={d}: {len(lab)} labelled polyhedral graphs; CNF exact: {agree}; rotation DNF exact: {rot_ok}')
        assert agree and rot_ok
    print('  CHECK verify PASSED')


def exp_sanity():
    print('== sanity: known 4-polytopes must pass, known non-polytopal graphs must fail')
    rel = lambda G: nx.convert_node_labels_to_integers(G, ordering='sorted')
    cases = [
        ('C5 x C5 (product of pentagons)', rel(nx.cartesian_product(nx.cycle_graph(5), nx.cycle_graph(5))), True),
        ('J(5,2) (rectified 5-cell)', rel(nx.complement(petersen())), True),
        ('icosahedron x K2 (icosahedral prism)', rel(nx.cartesian_product(nx.icosahedral_graph(), nx.path_graph(2))), True),
        ('K33 x K2 (PPS Prop. 2.11: not polytopal)', rel(nx.cartesian_product(nx.complete_bipartite_graph(3, 3), nx.path_graph(2))), False),
        ('P x K2 (PPS Thm 2.3: not polytopal)', product(petersen(), nx.path_graph(2)), False),
    ]
    for name, G, expect in cases:
        print(' ', name, G.number_of_nodes(), 'vertices')
        m = Model(G, induced_cycles(G, 8))
        r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls))
        got = r is not None
        if got:
            ch, fac = r
            print(f'    2-faces {sorted(Counter(len(f) for f in ch).items())}, facets {sorted(Counter(len(F["verts"]) for F in fac).items())}')
        assert got == expect, name
    print('  CHECK sanity PASSED')


def exp_structure():
    print('== structure of P x P')
    P = petersen(); G = pp_graph()
    cyc = induced_cycles(G, 10)
    print('  induced cycles of P x P up to 10 vertices:', sorted(Counter(len(c) for c in cyc).items()))
    # vertex-figure census (m mixed, pH, pV) over 1227 labelled graphs; local labels 0,1,2 = H ; 3,4,5 = V
    lab = polyhedral_labelled(6)
    cen = Counter()
    for _, E, _ in lab:
        mixed = sum(1 for a, b in E if (a < 3) != (b < 3)); pH = sum(1 for a, b in E if b < 3)
        cen[(mixed, pH, len(E) - mixed - pH)] += 1
    print('  vertex figures (mixed, pureH, pureV):', sorted(cen.items()))
    assert all(k[1] >= 1 and k[2] >= 1 and k[0] <= 8 for k in cen)
    print('  -> every vertex figure has a pure H edge, a pure V edge and at most 8 of the 9 mixed edges')
    # fibre families
    pc = induced_cycles(P, 10)
    pents = [c for c in pc if len(c) == 5]; hexes = [c for c in pc if len(c) == 6]
    print('  induced cycles of P:', Counter(len(c) for c in pc))

    def proper(c1, c2):
        s = set(c1) & set(c2)
        return len(s) <= 1 or (len(s) == 2 and P.has_edge(*s))
    allc = pents + hexes
    GA = nx.Graph(); GA.add_nodes_from(range(len(allc)))
    for i, j in itertools.combinations(range(len(allc)), 2):
        if proper(allc[i], allc[j]): GA.add_edge(i, j)
    fams = Counter(tuple(sorted(Counter(len(allc[k]) for k in f).items())) for f in nx.find_cliques(GA))
    print('  maximal families of pairwise properly meeting 5/6-cycles of P:', dict(fams))
    # Phi on symmetric forms
    E = sorted(tuple(sorted(e)) for e in P.edges()); eid = {e: i for i, e in enumerate(E)}

    def vec(c):
        v = 0
        for i in range(len(c)):
            v ^= 1 << eid[tuple(sorted((c[i], c[(i + 1) % len(c)])))]
        return v
    pv = [vec(c) for c in pents]
    basis, piv = [], {}
    for v0 in pv:
        v = v0
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; basis.append(v0); break
    assert len(basis) == 6

    def coords(v):
        for mask in range(64):
            s = 0
            for k in range(6):
                if mask >> k & 1: s ^= basis[k]
            if s == v: return [mask >> k & 1 for k in range(6)]
    co = [coords(v) for v in pv]
    pairs = [(i, j) for i, j in itertools.combinations(range(12), 2) if not set(pents[i]) & set(pents[j])]
    assert len(pairs) == 6
    vals = []
    for k in range(6):
        for l in range(k, 6):
            if k == l: tot = sum(co[i][k] * co[j][k] for i, j in pairs)
            else: tot = sum(co[i][k] * co[j][l] + co[i][l] * co[j][k] for i, j in pairs)
            vals.append(tot % 2)
    print('  Phi(beta) = sum over the 6 complementary pentagon pairs of beta(C,C\') vanishes on all symmetric '
          'bilinear forms on Z_1(P;Z/2):', all(v == 0 for v in vals))
    assert all(v == 0 for v in vals)
    print('  CHECK structure PASSED')


def nobent_candidates(G, L=6):
    return [f for f in induced_cycles(G, L) if kind(f) != 'bent']


def exp_pp_packing():
    print('== P x P, sub-case: every 2-face is a square or lies in a fibre, and every vertex lies in exactly 8 square '
          '2-faces')
    G = pp_graph()
    m = Model(G, nobent_candidates(G))
    extra = []
    for x in G.nodes():
        loc = m.N[x]
        mixed = [m.A(x, a, b) for a, b in itertools.combinations(loc, 2) if typ(x, a) != typ(x, b)]
        extra += CardEnc.equals(mixed, bound=8, vpool=m.pool, encoding=EncType.seqcounter).clauses
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls + extra))
    print('  result:', 'UNSAT (excluded)' if r is None else 'SAT')
    assert r is None


def exp_pp_prismgerm():
    print('== P x P, sub-case: no bent 2-faces and every vertex figure a triangulation whose faces all have both H and '
          'V vertices')
    G = pp_graph()
    m = Model(G, nobent_candidates(G))
    pairs = list(itertools.combinations(range(6), 2))
    extra, cache = [], {}
    for x in G.nodes():
        loc = m.N[x]
        T = tuple(typ(x, a) for a in loc)
        if T not in cache:
            allowed = []
            for mask, E, _ in polyhedral_labelled(6):
                if len(E) != 12: continue
                H = nx.Graph(); H.add_nodes_from(range(6)); H.add_edges_from(E)
                if all(len(f) == 3 and len(set(T[a] for a in f)) == 2 for f in faces_of_planar(H)):
                    allowed.append(E)
            cache[T] = allowed
        sels = []
        for t, E in enumerate(cache[T]):
            s = m.pool.id(('PS', x, t)); sels.append(s)
            for a, b in pairs:
                v = m.A(x, loc[a], loc[b]); extra.append([-s, v if (a, b) in E else -v])
        extra.append(sels)
    print('  allowed labelled vertex figures per vertex:', sorted(set(len(v) for v in cache.values())))
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls + extra))
    print('  result:', 'UNSAT (excluded)' if r is None else 'SAT')
    assert r is None


def exp_pp_nobent_norot():
    print('== P x P, no bent faces, vertex figures + proper intersections only (no edge-figure orientation): '
          'the purely local 2-face conditions are satisfiable')
    G = pp_graph()
    m = Model(G, nobent_candidates(G), rotation=False)
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls), facets_on=False, homology_on=False)
    ch, _ = r
    print('  a 2-face system with polyhedral vertex figures:', sorted(Counter((len(f), kind(f)) for f in ch).items()))


def exp_pp_levels(k, extra=0):
    print(f'== P x P, no bent faces, constraints (V)+(R) imposed only at the vertices of {k} of the 10 H-fibres'
          + (f' and at {extra} vertices of fibre {k}' if extra else ''))
    G = pp_graph()
    X = set(10 * u + v for u in range(10) for v in range(k)) | set(10 * u + k for u in range(extra))
    m = Model(G, nobent_candidates(G), X=X)
    t0 = time.time()
    ok = Solver(name=SOLVER, bootstrap_with=m.cls).solve()
    print(f'  {"SAT" if ok else "UNSAT"} ({time.time()-t0:.1f}s)')


def exp_pc3(L, crosscheck=False):
    print(f'== P x C3 (Petersen x triangle), candidate 2-faces = all induced cycles with <= {L} vertices')
    G = product(petersen(), nx.cycle_graph(3))
    cyc = induced_cycles(G, L)
    print('  candidates:', len(cyc), sorted(Counter(len(c) for c in cyc).items()))
    m = Model(G, cyc)
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls), report=100,
                   crosscheck=('glucose4' if (crosscheck or '--crosscheck' in sys.argv) else None))
    print('  result:', 'UNSAT (excluded)' if r is None else ('INCONCLUSIVE (iteration cap)' if r == 'maxit' else 'SAT'))
    return r


def exp_pp_packing_norot():
    print('== P x P, packing sub-case WITHOUT the edge-figure orientation (R): satisfiable (so (R) is what excludes it)')
    G = pp_graph()
    m = Model(G, nobent_candidates(G), rotation=False)
    extra = []
    for x in G.nodes():
        mixed = [m.A(x, a, b) for a, b in itertools.combinations(m.N[x], 2) if typ(x, a) != typ(x, b)]
        extra += CardEnc.equals(mixed, bound=8, vpool=m.pool, encoding=EncType.seqcounter).clauses
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls + extra), facets_on=False, homology_on=False)
    ch, _ = r
    print('  a 2-face system:', sorted(Counter((len(f), kind(f)) for f in ch).items()))
    assert all(sum(1 for f in ch if kind(f) == 'square' and x in f) == 8 for x in G.nodes())
    print('  every vertex lies in exactly 8 square 2-faces: True')


def exp_orientation():
    print('== orientability layer (O): validation')
    rel = lambda G: nx.convert_node_labels_to_integers(G, ordering='sorted')
    for name, G in [('C5 x C5', rel(nx.cartesian_product(nx.cycle_graph(5), nx.cycle_graph(5)))),
                    ('J(5,2)', rel(nx.complement(petersen()))),
                    ('icosahedron x K2', rel(nx.cartesian_product(nx.icosahedral_graph(), nx.path_graph(2))))]:
        m = Model(G, induced_cycles(G, 8), verbose=False)
        r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls + add_orientation(m)))
        ch, fac = r
        c = orientation_conflicts(G, ch)
        print(f'  {name}: SAT with (O); independent orientation check: {c} conflicts')
        assert c == 0
    K = nx.convert_node_labels_to_integers(nx.complete_bipartite_graph(3, 3))
    G = product(K, nx.cycle_graph(3))
    m = Model(G, induced_cycles(G, 19), verbose=False)
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls), homology_on=False)
    print('  K33 x K3, (H) off: the cellular S^1 x RP^2 found above has', orientation_conflicts(G, r[0]),
          'orientation conflicts (non-orientable)')
    m = Model(G, induced_cycles(G, 19), verbose=False)
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls + add_orientation(m)), homology_on=False)
    assert r is None
    print('  K33 x K3, (H) off, (O) on: UNSAT')


def exp_pp_levels_orient(k):
    print(f'== P x P, no bent faces, (V)+(R)+(O) imposed only at the vertices of {k} of the 10 H-fibres')
    G = pp_graph()
    X = set(10 * u + v for u in range(10) for v in range(k))
    m = Model(G, nobent_candidates(G), X=X)
    t0 = time.time()
    s = Solver(name=SOLVER, bootstrap_with=m.cls + add_orientation(m, X=X))
    ok = s.solve()
    print(f'  {"SAT" if ok else "UNSAT"} ({time.time()-t0:.1f}s)')
    if ok:
        vs = set(l for l in s.get_model() if l > 0)
        ch = [m.faces[i] for i in range(len(m.faces)) if m.F(i) in vs]
        print('  independent orientation check of the witness on the constrained fibres:',
              orientation_conflicts(G, ch, X), 'conflicts')


def exp_k33c3():
    print('== validation: K33 x K3.  PPS (end of their paper) announce that the unique strongly regular combinatorial')
    print('   manifold with this graph is a cellular S^1 x RP^2 (6 triangular prisms + 9 cubes).')
    K = nx.convert_node_labels_to_integers(nx.complete_bipartite_graph(3, 3))
    G = product(K, nx.cycle_graph(3))
    cyc = induced_cycles(G, 19)
    print('  all induced cycles:', len(cyc), sorted(Counter(len(c) for c in cyc).items()))
    m = Model(G, cyc)
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls), homology_on=False)
    ch, fac = r
    print('  without (H): 2-faces', sorted(Counter(len(f) for f in ch).items()), 'facets',
          sorted(Counter(len(F['verts']) for F in fac).items()), '; (H) fails on it:',
          m.homology_cut([m.fidx[f] for f in ch]) is not None)
    assert sorted(Counter(len(F['verts']) for F in fac).items()) == [(6, 6), (8, 9)]
    m = Model(G, cyc)
    r = solve_full(m, Solver(name=SOLVER, bootstrap_with=m.cls))
    assert r is None
    print('  with (H): UNSAT -- K33 x K3 is not 4-polytopal (as PPS announce)')


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    T0 = time.time()
    if what in ('verify', 'all'): exp_verify()
    if what in ('sanity', 'all'): exp_sanity()
    if what in ('structure', 'all'): exp_structure()
    if what in ('pp-nobent-norot', 'all'): exp_pp_nobent_norot()
    if what in ('pp-prismgerm', 'all'): exp_pp_prismgerm()
    if what in ('pp-packing', 'all'): exp_pp_packing()
    if what in ('pp-levels',): exp_pp_levels(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else 0)
    if what in ('pc3',): exp_pc3(int(sys.argv[2]))
    if what in ('k33c3', 'all'): exp_k33c3()
    if what in ('orientation', 'all'): exp_orientation()
    if what in ('pp-packing-norot', 'all'): exp_pp_packing_norot()
    if what in ('pp-levels-orient',): exp_pp_levels_orient(int(sys.argv[2]))
    if what == 'all':
        exp_pp_levels(8)
        exp_pc3(10, crosscheck=True)
    if what in ('all-long',):
        exp_pp_levels(9)
        exp_pp_levels_orient(8)
        exp_pc3(31)
    print(f'ALL REQUESTED CHECKS DONE ({time.time()-T0:.0f}s)')
