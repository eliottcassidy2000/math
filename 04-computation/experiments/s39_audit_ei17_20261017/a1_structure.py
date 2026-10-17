#!/usr/bin/env python3
"""Audit of THM-4665, item 1: direction structure of EI17 derived independently.

Independent derivation (does not use the VERT tables of the S39 scripts until the final comparison):
  (a) graph facts of PE (17 vertices, 31 edges, triangle-free, |Aut| = 2, chromatic number 4, Laman count,
      isomorphic to the House of Graphs adjacency of the S24 script);
  (b) the unit 4-cycles and the signed union-find closure of the parallelogram identifications
      (vec(a->b) = vec(d->c), vec(b->c) = vec(a->d) for every 4-cycle a-b-c-d);
  (c) the edge directions of the certified realization (24-digit TABLE, refined by Newton to 80 digits), clustered
      up to sign; comparison with (b);
  (d) a BFS spanning tree from vertex 17 gives every vertex as an integer combination of class variables, the
      non-tree edges give the cycle relations; rank, and the solution for the six dependent classes in terms of the
      script's free classes (sympy, exact rationals: integrality is checked, not assumed);
  (e) only now: comparison with the VERT table and the six displayed formulas of THM-4665, numerically at the
      certified realization, and the degenerate-branch collision identities.
Run: python -B a1_structure.py      (from anywhere; paths are resolved from this file)
"""
import importlib.util
import itertools
import os
os.environ.setdefault('NUMBA_CACHE_DIR', os.path.join(os.environ.get('TEMP', os.environ.get('TMPDIR', '/tmp')), 'nbc_s39_audit_ei17'))  # short path outside the repo (Windows MAX_PATH)
import sys

import mpmath as mp
import networkx as nx
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(EXP, name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S24 = load('ei17_placement_q_zeta55_20261012')
SPL = load('ei17_finite_field_realizations_20261017')
TW = load('ei17_twisted_unit_graph_embeddings_20261017')
PE, TABLE, HOG = S24.PE, S24.TABLE, S24.HOG
assert PE == SPL.PE == TW.PE, 'edge lists differ between the scripts'
assert SPL.VERT == TW.VERT, 'VERT tables differ between the two S39 scripts'
print('edge list identical in the three scripts; VERT identical in the two S39 scripts')

# ------------------------------------------------------------------ (a) graph facts
G = nx.Graph(PE)
assert G.number_of_nodes() == 17 and G.number_of_edges() == 31
tri = sum(nx.triangles(G).values()) // 3
H = nx.Graph()
for i, nb in enumerate(HOG):
    for j in nb:
        H.add_edge(i, j)
iso = nx.algorithms.isomorphism.GraphMatcher(G, H)
n_iso = sum(1 for _ in iso.isomorphisms_iter())
n_aut = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())
autos = [m for m in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter()]
nontriv = [m for m in autos if any(m[v] != v for v in m)][0]
moved = {v: nontriv[v] for v in sorted(nontriv) if nontriv[v] != v}


def colourable(k):
    from pysat.solvers import Glucose3
    var = lambda v, c: (v - 1) * k + c + 1
    s = Glucose3()
    for v in range(1, 18):
        s.add_clause([var(v, c) for c in range(k)])
    for i, j in PE:
        for c in range(k):
            s.add_clause([-var(i, c), -var(j, c)])
    ok = s.solve()
    s.delete()
    return ok


laman_ok = True
for r in range(2, 17):
    for S in itertools.combinations(range(1, 18), r):
        if G.subgraph(S).number_of_edges() > 2 * r - 3:
            laman_ok = False
            break
    if not laman_ok:
        break
print('(a) triangles %d, isomorphisms PE -> HoG %d, |Aut| %d, nontrivial automorphism moves %s' % (tri, n_iso, n_aut, moved))
print('    3-colourable %s, 4-colourable %s, Laman (every r-subset spans <= 2r-3 edges, 31 = 2*17-3) %s'
      % (colourable(3), colourable(4), laman_ok))

# ------------------------------------------------------------------ (b) 4-cycles and signed closure
adj = {v: set(G[v]) for v in G}
cyc = set()
for a in adj:
    for b in adj[a]:
        for c in adj[b]:
            if c == a:
                continue
            for d in adj[c]:
                if d in (a, b) or a not in adj[d]:
                    continue
                t = (a, b, c, d)
                rots = [t[i:] + t[:i] for i in range(4)]
                rots += [tuple(reversed(x)) for x in rots]
                cyc.add(min(rots))
cyc = sorted(cyc)
print('(b) unit 4-cycles: %d  %s' % (len(cyc), cyc))

EIDX = {e: k for k, e in enumerate(PE)}


def oriented(x, y):
    """(edge index, sign) with vec(x->y) = sign * D(edge), D(i,j) = v_j - v_i for i < j."""
    return (EIDX[(x, y)], 1) if (x, y) in EIDX else (EIDX[(y, x)], -1)


parent = list(range(31)); psign = [1] * 31   # D(e) = psign[e] * D(parent[e])


def find(e):
    if parent[e] == e:
        return e, 1
    r, s = find(parent[e])
    parent[e] = r; psign[e] = psign[e] * s
    return r, psign[e]


conflict = False


def union(e1, s1, e2, s2):
    """impose s1 D(e1) = s2 D(e2)."""
    global conflict
    r1, t1 = find(e1); r2, t2 = find(e2)
    # D(e1) = t1 D(r1), D(e2) = t2 D(r2): s1 t1 D(r1) = s2 t2 D(r2)
    if r1 == r2:
        if s1 * t1 != s2 * t2:
            conflict = True
        return
    parent[r1] = r2; psign[r1] = s1 * t1 * s2 * t2


for a, b, c, d in cyc:
    e_ab, s_ab = oriented(a, b); e_dc, s_dc = oriented(d, c)
    e_bc, s_bc = oriented(b, c); e_ad, s_ad = oriented(a, d)
    union(e_ab, s_ab, e_dc, s_dc)
    union(e_bc, s_bc, e_ad, s_ad)
roots = sorted({find(e)[0] for e in range(31)})
UF = {e: find(e) for e in range(31)}
print('    parallelogram closure: %d classes, sign conflict: %s' % (len(roots), conflict))

# ------------------------------------------------------------------ (c) certified realization
mp.mp.dps = 100
FREEV = [v for v in range(1, 18) if v not in (11, 17)]
E30 = [e for e in PE if e != (11, 17)]
X = mp.matrix([mp.mpf(TABLE[v][c]) for v in FREEV for c in (0, 1)])


def unpack(X):
    z = {11: (mp.mpf(-1), mp.mpf(0)), 17: (mp.mpf(0), mp.mpf(0))}
    for k, v in enumerate(FREEV):
        z[v] = (X[2 * k], X[2 * k + 1])
    return z


for it in range(12):
    z = unpack(X)
    r = mp.matrix([(z[i][0] - z[j][0]) ** 2 + (z[i][1] - z[j][1]) ** 2 - 1 for i, j in E30])
    J = mp.matrix(30, 30); col = {v: k for k, v in enumerate(FREEV)}
    for rr, (i, j) in enumerate(E30):
        dx = 2 * (z[i][0] - z[j][0]); dy = 2 * (z[i][1] - z[j][1])
        if i in col:
            J[rr, 2 * col[i]] += dx; J[rr, 2 * col[i] + 1] += dy
        if j in col:
            J[rr, 2 * col[j]] -= dx; J[rr, 2 * col[j] + 1] -= dy
    X = X - mp.lu_solve(J, r)
z = unpack(X)
res = max(abs((z[i][0] - z[j][0]) ** 2 + (z[i][1] - z[j][1]) ** 2 - 1) for i, j in PE)
shift = max(abs(z[v][c] - mp.mpf(TABLE[v][c])) for v in TABLE for c in (0, 1))
P = {v: mp.mpc(z[v][0], z[v][1]) for v in z}
mind = min(abs(P[v] - P[w]) for v in P for w in P if v < w)
print('(c) Newton from the 24-digit table: max edge residual %s, max shift from the table %s, min vertex distance %s'
      % (mp.nstr(res, 3), mp.nstr(shift, 3), mp.nstr(mind, 6)))
D = [P[j] - P[i] for i, j in PE]
reps = []; ncls = {}
for k, u in enumerate(D):
    for c, r_ in enumerate(reps):
        if abs(u - r_) < mp.mpf(10) ** -60:
            ncls[k] = (c, 1); break
        if abs(u + r_) < mp.mpf(10) ** -60:
            ncls[k] = (c, -1); break
    else:
        reps.append(u); ncls[k] = (len(reps) - 1, 1)
sep = min(min(abs(reps[a] - reps[b]), abs(reps[a] + reps[b])) for a in range(len(reps)) for b in range(a + 1, len(reps)))
print('    numerical classes up to sign: %d; least separation between distinct classes (min |u - w|, |u + w|) %s'
      % (len(reps), mp.nstr(sep, 6)))
# same partition and same relative signs as the union-find?
same = True
for e1 in range(31):
    for e2 in range(31):
        r1, s1 = UF[e1]; r2, s2 = UF[e2]
        c1, t1 = ncls[e1]; c2, t2 = ncls[e2]
        if (r1 == r2) != (c1 == c2):
            same = False
        elif r1 == r2 and s1 * s2 != t1 * t2:
            same = False
print('    numerical partition and signs = parallelogram closure: %s' % same)
# no other near-coincidence: least |u -+ w| between classes at 100 digits is far from 0 (printed above)

# ------------------------------------------------------------------ (d) relations from a spanning tree
NC = len(reps)
Xs = sp.symbols('X0:%d' % NC)
gauge_cls = ncls[EIDX[(11, 17)]][0]
pos = {17: sp.Integer(0)}
order = [17]; tree = set()
while len(pos) < 17:
    for v in list(order):
        for w in sorted(adj[v]):
            if w in pos:
                continue
            e, s = oriented(v, w)
            c, t = ncls[e]
            pos[w] = sp.expand(pos[v] + s * t * Xs[c]); tree.add(e); order.append(w)
rels = []
for k, (i, j) in enumerate(PE):
    if k in tree:
        continue
    c, t = ncls[k]
    r_ = sp.expand(pos[j] - pos[i] - t * Xs[c])
    if r_ != 0:
        rels.append(r_)
M = sp.Matrix([[r_.coeff(x) for x in Xs] for r_ in rels])
print('(d) spanning tree: %d non-tree edges, %d nontrivial cycle relations, rank %d (cycle space dim %d, minus 9 parallelograms)'
      % (31 - len(tree), len(rels), M.rank(), 31 - 17 + 1))

# ------------------------------------------------------------------ (e) comparison with the S39 scripts
VERT = SPL.VERT
BASIS = ['1', 'U0', 'U1', 'U2', 'U3', 'U4', 'U6']
FORM = {'U5': (1, 0, 0, 0, 0, 0, 0), 'U0': (0, 1, 0, 0, 0, 0, 0), 'U1': (0, 0, 1, 0, 0, 0, 0), 'U2': (0, 0, 0, 1, 0, 0, 0),
        'U3': (0, 0, 0, 0, 1, 0, 0), 'U4': (0, 0, 0, 0, 0, 1, 0), 'U6': (0, 0, 0, 0, 0, 0, 1),
        # the six displayed formulas of THM-4665 in the basis (1, U0, U1, U2, U3, U4, U6)
        'U7': (-1, 0, 1, 1, -1, 0, 0), 'U8': (1, -1, 0, -1, 0, 1, 0), 'U9': (1, -1, 1, -1, 1, 0, 1),
        'U10': (-1, 0, 0, 2, -1, 0, -1), 'U11': (-1, 0, 1, 1, 1, -1, 1), 'U12': (2, 0, -1, -1, 0, 0, -1)}
lab = {}
for k, (i, j) in enumerate(PE):
    d = tuple(VERT[j][t] - VERT[i][t] for t in range(7))
    hit = [(nm, 1) for nm, f in FORM.items() if f == d] + [(nm, -1) for nm, f in FORM.items() if tuple(-x for x in f) == d]
    assert len(hit) == 1, ('edge', (i, j), 'VERT difference', d, 'is not +-one of the 13 forms')
    lab[k] = hit[0]
print('(e) every VERT edge difference is +- exactly one of the 13 forms (U0..U6 free/gauge, U7..U12 displayed): True')
# labels vs numerical classes: a bijection of classes with consistent signs
cmap = {}
ok_map = True
for k in range(31):
    c, t = ncls[k]; nm, s = lab[k]
    if c in cmap:
        if cmap[c][0] != nm or cmap[c][1] != s * t:
            ok_map = False
    else:
        cmap[c] = (nm, s * t)
print('    class -> label map consistent and bijective: %s (%d labels)' % (ok_map and len(set(v[0] for v in cmap.values())) == 13, len(set(v[0] for v in cmap.values()))))
classes_by_label = {}
for k in range(31):
    classes_by_label.setdefault(lab[k][0], []).append(PE[k])
for nm in ['U0', 'U1', 'U2', 'U3', 'U4', 'U5', 'U6', 'U7', 'U8', 'U9', 'U10', 'U11', 'U12']:
    print('      %-4s edges %s' % (nm, classes_by_label[nm]))
# numerical values of the labelled directions (sign: D(edge) = s * U)
Uval = {}
for k in range(31):
    nm, s = lab[k]
    val = s * D[k]
    if nm in Uval:
        assert abs(Uval[nm] - val) < mp.mpf(10) ** -60
    Uval[nm] = val
print('    U5 (gauge) = %s' % mp.nstr(Uval['U5'], 8))
vals = [mp.mpf(1), Uval['U0'], Uval['U1'], Uval['U2'], Uval['U3'], Uval['U4'], Uval['U6']]
err_form = max(abs(sum(FORM[nm][t] * vals[t] for t in range(7)) - Uval[nm]) for nm in FORM)
err_vert = max(abs(sum(VERT[v][t] * vals[t] for t in range(7)) - P[v]) for v in range(1, 18))
err_unit = max(abs(abs(Uval[nm]) - 1) for nm in Uval)
print('    at the certified realization: max |formula - U| %s, max |VERT - vertex| %s, max ||U| - 1| %s'
      % (mp.nstr(err_form, 3), mp.nstr(err_vert, 3), mp.nstr(err_unit, 3)))
# exact: solve the cycle relations for the dependent classes in terms of the script's free ones
lab_of_cls = {c: cmap[c] for c in cmap}      # class c: X_c = sgn * U_label
Usym = {nm: sp.Symbol(nm) for nm in FORM}
subs = {Xs[c]: lab_of_cls[c][1] * Usym[lab_of_cls[c][0]] for c in range(NC)}
eqs = [sp.expand(r_.subs(subs)) for r_ in rels]
dep = [Usym['U%d' % k] for k in range(7, 13)]
sol = sp.solve(eqs, dep, dict=True)
assert len(sol) == 1
sol = sol[0]
ok_formulas = True
for k in range(7, 13):
    nm = 'U%d' % k
    expr = sp.expand(sol[Usym[nm]].subs(Usym['U5'], 1))
    want = sum(FORM[nm][t] * (sp.Integer(1) if t == 0 else Usym[BASIS[t]]) for t in range(7))
    integral = all(sp.Rational(cf).q == 1 for cf in sp.Poly(expr, *[Usym[b] for b in BASIS[1:]]).coeffs())
    print('    exact solve: %-3s = %-36s integral %s  matches displayed formula %s' % (nm, expr, integral, sp.expand(expr - want) == 0))
    ok_formulas &= integral and sp.expand(expr - want) == 0
Mdep = sp.Matrix([[sp.expand(e).coeff(d) for d in dep] for e in eqs])
print('    coefficient matrix of the six dependent classes in the six relations: rank %d, det %s'
      % (Mdep.rank(), Mdep.det() if Mdep.shape == (6, 6) else 'n/a (non-square: %s)' % (Mdep.shape,)))
# vertex formulas from the tree, in script labels and gauge, vs VERT
ok_vert = True
for v in range(1, 18):
    expr = sp.expand(pos[v].subs(subs).subs(Usym['U5'], 1))
    for k in range(7, 13):
        expr = sp.expand(expr.subs(Usym['U%d' % k], sol[Usym['U%d' % k]].subs(Usym['U5'], 1)))
    want = sum(VERT[v][t] * (sp.Integer(1) if t == 0 else Usym[BASIS[t]]) for t in range(7))
    if sp.expand(expr - want) != 0:
        ok_vert = False
        print('      vertex %d: tree gives %s, VERT gives %s' % (v, expr, want))
print('    the 17 vertex formulas derived from the tree equal the VERT table: %s' % ok_vert)
print('    all six displayed dependent formulas reproduced exactly: %s' % ok_formulas)

# degenerate branches: c = 0 of each quadratic forces the stated collision
U0, U1, U2, U3, U4, U6 = [Usym[b] for b in BASIS[1:]]
vt = {v: sum(VERT[v][t] * (sp.Integer(1) if t == 0 else Usym[BASIS[t]]) for t in range(7)) for v in VERT}
checks = [('U1 + U2 - 1 (U3 from |U7| = 1)', U1 + U2 - 1, (5, 11)),
          ('2 - U1 - U2 (U6 from |U12| = 1)', 2 - U1 - U2, (7, 10)),
          ('U1 + U2 + U3 + U6 - 1 (U4 from |U11| = 1)', U1 + U2 + U3 + U6 - 1, (2, 11)),
          ('1 - U2 + U4 (U0 from |U8| = 1)', 1 - U2 + U4, (1, 6))]
for txt, cvec, (a, b) in checks:
    diff = sp.expand(vt[a] - vt[b])
    prop = sp.simplify(diff / cvec)
    print('    constant vector %-44s v%d - v%d = %-24s = (%s) * c' % (txt, a, b, diff, prop))
# which vertex pairs differ by a constant multiple of each c (all collisions forced on that branch)
for txt, cvec, _ in checks:
    forced = [(a, b) for a in range(1, 18) for b in range(a + 1, 18)
              if sp.simplify(sp.expand(vt[a] - vt[b]) / cvec).is_number and sp.expand(vt[a] - vt[b]) != 0]
    print('      pairs proportional to that c: %s' % forced)
# collision forms: pairs of vertices whose difference is a single unit or zero identically (never happens)
zero_pairs = [(a, b) for a in range(1, 18) for b in range(a + 1, 18) if sp.expand(vt[a] - vt[b]) == 0]
print('    identically coincident vertex pairs: %s' % zero_pairs)
