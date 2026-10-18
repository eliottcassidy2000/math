#!/usr/bin/env python3
"""S42 audit a7: six-colour coupled-family and reflection-forcing checks.

1. Sandbox re-runs of the session scripts golden_forcing_coupled_20261017.py (base base, base MM) and
   golden_reflection_forcing_20261017.py (C, Cd): module __file__ redirected so that their JSON goes to
   s42_audit_20261017/sb/05-knowledge/results (no session file touched); compared with the committed JSON.
2. Independent re-encoding (base x base, all three families): families rebuilt here from the field arithmetic of the
   S41 Geo class; additive constraints encoded DIRECTLY (one clause per forbidden 6-tuple of colours, no auxiliary
   difference variables); blocked count over the 570 patterns; for 20 patterns the SAT model is decoded and every
   axis and mixed constraint re-checked by plain evaluation.
3. Positive control: plain 6-colouring of the golden two-distance graph on A (unit edges + d q steps) should block
   465 of 570 patterns (S41).
"""
import os, sys, json, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
sys.dont_write_bytecode = True
sys.path.insert(0, EXP)
from pysat.solvers import Solver
from golden_boundary_20261017c import PORT_EDGES, patterns
from golden_forcing_search_20261017 import Geo, recipe

SB = os.path.join(HERE, 'sb', 'x', 'y')
RES_SB = os.path.join(HERE, 'sb', '05-knowledge', 'results')
os.makedirs(SB, exist_ok=True); os.makedirs(RES_SB, exist_ok=True)
report = {}


def sandbox_run(modname, argv, outname):
    mod = __import__(modname)
    mod.__file__ = os.path.join(SB, modname + '.py')
    old = sys.argv
    sys.argv = [modname + '.py'] + argv
    try:
        mod.main()
    finally:
        sys.argv = old
    new = json.load(open(os.path.join(RES_SB, outname)))
    ref = json.load(open(os.path.join(EXP, '..', '..', '05-knowledge', 'results', outname)))
    return new, ref


t0 = time.time()
cmp = {}
for a, c in (('base', 'base'), ('base', 'MM')):
    new, ref = sandbox_run('golden_forcing_coupled_20261017', [a, c], 'golden_forcing_coupled_%s_%s_20261017.json' % (a, c))
    keys = ('A', 'C', 'edges_A', 'edges_C', 'directions', 'v_contact', 'reflected', 'e2_contact')
    same = all(new[k] == ref[k] for k in keys) and all(new['results'][s]['blocked'] == ref['results'][s]['blocked']
                                                        for s in ('shared', 'independent'))
    cmp['coupled_%s_%s' % (a, c)] = dict(rerun={k: new[k] for k in keys}, blocked={s: new['results'][s]['blocked'] for s in new['results']},
                                         identical_to_committed=same)
for w in ('C', 'Cd'):
    new, ref = sandbox_run('golden_reflection_forcing_20261017', [w], 'golden_reflection_forcing_%s_20261017.json' % w)
    same = (new['W'] == ref['W'] and new['edges_W'] == ref['edges_W'] and new['reflected'] == ref['reflected'] and
            all(new['results'][s]['blocked'] == ref['results'][s]['blocked'] for s in ('shared', 'independent')))
    cmp['reflection_%s' % w] = dict(W=new['W'], edges_W=new['edges_W'], reflected=new['reflected'],
                                    blocked={s: new['results'][s]['blocked'] for s in new['results']},
                                    identical_to_committed=same,
                                    difference_sets_identical=all(new['results'][s]['w_difference_sets'] == ref['results'][s]['w_difference_sets']
                                                                  for s in ('shared', 'independent')))
report['sandbox_reruns'] = cmp
print('sandbox re-runs: %s [%.0f s]' % (json.dumps(cmp), time.time() - t0), flush=True)

# ---- independent re-encoding, base x base ----
geo = Geo(); f = geo.f; d = geo.d; dbar = geo.bar(d)
A, C = recipe(geo, 'base')
assert C[0] == f.zero
EA = geo.unit_edges(A); EC = geo.unit_edges(C)
IA = {x: i for i, x in enumerate(A)}; IC = {x: i for i, x in enumerate(C)}
qdirs = {}
for i in range(len(C)):
    for j in range(len(C)):
        if i != j:
            q = f.sub(C[j], C[i])
            if geo.is_unit(q):
                qdirs.setdefault(q, []).append((i, j))
famV, famR, famE = [], [], []
for q, cp in qdirs.items():
    dq = f.mul(d, q); dbq = f.mul(dbar, q)
    ap = [(i, IA[f.add(A[i], dq)]) for i in range(len(A)) if f.add(A[i], dq) in IA]
    lp = [(i, IC[f.add(C[i], dbq)]) for i in range(len(C)) if f.add(C[i], dbq) in IC]
    for (b, b2) in cp:
        for (x, x2) in ap:
            for (c, c2) in cp:
                famV.append(((x, x2), (b, b2), (c, c2)))        # (a,b,c) -> (a+dq, b+q, c+q); w, l fixed
            for (l, l2) in lp:
                famE.append(((x, x2), (b, b2), (l, l2)))        # (a,b,l) -> (a+dq, b+q, l+conj(d)q); c, w fixed
        for (c, c2) in cp:
            for (l, l2) in lp:
                famR.append(((b, b2), (c, c2), (l, l2)))        # (b,c,w) -> (b+q, c+q, w+conj(d)q); a, l fixed
print('independent families (directed): v-contact %d, reflected %d, e^2 %d; unit directions %d' %
      (len(famV), len(famR), len(famE), len(qdirs)), flush=True)
nA, nC = len(A), len(C)
pats = patterns(8, PORT_EDGES)


def encode(shared):
    var = {}
    def v(layer, i, k):
        key = (layer, i, k)
        if key not in var:
            var[key] = len(var) + 1
        return var[key]
    layers = ['A', '1', '3', '4'] if shared else ['A', '1', '2', '3', '4']
    size = {'A': nA, '1': nC, '2': nC, '3': nC, '4': nC}
    cl = []
    for L in layers:
        for i in range(size[L]):
            cl.append([v(L, i, k) for k in range(6)])
            for k1, k2 in itertools.combinations(range(6), 2):
                cl.append([-v(L, i, k1), -v(L, i, k2)])
    for (i, j) in EA:
        for k in range(6):
            cl.append([-v('A', i, k), -v('A', j, k)])
    for L in layers[1:]:
        for (i, j) in EC:
            for k in range(6):
                cl.append([-v(L, i, k), -v(L, j, k)])
        cl.append([v(L, 0, 0)])
    L2 = '1' if shared else '2'
    def forbid(t1, t2, t3):
        # t = (layer, i, i2): forbid every colouring with total shift 0 mod 6
        (La, i1, j1), (Lb, i2, j2), (Lc, i3, j3) = t1, t2, t3
        for x1, y1, x2, y2, x3 in itertools.product(range(6), repeat=5):
            y3 = (x3 - ((y1 - x1) + (y2 - x2))) % 6        # (y1-x1) + (y2-x2) + (y3-x3) = 0
            lits = [-v(La, i1, x1), -v(La, j1, y1), -v(Lb, i2, x2), -v(Lb, j2, y2), -v(Lc, i3, x3), -v(Lc, j3, y3)]
            cl.append(sorted(set(lits)))
    for (aa, bb, cc) in famV:
        forbid(('A',) + aa, ('1',) + bb, (L2,) + cc)
    for (bb, cc, ll) in famR:
        forbid(('1',) + bb, (L2,) + cc, ('3',) + ll)
    for (aa, bb, ll) in famE:
        forbid(('A',) + aa, ('1',) + bb, ('4',) + ll)
    return var, cl, layers


def decode(model, var, layers):
    ms = set(x for x in model if x > 0)
    col = {}
    for (L, i, k), x in var.items():
        if x in ms:
            col[(L, i)] = k
    return col


def check(col, shared, pat):
    L2 = '1' if shared else '2'
    ok = all(col[('A', i)] == pat[i] for i in range(8))
    ok &= all(col[('A', i)] != col[('A', j)] for i, j in EA)
    for L in (['1', '3', '4'] if shared else ['1', '2', '3', '4']):
        ok &= all(col[(L, i)] != col[(L, j)] for i, j in EC) and col[(L, 0)] == 0
    sh = lambda L, i, j: col[(L, j)] - col[(L, i)]
    ok &= all((sh('A', *aa) + sh('1', *bb) + sh(L2, *cc)) % 6 != 0 for aa, bb, cc in famV)
    ok &= all((sh('1', *bb) + sh(L2, *cc) + sh('3', *ll)) % 6 != 0 for bb, cc, ll in famR)
    ok &= all((sh('A', *aa) + sh('1', *bb) + sh('4', *ll)) % 6 != 0 for aa, bb, ll in famE)
    return ok


indep = {}
for shared in (True, False):
    var, cl, layers = encode(shared)
    blocked = 0; checked = 0; check_ok = True
    with Solver(name='cadical195', bootstrap_with=cl) as s:
        for n, pat in enumerate(pats):
            if not s.solve(assumptions=[var[('A', i, c)] for i, c in enumerate(pat)]):
                blocked += 1
            elif n % 29 == 0:
                col = decode(s.get_model(), var, layers)
                check_ok &= check(col, shared, pat)
                checked += 1
    key = 'shared' if shared else 'independent'
    indep[key] = dict(clauses=len(cl), blocked=blocked, models_checked=checked, all_checked_models_valid=bool(check_ok))
    print('independent direct encoding (%s): %d clauses; blocked %d of %d; %d decoded models re-checked by evaluation: valid %s [%.0f s]'
          % (key, len(cl), blocked, len(pats), checked, check_ok, time.time() - t0), flush=True)
report['independent_base_base'] = dict(families=dict(v=len(famV), reflected=len(famR), e2=len(famE)), results=indep)

# ---- positive control: two-distance graph on A ----
dsteps = set()
for q in qdirs:
    dq = f.mul(d, q)
    for i in range(nA):
        y = f.add(A[i], dq)
        if y in IA:
            dsteps.add(tuple(sorted((i, IA[y]))))
edges2 = set(tuple(sorted(e)) for e in EA) | dsteps
blocked2 = 0
with Solver(name='cadical195') as s:
    for i in range(nA):
        s.add_clause([6 * i + k + 1 for k in range(6)])
    for i, j in edges2:
        for k in range(6):
            s.add_clause([-(6 * i + k + 1), -(6 * j + k + 1)])
    for pat in pats:
        if not s.solve(assumptions=[6 * i + c + 1 for i, c in enumerate(pat)]):
            blocked2 += 1
report['two_distance_control'] = dict(unit_edges=len(EA), d_step_edges=len(dsteps), blocked=blocked2)
print('positive control: two-distance graph on A (%d unit + %d d-step edges): blocked %d of 570 (S41: 465)'
      % (len(EA), len(dsteps), blocked2), flush=True)
json.dump(report, open(os.path.join(HERE, 'a7_sixcolour_checks.json'), 'w'), indent=1)
