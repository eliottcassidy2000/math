#!/usr/bin/env python3
"""port_probe_334.py -- sample (3,3,4;29) colourings by SAT (pysat) and record their colour-class edge counts and
degree data.  Colours 0,1 are triangle-free, colour 2 is K4-free.  Each found colouring is blocked (on its colour-2
edge set) and the solver is re-run, so successive solutions are pairwise distinct as labelled colourings.
These are the objects that appear as the 29-vertex ports of a hypothetical (3,3,5;56) or (3,4,4;76) colouring; if
their edge counts vary, no Boza-type forced count is available there.  Usage: port_probe_334.py NSOL TIMEOUT_S"""
import sys, time, itertools, json
from pysat.solvers import Cadical153, Glucose4
from pysat.card import CardEnc, EncType
n = 29; NSOL = int(sys.argv[1]) if len(sys.argv) > 1 else 6; TL = float(sys.argv[2]) if len(sys.argv) > 2 else 900
pairs = list(itertools.combinations(range(n), 2)); pid = {p: i for i, p in enumerate(pairs)}
def var(i, j, c):
    if i > j: i, j = j, i
    return 3 * pid[(i, j)] + c + 1
clauses = []
for p in pairs:
    i, j = p
    clauses.append([var(i, j, 0), var(i, j, 1), var(i, j, 2)])
    for c in range(3):
        for d in range(c + 1, 3): clauses.append([-var(i, j, c), -var(i, j, d)])
for a, b, c in itertools.combinations(range(n), 3):
    for col in (0, 1): clauses.append([-var(a, b, col), -var(b, c, col), -var(a, c, col)])
for q in itertools.combinations(range(n), 4):
    clauses.append([-var(x, y, 2) for x, y in itertools.combinations(q, 2)])
# symmetry breaking: vertex 0 has its colour-2 neighbours first, then colour-0, then colour-1 (a fixed ordering of its port types)
# (we only need samples, so keep it light: force 0-1 colour 2 and 0-28 colour 1)
clauses.append([var(0, 1, 2)]); clauses.append([var(0, 28, 1)])
MODE = sys.argv[3] if len(sys.argv) > 3 else 'plain'
top = 3 * len(pairs)
if MODE == 'deg7':   # ask for a colouring in which vertex 0 has colour-0 degree <= 7 (non-rigid degrees)
    lits = [var(0, j, 0) for j in range(1, n)]
    cnf = CardEnc.atmost(lits=lits, bound=7, top_id=top, encoding=EncType.seqcounter)
    clauses.extend(cnf.clauses); top = cnf.nv
if MODE == 'edges115':   # colour-0 edge count <= 115
    lits = [var(i, j, 0) for (i, j) in pairs]
    cnf = CardEnc.atmost(lits=lits, bound=115, top_id=top, encoding=EncType.seqcounter)
    clauses.extend(cnf.clauses); top = cnf.nv
print('mode', MODE, 'clauses', len(clauses), flush=True)
t0 = time.time(); results = []
s = (Glucose4 if MODE == 'glucose' else Cadical153)(bootstrap_with=clauses)
while len(results) < NSOL and time.time() - t0 < TL:
    ok = s.solve()
    if not ok: print('UNSAT after', len(results), 'solutions'); break
    m = set(l for l in s.get_model() if l > 0)
    col = {}
    for (i, j) in pairs:
        for c in range(3):
            if var(i, j, c) in m: col[(i, j)] = c
    e = [sum(1 for p in pairs if col[p] == c) for c in range(3)]
    deg = [[0]*n for _ in range(3)]
    for (i, j), c in col.items(): deg[c][i] += 1; deg[c][j] += 1
    degseq = [sorted(deg[c], reverse=True) for c in range(3)]
    results.append({'edges': e, 'degmax': [max(d) for d in degseq], 'degmin': [min(d) for d in degseq], 'deg2_hist': sorted(deg[2]), 't': round(time.time() - t0, 1)})
    print('solution %d at %.1fs: edges (c0,c1,c2) = %s, degree ranges %s..%s' % (len(results), time.time() - t0, e, [min(d) for d in degseq], [max(d) for d in degseq]), flush=True)
    # block this colour-2 edge set
    s.add_clause([-var(i, j, 2) if col[(i, j)] == 2 else var(i, j, 2) for (i, j) in pairs])
json.dump(results, open('port_probe_334_%s.json' % MODE, 'w'), indent=1)
print('distinct edge-count triples:', sorted({tuple(r['edges']) for r in results}))
