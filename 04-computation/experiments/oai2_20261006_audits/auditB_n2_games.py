#!/usr/bin/env python3
"""Independent audit (B) of the n = 2 tree-grid games (THM-453 E/F/G, THM-470, MISTAKE-577).

Three games on the grid [t]^2 (points (a, c), lex order):
  free   : one variable per unordered pair of points (THM-453 D, Q(2,t))
  row    : row-invariant (THM-453 F): within-row edge depends on the unordered column pair {c, c'},
           cross-row edge (a,c)-(a',c'), a < a', depends on (a'-a, c, c')
  gap    : fully gap-determined (THM-470 Finv): edge x <lex y depends only on y - x
Constraints, identical for all three games and generated from the POINTS (not from a
sum-free reformulation):
  * triangle-free: for every 3 points, not all three pair-variables true;
  * hitting: every binary subgrid {a<b} x ({c1<c2} for a, {c3<c4} for b) contains a true pair.
Two different SAT solvers (CaDiCaL 1.5.3 and Glucose 4) are run on each instance; SAT models are
re-verified by a direct brute-force check written separately (no use of the clause list).
"""
import itertools, sys, time
from pysat.solvers import Solver


def points(t):
    return [(a, c) for a in range(t) for c in range(t)]


def subgrids(t):
    pairs = list(itertools.combinations(range(t), 2))
    for a, b in itertools.combinations(range(t), 2):
        for p in pairs:
            for q in pairs:
                yield [(a, p[0]), (a, p[1]), (b, q[0]), (b, q[1])]


def key(game, x, y):
    if y < x:
        x, y = y, x
    if game == 'free':
        return ('p', x, y)
    if game == 'gap':
        return ('d', y[0] - x[0], y[1] - x[1])
    if game == 'row':
        if x[0] == y[0]:
            return ('R', x[1], y[1])           # x[1] < y[1] since x <lex y
        return ('B', y[0] - x[0], x[1], y[1])
    raise ValueError


def build(game, t):
    vid = {}
    def v(x, y):
        k = key(game, x, y)
        if k not in vid:
            vid[k] = len(vid) + 1
        return vid[k]
    P = points(t)
    clauses = set()
    for x, y, z in itertools.combinations(P, 3):
        cl = tuple(sorted({-v(x, y), -v(x, z), -v(y, z)}))
        clauses.add(cl)
    for g in subgrids(t):
        cl = tuple(sorted({v(x, y) for x, y in itertools.combinations(g, 2)}))
        clauses.add(cl)
    return vid, [list(c) for c in clauses]


def brute_verify(game, t, true_keys):
    """Direct check of a model: build the graph on points, test triangles and subgrids."""
    P = points(t)
    adj = {}
    for x, y in itertools.combinations(P, 2):
        adj[(x, y)] = adj[(y, x)] = key(game, x, y) in true_keys
    for x, y, z in itertools.combinations(P, 3):
        if adj[(x, y)] and adj[(x, z)] and adj[(y, z)]:
            return False, ('triangle', x, y, z)
    for g in subgrids(t):
        if not any(adj[(x, y)] for x, y in itertools.combinations(g, 2)):
            return False, ('independent subgrid', g)
    return True, None


def solve(game, t, solver_name):
    vid, cls = build(game, t)
    inv = {i: k for k, i in vid.items()}
    with Solver(name=solver_name, bootstrap_with=cls) as s:
        t0 = time.time()
        r = s.solve()
        dt = time.time() - t0
        model = s.get_model() if r else None
    true_keys = None
    if r:
        true_keys = {inv[l] for l in model if l > 0 and l in inv}
    return r, len(vid), len(cls), dt, true_keys


if __name__ == '__main__':
    out = []
    for game, ts in (('free', (3, 4, 5)), ('row', (3, 4, 5)), ('gap', (2, 3, 4, 5))):
        for t in ts:
            res = []
            for sn in ('cadical153', 'glucose4'):
                r, nv, nc, dt, tk = solve(game, t, sn)
                ver = ''
                if r:
                    ok, why = brute_verify(game, t, tk)
                    ver = f' brute-verified={ok}' + ('' if ok else f' {why}')
                    # also report graph edge count
                    P = points(t)
                    ne = sum(1 for x, y in itertools.combinations(P, 2) if key(game, x, y) in tk)
                    ver += f' |true classes|={len(tk)} graph-edges={ne}'
                res.append(f'{sn}: {"SAT" if r else "UNSAT"} ({dt:.2f}s){ver}')
            line = f'{game:4s} t={t}: vars={nv} clauses={nc} | ' + ' | '.join(res)
            print(line, flush=True)
