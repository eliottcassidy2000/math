#!/usr/bin/env python3
"""Audit B: two more n = 2 games, same point-level encoding as auditB_n2_games.py:
  F2   : gap-determined (sign, v_2) per coordinate (THM-465/469 bare 2-adic valuation grading)
  dyad : THM-453 G 'dyadic' row-invariant family: within-row R(c,c') arbitrary, cross-row B depends on (v_2(a'-a), c, c')
"""
import itertools
from pysat.solvers import Solver
import auditB_n2_games as G

def v2(x):
    x = abs(x); v = 0
    while x % 2 == 0:
        x //= 2; v += 1
    return v

def feat(d):
    return ('=',) if d == 0 else ((1 if d > 0 else -1), v2(d))

def key(game, x, y):
    if y < x:
        x, y = y, x
    if game == 'F2':
        return ('f', feat(y[0] - x[0]), feat(y[1] - x[1]))
    if game == 'dyad':
        if x[0] == y[0]:
            return ('R', x[1], y[1])
        return ('B', v2(y[0] - x[0]), x[1], y[1])

for game in ('F2', 'dyad'):
    for t in (3, 4, 5, 6):
        vid = {}
        def v(x, y):
            k = key(game, x, y)
            if k not in vid:
                vid[k] = len(vid) + 1
            return vid[k]
        P = G.points(t)
        cls = set()
        for x, y, z in itertools.combinations(P, 3):
            cls.add(tuple(sorted({-v(x, y), -v(x, z), -v(y, z)})))
        for g in G.subgrids(t):
            cls.add(tuple(sorted({v(x, y) for x, y in itertools.combinations(g, 2)})))
        with Solver(name='cadical153', bootstrap_with=[list(c) for c in cls]) as s:
            r = s.solve()
        print(f'{game} n=2 t={t}: {"SAT" if r else "UNSAT"} ({len(vid)} classes)', flush=True)
        if not r:
            break
