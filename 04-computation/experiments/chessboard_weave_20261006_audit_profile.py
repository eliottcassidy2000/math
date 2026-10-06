#!/usr/bin/env python3
"""Independent audit (written from scratch) of the 'solver-optimal ring profile of closed 8x8 tours' table in
Section 5 of 05-knowledge/results/chessboard_weave_20261006.md.
CP-SAT (2 workers) min/max of each quantity over closed knight's tours; every optimum tour is re-validated.
Status OPTIMAL = optimality proved by the solver.
"""
from collections import Counter
from ortools.sat.python import cp_model

KN = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
N = 8


def ring(i, j):
    return (max(abs(2 * i - 7), abs(2 * j - 7)) - 1) // 2


V = [(i, j) for i in range(N) for j in range(N)]
idx = {v: k for k, v in enumerate(V)}
E = sorted(set(tuple(sorted([(i, j), (i + dx, j + dy)])) for i in range(N) for j in range(N) for dx, dy in KN
               if 0 <= i + dx < N and 0 <= j + dy < N))


def solve(expr_fn, sense, tl=180.0):
    m = cp_model.CpModel()
    x, arcs = {}, []
    for (u, v) in E:
        for a, b in ((u, v), (v, u)):
            lit = m.NewBoolVar("")
            x[(a, b)] = lit
            arcs.append((idx[a], idx[b], lit))
    m.AddCircuit(arcs)
    ms = {}
    for (u, v) in E:
        key = tuple(sorted((ring(*u), ring(*v))))
        ms.setdefault(key, []).extend([x[(u, v)], x[(v, u)]])
    ms = {k: sum(vs) for k, vs in ms.items()}
    expr = expr_fn(ms)
    (m.Minimize if sense == 'min' else m.Maximize)(expr)
    s = cp_model.CpSolver()
    s.parameters.num_workers = 2
    s.parameters.max_time_in_seconds = tl
    st = s.Solve(m)
    tour = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        succ = {a: b for (a, b), lit in x.items() if s.Value(lit)}
        tour = [V[0]]
        while len(tour) < 64:
            tour.append(succ[tour[-1]])
        assert len(set(tour)) == 64 and all(sorted([abs(tour[t][0] - tour[(t + 1) % 64][0]), abs(tour[t][1] - tour[(t + 1) % 64][1])]) == [1, 2] for t in range(64))
        prof = Counter(tuple(sorted((ring(*tour[t]), ring(*tour[(t + 1) % 64])))) for t in range(64))
        assert prof[(3, 3)] == prof[(1, 1)] + prof[(1, 2)] + prof[(2, 2)]
    return s.StatusName(st), s.ObjectiveValue(), s.BestObjectiveBound()


Q = {
    'within-ring': (lambda ms: ms[(1, 1)] + ms[(2, 2)] + ms[(3, 3)], (1, 16)),
    'two-ring': (lambda ms: ms[(0, 2)] + ms[(1, 3)], (2, 32)),
    'm01': (lambda ms: ms[(0, 1)], (0, 8)),
    'm02': (lambda ms: ms[(0, 2)], (0, 8)),
    'm11': (lambda ms: ms[(1, 1)], (0, 6)),
    'm22': (lambda ms: ms[(2, 2)], (0, 8)),
    'm33': (lambda ms: ms[(3, 3)], (1, 8)),
    'm12': (lambda ms: ms[(1, 2)], (0, 8)),
    'm13': (lambda ms: ms[(1, 3)], (2, 24)),
    'm23': (lambda ms: ms[(2, 3)], (17, 40)),
}
for name, (fn, (cmin, cmax)) in Q.items():
    out = []
    for sense, claim in (('min', cmin), ('max', cmax)):
        st, val, bnd = solve(fn, sense)
        out.append(f"{sense}={val:g} ({st}, bound {bnd:g}, note {claim}, {'AGREE' if st == 'OPTIMAL' and val == claim else 'CHECK'})")
    print(f"{name}: " + "; ".join(out), flush=True)
