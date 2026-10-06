#!/usr/bin/env python3
"""Which weave numbers k = #(tour moves inside ring 3) occur among closed 8x8 knight tours? (CP-SAT feasibility per k)"""
from ortools.sat.python import cp_model
N = 8
def ring(i, j):
    c = (N - 1) / 2
    return int(max(abs(i - c), abs(j - c)) - 0.5)
V = [(i, j) for i in range(N) for j in range(N)]
idx = {v: k for k, v in enumerate(V)}
arcs = []
for (i, j) in V:
    for a, b in ((1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)):
        u = (i + a, j + b)
        if 0 <= u[0] < N and 0 <= u[1] < N:
            arcs.append((idx[(i, j)], idx[u]))
for k in range(0, 9):
    m = cp_model.CpModel()
    lit = {a: m.NewBoolVar("") for a in arcs}
    m.AddCircuit([(u, v, lit[(u, v)]) for (u, v) in arcs])
    inner3 = [lit[(u, v)] for (u, v) in arcs if ring(*V[u]) == 3 and ring(*V[v]) == 3]
    mid = [lit[(u, v)] for (u, v) in arcs if ring(*V[u]) in (1, 2) and ring(*V[v]) in (1, 2)]
    m.Add(sum(inner3) == k)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 60
    st = s.Solve(m)
    extra = ""
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        extra = f"; middle-annulus moves in the witness = {sum(s.Value(x) for x in mid)}"
    print(f"k = {k}: {s.StatusName(st)}{extra}")
