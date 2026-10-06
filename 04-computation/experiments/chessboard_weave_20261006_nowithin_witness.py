#!/usr/bin/env python3
"""Look for a short certificate that the 8x8 knight graph minus its 24 within-ring edges (G')
has no Hamiltonian cycle: (1) forced edges from degree-2 squares, (2) a bipartite Hall-surplus
violation |N(X)| <= |X| for a nonempty proper subset X of one colour class (a Hamiltonian cycle
in a balanced bipartite graph forces |N(X)| >= |X| + 1), searched exactly with CP-SAT."""
from ortools.sat.python import cp_model

N = 8
F = "abcdefgh"
def ring(i, j):
    c = (N - 1) / 2
    return int(max(abs(i - c), abs(j - c)) - 0.5)
def nm(s): return F[s[0]] + str(s[1] + 1)
V = [(i, j) for i in range(N) for j in range(N)]
adj = {v: set() for v in V}
for (i, j) in V:
    for a, b in ((1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)):
        u = (i + a, j + b)
        if 0 <= u[0] < N and 0 <= u[1] < N and ring(i, j) != ring(*u):
            adj[(i, j)].add(u)
deg2 = [v for v in V if len(adj[v]) == 2]
print("G': edges", sum(len(a) for a in adj.values()) // 2, "; degree-2 squares:", sorted(nm(v) for v in deg2))
forced = {}
for v in deg2:
    for u in adj[v]:
        forced.setdefault(u, set()).add(v)
over = {nm(u): sorted(nm(w) for w in ws) for u, ws in forced.items() if len(ws) >= 2}
print("squares receiving two forced edges (saturated):", over)
for colour in (0, 1):
    W = [v for v in V if (v[0] + v[1]) % 2 == colour]
    m = cp_model.CpModel()
    x = {v: m.NewBoolVar(nm(v)) for v in W}
    nb = sorted({u for v in W for u in adj[v]})
    y = {u: m.NewBoolVar("n" + nm(u)) for u in nb}
    for v in W:
        for u in adj[v]:
            m.AddImplication(x[v], y[u])
    m.Add(sum(x.values()) >= 1)
    m.Add(sum(x.values()) <= len(W) - 1)
    m.Minimize(sum(y.values()) - sum(x.values()))
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 120
    st = s.Solve(m)
    X = [v for v in W if s.Value(x[v])]
    NX = [u for u in nb if s.Value(y[u])]
    print(f"colour {colour}: min |N(X)| - |X| = {int(s.ObjectiveValue())} ({s.StatusName(st)}); witness |X|={len(X)} |N(X)|={len(NX)}")
    if s.ObjectiveValue() <= 0:
        print("   X =", " ".join(sorted(nm(v) for v in X)))
        print("   N(X) =", " ".join(sorted(nm(u) for u in NX)))
