#!/usr/bin/env python3
"""The weave law (chessboard-weave session 2026-10-06).
Theorem: in every closed knight's tour of the 8x8 board, (#moves inside ring 3) = (#moves inside rings 1-2) =: k >= 1.
Proof is counting (|rings 0,3| = |rings 1,2| = 32) plus the colour lock; this script checks it on
random closed tours (CP-SAT with random edge weights) and checks the 6x6 analogue
(#moves inside ring 2) = 4 + (#moves inside rings 0-1) on random 6x6 tours."""
import random
from ortools.sat.python import cp_model

def ring(i, j, n):
    c = (n - 1) / 2
    return int(max(abs(i - c), abs(j - c)) - 0.5)

def knight_arcs(n):
    V = [(i, j) for i in range(n) for j in range(n)]
    idx = {v: k for k, v in enumerate(V)}
    arcs = []
    for (i, j) in V:
        for a, b in ((1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)):
            u = (i + a, j + b)
            if 0 <= u[0] < n and 0 <= u[1] < n:
                arcs.append((idx[(i, j)], idx[u]))
    return V, arcs

def random_tour(n, seed):
    V, arcs = knight_arcs(n)
    rnd = random.Random(seed)
    m = cp_model.CpModel()
    lit = {}
    circ = []
    for (u, v) in arcs:
        lit[(u, v)] = m.NewBoolVar("")
        circ.append((u, v, lit[(u, v)]))
    m.AddCircuit(circ)
    m.Maximize(sum(rnd.randint(0, 1000) * lit[a] for a in arcs))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 2
    s.parameters.max_time_in_seconds = 5
    s.parameters.random_seed = seed
    st = s.Solve(m)
    assert st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    succ = {u: v for (u, v) in arcs if s.Value(lit[(u, v)])}
    tour, x = [0], succ[0]
    while x != 0:
        tour.append(x); x = succ[x]
    assert len(tour) == n * n
    return [V[k] for k in tour]

def stays(tour, n, ringset):
    c = 0
    for k in range(len(tour)):
        a, b = tour[k], tour[(k + 1) % len(tour)]
        if ring(*a, n) in ringset and ring(*b, n) in ringset:
            c += 1
    return c

ok8, ks = True, []
for seed in range(60):
    t = random_tour(8, seed)
    k_out, k_mid = stays(t, 8, {3}), stays(t, 8, {1, 2})
    ok8 &= (k_out == k_mid and k_out >= 1)
    ks.append(k_out)
print("8x8: 60 random closed tours; law (#inside ring 3) == (#inside rings 1-2) >= 1 holds:", ok8,
      "; k values seen:", sorted(set(ks)))
ok6, ds = True, []
for seed in range(60):
    t = random_tour(6, seed)
    o, i = stays(t, 6, {2}), stays(t, 6, {0, 1})
    ok6 &= (o == i + 4)
    ds.append(o)
print("6x6: 60 random closed tours; law (#inside ring 2) == 4 + (#inside rings 0-1) holds:", ok6,
      "; outer-ring stays seen:", sorted(set(ds)))
