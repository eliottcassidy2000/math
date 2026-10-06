#!/usr/bin/env python3
"""Closed knight's tours on 8x8 vs the Chebyshev rings: exact min/max of
ring-pair move counts via CP-SAT AddCircuit (chessboard-weave, 2026-10-06).

Model: one Boolean per directed knight arc (336 arcs), AddCircuit over the 64
squares (no self-loops => Hamiltonian cycle). Orientation symmetry broken by
fixing the arc c2 -> a1 (a1 has degree 2, so every closed tour passes
c2-a1-b3; fixing the direction just picks one of the two orientations).
Undirected edge usage y_e = x_uv + x_vu.  Objectives are linear in y.
Every reported optimum is re-verified on the witness tour by independent code.

Run: python3 chessboard_weave_20261006_knight_cpsat8.py [time_limit_s] [workers]
"""
import sys, time
from collections import Counter
from itertools import product
from ortools.sat.python import cp_model

N = 8
TL = float(sys.argv[1]) if len(sys.argv) > 1 else 600.0
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
ONLY = sys.argv[3].split(",") if len(sys.argv) > 3 else None

def ring(s):
    c = (N - 1) / 2
    return int(max(abs(s[0] - c), abs(s[1] - c)) - 0.5)

SQ = [(i, j) for i in range(N) for j in range(N)]
IDX = {s: k for k, s in enumerate(SQ)}
E = set()
for (i, j) in SQ:
    for di, dj in [(1, 2), (2, 1), (-1, 2), (-2, 1)]:
        a, b = i + di, j + dj
        if 0 <= a < N and 0 <= b < N:
            E.add(tuple(sorted((IDX[(i, j)], IDX[(a, b)]))))
E = sorted(E)
assert len(E) == 168
RG = [ring(s) for s in SQ]
def epair(e):
    return tuple(sorted((RG[e[0]], RG[e[1]])))
def alg(k):
    s = SQ[k]
    return "abcdefgh"[s[1]] + str(s[0] + 1)

PAIRS = sorted(set(epair(e) for e in E))

def build():
    m = cp_model.CpModel()
    x = {}
    arcs = []
    for (u, v) in E:
        x[(u, v)] = m.NewBoolVar(f"x{u}_{v}")
        x[(v, u)] = m.NewBoolVar(f"x{v}_{u}")
        arcs.append((u, v, x[(u, v)]))
        arcs.append((v, u, x[(v, u)]))
    m.AddCircuit(arcs)
    a1, c2 = IDX[(0, 0)], IDX[(1, 2)]
    m.Add(x[(c2, a1)] == 1)
    y = {e: x[e] + x[(e[1], e[0])] for e in E}
    return m, x, y

def tour_from(solver, x):
    succ = {}
    for (u, v), var in x.items():
        if solver.Value(var):
            succ[u] = v
    t = [0]
    while len(t) < 64:
        t.append(succ[t[-1]])
    assert succ[t[-1]] == 0
    return t

def verify(t):
    """independent check: Hamiltonian cycle in the knight graph; returns stats."""
    assert sorted(t) == list(range(64))
    Es = set(E)
    moves = [tuple(sorted((t[k], t[(k + 1) % 64]))) for k in range(64)]
    assert all(mv in Es for mv in moves)
    cnt = Counter(epair(mv) for mv in moves)
    within = sum(c for (r, s), c in cnt.items() if r == s)
    dr2 = sum(c for (r, s), c in cnt.items() if s - r == 2)
    blocks = sum(1 for k in range(64) if RG[t[k]] != RG[t[(k + 1) % 64]])
    return cnt, within, dr2, blocks

objectives = [("within", [e for e in E if RG[e[0]] == RG[e[1]]]),
              ("ringchange2", [e for e in E if abs(RG[e[0]] - RG[e[1]]) == 2])]
for p in PAIRS:
    objectives.append((f"pair{p[0]}{p[1]}", [e for e in E if epair(e) == p]))

print(f"8x8 knight graph: {len(E)} edges; ring pairs with edges: {PAIRS}")
print(f"CP-SAT {cp_model.__name__}, time limit {TL}s per solve, {WORKERS} workers\n")
results = {}
for name, cls in objectives:
    if ONLY and name not in ONLY:
        continue
    for sense in ("min", "max"):
        m, x, y = build()
        obj = sum(y[e] for e in cls)
        (m.Minimize if sense == "min" else m.Maximize)(obj)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = TL
        solver.parameters.num_search_workers = WORKERS
        solver.parameters.random_seed = 1
        t0 = time.time()
        st = solver.Solve(m)
        dt = time.time() - t0
        sname = solver.StatusName(st)
        if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            t = tour_from(solver, x)
            cnt, within, dr2, blocks = verify(t)
            val = round(solver.ObjectiveValue())
            # recompute objective independently from the tour
            ind = sum(1 for k in range(64) if tuple(sorted((t[k], t[(k + 1) % 64]))) in set(cls))
            assert ind == val, (ind, val)
            bnd = solver.BestObjectiveBound()
            results[(name, sense)] = (sname, val, bnd)
            print(f"[{name:12s} {sense}] status={sname:8s} value={val:3d} bound={bnd:7.2f} time={dt:7.1f}s"
                  f"  witness: within={within} ringchange2={dr2} blocks={blocks} pairs={dict(sorted(cnt.items()))}")
            print("    tour:", " ".join(alg(k) for k in t))
        else:
            results[(name, sense)] = (sname, None, solver.BestObjectiveBound())
            print(f"[{name:12s} {sense}] status={sname} time={dt:.1f}s bound={solver.BestObjectiveBound()}")
        sys.stdout.flush()

print("\nSUMMARY (value, status, bound):")
for name, _ in objectives:
    if ONLY and name not in ONLY:
        continue
    lo = results.get((name, "min")); hi = results.get((name, "max"))
    print(f"  {name:12s} min={lo[1]} ({lo[0]}, bound {lo[2]:.2f})   max={hi[1]} ({hi[0]}, bound {hi[2]:.2f})")
