#!/usr/bin/env python3
"""Census of ALL ring-transition vectors (m01,m02,m11,m12,m13,m22,m23,m33) realised
by closed knight's tours of the 8x8 board, by CP-SAT feasibility of every
candidate allowed by the half-edge identities (chessboard-weave, 2026-10-06).

Candidates: the 4 half-edge equations leave 4 free parameters (m02,m11,m22,m33);
m01 = 8-m02, m12 = m33-m11-m22, m13 = 24-2m11-m01-m12, m23 = 56-2m33-m13, all
within [0, #edges of that pair].  Each candidate: AddCircuit model (as in
..._knight_cpsat8.py) + equality constraints on all eight pair counts.
Second part: independent cross-check of the extremal bounds with a DIFFERENT
formulation (undirected degree-2 edge variables + single-commodity flow
connectivity, no AddCircuit).
Run: python3 chessboard_weave_20261006_knight_cpsat8_vectors.py [workers] [time_limit]
"""
import sys, time
from collections import Counter
from itertools import product
from ortools.sat.python import cp_model

WORKERS = int(sys.argv[1]) if len(sys.argv) > 1 else 2
TL = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
N = 8
def ring(s):
    return int(max(abs(s[0] - 3.5), abs(s[1] - 3.5)) - 0.5)
SQ = [(i, j) for i in range(N) for j in range(N)]
RG = [ring(s) for s in SQ]
E = set()
for k, (i, j) in enumerate(SQ):
    for di, dj in [(1, 2), (2, 1), (-1, 2), (-2, 1)]:
        a, b = i + di, j + dj
        if 0 <= a < N and 0 <= b < N:
            E.add(tuple(sorted((k, a * N + b))))
E = sorted(E)
PAIRS = [(0, 1), (0, 2), (1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]
def ep(e):
    return tuple(sorted((RG[e[0]], RG[e[1]])))
NEDGE = Counter(ep(e) for e in E)

def circuit_model():
    m = cp_model.CpModel()
    x, arcs = {}, []
    for (u, v) in E:
        for a, b in ((u, v), (v, u)):
            x[(a, b)] = m.NewBoolVar("")
            arcs.append((a, b, x[(a, b)]))
    m.AddCircuit(arcs)
    m.Add(x[(SQ.index((1, 2)), 0)] == 1)
    y = {e: x[e] + x[(e[1], e[0])] for e in E}
    return m, x, y

def check_tour(succ):
    t = [0]
    while len(t) < 64:
        t.append(succ[t[-1]])
    assert succ[t[-1]] == 0 and sorted(t) == list(range(64))
    mv = [tuple(sorted((t[k], t[(k + 1) % 64]))) for k in range(64)]
    assert all(e in set(E) for e in mv)
    c = Counter(ep(e) for e in mv)
    return tuple(c[p] for p in PAIRS)

cands = []
for m02, m11, m22, m33 in product(range(9), range(9), range(9), range(9)):
    m01 = 8 - m02
    m12 = m33 - m11 - m22
    m13 = 24 - 2 * m11 - m01 - m12
    m23 = 56 - 2 * m33 - m13
    v = (m01, m02, m11, m12, m13, m22, m23, m33)
    if all(0 <= val <= NEDGE[p] for val, p in zip(v, PAIRS)):
        cands.append(v)
print(f"candidates allowed by the half-edge identities and edge supplies: {len(cands)}")
feas, infeas, unk = [], [], []
t0 = time.time()
for v in cands:
    m, x, y = circuit_model()
    for val, p in zip(v, PAIRS):
        m.Add(sum(y[e] for e in E if ep(e) == p) == val)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = TL
    s.parameters.num_search_workers = WORKERS
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        succ = {a: b for (a, b), l in x.items() if s.Value(l)}
        assert check_tour(succ) == v
        feas.append(v)
    elif st == cp_model.INFEASIBLE:
        infeas.append(v)
    else:
        unk.append(v)
print(f"feasible {len(feas)}, infeasible {len(infeas)}, unresolved {len(unk)}  ({time.time()-t0:.0f}s)")
print("vector order:", " ".join(f"m{a}{b}" for a, b in PAIRS))
print("FEASIBLE vectors:")
for v in feas:
    print("   ", v)
print("INFEASIBLE vectors (allowed by counting, no closed tour):")
for v in infeas:
    print("   ", v)
if unk:
    print("UNRESOLVED:", unk)
for k, p in enumerate(PAIRS):
    vals = sorted({v[k] for v in feas})
    print(f"  m{p[0]}{p[1]} realised values: {vals}")
W = sorted({v[2] + v[5] + v[7] for v in feas})
print("within-ring moves W realised:", W, " -> ring blocks 64-W:", sorted(64 - w for w in W))
print("|ring change|=2 moves realised:", sorted({v[1] + v[4] for v in feas}))

# ---- independent formulation cross-check of the extremal values
print("\nCross-check with flow formulation (no AddCircuit):")
def flow_model():
    m = cp_model.CpModel()
    y = {e: m.NewBoolVar("") for e in E}
    for v in range(64):
        m.Add(sum(y[e] for e in E if v in e) == 2)
    # single-commodity flow from root 0: sends 63 units, every other square absorbs 1
    f = {}
    for (u, v) in E:
        for a, b in ((u, v), (v, u)):
            f[(a, b)] = m.NewIntVar(0, 63, "")
            m.Add(f[(a, b)] <= 63 * y[(u, v)])
    for v in range(64):
        out = sum(f[(v, w)] for (a, w) in f if a == v)
        inn = sum(f[(w, v)] for (w, b) in f if b == v)
        m.Add(out - inn == (63 if v == 0 else -1))
    return m, y
def solve_flow(obj_edges, sense):
    m, y = flow_model()
    o = sum(y[e] for e in obj_edges)
    (m.Minimize if sense == "min" else m.Maximize)(o)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = 600
    s.parameters.num_search_workers = WORKERS
    st = s.Solve(m)
    val = round(s.ObjectiveValue()) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    if val is not None:
        # verify the solution is one Hamiltonian cycle
        used = [e for e in E if s.Value(y[e])]
        adj = {v: [] for v in range(64)}
        for a, b in used:
            adj[a].append(b); adj[b].append(a)
        seen, st_ = {0}, [0]
        while st_:
            a = st_.pop()
            for b in adj[a]:
                if b not in seen:
                    seen.add(b); st_.append(b)
        assert len(used) == 64 and len(seen) == 64 and all(len(adj[v]) == 2 for v in range(64))
    return s.StatusName(st), val, s.BestObjectiveBound()
checks = [("within (m11+m22+m33)", [e for e in E if RG[e[0]] == RG[e[1]]], "min"),
          ("within (m11+m22+m33)", [e for e in E if RG[e[0]] == RG[e[1]]], "max"),
          ("|ring change|=2", [e for e in E if abs(RG[e[0]] - RG[e[1]]) == 2], "min"),
          ("m13", [e for e in E if ep(e) == (1, 3)], "min"),
          ("m23", [e for e in E if ep(e) == (2, 3)], "min"),
          ("m33", [e for e in E if ep(e) == (3, 3)], "min")]
for nm, cls, sense in checks:
    t0 = time.time()
    st, val, bnd = solve_flow(cls, sense)
    print(f"  {nm:22s} {sense}: {st} value {val} bound {bnd:.2f} ({time.time()-t0:.1f}s)")
