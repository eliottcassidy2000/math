#!/usr/bin/env python3
"""Knight's tours on the 8x8 with inner rings removed (chessboard-weave, 2026-10-06).

 (a) frame = rings 2,3 (2-wide border, 48 squares): closed tour?  count (exact DFS)
 (b) rings 1,2,3 (8x8 minus the central 2x2, 60 squares): closed tour? (CP-SAT)
 (c) central 4x4 = rings 0,1: Hamiltonian path? (exhaustive)
 plus single-ring / two-ring variants.  Positive controls for the DFS counter:
 the 6x6 and 5x8 sub-boards embedded in 8x8 (9862 / 44202 cycles, 6637920 paths).
Run: python3 chessboard_weave_20261006_knight_holes.py
"""
import os, subprocess, tempfile, time
from itertools import product
from ortools.sat.python import cp_model

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="knightholes_")
EXE = os.path.join(TMP, "subham")
subprocess.run(["cc", "-O2", "-o", EXE, os.path.join(HERE, "chessboard_weave_20261006_knight_subham.c")], check=True)

def ring(i, j):
    return int(max(abs(i - 3.5), abs(j - 3.5)) - 0.5)
def mask(pred):
    return sum(1 << (8 * i + j) for i, j in product(range(8), repeat=2) if pred(i, j))
def alg(k):
    return "abcdefgh"[k % 8] + str(k // 8 + 1)
def run(m, mode):
    return subprocess.run([EXE, f"{m:016x}", mode], capture_output=True, text=True, check=True).stdout.strip()

print("positive controls (sub-boards of the 8x8):")
print(" ", run(mask(lambda i, j: i < 6 and j < 6), "c"), " [6x6 known 9862]")
print(" ", run(mask(lambda i, j: i < 5), "c"), " [5x8 known 44202]")
print(" ", run(mask(lambda i, j: i < 6 and j < 6), "p"), " [6x6 known 6637920]")
print(" ", run(mask(lambda i, j: i < 5 and j < 6), "p"), " [5x6 known 37568]")

def cpsat_cycle(m, tl=600):
    V = [k for k in range(64) if m >> k & 1]
    arcs, x = [], {}
    mdl = cp_model.CpModel()
    for k in V:
        i, j = divmod(k, 8)
        for di, dj in [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]:
            a, b = i + di, j + dj
            if 0 <= a < 8 and 0 <= b < 8 and (m >> (8 * a + b) & 1):
                x[(k, 8 * a + b)] = mdl.NewBoolVar("")
                arcs.append((k, 8 * a + b, x[(k, 8 * a + b)]))
    # AddCircuit needs node indices 0..n-1: relabel
    lab = {k: t for t, k in enumerate(V)}
    mdl.AddCircuit([(lab[u], lab[v], l) for u, v, l in arcs])
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl
    s.parameters.num_search_workers = 2
    t0 = time.time()
    st = s.Solve(mdl)
    res = s.StatusName(st)
    tour = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        succ = {u: v for (u, v), l in x.items() if s.Value(l)}
        tour = [V[0]]
        while len(tour) < len(V):
            tour.append(succ[tour[-1]])
        assert succ[tour[-1]] == V[0] and len(set(tour)) == len(V)
        for a, b in zip(tour, tour[1:] + tour[:1]):
            da, db = abs(a // 8 - b // 8), abs(a % 8 - b % 8)
            assert {da, db} == {1, 2}
    return res, tour, time.time() - t0

cases = [("(a) frame = rings 2,3 (48 sq)", lambda i, j: ring(i, j) >= 2, True),
         ("(b) rings 1,2,3 = 8x8 minus centre 2x2 (60 sq)", lambda i, j: ring(i, j) >= 1, False),
         ("(c) central 4x4 = rings 0,1 (16 sq)", lambda i, j: ring(i, j) <= 1, True),
         ("    rings 1,2 (32 sq)", lambda i, j: ring(i, j) in (1, 2), True),
         ("    rings 0,1,2 = central 6x6 (36 sq)", lambda i, j: ring(i, j) <= 2, False),
         ("    ring 3 alone (28 sq)", lambda i, j: ring(i, j) == 3, True),
         ("    rings 0,2,3 (52 sq)", lambda i, j: ring(i, j) != 1, True),
         ("    rings 0,1,3 (44 sq)", lambda i, j: ring(i, j) != 2, True)]
for name, pred, count in cases:
    m = mask(pred)
    res, tour, dt = cpsat_cycle(m)
    print(f"\n{name}: CP-SAT Hamiltonian cycle -> {res} ({dt:.2f}s)")
    if tour:
        print("   witness tour:", " ".join(alg(k) for k in tour))
    if count:
        t0 = time.time()
        print("   exact DFS:", run(m, "c"), f"({time.time()-t0:.1f}s)")
        t0 = time.time()
        print("   exact DFS:", run(m, "p"), f"({time.time()-t0:.1f}s)")
