#!/usr/bin/env python3
"""
sens_cpsat.py — CP-SAT decision of: exists Boolean f on {0,1}^n, real degree <= d,
f(0)=0, f(e_i)=1 for all i.  (Same question as sens_search.c, second engine.)

Encoding: Boolean vars f[x]; for every S with |S| > d the Moebius coefficient
vanishes:  sum_{T subseteq S} (-1)^{|S|-|T|} f[T] == 0.
Symmetry breaking: coordinates sorted by level-2 degree sequence
(deg_i = sum_j f[e_i+e_j] non-increasing), plus lexicographic tie-break on
level-2 rows.  Valid: some coordinate permutation achieves this normal form.
"""
import sys
from ortools.sat.python import cp_model

def decide(n, d, workers=8, log=True):
    N = 1 << n
    m = cp_model.CpModel()
    f = [m.NewBoolVar(f"f{x}") for x in range(N)]
    m.Add(f[0] == 0)
    for i in range(n):
        m.Add(f[1 << i] == 1)
    pc = [bin(x).count("1") for x in range(N)]
    for S in range(N):
        if pc[S] <= d:
            continue
        terms, coefs = [], []
        T = S
        while True:
            sign = -1 if ((pc[S] - pc[T]) & 1) else 1
            terms.append(f[T]); coefs.append(sign)
            if T == 0:
                break
            T = (T - 1) & S
        m.Add(cp_model.LinearExpr.WeightedSum(terms, coefs) == 0)
    # symmetry breaking on pair-degrees
    degs = []
    for i in range(n):
        di = m.NewIntVar(0, n - 1, f"d{i}")
        m.Add(di == sum(f[(1 << i) | (1 << j)] for j in range(n) if j != i))
        degs.append(di)
    for i in range(n - 1):
        m.Add(degs[i] >= degs[i + 1])
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_memory_in_mb = 2800
    solver.parameters.log_search_progress = log
    status = solver.Solve(m)
    name = solver.StatusName(status)
    print(f"CPSAT n={n} d={d}: {name}")
    if name == "FEASIBLE" or name == "OPTIMAL":
        tt = [x for x in range(N) if solver.Value(f[x])]
        print("truth table:", tt)
        # Moebius coefficients
        vals = [solver.Value(f[x]) for x in range(N)]
        poly = []
        for S in range(N):
            c = 0
            T = S
            while True:
                c += (-1 if ((pc[S] - pc[T]) & 1) else 1) * vals[T]
                if T == 0:
                    break
                T = (T - 1) & S
            if c:
                mono = "*".join(f"x{i+1}" for i in range(n) if S & (1 << i)) or "1"
                poly.append(f"{'+' if c > 0 else '-'} {abs(c) if abs(c)!=1 or S==0 else ''}{mono if S else ''}")
        print("poly:", " ".join(poly))
    return name

if __name__ == "__main__":
    n, d = int(sys.argv[1]), int(sys.argv[2])
    w = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    decide(n, d, workers=w)
