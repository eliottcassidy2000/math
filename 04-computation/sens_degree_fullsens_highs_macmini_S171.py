#!/usr/bin/env python3
"""Second engine for THM-2832/2835 (max point-sensitivity of low-degree
Boolean functions): HiGHS MIP feasibility of the Moebius-vanishing system.
Usage: python3 sens_degree_fullsens_highs_macmini_S171.py n d
Verdicts obtained at S171: (10,4) kInfeasible; (9,4) feasible; (7,3)
kInfeasible; LP relaxation of (10,4) feasible (integral obstruction)."""
import highspy, numpy as np, sys

def decide(n, d, integer=True):
    N = 1 << n
    pc = [bin(x).count("1") for x in range(N)]
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    lower = np.zeros(N); upper = np.ones(N)
    lower[0] = upper[0] = 0.0
    for i in range(n):
        lower[1 << i] = upper[1 << i] = 1.0
    h.addVars(N, lower, upper)
    if integer:
        h.changeColsIntegrality(N, np.arange(N, dtype=np.int32),
                                np.array([highspy.HighsVarType.kInteger] * N))
    starts, idxs, vals = [0], [], []
    rows = 0
    for S in range(N):
        if pc[S] <= d:
            continue
        T = S
        while True:
            idxs.append(T)
            vals.append(-1.0 if ((pc[S] - pc[T]) & 1) else 1.0)
            if T == 0:
                break
            T = (T - 1) & S
        starts.append(len(idxs))
        rows += 1
    h.addRows(rows, np.zeros(rows), np.zeros(rows), len(idxs),
              np.array(starts[:-1], dtype=np.int32),
              np.array(idxs, dtype=np.int32), np.array(vals))
    h.run()
    return h.getModelStatus()

if __name__ == "__main__":
    n, d = int(sys.argv[1]), int(sys.argv[2])
    print(f"HiGHS n={n} d={d}: {decide(n, d)}")
