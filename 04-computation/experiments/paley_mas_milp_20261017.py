#!/usr/bin/env python3
"""Exact MAS of a Paley tournament as a linear-ordering MIP (HiGHS via scipy.optimize.milp), S40.

Variables x_ij (i < j) = 1 if i precedes j; 3-dicycle inequalities 0 <= x_ij + x_jk - x_ik <= 1 for i < j < k; vertex 0
first (vertex-transitivity). Objective: number of forward arcs. Intended to test HYP-9287 beyond the subset-DP range
(p <= 31), after validation against the DP values at p = 23, 31 and the certified MAS(P_43) = 543; the validation
itself did not finish (below).
Usage: python paley_mas_milp_20261017.py P [TIME_LIMIT_SECONDS]
Result (S40, negative): HiGHS did not finish p = 23 within 600 s (the subset DP of paley_mas_exact_20261017.py takes
under a second), so this plain formulation cannot certify p = 47; a symmetry-aware branch and bound is needed.
Session opus-2026-10-17-S40.
"""
import sys, os, json, time, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import kronecker_table, class_number


def main():
    p = int(sys.argv[1])
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 3600.0
    chi = kronecker_table(-p)
    idx = {}
    for i in range(p):
        for j in range(i + 1, p):
            idx[(i, j)] = len(idx)
    nv = len(idx)
    # objective: forward arcs = sum_{i<j} [a_ij x_ij + a_ji (1 - x_ij)]; a_ij = 1 iff chi(j - i) = 1
    c = np.zeros(nv)
    const = 0
    for (i, j), k in idx.items():
        if chi[(j - i) % p] == 1:
            c[k] = 1.0
        else:
            c[k] = -1.0
            const += 1
    rows, cols, vals = [], [], []
    r = 0
    for i, j, k in itertools.combinations(range(p), 3):
        rows += [r, r, r]; cols += [idx[(i, j)], idx[(j, k)], idx[(i, k)]]; vals += [1, 1, -1]
        r += 1
    A = coo_matrix((vals, (rows, cols)), shape=(r, nv)).tocsr()
    lb = np.zeros(nv); ub = np.ones(nv)
    for j in range(1, p):
        lb[idx[(0, j)]] = 1.0            # vertex 0 first
    h = class_number(-p)
    natural = (p * (p - 1) // 2 + 2 * p * h // (6 if p == 3 else 2)) // 2
    t0 = time.time()
    res = milp(-c, constraints=[LinearConstraint(A, 0, 1)], integrality=np.ones(nv), bounds=Bounds(lb, ub),
               options=dict(time_limit=tl, disp=False, mip_rel_gap=0.0))
    dt = time.time() - t0
    out = dict(p=p, h=h, natural=natural, status=int(res.status), message=str(res.message), seconds=dt)
    if res.x is not None:
        val = int(round(-res.fun)) + const
        out['best_found'] = val
        bound = getattr(res, 'mip_dual_bound', None)
        if bound is not None:
            out['dual_bound'] = float(-bound) + const
    print(out, flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_mas_milp_p%d_20261017.json' % p), 'w'), indent=1)


if __name__ == '__main__':
    main()
