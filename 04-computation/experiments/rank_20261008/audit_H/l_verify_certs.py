#!/usr/bin/env python3
"""audit H: verify dual_certificates.json exactly.  For each (map, k): recompute the set of k-step blocks with hdp.Tables
(independently of the LP code), check every certificate block is one of them, and check that
N = sum c_l [A_l - 2 (A_l z_l)(A_l z_l)^T/(z_l^T A_l z_l)] is negative definite (Fractions; principal-minor sums of -N).
A negative definite N proves: no positive definite Q balances all k-step blocks of that map (no fixed-length-k certificate)."""
import json, itertools
from fractions import Fraction as Fr
import numpy as np
from hdp import Tables

def det(M):
    n = len(M)
    if n == 1: return M[0][0]
    return sum((-1) ** c * M[0][c] * det([r[:c] + r[c + 1:] for r in M[1:]]) for c in range(n))

for cert in json.load(open('dual_certificates.json')):
    m, k = cert['m'], cert['k']
    tb = Tables(5, m)
    T = tb.table(k).reshape(-1, tb.rho, tb.rho)
    S = {tuple(x.ravel()) for x in T}
    if cert.get('rank_obstruction'):
        rk = min(np.linalg.matrix_rank(np.array(x, dtype=float).reshape(3, 3)) for x in S if any(x))
        print(f"m = {m}, k = {k}: rank obstruction, min block rank {rk}"); continue
    N = [[Fr(0)] * 3 for _ in range(3)]
    member = True
    for c in cert['cuts']:
        A = c['A']; z = [Fr(x) for x in c['z']]; w = Fr(c['c'])
        member &= tuple(np.array(A).ravel()) in S
        assert w >= 0
        Az = [sum(A[p][q] * z[q] for q in range(3)) for p in range(3)]
        zAz = sum(z[p] * Az[p] for p in range(3)); assert zAz > 0
        for p in range(3):
            for q in range(3):
                N[p][q] += w * (A[p][q] - 2 * Az[p] * Az[q] / zAz)
    negN = [[-x for x in r] for r in N]
    nd = all(sum(det([[negN[a][b] for b in idx] for a in idx]) for idx in itertools.combinations(range(3), s)) > 0 for s in (1, 2, 3))
    print(f"m = {m}, k = {k}: {len(cert['cuts'])} cuts, all blocks genuine: {member}, sum negative definite: {nd} "
          f"-> no fixed-length-{k} certificate: {member and nd}")
