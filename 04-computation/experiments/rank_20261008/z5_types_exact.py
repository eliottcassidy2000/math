#!/usr/bin/env python3
"""Exact balanced forms for the two rank-3 types of contracting translation-only maps on Z_5.
Type A: exponent vectors (0, 0, a, b, c) (unit multipliers adjacent);  type B: (0, a, 0, b, c) (distance 2).
In the basis a, b, c (any three independent vectors are GL(3)-equivalent, and the balanced property is
GL-invariant), the lag covariances are C_b = (1/5) sum_j (v_(j+b) - v_j)(v_(j+b) - v_j)^T, b = 1..4.
We look for a rational positive definite Q with  S_b := (1/2) tr(C_b Q^-1) Q - C_b  positive definite for all b
(Lamperti: V(x) = (x^T Q^-1 x)^(-alpha/2) is then a supermartingale far out, for 0 < alpha < alpha_max).
Also: alpha_max = min_b tr(C_b Q^-1)/lambda_max(C_b Q^-1) - 2, and lag-covariance spectra in the standard form."""
import itertools, random, math
from fractions import Fraction as Fr
from balanced_lamperti import matinv, posdef_exact, jacobi_eigs
def lagcovs(v):
    d = len(v); out = []
    for b in range(1, d):
        C = [[Fr(0)] * 3 for _ in range(3)]
        for j in range(d):
            w = [x - y for x, y in zip(v[(j + b) % d], v[j])]
            for p in range(3):
                for q in range(3): C[p][q] += w[p] * w[q] / d
        out.append(C)
    return out
E = [[Fr(1), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0)], [Fr(0), Fr(0), Fr(1)]]; Z = [Fr(0)] * 3
types = {'A (units adjacent)': [Z, Z, E[0], E[1], E[2]], 'B (units at distance 2)': [Z, E[0], Z, E[1], E[2]]}
def verify(Q, Cs):
    if not posdef_exact(Q): return None
    Qi = matinv(Q); amax = 1e9
    for C in Cs:
        t = sum(C[a][c] * Qi[c][a] for a in range(3) for c in range(3))
        S = [[t / 2 * Q[a][c] - C[a][c] for c in range(3)] for a in range(3)]
        if not posdef_exact(S): return None
        # lambda_max of C Q^-1 (float) for alpha_max
        CQ = [[float(sum(C[a][k] * Qi[k][c] for k in range(3))) for c in range(3)] for a in range(3)]
        # symmetrize via Q^-1/2 is awkward; use power iteration on CQ (real positive spectrum)
        x = [1.0, 0.7, 0.3]
        for _ in range(500):
            y = [sum(CQ[a][c] * x[c] for c in range(3)) for a in range(3)]
            nrm = max(abs(z) for z in y); x = [z / nrm for z in y]
        lam = nrm; amax = min(amax, float(t) / lam - 2)
    return amax
rnd = random.Random(1)
for name, v in types.items():
    Cs = lagcovs(v)
    spec = [sorted(round(x, 4) for x in jacobi_eigs([[float(x) for x in row] for row in C])) for C in Cs]
    print(f"type {name}: lag covariance spectra (standard form, x5): {[[round(5*s, 4) for s in sp] for sp in spec]}")
    best = None
    # small-integer search for Q (entries in -3..6, symmetric), then report the one with the largest alpha_max
    cands = []
    for _ in range(60000):
        a, b_, c = rnd.randint(1, 8), rnd.randint(1, 8), rnd.randint(1, 8)
        x, y, z = rnd.randint(-4, 4), rnd.randint(-4, 4), rnd.randint(-4, 4)
        Q = [[Fr(a), Fr(x), Fr(y)], [Fr(x), Fr(b_), Fr(z)], [Fr(y), Fr(z), Fr(c)]]
        am = verify(Q, Cs)
        if am is not None and (best is None or am > best[0]): best = (am, Q)
    if best:
        am, Q = best
        print(f"   exact balanced integer form found: Q = {[[int(x) for x in row] for row in Q]};  alpha_max = {am:.4f}")
    else:
        print("   no integer form found in the box")
