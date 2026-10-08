#!/usr/bin/env python3
"""Beyond translation-only: for a map whose multiplier residues generate G = <m_i mod d> in (Z/d)^*, the possible
couplings are the affine permutations j -> a j + b (a in G, b in Z/d).  One-step Lamperti transience holds if a single
form Q balances the covariance of EVERY non-identity coupling (any of them may occur, in any order).
Covariance of coupling pi: C_pi = (1/d) sum_j (v_pi(j) - v_j)(v_pi(j) - v_j)^T (= transported Laplacian of the cycle
graph of pi).  Test independent rank-3 and rank-4 configurations on Z_7 and Z_11 for every subgroup G."""
import itertools, random, math
from fractions import Fraction as Fr
from balanced_lamperti import matinv, posdef_exact, jacobi_eigs
from balanced_whitened import cholesky, inv_lower, mm, tr_, margin_f, search
from balanced_types import exact_verify, orbit_rep
def subgroups(d):
    units = [a for a in range(1, d)]
    subs = set()
    for g in units:
        H = {1}; x = g
        while x != 1: H.add(x); x = x * g % d
        subs.add(frozenset(H))
    return sorted(subs, key=len)
def config(d, Z):
    rho = d - len(Z); idx = 0; v = []
    for j in range(d):
        if j in Z: v.append([Fr(0)] * rho)
        else:
            e = [Fr(0)] * rho; e[idx] = Fr(1); idx += 1; v.append(e)
    return rho, v
def cov(v, pi, d, rho):
    C = [[Fr(0)] * rho for _ in range(rho)]
    for j in range(d):
        w = [x - y for x, y in zip(v[pi[j]], v[j])]
        for p in range(rho):
            for q in range(rho): C[p][q] += w[p] * w[q] / d
    return C
rnd = random.Random(9)
for d in (7, 11):
    subs = subgroups(d)
    reps = sorted({orbit_rep(d, Z) for zs in (d - 3, d - 4) for Z in itertools.combinations(range(d), zs)}, key=lambda z: (-len(z), z))
    print(f"d = {d}: subgroups of (Z/{d})^* by order: {[sorted(H) for H in subs]}")
    for Z in reps:
        rho, v = config(d, set(Z))
        line = []
        for H in subs:
            covs = []
            for a in sorted(H):
                for b in range(d):
                    if a == 1 and b == 0: continue
                    pi = [(a * j + b) % d for j in range(d)]
                    C = cov(v, pi, d, rho)
                    if any(any(x != 0 for x in row) for row in C): covs.append(C)
            ranks = sorted({sum(1 for ev in jacobi_eigs([[float(x) for x in row] for row in C]) if ev > 1e-9) for C in covs})
            S = [[float(sum(C[p][q] for C in covs)) for q in range(rho)] for p in range(rho)]
            L = cholesky(S); Li = inv_lower(L)
            Csw = [mm(mm(Li, [[float(x) for x in row] for row in C]), tr_(Li)) for C in covs]
            Qw, m = search(Csw, rho, rnd, iters=1500)
            ver = exact_verify(mm(mm(L, Qw), tr_(L)), covs, rho) if m > 0 else None
            line.append(f"|G|={len(H)}: {'BAL' if ver else ('no' if m <= 0 else '?')} (min cov rank {ranks[0]})")
        print(f"   units at {Z} (rank {rho}): " + "; ".join(line), flush=True)
