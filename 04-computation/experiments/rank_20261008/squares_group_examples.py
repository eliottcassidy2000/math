#!/usr/bin/env python3
"""Non-translation examples of one-step Lamperti transience: coupling group G = squares mod d (odd order, -1 not in G
for d = 3 mod 4).  Configurations: units at Z, independent vectors elsewhere.  Find and EXACTLY verify one form Q that
balances the covariance of every non-identity coupling j -> a j + b, a in G, b in Z/d.  Print Q.
Contracting example on Z_7: multipliers (1, 1, 1, 2, 11, 23, 29): residues 2, 4, 2, 1 lie in {1, 2, 4}; 2, 11, 23, 29
independent; product 14674 < 7^7; rank 4."""
import random, math
from fractions import Fraction as Fr
from balanced_lamperti import matinv, posdef_exact, jacobi_eigs
from balanced_whitened import cholesky, inv_lower, mm, tr_, search
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
def verify(Qr, covs, rho):
    if not posdef_exact(Qr): return False
    Qi = matinv(Qr)
    for C in covs:
        t = sum(C[a][c] * Qi[c][a] for a in range(rho) for c in range(rho))
        S = [[t / 2 * Qr[a][c] - C[a][c] for c in range(rho)] for a in range(rho)]
        if not posdef_exact(S): return False
    return True
rnd = random.Random(12)
for d, Z in ((7, {0, 1, 2}), (11, {0, 1, 2, 3, 4, 5, 6}), (11, {0, 1, 2, 3, 4, 5, 7})):
    G = sorted({(x * x) % d for x in range(1, d)})
    rho, v = config(d, Z)
    covs = [cov(v, [(a * j + b) % d for j in range(d)], d, rho) for a in G for b in range(d) if not (a == 1 and b == 0)]
    S = [[float(sum(C[p][q] for C in covs)) for q in range(rho)] for p in range(rho)]
    L = cholesky(S); Li = inv_lower(L)
    Csw = [mm(mm(Li, [[float(x) for x in row] for row in C]), tr_(Li)) for C in covs]
    best = None
    for trial in range(6):
        Qw, m = search(Csw, rho, rnd, iters=3000)
        if best is None or m > best[1]: best = (Qw, m)
    Q = mm(mm(L, best[0]), tr_(L)); found = None
    for den in (1, 2, 4, 10, 20, 100, 1000):
        Qr = [[Fr(round(Q[i][j] * den), den) for j in range(rho)] for i in range(rho)]
        Qr = [[(Qr[i][j] + Qr[j][i]) / 2 for j in range(rho)] for i in range(rho)]
        if verify(Qr, covs, rho): found = Qr; break
    print(f"d = {d}, G = squares {G}, units at {sorted(Z)} (rank {rho}): {len(covs)} couplings; whitened margin {best[1]:+.4f}; "
          f"exact balanced form: {[[str(x) for x in row] for row in found] if found else 'NOT FOUND'}", flush=True)
m = [1, 1, 1, 2, 11, 23, 29]
print("Z_7 example m =", m, "residues", [x % 7 for x in m], "Lambda =", round(sum(math.log(x / 7) for x in m) / 7, 4), "r_i =", [(-m[i] * i) % 7 for i in range(7)])
