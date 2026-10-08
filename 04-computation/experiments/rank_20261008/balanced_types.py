#!/usr/bin/env python3
"""Balanced-configuration types for translation-only maps with INDEPENDENT non-unit multipliers.
Configuration: v_j = 0 for j in Z (unit multipliers), v_j = e_(rank index) otherwise (standard basis of R^rho).
Balanced: exists Q > 0 with (1/2) tr(C_b Q^-1) Q - C_b > 0 for every lag b != 0.  The property is invariant under
GL(rho) and under the affine group AGL(1, d) acting on positions (which permutes the lags), so it depends only on
the AGL(1, d)-orbit of the zero set Z.  For each orbit (|Z| = d - rho, rho >= 3) we search Q (whitened local search)
and verify exactly.  Contraction needs prod(m_i) < d^d with m_i = 1 mod d independent: for d = 5 rho <= 3,
d = 7 rho <= 4, d = 11 rho <= 6 or so (the smallest independent m = 1 mod d multiply too fast)."""
import itertools, random, math, sys
from fractions import Fraction as Fr
from balanced_lamperti import matinv, posdef_exact, jacobi_eigs
from balanced_whitened import cholesky, inv_lower, mm, tr_, margin_f, search
def config(d, Z):
    rho = d - len(Z); idx = 0; v = []
    for j in range(d):
        if j in Z: v.append([Fr(0)] * rho)
        else:
            e = [Fr(0)] * rho; e[idx] = Fr(1); idx += 1; v.append(e)
    return rho, v
def lagcovs(v, d, rho):
    out = []
    for b in range(1, d):
        C = [[Fr(0)] * rho for _ in range(rho)]
        for j in range(d):
            w = [x - y for x, y in zip(v[(j + b) % d], v[j])]
            for p in range(rho):
                for q in range(rho): C[p][q] += w[p] * w[q] / d
        out.append(C)
    return out
def orbit_rep(d, Z):
    best = None
    for a in range(1, d):
        for c in range(d):
            img = tuple(sorted((a * z + c) % d for z in Z))
            if best is None or img < best: best = img
    return best
def exact_verify(Q, Cs, rho):
    for den in (1, 10, 100, 1000, 10000):
        Qr = [[Fr(round(Q[i][j] * den), den) for j in range(rho)] for i in range(rho)]
        Qr = [[(Qr[i][j] + Qr[j][i]) / 2 for j in range(rho)] for i in range(rho)]
        if not posdef_exact(Qr): continue
        Qi = matinv(Qr); ok = True
        for C in Cs:
            t = sum(C[a][c] * Qi[c][a] for a in range(rho) for c in range(rho))
            S = [[t / 2 * Qr[a][c] - C[a][c] for c in range(rho)] for a in range(rho)]
            if not posdef_exact(S): ok = False; break
        if ok: return den, Qr
    return None
rnd = random.Random(5)
for d in (5, 7, 11):
    reps = set()
    for zsize in range(1, d - 2):   # zsize 0 is the zsize-1 type shifted (only differences v_(j+b) - v_j matter)
        for Z in itertools.combinations(range(d), zsize):
            reps.add(orbit_rep(d, Z))
    print(f"d = {d}: {len(reps)} AGL(1,{d})-orbits of unit-position sets with rank >= 3", flush=True)
    for Z in sorted(reps, key=lambda z: (-len(z), z)):
        rho, v = config(d, set(Z))
        if rho > 7: continue
        Cs = lagcovs(v, d, rho)
        S = [[float(sum(C[a][c] for C in Cs)) for c in range(rho)] for a in range(rho)]
        L = cholesky(S); Li = inv_lower(L)
        Csw = [mm(mm(Li, [[float(x) for x in row] for row in C]), tr_(Li)) for C in Cs]
        bestm = -1e9; bestQ = None
        for trial in range(3):
            Qw, m = search(Csw, rho, rnd, iters=2500)
            if m > bestm: bestm, bestQ = m, Qw
        ver = None
        if bestm > 0:
            Q = mm(mm(L, bestQ), tr_(L)); ver = exact_verify(Q, Cs, rho)
        lagspec = sorted({tuple(sorted(round(x * d, 3) for x in jacobi_eigs([[float(x) for x in row] for row in C]))) for C in Cs})
        print(f"   units at {Z} (rank {rho}): whitened margin {bestm:+.4f}  exact: {'BALANCED (den %d)' % ver[0] if ver else '-'}  "
              f"distinct lag spectra x{d}: {len(lagspec)}", flush=True)
