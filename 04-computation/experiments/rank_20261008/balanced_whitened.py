#!/usr/bin/env python3
"""Basis-independent balanced-form search.  For a translation-only map on Z_d with multiplier exponent vectors v_j,
whiten by S = sum_b C_b (= 2 * scatter of the v_j): in coordinates with S = I search Q near I for the balanced
condition  (1/2) tr(C_b Q^-1) Q - C_b > 0  for all lags b != 0; map back, rationalize, verify EXACTLY (Sylvester).
The balanced property depends only on the configuration up to GL(rho), so equivalent maps get the same answer."""
import itertools, math, random, sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
from balanced_lamperti import rank_and_coords, lag_cov, matinv, posdef_exact, jacobi_eigs, det
def cholesky(A):
    n = len(A); L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = A[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(s) if i == j else s / L[j][j]
    return L
def inv_lower(L):
    n = len(L); X = [[0.0] * n for _ in range(n)]
    for i in range(n):
        X[i][i] = 1 / L[i][i]
        for j in range(i):
            X[i][j] = -sum(L[i][k] * X[k][j] for k in range(j, i)) / L[i][i]
    return X
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def tr_(A): return [list(r) for r in zip(*A)]
def margin_f(Q, Cs):
    n = len(Q)
    # Q^-1 via float Gauss
    A = [Q[i][:] + [float(i == j) for j in range(n)] for i in range(n)]
    for c in range(n):
        pr = max(range(c, n), key=lambda r: abs(A[r][c])); A[c], A[pr] = A[pr], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c: f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    Qi = [row[n:] for row in A]
    worst = 1e9
    for C in Cs:
        t = sum(C[a][c] * Qi[c][a] for a in range(n) for c in range(n))
        S = [[t / 2 * Q[a][c] - C[a][c] for c in range(n)] for a in range(n)]
        worst = min(worst, min(jacobi_eigs(S)) / t)
    return worst
def search(Csf, rho, rnd, iters=3000):
    Q = [[float(i == j) for j in range(rho)] for i in range(rho)]
    best, bm = Q, margin_f(Q, Csf); step = 0.3
    for it in range(iters):
        X = [[rnd.gauss(0, step) for _ in range(rho)] for _ in range(rho)]
        X = [[(X[i][j] + X[j][i]) / 2 for j in range(rho)] for i in range(rho)]
        Qn = [[best[i][j] + X[i][j] for j in range(rho)] for i in range(rho)]
        if min(jacobi_eigs(Qn)) <= 1e-3: continue
        m = margin_f(Qn, Csf)
        if m > bm: best, bm = Qn, m
        if it % 500 == 499: step *= 0.6
    return best, bm
def analyse(ms, d, rnd):
    rho, coords = rank_and_coords(list(ms))
    if rho < 3: return rho, None, None
    Cs = [lag_cov(coords, b) for b in range(1, d)]
    S = [[float(sum(C[a][c] for C in Cs)) for c in range(rho)] for a in range(rho)]
    L = cholesky(S); Li = inv_lower(L)
    Csw = [mm(mm(Li, [[float(x) for x in row] for row in C]), tr_(Li)) for C in Cs]
    Qw, m = search(Csw, rho, rnd)
    if m <= 0: return rho, m, None
    Q = mm(mm(L, Qw), tr_(L))                 # back to lattice coordinates
    for den in (10, 100, 1000, 10000, 100000):
        Qr = [[Fr(round(x * den), den) for x in row] for row in Q]
        Qr = [[(Qr[i][j] + Qr[j][i]) / 2 for j in range(rho)] for i in range(rho)]
        if not posdef_exact(Qr): continue
        Qi = matinv(Qr); ok = True
        for C in Cs:
            t = sum(C[a][c] * Qi[c][a] for a in range(rho) for c in range(rho))
            Sx = [[t / 2 * Qr[a][c] - C[a][c] for c in range(rho)] for a in range(rho)]
            if not posdef_exact(Sx): ok = False; break
        if ok: return rho, m, (den, Qr)
    return rho, m, None
if __name__ == '__main__':
    rnd = random.Random(3)
    for d in (5, 7):
        cand = [1] + [k * d + 1 for k in range(1, 14)]
        seen = set(); maps = []
        for ms in itertools.product(cand, repeat=d):
            if ms[0] != 1 or math.prod(ms) >= d ** d: continue
            key = min(min(tuple(ms[(j + c) % d] for j in range(d)) for c in range(d)),
                      min(tuple(ms[(-j + c) % d] for j in range(d)) for c in range(d)))
            if key in seen: continue
            seen.add(key); maps.append(key)
        res = {}
        for ms in maps:
            rho, m, ver = analyse(ms, d, rnd)
            if rho < 3: res.setdefault(('rank<3', rho), 0); res[('rank<3', rho)] += 1; continue
            tag = 'balanced' if ver else ('unbalanced?' if m is not None and m <= 0 else 'search-fail')
            res.setdefault((tag, rho), []).append((ms, m))
        print(f"d = {d}: {len(maps)} contracting translation-only maps (up to dihedral symmetry)")
        for k, v in sorted(res.items(), key=str):
            if isinstance(v, int): print(f"   {k}: {v}")
            else: print(f"   {k}: {len(v)} maps; examples {[x[0] for x in v[:4]]}; margins {[round(x[1], 4) for x in v[:4]]}", flush=True)
