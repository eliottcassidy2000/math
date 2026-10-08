#!/usr/bin/env python3
"""Balanced lag covariances for translation-only Matthews-Watts maps (all m_i = 1 mod d).
For such maps the coupling is always a translation j -> j + b (b = e mod d), so the debt step given the past is
the uniform law on the d 'cyclic roots' v_(j+b) - v_j (v_j = prime-exponent vector of m_j).  Its covariance is
C_b = (1/d) sum_j (v_(j+b) - v_j)(v_(j+b) - v_j)^T.  Lamperti: if one inner product Q has
   (1/2) tr(C_b Q^-1) Q - C_b  positive definite on the debt lattice span for every b != 0,
then V = |x|_Q^-alpha is a supermartingale far out and the debt walk is transient (rank >= 3 needed).
This script enumerates contracting translation-only maps on Z_d (d = 5, 7), searches Q, verifies exactly."""
import itertools, random, math
from fractions import Fraction as Fr
def factor(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def rank_and_coords(ms):
    primes = sorted({p for m in ms for p in factor(m)})
    vecs = [[Fr(factor(m).get(p, 0)) for p in primes] for m in ms]
    diffs = [[a - b for a, b in zip(vecs[i], vecs[0])] for i in range(len(ms))]
    # row-reduce to find a basis of span{v_i - v_0}
    basis = []
    for v in diffs:
        w = v[:]
        for bvec, piv in basis:
            if w[piv] != 0:
                c = w[piv] / bvec[piv]; w = [x - c * y for x, y in zip(w, bvec)]
        nz = [k for k, x in enumerate(w) if x != 0]
        if nz: basis.append((w, nz[0]))
    rho = len(basis)
    # coordinates of each v_i - v_0 in the (non-orthogonal) basis: solve by least squares over Q (exact via elimination)
    B = [b for b, _ in basis]
    def coords(v):
        # solve sum c_k B_k = v
        n = len(B); A = [[B[k][r] for k in range(n)] + [v[r]] for r in range(len(v))]
        # gaussian elimination on overdetermined consistent system
        rows = [row[:] for row in A]; piv_cols = []; r = 0
        for c in range(n):
            pr = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
            if pr is None: continue
            rows[r], rows[pr] = rows[pr], rows[r]
            pv = rows[r][c]; rows[r] = [x / pv for x in rows[r]]
            for i in range(len(rows)):
                if i != r and rows[i][c] != 0:
                    f = rows[i][c]; rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
            piv_cols.append(c); r += 1
        sol = [Fr(0)] * n
        for i, c in enumerate(piv_cols): sol[c] = rows[i][n]
        return sol
    return rho, [coords(d_) for d_ in diffs]
def lag_cov(coords, b):
    d = len(coords); rho = len(coords[0])
    C = [[Fr(0)] * rho for _ in range(rho)]
    for j in range(d):
        w = [x - y for x, y in zip(coords[(j + b) % d], coords[j])]
        for a in range(rho):
            for c in range(rho): C[a][c] += w[a] * w[c] / d
    return C
def matinv(M):
    n = len(M); A = [list(map(Fr, row)) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        pr = next(i for i in range(c, n) if A[i][c] != 0); A[c], A[pr] = A[pr], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for i in range(n):
            if i != c and A[i][c] != 0:
                f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return [row[n:] for row in A]
def det(M):
    n = len(M); A = [row[:] for row in M]; dv = Fr(1)
    for c in range(n):
        pr = next((i for i in range(c, n) if A[i][c] != 0), None)
        if pr is None: return Fr(0)
        if pr != c: A[c], A[pr] = A[pr], A[c]; dv = -dv
        dv *= A[c][c]
        for i in range(c + 1, n):
            f = A[i][c] / A[c][c]; A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return dv
def posdef_exact(S):
    return all(det([row[:k] for row in S[:k]]) > 0 for k in range(1, len(S) + 1))
def margin(Q, Cs):
    # float margin: min over b of lambda_min( (1/2)tr(C Q^-1) Q - C ) / tr-scale, via Jacobi
    Qi = matinv([[Fr(x) for x in row] for row in Q])
    worst = 1e9
    for C in Cs:
        t = sum(sum(C[a][c] * Qi[c][a] for c in range(len(C))) for a in range(len(C)))
        S = [[float(t / 2 * Q[a][c] - C[a][c]) for c in range(len(C))] for a in range(len(C))]
        worst = min(worst, min(jacobi_eigs(S)) / float(t))
    return worst
def jacobi_eigs(A):
    n = len(A); A = [row[:] for row in A]
    for _ in range(100):
        off = max((abs(A[i][j]), i, j) for i in range(n) for j in range(n) if i != j) if n > 1 else (0, 0, 0)
        if off[0] < 1e-13: break
        _, p, q = off
        th = 0.5 * math.atan2(2 * A[p][q], A[q][q] - A[p][p])
        c, s = math.cos(th), math.sin(th)
        for k in range(n):
            akp, akq = A[k][p], A[k][q]
            A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq
        for k in range(n):
            apk, aqk = A[p][k], A[q][k]
            A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk
    return [A[i][i] for i in range(n)]
def search_Q(Cs, rho, rnd, iters=4000):
    best = None; bestm = -1e9
    Q = [[float(i == j) for j in range(rho)] for i in range(rho)]
    for it in range(iters):
        if best is None or rnd.random() < 0.3:
            L = [[rnd.uniform(-1, 1) if j <= i else 0.0 for j in range(rho)] for i in range(rho)]
            for i in range(rho): L[i][i] = abs(L[i][i]) + 0.3
        else:
            L = [[best_L[i][j] + rnd.gauss(0, 0.05) if j <= i else 0.0 for j in range(rho)] for i in range(rho)]
        Q = [[sum(L[i][k] * L[j][k] for k in range(rho)) for j in range(rho)] for i in range(rho)]
        try: m = margin(Q, Cs)
        except Exception: continue
        if m > bestm: bestm, best, best_L = m, Q, L
    return best, bestm
def rationalize(Q, den=20):
    return [[Fr(round(x * den), den) for x in row] for row in Q]
def check_exact(Qr, Cs):
    if not posdef_exact(Qr): return False
    Qi = matinv(Qr); out = []
    for C in Cs:
        t = sum(sum(C[a][c] * Qi[c][a] for c in range(len(C))) for a in range(len(C)))
        S = [[t / 2 * Qr[a][c] - C[a][c] for c in range(len(C))] for a in range(len(C))]
        out.append(posdef_exact(S))
    return all(out)
if __name__ == '__main__':
    rnd = random.Random(7)
    for d in (5, 7):
        cand = [1] + [k * d + 1 for k in range(1, 12)]
        seen = set(); maps = []
        for ms in itertools.product(cand, repeat=d):
            if ms[0] != 1: continue                      # normalize: put a 1 at position 0 (rotation)
            if math.prod(ms) >= d ** d: continue
            key = min(tuple(ms[(j + c) % d] for j in range(d)) for c in range(d))
            key = min(key, min(tuple(ms[(-j + c) % d] for j in range(d)) for c in range(d)))
            if key in seen: continue
            seen.add(key); maps.append(key)
        print(f"d = {d}: {len(maps)} contracting translation-only maps up to rotation/reflection")
        stats = {}
        for ms in maps:
            rho, coords = rank_and_coords(list(ms))
            stats.setdefault(rho, []).append(ms)
            if rho < 3: continue
            Cs = [lag_cov(coords, b) for b in range(1, d)]
            Q, m = search_Q(Cs, rho, rnd, iters=1500 if d == 5 else 600)
            ok = False
            if m > 0:
                for den in (4, 10, 20, 50, 100):
                    Qr = rationalize(Q, den)
                    if check_exact(Qr, Cs): ok = True; break
            Lam = sum(math.log(x / d) for x in ms) / d
            lagranks = []
            for C in Cs:
                ev = sorted(jacobi_eigs([[float(x) for x in row] for row in C]))
                lagranks.append(sum(1 for x in ev if x > 1e-9))
            print(f"   m = {ms}  rank {rho}  Lambda {Lam:+.3f}  lag-cov ranks {lagranks}  best float margin {m:+.4f}  exact balanced Q: {'YES (den %d) Q=%s' % (den, [[str(x) for x in row] for row in Qr]) if ok else 'no'}", flush=True)
        print(f"   rank histogram: { {r: len(v) for r, v in sorted(stats.items())} }")
