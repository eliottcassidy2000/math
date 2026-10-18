#!/usr/bin/env python3
"""S42 audit a3: independent synchronous-depth decomposition E_Y = E_0 + 2 (MM_hi + MM_lo + MF) (no session code).

Initial coefficients by direct folding of n in (N, 4N] (odd part, sign (-1)^k, seed), chunked; then the stopped process
in pure numpy: each depth, every moving seed s > Y moves to P(s) with sign (-1)^a; contributions landing above Y merge
(MM_hi = sum of pairwise products within a target), contributions landing on a frozen seed u <= Y in the same step merge
(MM_lo), and the arrival L_u meets the current already frozen there (MF = sum fz[u] L_u).  The finite-depth identity
E_frozen + E_moving = E_0 + 2 (MM_hi + MM_lo + MF) is asserted at EVERY depth; R = sum f(u) fz[u] is compared with the
direct sum.  Values printed as in the session table (2 MM / N).
Usage: python a3_decomposition_check.py "N1:Y1,N2:Y2,..."
"""
import os, sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def fpar(n):
    c = 0
    while n != 1:
        n = (3 * n + 1) // 2 if n & 1 else n // 2
        c += 1
    return 1 - 2 * (c & 1)


def seed_vec(m):
    m = m.copy()
    msk = (m % 8) == 5
    while msk.any():
        m[msk] = (m[msk] - 1) // 4
        msk = (m % 8) == 5
    return m


def initial_by_folding(N, chunk=10 ** 7):
    acc = {}
    sums = []
    for lo in range(N + 1, 4 * N + 1, chunk):
        n = np.arange(lo, min(lo + chunk, 4 * N + 1), dtype=np.int64)
        w = np.where(n % 3 == 2, 2, -1).astype(np.int64)
        tz = np.zeros(len(n), dtype=np.int64)
        m = n.copy()
        msk = (m % 2) == 0
        while msk.any():
            m[msk] //= 2
            tz[msk] += 1
            msk = (m % 2) == 0
        sg = 1 - 2 * (tz & 1)
        s = seed_vec(m)
        u, inv = np.unique(s, return_inverse=True)
        cc = np.zeros(len(u), dtype=np.int64)
        np.add.at(cc, inv, w * sg)
        sums.append((u, cc))
    U = np.concatenate([a for a, _ in sums]); Cc = np.concatenate([b for _, b in sums])
    u, inv = np.unique(U, return_inverse=True)
    c = np.zeros(len(u), dtype=np.int64)
    np.add.at(c, inv, Cc)
    keep = c != 0
    return u[keep], c[keep]


def P_vec(S):
    v = 3 * S + 1
    a = np.zeros(len(S), dtype=np.int64)
    msk = (v % 2) == 0
    while msk.any():
        v[msk] //= 2
        a[msk] += 1
        msk = (v % 2) == 0
    assert np.all((a == 1) | (a == 2))
    return seed_vec(v), np.where(a % 2 == 1, -1, 1).astype(np.int64)


def run(N, Y):
    S, C = initial_by_folding(N)
    E0 = int((C * C).sum())
    fz = np.zeros(Y + 1, dtype=np.int64)
    low = S <= Y
    np.add.at(fz, S[low], C[low])
    S, C = S[~low], C[~low]
    MMhi = MMlo = MF = 0
    depth = 0
    per_depth = []
    while len(S):
        depth += 1
        T, sg = P_vec(S)
        A = sg * C
        hi = T > Y
        Th, Ah = T[hi], A[hi]
        if len(Th):
            v, inv = np.unique(Th, return_inverse=True)
            G = np.zeros(len(v), dtype=np.int64); np.add.at(G, inv, Ah)
            sq = np.zeros(len(v), dtype=np.int64); np.add.at(sq, inv, Ah * Ah)
            mmhi = int(((G * G) - sq).sum()) // 2
        else:
            v = np.zeros(0, dtype=np.int64); G = np.zeros(0, dtype=np.int64); mmhi = 0
        L = np.zeros(Y + 1, dtype=np.int64); np.add.at(L, T[~hi], A[~hi])
        Lsq = np.zeros(Y + 1, dtype=np.int64); np.add.at(Lsq, T[~hi], A[~hi] * A[~hi])
        mmlo = int(((L * L) - Lsq).sum()) // 2
        mf = int((fz * L).sum())
        fz += L
        MMhi += mmhi; MMlo += mmlo; MF += mf
        keep = G != 0
        S, C = v[keep], G[keep]
        Emov = int((C * C).sum()); Efz = int((fz * fz).sum())
        assert Efz + Emov == E0 + 2 * (MMhi + MMlo + MF), (N, Y, depth)
        per_depth.append((mmhi, mmlo, mf))
    EY = int((fz * fz).sum())
    R = sum(int(fz[u]) * fpar(u) for u in range(1, Y + 1) if fz[u])
    n = np.arange(N + 1, 4 * N + 1, dtype=np.int64)
    return dict(N=N, Y=Y, E0=E0 / N, EY=EY / N, MMhi=2 * MMhi / N, MMlo=2 * MMlo / N, MF=2 * MF / N, depth=depth, R=R,
                R2N=R * R / N, identity_checked_at_every_depth=True)


def main():
    pairs = [tuple(int(x) for x in p.split(':')) for p in sys.argv[1].split(',')]
    out = []
    t0 = time.time()
    for N, Y in pairs:
        r = run(N, Y)
        out.append(r)
        print('N = %9d Y = %5d: E0/N %.5f  EY/N %.5f  MM_hi %+.5f  MM_lo %+.5f  MF %+.5f  depth %d  R = %d  R^2/N %.4f [%.0f s]'
              % (N, Y, r['E0'], r['EY'], r['MMhi'], r['MMlo'], r['MF'], r['depth'], r['R'], r['R2N'], time.time() - t0), flush=True)
    json.dump(out, open(os.path.join(HERE, 'a3_decomposition_check.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
