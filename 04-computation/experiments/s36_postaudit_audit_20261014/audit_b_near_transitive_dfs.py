"""Referee check (B): independent re-implementation of the near-transitive sweep (explicit DFS over Hamiltonian paths,
Landau score test + Warshall closure for strongness, triangle test for transitivity of M).

For each (n, k): all k-sets of arcs lo -> hi with hi - lo >= 2 reversed in the transitive order (higher beats lower).
Per strong tournament: H; M(T) transitive?; F(T) strong?; vertex-lemma margin (2p - 1, best in-neighbour) overall and
by min(indeg, outdeg) = 1, 2, >= 3; majority-source margin; D_in; (SS) failures at k = 0..4 (only 2k + 1 <= d, the
other k being equivalent to k' = d - 1 - k); min E[N] - d/2 by in-degree d <= 6; in-degree of every vertex with D_in <= 0;
min (SS) k = 1 slack; ALSO quantities the session did not compute: the 3-cycle certificate (max over directed 3-cycles
of the least arc probability), max Pr[w first], min D_out.  The T_n layer (k = 1: only the extreme arc gives a strong
tournament) is run separately for n = 4..20.
Usage: python audit_b_near_transitive_dfs.py [n k ...]
"""
import sys, time, itertools
from fractions import Fraction as Fr
import numpy as np
import numba as nb
from pa_core import (enumerate_paths_dfs_first, strong_landau, strong_warshall, transitive_by_triangles, _popc)

BIG = 1 << 62
NS = 40


@nb.njit(cache=False)
def analyze(n, A, C, Nd, L, first, res):
    for i in range(NS):
        res[i] = 0
    sl = strong_landau(n, A)
    sw = strong_warshall(n, A)
    res[0] = 1 if sl else 0
    res[1] = 1 if sw else 0
    if not sl:
        return
    H = enumerate_paths_dfs_first(n, A, C, Nd, L, first)
    res[2] = H
    M = np.zeros((n, n), dtype=np.int8)
    Fm = np.zeros((n, n), dtype=np.int8)
    for a in range(n):
        for b in range(n):
            if a != b and 2 * C[a, b] > H:
                M[a, b] = 1
                if A[a, b]:
                    Fm[a, b] = 1
    res[3] = 1 if transitive_by_triangles(n, M) else 0
    res[4] = 1 if strong_warshall(n, Fm) else 0
    vl = BIG; ms = BIG; dm = BIG; c1 = BIG; c2 = BIG; c3 = BIG; dout = BIG
    for k in range(16, 22):
        res[k] = BIG
    res[24] = BIG
    res[25] = BIG
    for w in range(n):
        best = -BIG; anyb = -BIG; dsum = 0; d = 0; od = 0; dosum = 0
        for u in range(n):
            if u == w:
                continue
            mg = 2 * C[u, w] - H
            if mg > anyb:
                anyb = mg
            if A[u, w]:
                d += 1
                dsum += mg
                if mg > best:
                    best = mg
            else:
                od += 1
                dosum += 2 * C[w, u] - H
        cls = min(d, od)
        if best < vl:
            vl = best
        if cls == 1 and best < c1:
            c1 = best
        if cls == 2 and best < c2:
            c2 = best
        if cls >= 3 and best < c3:
            c3 = best
        if anyb < ms:
            ms = anyb
        if dsum < dm:
            dm = dsum
        if dosum < dout:
            dout = dosum
        if best <= 0:
            res[27] += 1
        if anyb <= 0:
            res[28] += 1
        if dsum <= 0:
            res[23] += 1
            if d < res[24]:
                res[24] = d
        # (SS)
        for kk in range(5):
            if 2 * kk + 1 > d:
                break
            hi_ = 0; lo_ = 0
            for j in range(d - kk, d + 1):
                hi_ += Nd[w, j]
            for j in range(0, kk + 1):
                lo_ += Nd[w, j]
            if hi_ < lo_:
                res[11 + kk] += 1
            if kk == 1 and hi_ - lo_ < res[25]:
                res[25] = hi_ - lo_
        if 1 <= d <= 6:
            s1 = 0
            for j in range(d + 1):
                s1 += j * Nd[w, j]
            val = 2 * s1 - d * H
            if val < res[15 + d]:
                res[15 + d] = val
        if first[w] > res[29]:
            res[29] = first[w]
    res[5] = vl; res[6] = ms; res[7] = dm; res[8] = c1; res[9] = c2; res[10] = c3; res[30] = dout
    # 3-cycle certificate: max over directed 3-cycles x->y->z->x of min(C[x,y], C[y,z], C[z,x])
    cert = -1
    for x in range(n):
        for y in range(n):
            if y == x or A[x, y] == 0:
                continue
            for z in range(n):
                if z == x or z == y:
                    continue
                if A[y, z] and A[z, x]:
                    m = C[x, y]
                    if C[y, z] < m:
                        m = C[y, z]
                    if C[z, x] < m:
                        m = C[z, x]
                    if m > cert:
                        cert = m
    res[22] = cert


@nb.njit(parallel=True, cache=False)
def sweep(n, arcsets, nthreads):
    m = arcsets.shape[0]; k = arcsets.shape[1]
    out = np.zeros((m, NS), dtype=np.int64)
    chunk = (m + nthreads - 1) // nthreads
    for t in nb.prange(nthreads):
        A = np.zeros((n, n), dtype=np.int8)
        C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64)
        L = np.zeros((n, n), dtype=np.int64); first = np.zeros(n, dtype=np.int64)
        res = np.zeros(NS, dtype=np.int64)
        for idx in range(t * chunk, min(m, (t + 1) * chunk)):
            for a in range(n):
                for b in range(n):
                    A[a, b] = 1 if b < a else 0
            for e in range(k):
                lo = arcsets[idx, e, 0]; hi = arcsets[idx, e, 1]
                A[hi, lo] = 0; A[lo, hi] = 1
            analyze(n, A, C, Nd, L, first, res)
            for i in range(NS):
                out[idx, i] = res[i]
    return out


def fmt_min(o, col, H):
    rel = o[:, col] / H
    i = int(np.argmin(rel))
    return i, f'{int(o[i, col])}/{int(o[i, 2])} = {rel[i]:.6f}'


def run(n, k):
    cand = [(lo, hi) for lo in range(n) for hi in range(lo + 2, n)]
    combos = np.array(list(itertools.combinations(cand, k)), dtype=np.int64).reshape(-1, k, 2)
    t0 = time.time()
    out = sweep(n, combos, 12)
    s = out[:, 0] == 1
    mism = int((out[:, 0] != out[:, 1]).sum())
    o = out[s]; cs = combos[s]
    print(f'n = {n}, k = {k}: {len(combos)} arc sets, {int(s.sum())} strong (Landau; Landau/Warshall mismatches {mism})  [{time.time() - t0:.1f} s]')
    H = o[:, 2].astype(float)
    print(f'  H odd for all: {bool((o[:, 2] % 2 == 1).all())};  M(T) transitive: {int(o[:, 3].sum())};  F(T) not strong: {int((o[:, 4] == 0).sum())}')
    for col, name in ((5, 'vertex-lemma margin 2p-1'), (6, 'majority-source margin'), (7, 'D_in'), (30, 'D_out')):
        i, txt = fmt_min(o, col, H)
        print(f'  min {name} = {txt} at {[tuple(map(int, a)) for a in cs[i]]}')
    print(f'  (T, w) vertex-lemma failures: {int(o[:, 27].sum())}; majority-source failures: {int(o[:, 28].sum())}; '
          f'tournaments with some D_in <= 0: {int((o[:, 23] > 0).sum())} ((T, w) pairs: {int(o[:, 23].sum())}; '
          f'least in-degree at such a vertex: {int(o[:, 24].min()) if (o[:, 23] > 0).any() else "-"})')
    cls = []
    for col, name in ((8, '1'), (9, '2'), (10, '>= 3')):
        mask = o[:, col] < BIG // 2
        if mask.any():
            rel = o[mask, col] / H[mask]
            i = int(np.argmin(rel))
            cls.append(f'{name}: {int(o[mask][i, col])}/{int(o[mask][i, 2])} = {rel[i]:.6f}')
        else:
            cls.append(f'{name}: -')
    print('  min vertex-lemma margin by min(indeg, outdeg): ' + '; '.join(cls))
    ss = o[:, 11:16].sum(axis=0)
    dd = []
    for d in range(1, 7):
        col = 15 + d
        mask = o[:, col] < BIG // 2
        if mask.any():
            rel = o[mask, col] / (2 * H[mask])
            dd.append(f'd={d}: {rel.min():.4f}')
    print('  (SS) failures (T, w) by k: ' + ', '.join(f'k={kk}: {int(ss[kk])}' for kk in range(5)) + ';  min E[N] - d/2 by d: ' + ', '.join(dd))
    rel = o[:, 25] / H
    i = int(np.argmin(rel))
    print(f'  min (SS) k=1 slack [Pr(N >= d-1) - Pr(N <= 1)] over d >= 3: {int(o[i, 25])}/{int(o[i, 2])} = {rel[i]:.6f}')
    Tn = (2 ** (n - 3) + 1) / (2 ** (n - 2) + 1)
    TnF = Fr(2 ** (n - 3) + 1, 2 ** (n - 2) + 1)
    crel = o[:, 22] / H
    i = int(np.argmin(crel))
    below = sum(1 for r in o if Fr(int(r[22]), int(r[2])) < TnF)
    eq = sum(1 for r in o if Fr(int(r[22]), int(r[2])) == TnF)
    print(f'  3-cycle certificate: min {int(o[i, 22])}/{int(o[i, 2])} = {crel[i]:.6f} at {[tuple(map(int, a)) for a in cs[i]]}; '
          f'T_n value {TnF} = {Tn:.6f}; tournaments below it: {below}, equal: {eq}')
    frel = o[:, 29] / H
    i = int(np.argmax(frel))
    Tf = Fr(2 ** (n - 3), 2 ** (n - 2) + 1)
    above = sum(1 for r in o if Fr(int(r[29]), int(r[2])) >= Tf)
    print(f'  max Pr[w first] = {int(o[i, 29])}/{int(o[i, 2])} = {frel[i]:.6f} (T_n value {Tf} = {float(Tf):.6f}; tournaments >= it: {above})', flush=True)
    return o, cs


def run_Tn(nmax):
    print('T_n layer (k = 1; transitive + extreme arc reversed):')
    for n in range(4, nmax + 1):
        t0 = time.time()
        combos = np.array([[[0, n - 1]]], dtype=np.int64)
        o = sweep(n, combos, 1)[0]
        H = int(o[2])
        TnF = Fr(2 ** (n - 3) + 1, 2 ** (n - 2) + 1)
        ss = o[11:16]
        print(f'  n = {n}: H = {H} (= 2^(n-2)+1: {H == 2 ** (n - 2) + 1}); M transitive {int(o[3])}; F strong {int(o[4])}; '
              f'VL margin {int(o[5])}/{H}; MS margin {int(o[6])}/{H}; D_in min {int(o[7])}/{H}; D_out min {int(o[30])}/{H}; '
              f'(SS) failures k=0..4: {list(map(int, ss))}; min (SS) k=1 slack {int(o[25])}/{H}; '
              f'cert {Fr(int(o[22]), H)} (= T_n value: {Fr(int(o[22]), H) == TnF}); max first {Fr(int(o[29]), H)}  [{time.time() - t0:.1f} s]', flush=True)


if __name__ == '__main__':
    args = [int(a) for a in sys.argv[1:]]
    if args:
        jobs = list(zip(args[0::2], args[1::2]))
        for n, k in jobs:
            if k == 1:
                run_Tn(n)
            else:
                run(n, k)
    else:
        run_Tn(20)
        for n, k in [(11, 2), (12, 2), (10, 3), (11, 3), (10, 2), (13, 2), (14, 2), (16, 2), (12, 3), (10, 4), (11, 4), (10, 5)]:
            run(n, k)
