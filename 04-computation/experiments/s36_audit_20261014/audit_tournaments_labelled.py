#!/usr/bin/env python3
"""S36 referee audit, part A: exhaustive check of the HYP-9280 claims over ALL labelled tournaments, n = 3..NMAX.

Independent of the session's code (see tourn_core.py).  For n <= PERMMAX every tournament's statistics are computed
twice, by DFS path enumeration and by a brute-force scan of all n! orders, and compared entry by entry.
Reports, over labelled strong tournaments:
  counts; M(T) transitive; F(T) strongly connected; inverted-arc distribution and degree types (outdeg y, indeg x) of
  inverted arcs x -> y; squeezed-triangle check; forward 3-cycle; extremal values (best 3-cycle certificate, vertex
  lemma min_w max_{u->w} p(u,w), min D_in, min D_out, max Pr[w first], min delta) with the number of labelled
  tournaments attaining each and how many of those are isomorphic to T_n; the stochastic symmetry (SS) at every
  vertex; violations of 'best in-neighbour + Pr[w first] >= 1'; minimum number of 3-cycles of M(T).
Usage: python audit_tournaments_labelled.py [NMAX=7] [PERMMAX=6]
"""
import sys, time, itertools
from fractions import Fraction
import numpy as np
import numba as nb
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tourn_core import decode, landau_strong, enum_paths, perm_paths, analyse, inmasks, is_Tn, NSTAT

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
PERMMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6

# extreme slots: (stat index, 'min' or 'max')
EXT = [(5, 'min'), (6, 'min'), (7, 'min'), (8, 'min'), (9, 'max'), (12, 'min')]
EXTNAME = ['best 3-cycle certificate', 'vertex lemma min_w max_(u->w) p', 'min_w D_in', 'min_w D_out',
           'max_w Pr[w first]', 'delta (best balance)']


@nb.njit(parallel=True, cache=True)
def pass1(n, pi, pj, nchunk, perms, useperm):
    E = pi.shape[0]
    total = 1 << E
    chunk = (total + nchunk - 1) // nchunk
    agg = np.zeros((nchunk, 16), dtype=np.int64)
    invh = np.zeros((nchunk, 64), dtype=np.int64)
    typh = np.zeros((nchunk, 256), dtype=np.int64)
    exn = np.zeros((nchunk, 6), dtype=np.int64)
    exd = np.zeros((nchunk, 6), dtype=np.int64)
    m3min = np.full(nchunk, 1 << 40, dtype=np.int64)
    for c in nb.prange(nchunk):
        adj = np.zeros(n, dtype=np.int64)
        inm = np.zeros(n, dtype=np.int64)
        C = np.zeros((n, n), dtype=np.int64)
        fi = np.zeros(n, dtype=np.int64)
        Nd = np.zeros((n, n), dtype=np.int64)
        C2 = np.zeros((n, n), dtype=np.int64)
        fi2 = np.zeros(n, dtype=np.int64)
        Nd2 = np.zeros((n, n), dtype=np.int64)
        out = np.zeros(NSTAT, dtype=np.int64)
        th = np.zeros(256, dtype=np.int64)
        init = np.zeros(6, dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            decode(n, code, pi, pj, adj)
            if not landau_strong(n, adj):
                continue
            inmasks(n, adj, inm)
            H = enum_paths(n, adj, C, fi, Nd, inm)
            if useperm:
                H2 = perm_paths(n, adj, perms, C2, fi2, Nd2, inm)
                bad = 0
                if H2 != H:
                    bad = 1
                for i in range(n):
                    if fi[i] != fi2[i]:
                        bad = 1
                    for j in range(n):
                        if C[i, j] != C2[i, j] or Nd[i, j] != Nd2[i, j]:
                            bad = 1
                agg[c, 12] += bad
            for k in range(256):
                th[k] = 0
            analyse(n, adj, C, fi, Nd, inm, H, out, th)
            for k in range(256):
                typh[c, k] += th[k]
            agg[c, 0] += 1
            agg[c, 1] += out[1]
            agg[c, 2] += 1 - out[2]
            agg[c, 3] += 1 - out[4]
            agg[c, 4] += out[10]
            if out[10] > 0:
                agg[c, 5] += 1
            agg[c, 6] += out[11]
            agg[c, 7] += out[15]
            agg[c, 8] += out[16]
            agg[c, 9] += out[17]
            agg[c, 10] += out[18]
            agg[c, 11] += out[19]
            if out[11] > 0:
                agg[c, 13] += 1
            invh[c, min(out[3], 63)] += 1
            if out[13] < m3min[c]:
                m3min[c] = out[13]
            for e in range(6):
                k = 5 + e if e < 4 else (9 if e == 4 else 12)
                num = out[k]
                den = H
                if init[e] == 0:
                    exn[c, e] = num
                    exd[c, e] = den
                    init[e] = 1
                else:
                    if e == 4:
                        better = num * exd[c, e] > exn[c, e] * den
                    else:
                        better = num * exd[c, e] < exn[c, e] * den
                    if better:
                        exn[c, e] = num
                        exd[c, e] = den
    return agg, invh, typh, exn, exd, m3min


@nb.njit(parallel=True, cache=True)
def pass2(n, pi, pj, nchunk, tn, td):
    """count labelled strong tournaments attaining each extreme and how many of them are T_n"""
    E = pi.shape[0]
    total = 1 << E
    chunk = (total + nchunk - 1) // nchunk
    att = np.zeros((nchunk, 6), dtype=np.int64)
    attT = np.zeros((nchunk, 6), dtype=np.int64)
    nTn = np.zeros(nchunk, dtype=np.int64)
    for c in nb.prange(nchunk):
        adj = np.zeros(n, dtype=np.int64)
        inm = np.zeros(n, dtype=np.int64)
        C = np.zeros((n, n), dtype=np.int64)
        fi = np.zeros(n, dtype=np.int64)
        Nd = np.zeros((n, n), dtype=np.int64)
        out = np.zeros(NSTAT, dtype=np.int64)
        th = np.zeros(256, dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            decode(n, code, pi, pj, adj)
            if not landau_strong(n, adj):
                continue
            tt = is_Tn(n, adj)
            if tt:
                nTn[c] += 1
            inmasks(n, adj, inm)
            H = enum_paths(n, adj, C, fi, Nd, inm)
            analyse(n, adj, C, fi, Nd, inm, H, out, th)
            for e in range(6):
                k = 5 + e if e < 4 else (9 if e == 4 else 12)
                if out[k] * td[e] == tn[e] * H:
                    att[c, e] += 1
                    if tt:
                        attT[c, e] += 1
    return att, attT, nTn


def main():
    for n in range(3, NMAX + 1):
        t0 = time.time()
        # pair order deliberately different from the session's: column-major (j outer)
        pairs = [(i, j) for j in range(n) for i in range(j)]
        pi = np.array([p[0] for p in pairs], dtype=np.int64)
        pj = np.array([p[1] for p in pairs], dtype=np.int64)
        useperm = n <= PERMMAX
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64) if useperm else np.zeros((1, n), dtype=np.int64)
        nchunk = 256
        agg, invh, typh, exn, exd, m3min = pass1(n, pi, pj, nchunk, perms, useperm)
        A = agg.sum(axis=0)
        ih = invh.sum(axis=0)
        th = typh.sum(axis=0)
        ext = []
        for e in range(6):
            vals = [Fraction(int(exn[c, e]), int(exd[c, e])) for c in range(nchunk) if exd[c, e] > 0]
            ext.append(max(vals) if e == 4 else min(vals))
        tn = np.array([x.numerator for x in ext], dtype=np.int64)
        td = np.array([x.denominator for x in ext], dtype=np.int64)
        att, attT, nTn = pass2(n, pi, pj, nchunk, tn, td)
        att = att.sum(axis=0); attT = attT.sum(axis=0)
        print(f'=== n = {n} ===  [{time.time() - t0:.1f} s]')
        print(f'  labelled strong tournaments: {A[0]}' + (f'; DFS vs brute-force-permutation mismatches: {A[12]}' if useperm else ' (DFS only)'))
        print(f'  internal consistency errors (Redei parity, C + C^T = H, N-distribution moments): {A[7]}')
        print(f'  M(T) transitive: {A[1]};  F(T) NOT strongly connected: {A[2]};  no all-forward 3-cycle of T: {A[3]}')
        print(f'  inverted-arc distribution: {dict((k, int(ih[k])) for k in range(64) if ih[k])}')
        types = {(k // 16, k % 16): int(th[k]) for k in range(256) if th[k]}
        print(f'  inverted-arc degree types (outdeg y, indeg x) -> #arcs: {types};  (1,1)-arcs whose y-out-neighbour is not the x-in-neighbour: {A[10]};  inverted arcs with p < 1/3: {A[11]}')
        print(f'  min #3-cycles of M(T): {int(m3min.min())} (n - 2 = {n - 2})')
        print(f'  (SS) violated at {A[4]} (T, w) pairs in {A[5]} tournaments;  vertex lemma fails at {A[9]} vertices;  D_in(w) <= 0 at {A[8]} vertices')
        print(f'  "max_(u->w) p(u,w) + Pr[w first] >= 1" violated at {A[6]} (T, w) pairs ({A[13]} tournaments)')
        print(f'  labelled copies of T_n among strong tournaments: {int(nTn.sum())}')
        tval = [Fraction(2 ** (n - 3) + 1, 2 ** (n - 2) + 1), Fraction(2 ** (n - 3) + 1, 2 ** (n - 2) + 1),
                Fraction(1, 2 ** (n - 2) + 1), Fraction(1, 2 ** (n - 2) + 1), Fraction(2 ** (n - 3), 2 ** (n - 2) + 1), None]
        for e in range(6):
            pred = f' (session formula {tval[e]}: {"MATCH" if tval[e] == ext[e] else "DIFFERENT"})' if tval[e] is not None else ''
            print(f'  {EXTNAME[e]}: extreme = {ext[e]}{pred}; attained by {int(att[e])} labelled tournaments, of which {int(attT[e])} are T_n')
        sys.stdout.flush()


if __name__ == '__main__':
    main()
