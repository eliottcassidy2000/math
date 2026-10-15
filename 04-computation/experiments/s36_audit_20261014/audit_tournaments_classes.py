#!/usr/bin/env python3
"""S36 referee audit, part C: isomorphism classes of tournaments up to n = NMAX (default 8), and the HYP-9280 statistics
on every strong class (an independent check of the pending [N8] sweep).

Generation: classes at n from classes at n - 1 by adding a vertex with every out-set; canonical code = minimum code
over all relabelings that list vertices by non-decreasing score.  Checks: class counts 1, 1, 2, 4, 12, 56, 456, 6880
(OEIS A000568), labelled total 2^C(n,2) via n!/|Aut|, strong classes 1, 0, 1, 1, 6, 35, 353, 6008 (A051337) and
labelled strong counts (A054946).
Statistics per strong class are those of tourn_core.analyse (DFS enumeration of all Hamiltonian paths).
Usage: python audit_tournaments_classes.py [NMAX=8]
"""
import sys, time, itertools, os
from fractions import Fraction
from math import factorial
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tourn_core import landau_strong, enum_paths, analyse, inmasks, is_Tn, NSTAT, popcount

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8


@nb.njit(cache=True)
def pair_index(n):
    idx = np.full((n, n), -1, dtype=np.int64)
    e = 0
    for i in range(n):
        for j in range(i + 1, n):
            idx[i, j] = e
            e += 1
    return idx


@nb.njit(parallel=True, cache=True)
def canon_all(n, cand, perms, idx):
    """cand: (m, n) adjacency rows. Returns canonical codes and automorphism counts."""
    m = cand.shape[0]
    codes = np.zeros(m, dtype=np.int64)
    auts = np.zeros(m, dtype=np.int64)
    for t in nb.prange(m):
        sc = np.empty(n, dtype=np.int64)
        for v in range(n):
            sc[v] = popcount(cand[t, v])
        best = -1
        cnt = 0
        for r in range(perms.shape[0]):
            ok = True
            for i in range(n - 1):
                if sc[perms[r, i]] > sc[perms[r, i + 1]]:
                    ok = False
                    break
            if not ok:
                continue
            code = 0
            for i in range(n):
                a = cand[t, perms[r, i]]
                for j in range(i + 1, n):
                    if (a >> perms[r, j]) & 1:
                        code |= 1 << idx[i, j]
            if best < 0 or code < best:
                best = code
                cnt = 1
            elif code == best:
                cnt += 1
        codes[t] = best
        auts[t] = cnt
    return codes, auts


def code_to_adj(n, code):
    adj = [0] * n
    e = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> e) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            e += 1
    return adj


@nb.njit(parallel=True, cache=True)
def stats_all(n, A):
    m = A.shape[0]
    res = np.zeros((m, NSTAT), dtype=np.int64)
    typ = np.zeros((m, 256), dtype=np.int64)
    strong = np.zeros(m, dtype=np.int64)
    tn = np.zeros(m, dtype=np.int64)
    for t in nb.prange(m):
        adj = A[t].copy()
        if not landau_strong(n, adj):
            continue
        strong[t] = 1
        if is_Tn(n, adj):
            tn[t] = 1
        inm = np.zeros(n, dtype=np.int64)
        inmasks(n, adj, inm)
        C = np.zeros((n, n), dtype=np.int64)
        fi = np.zeros(n, dtype=np.int64)
        Nd = np.zeros((n, n), dtype=np.int64)
        H = enum_paths(n, adj, C, fi, Nd, inm)
        out = np.zeros(NSTAT, dtype=np.int64)
        th = np.zeros(256, dtype=np.int64)
        analyse(n, adj, C, fi, Nd, inm, H, out, th)
        for k in range(NSTAT):
            res[t, k] = out[k]
        for k in range(256):
            typ[t, k] = th[k]
    return res, typ, strong, tn


def main():
    classes = {1: [[0]]}
    auts = {1: [1]}
    t0 = time.time()
    for n in range(2, NMAX + 1):
        idx = pair_index(n)
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
        cand = []
        for adj in classes[n - 1]:
            for S in range(1 << (n - 1)):
                new = list(adj) + [S]
                for v in range(n - 1):
                    if not (S >> v) & 1:
                        new[v] |= 1 << (n - 1)
                cand.append(new)
        cand = np.array(cand, dtype=np.int64)
        codes, au = canon_all(n, cand, perms, idx)
        uniq, first_idx = np.unique(codes, return_index=True)
        classes[n] = [code_to_adj(n, int(c)) for c in uniq]
        auts[n] = [int(au[i]) for i in first_idx]
        lab = sum(factorial(n) // a for a in auts[n])
        print(f'n = {n}: {len(uniq)} classes from {len(cand)} candidates; sum n!/|Aut| = {lab} (2^C(n,2) = {2 ** (n * (n - 1) // 2)})  [{time.time() - t0:.1f} s]', flush=True)
    for n in range(3, NMAX + 1):
        A = np.array(classes[n], dtype=np.int64)
        w = np.array([factorial(n) // a for a in auts[n]], dtype=np.int64)
        t1 = time.time()
        res, typ, strong, tn = stats_all(n, A)
        S = strong == 1
        print(f'\n=== n = {n} ===  [{time.time() - t1:.1f} s]')
        print(f'  strong classes: {int(S.sum())}; labelled strong: {int(w[S].sum())}; T_n classes among them: {int(tn.sum())}')
        print(f'  consistency errors: {int(res[S, 15].sum())}')
        print(f'  classes with M(T) transitive: {int(res[S, 1].sum())}; with F(T) NOT strongly connected: {int((res[S, 2] == 0).sum())}; '
              f'without an all-forward 3-cycle: {int((res[S, 4] == 0).sum())}')
        inv = res[S, 3]; ws = w[S]
        print(f'  inverted-arc distribution (classes): {dict((int(k), int((inv == k).sum())) for k in np.unique(inv))}; '
              f'(labelled): {dict((int(k), int(ws[inv == k].sum())) for k in np.unique(inv))}')
        th = (typ[S] * ws[:, None]).sum(axis=0)
        print(f'  inverted-arc degree types (outdeg y, indeg x) -> #labelled arcs: {dict(((k // 16, k % 16), int(th[k])) for k in range(256) if th[k])}')
        print(f'  squeezed-triangle failures: {int(res[S, 18].sum())}; inverted arcs with p < 1/3 (classes having one): {int((res[S, 19] > 0).sum())}')
        print(f'  min #3-cycles of M(T): {int(res[S, 13].min())} (n - 2 = {n - 2})')
        print(f'  (SS) violations: {int(res[S, 10].sum())} (T, w) pairs in {int((res[S, 10] > 0).sum())} classes '
              f'({int(ws[res[S, 10] > 0].sum())} labelled tournaments)')
        print(f'  vertex lemma failures: {int(res[S, 17].sum())}; vertices with D_in <= 0: {int(res[S, 16].sum())}')
        H = res[S, 0]
        for name, k, mode in (('best 3-cycle certificate', 5, 'min'), ('vertex lemma value', 6, 'min'), ('min D_in', 7, 'min'),
                              ('min D_out', 8, 'min'), ('max Pr[w first]', 9, 'max'), ('delta', 12, 'min')):
            vals = [Fraction(int(a), int(b)) for a, b in zip(res[S, k], H)]
            ext = min(vals) if mode == 'min' else max(vals)
            att = [i for i, v in enumerate(vals) if v == ext]
            ntn = sum(int(tn[S][i]) for i in att)
            print(f'  {name}: {ext} attained by {len(att)} classes ({ntn} of them T_n)')
        print(f'  max #inverted arcs: {int(inv.max())}; max H: {int(H.max())}', flush=True)
        if n == NMAX:
            # list SS violators and the classes with the most inverted arcs
            Sidx = np.nonzero(S)[0]
            for i in Sidx[res[Sidx, 10] > 0][:10]:
                print(f'    SS violator: adj = {list(map(int, A[i]))}, H = {int(res[i, 0])}')
            mx = int(res[Sidx, 3].max())
            for i in Sidx[res[Sidx, 3] == mx][:5]:
                print(f'    {mx} inverted arcs: adj = {list(map(int, A[i]))}, H = {int(res[i, 0])}, M 3-cycles = {int(res[i, 13])}')


if __name__ == '__main__':
    main()
