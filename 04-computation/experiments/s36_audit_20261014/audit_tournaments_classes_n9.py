#!/usr/bin/env python3
"""S36 referee audit, part D: all 191,536 isomorphism classes of 9-vertex tournaments and the HYP-9280 statistics on the
178,133 strong ones (beyond the session's n <= 7 exhaustive range and its pending [N8]).

Classes at n = 8 come from audit_tournaments_classes.canon-style generation (rebuilt here); classes at n = 9 by adding
a ninth vertex with each of the 256 out-sets to every 8-class.  Canonical code: minimum code over the relabelings that
sort vertices by a refined invariant (score, then the sorted multiset of out-neighbour scores, then in-neighbour
score sum) - an isomorphism invariant, so the minimum is a canonical form.  Checks: 191,536 classes (A000568), labelled
total 2^36, strong classes 178,133 (A051337).
Usage: python audit_tournaments_classes_n9.py
"""
import sys, time, itertools, os
from fractions import Fraction
from math import factorial
import numpy as np
import numba as nb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tourn_core import landau_strong, enum_paths, analyse, inmasks, is_Tn, NSTAT, popcount
from audit_tournaments_classes import canon_all, pair_index, code_to_adj, stats_all


@nb.njit(cache=True)
def vertex_invariants(n, adj, inv):
    sc = np.empty(n, dtype=np.int64)
    for v in range(n):
        sc[v] = popcount(adj[v])
    for v in range(n):
        # histogram of out-neighbour scores (base 16 digits) and in-neighbour score sum
        h = 0
        si = 0
        for u in range(n):
            if u == v:
                continue
            if (adj[v] >> u) & 1:
                h += 1 << (4 * sc[u])
            else:
                si += sc[u]
        inv[v] = (sc[v] << 52) + (h << 8) + si   # h < 2^40 for n <= 9 (scores <= 8, counts <= 8 -> 4 bits each)


@nb.njit(parallel=True, cache=True)
def canon_refined(n, cand, idx, permtab, permoff, fact):
    m = cand.shape[0]
    codes = np.zeros(m, dtype=np.int64)
    auts = np.zeros(m, dtype=np.int64)
    for t in nb.prange(m):
        adj = cand[t].copy()
        inv = np.empty(n, dtype=np.int64)
        vertex_invariants(n, adj, inv)
        order = np.argsort(inv)
        # blocks of equal invariant in sorted order
        bstart = np.empty(n, dtype=np.int64)
        bsize = np.empty(n, dtype=np.int64)
        nb_ = 0
        i = 0
        while i < n:
            j = i
            while j + 1 < n and inv[order[j + 1]] == inv[order[i]]:
                j += 1
            bstart[nb_] = i
            bsize[nb_] = j - i + 1
            nb_ += 1
            i = j + 1
        ctr = np.zeros(nb_, dtype=np.int64)
        perm = np.empty(n, dtype=np.int64)
        best = -1
        cnt = 0
        while True:
            # build relabeling: new position bstart[b] + q gets old vertex order[bstart[b] + permtab[...][q]]
            for b in range(nb_):
                k = bsize[b]
                base = permoff[k] + ctr[b] * k
                for q in range(k):
                    perm[bstart[b] + q] = order[bstart[b] + permtab[base + q]]
            code = 0
            for a in range(n):
                row = adj[perm[a]]
                for c in range(a + 1, n):
                    if (row >> perm[c]) & 1:
                        code |= 1 << idx[a, c]
            if best < 0 or code < best:
                best = code
                cnt = 1
            elif code == best:
                cnt += 1
            # advance mixed-radix counter
            b = 0
            while b < nb_:
                ctr[b] += 1
                if ctr[b] < fact[bsize[b]]:
                    break
                ctr[b] = 0
                b += 1
            if b == nb_:
                break
        codes[t] = best
        auts[t] = cnt
    return codes, auts


def main():
    t0 = time.time()
    # flat table of all permutations of 0..k-1 for k <= 9
    fact = np.array([factorial(k) for k in range(10)], dtype=np.int64)
    permoff = np.zeros(10, dtype=np.int64)
    chunks = []
    off = 0
    for k in range(10):
        permoff[k] = off
        P = np.array(list(itertools.permutations(range(k))), dtype=np.int64).reshape(-1) if k > 0 else np.zeros(0, dtype=np.int64)
        chunks.append(P)
        off += len(P)
    permtab = np.concatenate(chunks)
    # rebuild classes to n = 8 with the refined canonical form, checking counts
    classes = [[0]]
    for n in range(2, 10):
        idx = pair_index(n)
        cand = []
        for adj in classes:
            for S in range(1 << (n - 1)):
                new = list(adj) + [S]
                for v in range(n - 1):
                    if not (S >> v) & 1:
                        new[v] |= 1 << (n - 1)
                cand.append(new)
        cand = np.array(cand, dtype=np.int64)
        codes, au = canon_refined(n, cand, idx, permtab, permoff, fact)
        uniq, first_idx = np.unique(codes, return_index=True)
        classes = [code_to_adj(n, int(c)) for c in uniq]
        autn = au[first_idx]
        lab = sum(factorial(n) // int(a) for a in autn)
        print(f'n = {n}: {len(uniq)} classes from {len(cand)} candidates; sum n!/|Aut| = {lab} (2^C(n,2) = {2 ** (n * (n - 1) // 2)})  [{time.time() - t0:.1f} s]', flush=True)
    n = 9
    A = np.array(classes, dtype=np.int64)
    w = np.array([factorial(n) // int(a) for a in autn], dtype=np.int64)
    t1 = time.time()
    res, typ, strong, tn = stats_all(n, A)
    S = strong == 1
    ws = w[S]
    print(f'\n=== n = 9 (all classes) ===  [{time.time() - t1:.1f} s]')
    print(f'  strong classes: {int(S.sum())}; labelled strong: {int(ws.sum())}; T_n classes: {int(tn.sum())}; consistency errors: {int(res[S, 15].sum())}')
    print(f'  classes with M(T) transitive: {int(res[S, 1].sum())}; F(T) NOT strongly connected: {int((res[S, 2] == 0).sum())}; '
          f'no all-forward 3-cycle: {int((res[S, 4] == 0).sum())}')
    inv = res[S, 3]
    print(f'  inverted-arc distribution (classes): {dict((int(k), int((inv == k).sum())) for k in np.unique(inv))}; '
          f'(labelled): {dict((int(k), int(ws[inv == k].sum())) for k in np.unique(inv))}')
    print(f'  min #3-cycles of M(T): {int(res[S, 13].min())} (n - 2 = 7)')
    print(f'  (SS) violations: {int(res[S, 10].sum())} (T, w) pairs in {int((res[S, 10] > 0).sum())} classes ({int(ws[res[S, 10] > 0].sum())} labelled)')
    print(f'  vertex lemma failures: {int(res[S, 17].sum())}; vertices with D_in <= 0: {int(res[S, 16].sum())}')
    print(f'  classes with an inverted arc of p < 1/3: {int((res[S, 19] > 0).sum())}')
    H = res[S, 0]
    for name, k, mode in (('best 3-cycle certificate', 5, 'min'), ('vertex lemma value', 6, 'min'), ('min D_in', 7, 'min'),
                          ('min D_out', 8, 'min'), ('max Pr[w first]', 9, 'max'), ('delta', 12, 'min')):
        vals = [Fraction(int(a), int(b)) for a, b in zip(res[S, k], H)]
        ext = min(vals) if mode == 'min' else max(vals)
        att = [i for i, v in enumerate(vals) if v == ext]
        ntn = sum(int(tn[S][i]) for i in att)
        print(f'  {name}: {ext} attained by {len(att)} classes ({ntn} of them T_n)')
    print(f'  max #inverted arcs: {int(inv.max())}; max H: {int(H.max())}', flush=True)
    Sidx = np.nonzero(S)[0]
    for i in Sidx[res[Sidx, 10] > 0][:12]:
        print(f'    SS violator: adj = {list(map(int, A[i]))}, H = {int(res[i, 0])}, #SS-violating vertices = {int(res[i, 10])}')
    mx = int(res[Sidx, 3].max())
    for i in Sidx[res[Sidx, 3] == mx][:5]:
        print(f'    {mx} inverted arcs: adj = {list(map(int, A[i]))}, H = {int(res[i, 0])}, M 3-cycles = {int(res[i, 13])}')
    print(f'total {time.time() - t0:.1f} s')


if __name__ == '__main__':
    main()
