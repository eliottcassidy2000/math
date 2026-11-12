#!/usr/bin/env python3
"""Verifier: independent recomputation of M6 (arc cube n = 4..7), M3 (n = 5, 6 inputs), table sanity checks.
usage: python v1_tables_cubes.py"""
import sys, time, itertools, random, math
from collections import Counter
import numpy as np
from numba import njit, prange
from vlib import *


@njit(parallel=True, cache=True)
def pair_scan(H, S, E, F):
    npairs = E.shape[0]
    N = H.shape[0]
    out = np.zeros((npairs, 5), np.int64)
    for p in prange(npairs):
        be = 1 << E[p]; bf = 1 << F[p]
        free = (N - 1) ^ (be | bf)
        m = 0
        tot = 0; z = 0; mn = 1 << 30; nt = 0; nts = 0
        while True:
            d = H[m] - H[m | be] - H[m | bf] + H[m | be | bf]
            tot += 1
            a = d if d >= 0 else -d
            if a < mn:
                mn = a
            if d == 0:
                z += 1
                if H[m | be] != H[m] and H[m | bf] != H[m]:
                    nt += 1
                    if S[m] and S[m | be] and S[m | bf] and S[m | be | bf]:
                        nts += 1
            m = ((m | ~free) + 1) & free
            if m == 0:
                break
        out[p, 0] = tot; out[p, 1] = z; out[p, 2] = mn; out[p, 3] = nt; out[p, 4] = nts
    return out


@njit(cache=True)
def cube3_one(H, S, e, f, g, rec, cap):
    N = H.shape[0]
    be = 1 << e; bf = 1 << f; bg = 1 << g
    free = (N - 1) ^ (be | bf | bg)
    m = 0
    c0 = 0; c1 = 0; c2 = 0; k = 0
    done = False
    while not done:
        h0 = H[m]; de = H[m | be] - h0; df = H[m | bf] - h0; dg = H[m | bg] - h0
        ok = H[m | be | bf] == h0 + de + df
        ok = ok and H[m | be | bg] == h0 + de + dg
        ok = ok and H[m | bf | bg] == h0 + df + dg
        ok = ok and H[m | be | bf | bg] == h0 + de + df + dg
        if ok:
            c0 += 1
            if de != 0 and df != 0 and dg != 0:
                c1 += 1
                st = 1
                for a in range(8):
                    mm = m
                    if a & 1: mm |= be
                    if a & 2: mm |= bf
                    if a & 4: mm |= bg
                    if S[mm] == 0:
                        st = 0
                if st:
                    c2 += 1
                if k < cap:
                    rec[k, 0] = h0; rec[k, 1] = de; rec[k, 2] = df; rec[k, 3] = dg; rec[k, 4] = st
                k += 1
        m = ((m | ~free) + 1) & free
        if m == 0:
            done = True
    return c0, c1, c2, k


@njit(parallel=True, cache=True)
def cube3_scan(H, S, T, cap):
    nt = T.shape[0]
    out = np.zeros((nt, 3), np.int64)
    rec = np.zeros((nt, cap, 5), np.int32)
    nrec = np.zeros(nt, np.int64)
    for p in prange(nt):
        c0, c1, c2, k = cube3_one(H, S, T[p, 0], T[p, 1], T[p, 2], rec[p], cap)
        out[p, 0] = c0; out[p, 1] = c1; out[p, 2] = c2; nrec[p] = k
    return out, rec, nrec


def main():
    random.seed(20261111)
    tables = {}
    for n in range(3, 8):
        pi, pj = pair_arrays(n)
        t0 = time.time()
        H, S = build_table(n, pi, pj)
        tables[n] = (H, S)
        N = H.shape[0]
        vals = sorted(set(H.tolist()))
        tot = int(H.astype(np.int64).sum())
        expect = math.factorial(n) * 2 ** (n * (n - 1) // 2 - (n - 1))
        smin = int(H[S == 1].min()) if S.sum() else None
        svals = sorted(set(H[S == 1].tolist()))
        print(f"n = {n}: {N} labelled tournaments ({time.time()-t0:.1f}s); sum H = {tot} (n! 2^(C(n,2)-n+1) = {expect}, "
              f"match {tot == expect}); max H = {int(H.max())}; distinct values {len(vals)}; strong: {int(S.sum())}; "
              f"min H over strong = {smin}; all odd: {bool(np.all(H % 2 == 1))}")
        if n <= 6:
            print(f"   values: {vals};  strong values: {svals}")
        else:
            miss = [v for v in range(1, int(H.max()) + 1, 2) if v not in set(vals)]
            print(f"   odd numbers <= max missing: {miss}")
        # brute-force spot check of the table
        bad = 0
        ns = N if n <= 5 else 60
        for t in range(ns):
            m = t if n <= 5 else random.randrange(N)
            if brute_H(n, adjm_from_mask(n, m)) != int(H[m]):
                bad += 1
        print(f"   brute-force (permutations) check of {ns} table entries: mismatches {bad}")
    sys.stdout.flush()
    # ---- M6 (A), (B)
    for n in range(4, 8):
        H, S = tables[n]
        prs = pairs_of(n); na = len(prs)
        E = []; F = []; sh = []
        for e in range(na):
            for f in range(e + 1, na):
                E.append(e); F.append(f); sh.append(len(set(prs[e]) & set(prs[f])))
        out = pair_scan(H, S, np.array(E, np.int64), np.array(F, np.int64))
        sh = np.array(sh)
        adjm = out[sh == 1]; dis = out[sh == 0]
        print(f"n = {n} (A) adjacent: squares {int(adjm[:,0].sum())}, D = 0: {int(adjm[:,1].sum())}, min |D| = {int(adjm[:,2].min())}"
              f" | (B) disjoint: squares {int(dis[:,0].sum())}, D = 0: {int(dis[:,1].sum())}, both increments nonzero: "
              f"{int(dis[:,3].sum())}, of these all four strong: {int(dis[:,4].sum())}")
        sys.stdout.flush()
    # ---- M6 (C) affine 3-cubes on 3-matchings
    for n in (6, 7):
        H, S = tables[n]
        prs = pairs_of(n); na = len(prs)
        T = [c for c in itertools.combinations(range(na), 3)
             if len(set(prs[c[0]]) | set(prs[c[1]]) | set(prs[c[2]])) == 6]
        cap = 8192
        out, rec, nrec = cube3_scan(H, S, np.array(T, np.int64), cap)
        assert nrec.max() <= cap
        inc = Counter(); consec = 0; dist8 = 0; allstrong = 0
        ex = {}
        for p in range(len(T)):
            for k in range(int(nrec[p])):
                h0, de, df, dg, st = [int(x) for x in rec[p, k]]
                key = tuple(sorted((abs(de), abs(df), abs(dg))))
                inc[(key, st)] += 1
                vs = sorted(h0 + a * de + b * df + c * dg for a in (0, 1) for b in (0, 1) for c in (0, 1))
                if len(set(vs)) == 8:
                    dist8 += 1
                    if vs[-1] - vs[0] == 14:
                        consec += 1
                ex.setdefault(key, vs)
        print(f"n = {n} (C): 3-matchings {len(T)}; affine 3-cubes {int(out[:,0].sum())}; with three nonzero increments "
              f"{int(out[:,1].sum())}; of these all 8 strong {int(out[:,2].sum())}; 8 distinct values {dist8}; "
              f"8 consecutive odd {consec}")
        for key, cnt in sorted(inc.items()):
            print(f"      |increments| {key[0]} strong={key[1]}: {cnt}   example values {ex[key[0]]}")
        # also: affine 3-cubes on NON-matchings should not exist (Theorem A(a)): count them
        T2 = [c for c in itertools.combinations(range(na), 3)
              if len(set(prs[c[0]]) | set(prs[c[1]]) | set(prs[c[2]])) < 6]
        out2, rec2, nrec2 = cube3_scan(H, S, np.array(T2, np.int64), 16)
        print(f"      affine 3-cubes on arc triples that are NOT matchings ({len(T2)} triples): {int(out2[:,0].sum())}")
        sys.stdout.flush()
    # ---- M3: strong affine squares at n = 5 and 6 by value set; orbits of the (9,11,13,15) squares at n = 5
    for n in (5, 6):
        H, S = tables[n]
        prs = pairs_of(n); na = len(prs)
        N = 1 << na
        cnt = Counter(); cnt_any = Counter()
        squares = []
        for e in range(na):
            for f in range(e + 1, na):
                be = 1 << e; bf = 1 << f
                for m in range(N):
                    if m & (be | bf):
                        continue
                    a, b, c, d = int(H[m]), int(H[m | be]), int(H[m | bf]), int(H[m | be | bf])
                    if a - b - c + d != 0:
                        continue
                    vs = tuple(sorted(set((a, b, c, d))))
                    cnt_any[vs] += 1
                    if S[m] and S[m | be] and S[m | bf] and S[m | be | bf]:
                        cnt[vs] += 1
                        if n == 5 and vs == (9, 11, 13, 15):
                            squares.append((m, e, f))
        four = {k: v for k, v in cnt.items() if len(k) == 4}
        print(f"n = {n}: strong affine squares with 4 distinct values: {sorted(four.items())}")

        def covers(vs):
            s = set(vs); res = []
            for X in range(3, 200):
                odd = [v for v in range(X, 2 * X + 1) if v % 2 == 1]
                if all(v in s for v in odd):
                    res.append(X)
            return res
        print(f"   strong affine squares (any number of values) containing all odd numbers of some [X,2X], X>=3: "
              f"{[(k, v, covers(k)) for k, v in sorted(cnt.items()) if covers(k)]}")
        print(f"   (for comparison, not-necessarily-strong affine squares with that property: "
              f"{[(k, v, covers(k)) for k, v in sorted(cnt_any.items()) if covers(k)]})")
        if n == 5:
            # orbits under S_5 of (cube as set of 4 masks)
            def relabel(m, perm):
                A = adjm_from_mask(5, m)
                mm = 0
                for k, (i, j) in enumerate(prs):
                    # new tournament: perm[a] -> perm[b] iff a -> b
                    pass
                B = [[0] * 5 for _ in range(5)]
                for a in range(5):
                    for b in range(5):
                        if A[a][b]:
                            B[perm[a]][perm[b]] = 1
                for k, (i, j) in enumerate(prs):
                    if B[i][j]:
                        mm |= 1 << k
                return mm
            canon = Counter()
            for (m, e, f) in squares:
                members = [m, m | (1 << e), m | (1 << f), m | (1 << e) | (1 << f)]
                best = None
                for perm in itertools.permutations(range(5)):
                    key = tuple(sorted(relabel(x, perm) for x in members))
                    if best is None or key < best:
                        best = key
                canon[best] += 1
            print(f"   n = 5: labelled strong affine squares with values 9,11,13,15: {len(squares)}; isomorphism classes: "
                  f"{len(canon)} with sizes {sorted(canon.values())}")
            m, e, f = squares[0]
            print(f"   example: base mask {m}, arcs {prs[e]}, {prs[f]}, values "
                  f"{[int(H[m]), int(H[m|(1<<e)]), int(H[m|(1<<f)]), int(H[m|(1<<e)|(1<<f)])]}")
        sys.stdout.flush()


if __name__ == '__main__':
    main()
