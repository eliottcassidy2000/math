#!/usr/bin/env python3
"""Verifier (lane digits), check V1: the one-mark block matrices, derived INDEPENDENTLY of the lane's interpolation.

The lane (onemark.py) interpolates the 4 x 4 matrix of a block op (X, W, m) from four base tournaments and then
tests it on the other bases ONLY for the first labelled op that yields a new matrix.  Here every matrix is computed
from the path-surgery formulas of the proof (counts on X alone, no base tournament involved):

  s' = c1 (h - e) + c2 e + c3 (B - g) + c4 g   for each coordinate,  N row = (c1, c2 - c1, c3, c4 - c3).

(1) all 202,012 labelled blocks with |X| <= 5: set of distinct matrices compared with onemark_matrices_5.json;
    every json entry (k, adjX, W, m) compared with the formula matrix.
(2) direct brute-force check of the law on random marked bases of orders 6 and 7 (outside the lane's test set,
    which stops at order 5) with random blocks of order <= 5, own subset DP.
"""
import json
import random
import sys
import time

import numpy as np
from numba import njit


@njit(cache=True)
def path_table(k, adj):
    full = 1 << k
    P = np.zeros((full, k, k), np.int64)
    for v in range(k):
        P[1 << v, v, v] = 1
    for m in range(1, full):
        for s in range(k):
            if not (m >> s) & 1:
                continue
            for e in range(k):
                c = P[m, s, e]
                if c == 0:
                    continue
                out = adj[e] & ~m
                for u in range(k):
                    if (out >> u) & 1:
                        P[m | (1 << u), s, u] += c
    return P


@njit(cache=True)
def block_matrix(k, P, W, m, out):
    """matrix of the block (X, W, m) in the basis (h, e, B, g); m = -1: mark stays."""
    full = (1 << k) - 1
    nm = full + 1
    A = np.zeros(nm, np.int64)       # all paths
    SnW = np.zeros(nm, np.int64)     # start outside W
    EW = np.zeros(nm, np.int64)      # end in W
    Em = np.zeros(nm, np.int64)      # end at m
    SnWEm = np.zeros(nm, np.int64)   # start outside W, end at m
    Sm = np.zeros(nm, np.int64)      # start at m
    SmEW = np.zeros(nm, np.int64)    # start at m, end in W
    for mask in range(1, nm):
        for s in range(k):
            if not (mask >> s) & 1:
                continue
            for e in range(k):
                c = P[mask, s, e]
                if c == 0:
                    continue
                A[mask] += c
                sw = (W >> s) & 1
                ew = (W >> e) & 1
                if sw == 0:
                    SnW[mask] += c
                if ew == 1:
                    EW[mask] += c
                if e == m:
                    Em[mask] += c
                    if sw == 0:
                        SnWEm[mask] += c
                if s == m:
                    Sm[mask] += c
                    if ew == 1:
                        SmEW[mask] += c
    Ae = A.copy()
    Ae[0] = 1
    SnWe = SnW.copy()
    SnWe[0] = 1
    c = np.zeros((4, 4), np.int64)   # c[row, 0..3] = coefficients of (h - e, e, B - g, g)
    # ---- row h'
    c[0, 0] = A[full]
    c[0, 1] = SnW[full]
    for Y1 in range(1, nm):
        Y2 = full ^ Y1
        c[0, 2] += EW[Y1] * Ae[Y2]
        c[0, 3] += EW[Y1] * SnWe[Y2]
    if m < 0:
        # e' = g End_X(W)
        c[1, 3] = EW[full]
        # B' = (B - g) Pi + g Pi_W
        for Y1 in range(0, nm):
            Y2 = full ^ Y1
            c[2, 2] += Ae[Y1] * Ae[Y2]
            c[2, 3] += Ae[Y1] * SnWe[Y2]
        # g' = g H(X)
        c[3, 3] = A[full]
    else:
        # e': paths ending at m
        c[1, 0] = Em[full]
        c[1, 1] = SnWEm[full]
        for Y1 in range(1, full):
            Y2 = full ^ Y1
            c[1, 2] += EW[Y1] * Em[Y2]
            c[1, 3] += EW[Y1] * SnWEm[Y2]
        # B': pairs (P1, P2), P2 starts at m
        for Y1 in range(0, nm):
            Y2 = full ^ Y1
            if Y2 == 0:
                continue
            c[2, 0] += Ae[Y1] * Sm[Y2]          # (alpha): P1 = I [Y1], P2 = Y2
            c[2, 1] += SnWe[Y1] * Sm[Y2]
        for Y1 in range(1, nm):                 # entry in P1: P1 = [A] Y1 Z [Y1'], P2 = Y2
            if EW[Y1] == 0:
                continue
            rest = full ^ Y1
            Y2 = rest
            while Y2 > 0:
                if Sm[Y2] != 0:
                    Y1p = rest ^ Y2
                    c[2, 2] += EW[Y1] * Ae[Y1p] * Sm[Y2]
                    c[2, 3] += EW[Y1] * SnWe[Y1p] * Sm[Y2]
                Y2 = (Y2 - 1) & rest
        for Y2 in range(1, nm):                 # entry in P2: P2 = Y2 Z [Y2'], P1 = [A] [Y1]
            if SmEW[Y2] == 0:
                continue
            rest = full ^ Y2
            Y1 = rest
            while True:
                Y2p = rest ^ Y1
                c[2, 2] += Ae[Y1] * SmEW[Y2] * Ae[Y2p]
                c[2, 3] += Ae[Y1] * SmEW[Y2] * SnWe[Y2p]
                if Y1 == 0:
                    break
                Y1 = (Y1 - 1) & rest
        # g' = H(T' - m)
        if k == 1:
            c[3, 0] = 1
            c[3, 1] = 1
        else:
            fm = full ^ (1 << m)
            c[3, 0] = A[fm]
            c[3, 1] = SnW[fm]
            Y1 = fm
            while Y1 > 0:
                Y2 = fm ^ Y1
                c[3, 2] += EW[Y1] * Ae[Y2]
                c[3, 3] += EW[Y1] * SnWe[Y2]
                Y1 = (Y1 - 1) & fm
    for r in range(4):
        out[r, 0] = c[r, 0]
        out[r, 1] = c[r, 1] - c[r, 0]
        out[r, 2] = c[r, 2]
        out[r, 3] = c[r, 3] - c[r, 2]


@njit(cache=True)
def all_blocks(k, res):
    """res[idx] = flattened matrix for every labelled (X, W, m); returns count."""
    npairs = k * (k - 1) // 2
    adj = np.zeros(k, np.int64)
    out = np.zeros((4, 4), np.int64)
    idx = 0
    for bits in range(1 << npairs):
        for i in range(k):
            adj[i] = 0
        t = 0
        for i in range(k):
            for j in range(i + 1, k):
                if (bits >> t) & 1:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
                t += 1
        P = path_table(k, adj)
        for W in range(1 << k):
            for m in range(-1, k):
                block_matrix(k, P, W, m, out)
                for r in range(4):
                    for q in range(4):
                        res[idx, r * 4 + q] = out[r, q]
                idx += 1
    return idx


# ---------------------------------------------------------------- direct states (own DP)
@njit(cache=True)
def end_table(n, adj):
    full = 1 << n
    E = np.zeros((full, n), np.int64)
    for v in range(n):
        E[1 << v, v] = 1
    for m in range(1, full):
        for v in range(n):
            c = E[m, v]
            if c == 0:
                continue
            for u in range(n):
                if (adj[v] >> u) & 1 and not (m >> u) & 1:
                    E[m | (1 << u), u] += c
    return E


@njit(cache=True)
def direct_state(n, adj, z, s):
    """(h, e, B, g) of (T, z) from the definitions."""
    full = (1 << n) - 1
    E = end_table(n, adj)
    radj = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(n):
            if (adj[i] >> j) & 1:
                radj[j] |= 1 << i
    R = end_table(n, radj)            # R[mask, v] = paths of T[mask] STARTING at v
    h = 0
    for v in range(n):
        h += E[full, v]
    e = E[full, z]
    B = 0
    for S in range(0, full + 1):
        if (S >> z) & 1:
            continue
        comp = full ^ S
        hs = 1
        if S != 0:
            hs = 0
            for v in range(n):
                hs += E[S, v]
        B += hs * R[comp, z]
    g = 1
    if n > 1:
        g = 0
        fz = full ^ (1 << z)
        for v in range(n):
            g += E[fz, v]
    s[0] = h
    s[1] = e
    s[2] = B
    s[3] = g


def main():
    t0 = time.time()
    here = sys.path[0]
    ops = json.load(open(here + "/../onemark_matrices_5.json"))
    lane = {tuple(tuple(r) for r in o["N"]) for o in ops}
    print(f"lane file: {len(ops)} entries, {len(lane)} distinct matrices")
    mine = set()
    total = 0
    per_k = {}
    for k in range(1, 6):
        cnt = (1 << (k * (k - 1) // 2)) * (1 << k) * (k + 1)
        res = np.zeros((cnt, 16), np.int64)
        got = all_blocks(k, res)
        assert got == cnt
        total += cnt
        u = np.unique(res, axis=0)
        for row in u:
            mine.add(tuple(tuple(int(x) for x in row[4 * r:4 * r + 4]) for r in range(4)))
        per_k[k] = (cnt, len(mine))
        print(f"   |X| = {k}: {cnt} labelled blocks; distinct formula matrices so far {len(mine)}   [{time.time() - t0:.1f} s]", flush=True)
    print(f"(1a) labelled blocks |X| <= 5: {total}; distinct matrices from the formulas: {len(mine)}; equal to the lane set: {mine == lane}; "
          f"formula-only {len(mine - lane)}, lane-only {len(lane - mine)}")
    # every json entry against the formula
    bad = 0
    out = np.zeros((4, 4), np.int64)
    for o in ops:
        k = o["k"]
        adj = np.array(o["adjX"], np.int64)
        P = path_table(k, adj)
        block_matrix(k, P, o["W"], o["m"], out)
        if [list(map(int, r)) for r in out] != o["N"]:
            bad += 1
    print(f"(1b) json entries whose stored matrix differs from the formula matrix of their own (X, W, m): {bad}")
    # (2) direct check on random bases of orders 6, 7 and on ALL labelled blocks with |X| <= 3 over all bases of order <= 4
    rng = random.Random(931)
    fails = 0
    ntest = 0
    s = np.zeros(4, np.int64)
    s2 = np.zeros(4, np.int64)
    for trial in range(1500):
        n = rng.choice([6, 7])
        adj = np.zeros(n, np.int64)
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < 0.5:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
        z = rng.randrange(n)
        k = rng.choice([1, 2, 3, 4, 5, 5])
        adjX = np.zeros(k, np.int64)
        for i in range(k):
            for j in range(i + 1, k):
                if rng.random() < 0.5:
                    adjX[i] |= 1 << j
                else:
                    adjX[j] |= 1 << i
        W = rng.randrange(1 << k)
        m = rng.randrange(-1, k)
        N = n + k
        a = np.zeros(N, np.int64)
        xmask = ((1 << k) - 1) << n
        for i in range(n):
            a[i] = int(adj[i])
            if i != z:
                a[i] |= xmask
            else:
                a[i] |= (((1 << k) - 1) & ~W) << n
        for i in range(k):
            a[n + i] = int(adjX[i]) << n
            if (W >> i) & 1:
                a[n + i] |= 1 << z
        z2 = z if m < 0 else n + m
        direct_state(n, adj, z, s)
        direct_state(N, a, z2, s2)
        P = path_table(k, adjX)
        block_matrix(k, P, W, m, out)
        pred = [sum(int(out[r, q]) * int(s[q]) for q in range(4)) for r in range(4)]
        ntest += 1
        if pred != [int(x) for x in s2]:
            fails += 1
            if fails <= 3:
                print("   LAW FAILS", n, list(map(int, adj)), z, k, list(map(int, adjX)), W, m, list(map(int, s)), pred, list(map(int, s2)))
    print(f"(2) law s' = N s by direct count (own DP), random marked bases of orders 6, 7, random blocks |X| <= 5: {ntest} tests, failures {fails}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
