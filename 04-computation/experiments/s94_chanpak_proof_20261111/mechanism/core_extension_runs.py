#!/usr/bin/env python3
"""One-vertex extension spectrum of a single core (exhaustive, exact), core order c = 3..6.
For every labelled core C (order c) the set V(C) = {H(C + v_B) : B subset of V(C)} (2^c members) is read from the
table of order c+1.  Reported: the longest run of consecutive odd numbers in V(C), the largest A such that V(C)
contains ALL odd numbers of [A-1, 2(A-1)] (a dyadic interval), and the same restricted to strong extensions."""
import numpy as np
from arc_cube import table
for n in (4, 5, 6, 7):
    c = n - 1
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pi = np.array([p[0] for p in pairs], np.int64); pj = np.array([p[1] for p in pairs], np.int64)
    H, S = table(n, pi, pj)
    core_idx = [k for k, p in enumerate(pairs) if p[1] != n - 1]
    cut_idx = [k for k, p in enumerate(pairs) if p[1] == n - 1]
    ncore = 1 << len(core_idx); ncut = 1 << len(cut_idx)
    cm = np.zeros(ncore, np.int64)
    for t, k in enumerate(core_idx):
        cm |= ((np.arange(ncore) >> t) & 1) << k
    bm = np.zeros(ncut, np.int64)
    for t, k in enumerate(cut_idx):
        bm |= ((np.arange(ncut) >> t) & 1) << k
    M = cm[:, None] | bm[None, :]
    V = H[M]            # ncore x ncut
    St = S[M]
    best_run = (0, None); best_dy = (0, None); best_dy_strong = (0, None)
    ndy = 0
    for r in range(ncore):
        for strong_only in (False, True):
            vals = np.unique(V[r][St[r] == 1]) if strong_only else np.unique(V[r])
            if len(vals) == 0: continue
            # runs of consecutive odd numbers
            brk = np.where(np.diff(vals) != 2)[0]
            starts = np.concatenate(([0], brk + 1)); ends = np.concatenate((brk, [len(vals) - 1]))
            for s, e in zip(starts, ends):
                L = e - s + 1; A = int(vals[s])
                if not strong_only and L > best_run[0]: best_run = (int(L), (r, A))
                # all odd numbers of [A'-1, 2(A'-1)] for some A' in the run: need A' + 2(L'-1) >= 2A' - 3
                # best choice inside the run: largest A' with 2A' - 3 <= last
                last = int(vals[e])
                Ap = (last + 3) // 2
                if Ap % 2 == 0: Ap -= 1
                if Ap >= A and Ap >= 3:
                    if strong_only:
                        if Ap > best_dy_strong[0]: best_dy_strong = (Ap, (r, A, last))
                    else:
                        if Ap > best_dy[0]: best_dy = (Ap, (r, A, last))
    print(f"core order {c} ({ncore} labelled cores, {ncut} cuts each): max H of extensions {int(V.max())}")
    print(f"   longest run of consecutive odd values in one extension spectrum: {best_run[0]} (core mask, start) = {best_run[1]}")
    print(f"   largest A with all odd numbers of [A-1, 2A-2] in ONE extension spectrum: A = {best_dy[0]}  (core, run start, run end) = {best_dy[1]}")
    print(f"   same with strong extensions only: A = {best_dy_strong[0]}  {best_dy_strong[1]}")
