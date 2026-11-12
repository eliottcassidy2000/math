#!/usr/bin/env python3
"""Longest run of consecutive odd numbers in the one-vertex extension spectrum of a single core, via the cut law
   H(C + v_B) = sum_{y not in B} s(y) + sum_{x in B} e(x) + sum_{x in B, y not in B} Q(x, y)        (PROOFS.txt (0.1))
(exact int64; exhaustive over all labelled cores of order c).  Cross-check: c <= 6 must reproduce
core_extension_runs.out (table method).  usage: python core_runs_cutlaw.py c"""
import sys, time
import numpy as np
from numba import njit, prange

@njit(cache=True)
def core_stats(c, m, pi, pj, vals):
    adj = np.zeros(c, np.int64); radj = np.zeros(c, np.int64)
    na = c * (c - 1) // 2
    for k in range(na):
        if (m >> k) & 1:
            a = pi[k]; b = pj[k]
        else:
            a = pj[k]; b = pi[k]
        adj[a] |= 1 << b
        radj[b] |= 1 << a
    full = (1 << c) - 1
    f = np.zeros((1 << c, c), np.int64)   # paths covering S ending at x
    g = np.zeros((1 << c, c), np.int64)   # paths covering S starting at y
    for v in range(c):
        f[1 << v, v] = 1
        g[1 << v, v] = 1
    for S in range(1, 1 << c):
        for last in range(c):
            cf = f[S, last]
            if cf != 0:
                nxt = adj[last] & ~S
                for w in range(c):
                    if (nxt >> w) & 1:
                        f[S | (1 << w), w] += cf
            cg = g[S, last]
            if cg != 0:
                nxt = radj[last] & ~S
                for w in range(c):
                    if (nxt >> w) & 1:
                        g[S | (1 << w), w] += cg
    Q = np.zeros((c, c), np.int64)
    for S in range(1, full):
        R = full ^ S
        for x in range(c):
            fx = f[S, x]
            if fx == 0:
                continue
            for y in range(c):
                if (R >> y) & 1:
                    Q[x, y] += fx * g[R, y]
    h = 0
    for x in range(c):
        h += f[full, x]
    nv = 1 << c
    for B in range(nv):
        t = 0
        for x in range(c):
            if (B >> x) & 1:
                t += f[full, x]
                for y in range(c):
                    if not (B >> y) & 1:
                        t += Q[x, y]
            else:
                t += g[full, x]
        vals[B] = t
    vs = np.sort(vals[:nv])
    # longest run of consecutive odd numbers (distinct values)
    bestL = 1; bestA = vs[0]
    curL = 1; curA = vs[0]
    bestdy = 0
    for i in range(1, nv):
        if vs[i] == vs[i - 1]:
            continue
        if vs[i] == vs[i - 1] + 2:
            curL += 1
        else:
            curL = 1; curA = vs[i]
        if curL > bestL:
            bestL = curL; bestA = curA
        # dyadic: all odd numbers of [A'-1, 2A'-2] for some A' >= curA in this run:  2A' - 3 <= vs[i]
        Ap = (vs[i] + 3) // 2
        if Ap % 2 == 0:
            Ap -= 1
        if Ap >= curA and Ap > bestdy:
            bestdy = Ap
    return h, bestL, bestA, bestdy

@njit(parallel=True, cache=True)
def scan(c, pi, pj):
    na = c * (c - 1) // 2
    M = 1 << na
    L = np.zeros(M, np.int32); A = np.zeros(M, np.int64); Hc = np.zeros(M, np.int64); DY = np.zeros(M, np.int64)
    for m in prange(M):
        vals = np.zeros(1 << c, np.int64)
        h, l, a, dy = core_stats(c, m, pi, pj, vals)
        L[m] = l; A[m] = a; Hc[m] = h; DY[m] = dy
    return Hc, L, A, DY

def main():
    c = int(sys.argv[1])
    pairs = [(i, j) for i in range(c) for j in range(i + 1, c)]
    pi = np.array([p[0] for p in pairs], np.int64); pj = np.array([p[1] for p in pairs], np.int64)
    t0 = time.time()
    Hc, L, A, DY = scan(c, pi, pj)
    b = int(L.max()); idx = np.where(L == b)[0]
    print(f"core order {c}: {len(L)} labelled cores, {time.time()-t0:.1f}s")
    print(f"   longest run of consecutive odd values in one extension spectrum: {b}; attained by {len(idx)} labelled cores; "
          f"core H values {sorted(set(Hc[idx].tolist()))}; run starts {sorted(set(A[idx].tolist()))}; example mask {int(idx[0])}")
    print(f"   largest A with all odd numbers of [A-1, 2A-2] in one spectrum: {int(DY.max())}")
    # best run by core H
    print("   max run length by core H (H: run, start):", end=" ")
    for hv in sorted(set(Hc.tolist())):
        sel = Hc == hv
        j = np.argmax(L[sel]); print(f"{hv}: {int(L[sel][j])}@{int(A[sel][j])}", end="; ")
    print()

if __name__ == '__main__':
    main()
