#!/usr/bin/env python3
"""Verifier: cut law (M1), Lemma 1, and single-core extension spectra (M7), own code.
For each labelled core C of order c: f, g, Q by DP; vals[B] by the DIRECT decomposition (P1, v, P2) (first form);
second form H + sum_{x in B}(U'(x) - s(x)) - sum_{pairs of B} W_xy compared; for c <= 6 also compared with the table
of order c + 1; Lemma 1 W_xy >= 2 H(C - x - y); longest run of consecutive odd values; largest X (ANY integer X >= 2)
such that all odd numbers of [X, 2X] are in the spectrum.   usage: python v2_cores.py c"""
import sys, time, random
import numpy as np
from numba import njit, prange
from vlib import *

NB = 400


@njit(cache=True)
def core_one(c, m, pi, pj, Hnext, pidx, use_table, res):
    adj = np.zeros(c, np.int64); radj = np.zeros(c, np.int64)
    adj_from_mask(c, m, pi, pj, adj, radj)
    f = np.zeros((1 << c, c), np.int64); g = np.zeros((1 << c, c), np.int64); Q = np.zeros((c, c), np.int64)
    h = fgQ(c, adj, f, g, Q)
    full = (1 << c) - 1
    # base mask of the (c+1)-tournament: core arcs
    M0 = 0
    if use_table:
        for i in range(c):
            for j in range(i + 1, c):
                if (adj[i] >> j) & 1:
                    M0 |= 1 << pidx[i, j]
    U = np.zeros(c, np.int64)
    for x in range(c):
        t = f[full, x]
        for y in range(c):
            if y != x:
                t += Q[x, y]
        U[x] = t
    present = np.zeros(NB, np.uint8)
    bad_form = 0; bad_table = 0
    for B in range(1 << c):
        t = 0
        for x in range(c):
            if (B >> x) & 1:
                t += f[full, x]
                for y in range(c):
                    if not (B >> y) & 1:
                        t += Q[x, y]
            else:
                t += g[full, x]
        # second form
        t2 = h
        for x in range(c):
            if (B >> x) & 1:
                t2 += U[x] - g[full, x]
                for y in range(x + 1, c):
                    if (B >> y) & 1:
                        t2 -= Q[x, y] + Q[y, x]
        if t2 != t:
            bad_form += 1
        if use_table:
            M = M0
            for x in range(c):
                if (B >> x) & 1:
                    M |= 1 << pidx[x, c]
            if Hnext[M] != t:
                bad_table += 1
        present[(t - 1) // 2] = 1
    # Lemma 1
    lem_viol = 0; lem_eq = 0; npairs = 0
    for x in range(c):
        for y in range(x + 1, c):
            R = full ^ (1 << x) ^ (1 << y)
            hr = 1
            if R != 0:
                hr = 0
                for z in range(c):
                    hr += f[R, z]
            w = Q[x, y] + Q[y, x]
            npairs += 1
            if w < 2 * hr:
                lem_viol += 1
            if w == 2 * hr:
                lem_eq += 1
    # longest run
    best = 0; bstart = 0; cur = 0
    for b in range(NB):
        if present[b]:
            cur += 1
            if cur > best:
                best = cur; bstart = 2 * (b - cur + 1) + 1
        else:
            cur = 0
    # largest X with all odd numbers of [X, 2X] present; prefix sums
    pre = np.zeros(NB + 1, np.int64)
    for b in range(NB):
        pre[b + 1] = pre[b] + present[b]
    bestX = 0
    for X in range(2, NB):
        # odd v in [X, 2X]: bits from lo = X // 2 (v = X + 1 if X even, v = X if X odd) to X - 1
        lo = X // 2
        if pre[X] - pre[lo] == X - lo:
            bestX = X
    res[0] = h; res[1] = best; res[2] = bstart; res[3] = bestX; res[4] = bad_form; res[5] = bad_table
    res[6] = lem_viol; res[7] = lem_eq; res[8] = npairs


@njit(parallel=True, cache=True)
def scan(c, pi, pj, Hnext, pidx, use_table):
    na = c * (c - 1) // 2
    M = 1 << na
    R = np.zeros((M, 9), np.int64)
    for m in prange(M):
        core_one(c, m, pi, pj, Hnext, pidx, use_table, R[m])
    return R


def main():
    c = int(sys.argv[1])
    pi, pj = pair_arrays(c)
    use_table = c <= 6
    n1 = c + 1
    pidx = np.zeros((n1, n1), np.int64)
    for k, (i, j) in enumerate(pairs_of(n1)):
        pidx[i, j] = k; pidx[j, i] = k
    if use_table:
        p1, p2 = pair_arrays(n1)
        Hnext, _ = build_table(n1, p1, p2)
    else:
        Hnext = np.zeros(1, np.int32)
    t0 = time.time()
    R = scan(c, pi, pj, Hnext, pidx, use_table)
    print(f"core order {c}: {R.shape[0]} labelled cores x {1 << c} cuts ({time.time()-t0:.1f}s)")
    print(f"   second form vs direct form mismatches: {int(R[:,4].sum())}; "
          f"direct form vs table of order {c+1} mismatches: {int(R[:,5].sum()) if use_table else 'n/a (sampled below)'}")
    print(f"   Lemma 1 violations: {int(R[:,6].sum())}; equality in {int(R[:,7].sum())} of {int(R[:,8].sum())} pairs")
    b = int(R[:, 1].max()); idx = np.where(R[:, 1] == b)[0]
    print(f"   longest run of consecutive odd values in one spectrum: {b}; attained by {len(idx)} labelled cores; core H "
          f"{sorted(set(R[idx,0].tolist()))}; run starts {sorted(set(R[idx,2].tolist()))}; example core mask {int(idx[0])}")
    bx = int(R[:, 3].max()); idx2 = np.where(R[:, 3] == bx)[0]
    print(f"   largest X (any integer >= 2) with all odd numbers of [X, 2X] in one spectrum: {bx}; attained by "
          f"{len(idx2)} cores; core H {sorted(set(R[idx2,0].tolist()))}")
    # distribution of best X
    from collections import Counter
    print(f"   distribution of best X over cores: {sorted(Counter(R[:,3].tolist()).items())}")
    if not use_table:
        # direct check with the 8-vertex DP: the example core (all cuts) and random (core, cut) samples
        from numba import njit as _nj
        random.seed(7)
        adj = np.zeros(n1, np.int64); radj = np.zeros(n1, np.int64); dp = np.zeros((1 << n1, n1), np.int64)
        res = np.zeros(9, np.int64)

        def direct(mcore, B):
            A = adjm_from_mask(c, mcore)
            a = np.zeros(n1, np.int64)
            for i in range(c):
                for j in range(c):
                    if A[i][j]:
                        a[i] |= 1 << j
            for x in range(c):
                if (B >> x) & 1:
                    a[x] |= 1 << c
                else:
                    a[c] |= 1 << x
            return int(hp(n1, a, dp))
        ex = int(idx[0])
        spec = sorted(set(direct(ex, B) for B in range(1 << c)))
        runs = []; cur = [spec[0]]
        for v in spec[1:]:
            if v == cur[-1] + 2:
                cur.append(v)
            else:
                runs.append(cur); cur = [v]
        runs.append(cur)
        lr = max(runs, key=len)
        print(f"   direct 8-vertex DP on the example core {ex}: spectrum has {len(spec)} values, longest run "
              f"{len(lr)}: {lr[0]}..{lr[-1]}")
        # random samples: cut-law value (recomputed for single core) vs direct
        bad = 0; ns = 3000
        f = np.zeros((1 << c, c), np.int64); g = np.zeros((1 << c, c), np.int64); Q = np.zeros((c, c), np.int64)
        for t in range(ns):
            mcore = random.randrange(1 << (c * (c - 1) // 2)); B = random.randrange(1 << c)
            a = np.zeros(c, np.int64); ra = np.zeros(c, np.int64)
            adj_from_mask(c, mcore, pi, pj, a, ra)
            fgQ(c, a, f, g, Q)
            full = (1 << c) - 1
            val = 0
            for x in range(c):
                if (B >> x) & 1:
                    val += int(f[full, x]) + sum(int(Q[x, y]) for y in range(c) if not (B >> y) & 1)
                else:
                    val += int(g[full, x])
            if val != direct(mcore, B):
                bad += 1
        print(f"   {ns} random (core, cut) samples, cut-law value vs direct 8-vertex DP: mismatches {bad}")


if __name__ == '__main__':
    main()
