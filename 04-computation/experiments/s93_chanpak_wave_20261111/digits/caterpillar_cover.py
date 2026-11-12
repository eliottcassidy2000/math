#!/usr/bin/env python3
"""Lane "digits", part (3): value set of the one-mark one-vertex system (words in a, b, c over a point).

State (h, e, B, g) of a marked tournament (T, z); start (1, 1, 1, 1) (T = {z}).
    a: new v -> z only, mark stays      (h, e, B, g) -> (h - e + B, g, 2B - g, g)
    b: new v -> z only, mark := v       (h, e, B, g) -> (h - e + B, h - e, B + h, h)
    c: new v a sink, mark stays         (h, e, B, g) -> (h, 0, 2B, g)
(d: new sink, mark := v gives (h, h, h, h): the product with a fresh copy; omitted, or included with flag D = 1.)
Every word is an actual tournament ("caterpillar": each new vertex beats at most the current mark).
h never decreases and a state with B > LIM can never return below LIM, so the search below is exhaustive for
values <= LIM: depth-first over words, pruning h > LIM or B > LIM.

Output: the odd numbers <= LIM that are not attained; coverage by dyadic window; multiplicity statistics.

Usage: caterpillar_cover.py LIM [D]
"""
import sys
import time

import numpy as np
from numba import njit, prange


@njit(cache=True)
def dfs(h0, e0, B0, g0, lim, mult, use_d):
    """iterative DFS from one state; mult[(h-1)/2] += number of words reaching h (saturating)."""
    maxd = 400
    sh = np.zeros(maxd, np.int64)
    se = np.zeros(maxd, np.int64)
    sB = np.zeros(maxd, np.int64)
    sg = np.zeros(maxd, np.int64)
    so = np.zeros(maxd, np.int64)
    nodes = 0
    d = 0
    sh[0] = h0
    se[0] = e0
    sB[0] = B0
    sg[0] = g0
    so[0] = 0
    nops = 4 if use_d else 3
    # the start state itself is counted by the caller
    while d >= 0:
        o = so[d]
        if o == nops:
            d -= 1
            continue
        so[d] += 1
        h = sh[d]
        e = se[d]
        B = sB[d]
        g = sg[d]
        if o == 0:
            h2 = h - e + B
            e2 = g
            B2 = 2 * B - g
            g2 = g
            if h2 == h and e2 == e and B2 == B:
                continue          # z is a sink: the new vertex is a source-like padding, same state
        elif o == 1:
            h2 = h - e + B
            e2 = h - e
            B2 = B + h
            g2 = h
        elif o == 2:
            h2 = h
            e2 = 0
            B2 = 2 * B
            g2 = g
            if B2 > lim:
                continue          # same value h, already counted at the parent
        else:
            # d: reset; only useful once per value; the state (h, h, h, h)
            if e == h and B == h and g == h:
                continue
            h2 = h
            e2 = h
            B2 = h
            g2 = h
        if h2 > lim:
            continue
        nodes += 1
        if mult[(h2 - 1) >> 1] < 60000:
            mult[(h2 - 1) >> 1] += 1
        if B2 > lim:
            continue              # the word is counted; no descendant has a new value <= LIM
        if d + 2 >= maxd:
            continue              # never reached for LIM <= 10^9 (depth is O(log LIM)); guard only
        d += 1
        sh[d] = h2
        se[d] = e2
        sB[d] = B2
        sg[d] = g2
        so[d] = 0
    return nodes


@njit(cache=True)
def step(o, s):
    h, e, B, g = s
    if o == 0:
        return (h - e + B, g, 2 * B - g, g)
    if o == 1:
        return (h - e + B, h - e, B + h, h)
    return (h, 0, 2 * B, g)


@njit(parallel=True, cache=True)
def run(starts, lim, use_d, nthreads):
    n = starts.shape[0]
    half = (lim + 1) // 2
    mults = np.zeros((nthreads, half), np.uint16)
    nodes = np.zeros(nthreads, np.int64)
    for t in prange(nthreads):
        for i in range(t, n, nthreads):
            nodes[t] += dfs(starts[i, 0], starts[i, 1], starts[i, 2], starts[i, 3], lim, mults[t], use_d)
    return mults, nodes


def main():
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    use_d = bool(int(sys.argv[2])) if len(sys.argv) > 2 else False
    t0 = time.time()
    # prefixes of length PRE (no d inside the prefix), then a parallel DFS from each
    pre = 7
    states = [(1, 1, 1, 1)]
    seen_prefix_values = []
    frontier = [(1, 1, 1, 1)]
    for ln in range(pre):
        nxt = []
        for s in frontier:
            for o in range(3):
                h, e, B, g = s
                if o == 0:
                    s2 = (h - e + B, g, 2 * B - g, g)
                    if s2 == s:
                        continue
                elif o == 1:
                    s2 = (h - e + B, h - e, B + h, h)
                else:
                    s2 = (h, 0, 2 * B, g)
                if s2[0] <= lim:
                    seen_prefix_values.append(s2[0])
                    if s2[2] <= lim:
                        nxt.append(s2)
        frontier = nxt
    starts = np.array(frontier, dtype=np.int64)
    nthreads = 12
    mults, nodes = run(starts, lim, use_d, nthreads)
    mult = mults.astype(np.int64).sum(axis=0)
    for v in [1] + seen_prefix_values:
        mult[(v - 1) >> 1] += 1
    total_nodes = int(nodes.sum())
    half = (lim + 1) // 2
    missing = [2 * i + 1 for i in range(half) if mult[i] == 0]
    print(f"one-mark words in a, b, c{', d' if use_d else ''} over a point: values <= {lim}; word-nodes visited {total_nodes}  [{time.time() - t0:.1f} s]")
    print(f"odd numbers <= {lim} not attained: {len(missing)}")
    print("   first 80:", missing[:80])
    print("   last 20:", missing[-20:])
    k = 3
    while (1 << (k + 1)) <= lim + 1:
        lo, hi = (1 << k) >> 1, (1 << (k + 1)) >> 1
        seg = mult[lo:hi]
        print(f"   odd N in [2^{k}, 2^{k + 1}): attained {int((seg > 0).sum())} of {hi - lo} ({(seg > 0).mean():.4f}); min mult {int(seg.min())}, median {int(np.median(seg))}")
        k += 1


if __name__ == "__main__":
    main()
