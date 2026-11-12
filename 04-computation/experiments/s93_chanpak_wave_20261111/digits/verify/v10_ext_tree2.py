#!/usr/bin/env python3
"""Verifier (lane digits), check V10: the extension tree H -> 2H + 1, 2H - 1 when a step adds ONE or TWO vertices
with arbitrary arcs.  K_1, K_2 as in ext_tree2.py.  Own counting kernels (v6_ext_tree.ext_vals); children are
deduplicated with the lane's exact canonical form (tour.canon).

For every root T of order <= ROOTMAX in K_1: T is not in K_2 as soon as ALL children of one sign are outside K_1;
the sign with fewer children is examined first, the other one only if needed.

Usage: v10_ext_tree2.py ROOTMAX
"""
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, sys.path[0] + "/..")
from tour import canon, iso_classes  # noqa: E402

from v6_ext_tree import ext_vals, extend  # noqa: E402


@njit(cache=True)
def children(n, adj, sign, out1, out2, maxc):
    """labelled children with H = 2h + sign by one- or two-vertex steps.
    out1[i] = out-set of the first new vertex, out2[i] = out-set of the second one or -1.  Returns (h, count, has_other)."""
    buf = np.zeros(1 << n, np.int64)
    buf2 = np.zeros(1 << (n + 1), np.int64)
    a = np.zeros(n + 1, np.int64)
    h = ext_vals(n, adj, buf)
    tgt = 2 * h + sign
    oth = 2 * h - sign
    cnt = 0
    has_other = False
    for O in range(1 << n):
        if buf[O] == tgt:
            if cnt < maxc:
                out1[cnt] = O
                out2[cnt] = -1
            cnt += 1
        elif buf[O] == oth:
            has_other = True
    for O1 in range(1 << n):
        extend(n, adj, O1, a)
        ext_vals(n + 1, a, buf2)
        for O2 in range(1 << (n + 1)):
            if buf2[O2] == tgt:
                if cnt < maxc:
                    out1[cnt] = O1
                    out2[cnt] = O2
                cnt += 1
            elif buf2[O2] == oth:
                has_other = True
    return h, cnt, has_other


@njit(cache=True)
def in_K1_two(m, adj):
    buf = np.zeros(1 << m, np.int64)
    buf2 = np.zeros(1 << (m + 1), np.int64)
    a = np.zeros(m + 1, np.int64)
    h = ext_vals(m, adj, buf)
    p = False
    q = False
    for O in range(1 << m):
        if buf[O] == 2 * h + 1:
            p = True
        elif buf[O] == 2 * h - 1:
            q = True
    if p and q:
        return True
    for O1 in range(1 << m):
        extend(m, adj, O1, a)
        ext_vals(m + 1, a, buf2)
        for O2 in range(1 << (m + 1)):
            if buf2[O2] == 2 * h + 1:
                p = True
            elif buf2[O2] == 2 * h - 1:
                q = True
        if p and q:
            return True
    return False


def build(n, adj, O1, O2):
    a = np.zeros(n + 1, np.int64)
    extend(n, adj, O1, a)
    if O2 < 0:
        return n + 1, a
    b = np.zeros(n + 2, np.int64)
    extend(n + 1, a, O2, b)
    return n + 2, b


def main():
    rootmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    t0 = time.time()
    classes = iso_classes(rootmax)
    maxc = 200000
    o1 = np.zeros(maxc, np.int64)
    o2 = np.zeros(maxc, np.int64)
    memo = {}
    tested = 0

    def k1(m, a):
        nonlocal tested
        key = canon(m, a)
        if key not in memo:
            memo[key] = bool(in_K1_two(m, a))
            tested += 1
        return memo[key]

    for n in range(1, rootmax + 1):
        nk1 = 0
        nk2 = 0
        for key, adj in classes[n].items():
            adj = np.array(adj, np.int64)
            lists = {}
            for sign in (1, -1):
                h, cnt, has_other = children(n, adj, sign, o1, o2, maxc)
                assert cnt <= maxc
                lists[sign] = [(int(o1[i]), int(o2[i])) for i in range(cnt)]
            if not lists[1] or not lists[-1]:
                continue
            nk1 += 1
            order = sorted((1, -1), key=lambda s: len(lists[s]))
            ok = {}
            for s in order:
                found = False
                seen = set()
                for (O1, O2) in lists[s]:
                    m, a = build(n, adj, O1, O2)
                    ck = canon(m, a)
                    if ck in seen:
                        continue
                    seen.add(ck)
                    if k1(m, a):
                        found = True
                        break
                ok[s] = found
                if not found:
                    break
            if all(ok.get(s, False) for s in (1, -1)):
                nk2 += 1
                print("      IN K_2:", n, [int(x) for x in adj])
        print(f"   order {n}: classes {len(classes[n])}, in K_1 (one- or two-vertex steps): {nk1}, in K_2: {nk2}; "
              f"distinct children tested so far {tested}   [{time.time() - t0:.1f} s]", flush=True)
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
