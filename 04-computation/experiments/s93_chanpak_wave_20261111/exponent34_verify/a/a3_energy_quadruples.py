#!/usr/bin/env python3
"""Independent verifier, part (A)3 of THM-5141 (draft): additive energy of G_w = {(x, rev_w x) : 0 <= x < 2^w}.

E1  literal enumeration of all 2^(4w) quadruples (x1, x2, x3, x4), counting
       sum form:        p(x1) + p(x2) = p(x3) + p(x4)
       difference form: p(x1) - p(x2) = p(x3) - p(x4)          (p(x) = (x, rev_w x)),  w <= W4.
E2  three free points, the fourth forced (x4 = x3 - x1 + x2 must be in range and satisfy the reversed equation),
    w <= W3 (2^(3w) triples, 4 threads).
Both must equal 6^w.  No signed strings, no difference tables: only the definition of additive energy.
Usage: python a3_energy_quadruples.py [W4 (default 7)] [W3 (default 11)]
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
import sys
import time
import numpy as np
from numba import njit, prange


@njit(cache=False)
def rev_table(w):
    n = 1 << w
    rev = np.zeros(n, np.int64)
    for x in range(n):
        r = 0
        for i in range(w):
            if (x >> i) & 1:
                r |= 1 << (w - 1 - i)
        rev[x] = r
    return rev


@njit(parallel=True, cache=False)
def energy4(w):
    n = 1 << w
    rev = rev_table(w)
    es = np.zeros(n, np.int64)
    ed = np.zeros(n, np.int64)
    for x1 in prange(n):
        cs = 0
        cd = 0
        for x2 in range(n):
            for x3 in range(n):
                for x4 in range(n):
                    if x1 + x2 == x3 + x4 and rev[x1] + rev[x2] == rev[x3] + rev[x4]:
                        cs += 1
                    if x1 - x2 == x3 - x4 and rev[x1] - rev[x2] == rev[x3] - rev[x4]:
                        cd += 1
        es[x1] = cs
        ed[x1] = cd
    return es.sum(), ed.sum()


@njit(parallel=True, cache=False)
def energy3(w):
    n = 1 << w
    rev = rev_table(w)
    e = np.zeros(n, np.int64)
    for x1 in prange(n):
        c = 0
        r1 = rev[x1]
        for x2 in range(n):
            dx = x1 - x2
            dr = r1 - rev[x2]
            for x3 in range(n):
                x4 = x3 - dx
                if x4 >= 0 and x4 < n and rev[x3] - rev[x4] == dr:
                    c += 1
        e[x1] = c
    return e.sum()


def main():
    W4 = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    W3 = int(sys.argv[2]) if len(sys.argv) > 2 else 11
    ok = True
    print("=== E1: all 2^(4w) quadruples ===")
    for w in range(1, W4 + 1):
        t0 = time.time()
        es, ed = energy4(w)
        good = (int(es) == 6 ** w and int(ed) == 6 ** w)
        ok &= good
        print(f"  w = {w}: sum-form energy = {int(es)}, difference-form energy = {int(ed)}, 6^w = {6 ** w}: {good}   [{time.time()-t0:.1f}s]", flush=True)
    print("=== E2: three free points ===")
    for w in range(1, W3 + 1):
        t0 = time.time()
        e = energy3(w)
        good = (int(e) == 6 ** w)
        ok &= good
        print(f"  w = {w:2d}: energy = {int(e)}, 6^w = {6 ** w}: {good}   [{time.time()-t0:.1f}s]", flush=True)
    # control: a set that is NOT the reversal graph (identity graph {(x, x)}) has energy (2*8^w + 4^w)/3 != 6^w
    print("=== control (the code can tell 6^w from something else): graph of the identity ===")
    for w in range(1, 7):
        n = 1 << w
        cnt = 0
        for x1 in range(n):
            for x2 in range(n):
                for x3 in range(n):
                    x4 = x3 - (x1 - x2)
                    if 0 <= x4 < n:
                        cnt += 1
        print(f"  w = {w}: energy of {{(x, x)}} = {cnt}  (6^w = {6 ** w})")
    print("ALL CHECKS PASSED" if ok else "*** SOME CHECK FAILED ***")


if __name__ == "__main__":
    main()
