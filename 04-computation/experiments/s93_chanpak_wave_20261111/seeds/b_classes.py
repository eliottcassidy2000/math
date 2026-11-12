#!/usr/bin/env python3
"""S93 lane "seeds": isomorphism classes of tournaments of order <= NMAX (default 9), own generator.

Method: every tournament of order n has a vertex of minimum score; deleting it leaves a tournament of order n - 1.
So the classes of order n are the canonical forms of T + v, T a class representative of order n - 1 and v a new
vertex whose score is <= every score of T + v.  Canonical form: least code (column-major upper triangle, compared
chunk by chunk) over all vertex orders that respect the cells of an isomorphism-invariant vertex colouring.
Check: the class counts 1, 1, 2, 4, 12, 56, 456, 6880, 191536 (OEIS A000568).
Output: cls_<n>.txt, one column-major code per line (seedlib.out_from_colcode).
"""
import sys
import time

import numpy as np
from numba import njit

MASK = (1 << 61) - 1


@njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def canon(n, out):
    inv = np.zeros(n, np.int64)
    tmp = np.zeros(n, np.int64)
    for v in range(n):
        inv[v] = popcount(out[v])
    for rnd in range(3):
        for v in range(n):
            a = 0
            b = 0
            for u in range(n):
                if u != v:
                    x = inv[u]
                    h = (x * 2654435761 + 12345) & MASK
                    h = ((h ^ (h >> 13)) * 1000003) & MASK
                    if (out[v] >> u) & 1:
                        a = (a + h) & MASK
                    else:
                        b = (b + ((h * 7919 + 1) & MASK)) & MASK
            tmp[v] = ((inv[v] * 1315423911) + a * 31 + b * 131071) & MASK
        for v in range(n):
            inv[v] = tmp[v]
    order = np.argsort(inv)
    cv = np.zeros(n, np.int64)          # cell of a vertex
    cid = np.zeros(n, np.int64)         # cell of a position
    c = 0
    for k in range(n):
        if k > 0 and inv[order[k]] != inv[order[k - 1]]:
            c += 1
        cv[order[k]] = c
        cid[k] = c
    best = np.zeros(n, np.int64)
    chunk = np.zeros(n, np.int64)
    perm = np.zeros(n, np.int64)
    cidx = np.full(n, -1, np.int64)
    tied = np.zeros(n + 1, np.bool_)
    have = False
    used = 0
    depth = 0
    while depth >= 0:
        v = cidx[depth] + 1
        while v < n and (((used >> v) & 1) == 1 or cv[v] != cid[depth]):
            v += 1
        if v == n:
            cidx[depth] = -1
            depth -= 1
            if depth >= 0:
                used ^= 1 << perm[depth]
            continue
        cidx[depth] = v
        ch = 0
        for i in range(depth):
            if (out[perm[i]] >> v) & 1:
                ch |= 1 << i
        t = tied[depth]
        tnew = t
        if have and t:
            if ch > best[depth]:
                continue
            if ch < best[depth]:
                tnew = False
        perm[depth] = v
        chunk[depth] = ch
        if depth == n - 1:
            if (not have) or (not tnew):
                for k in range(n):
                    best[k] = chunk[k]
                    tied[k] = True
                have = True
            continue
        used |= 1 << v
        tied[depth + 1] = tnew and have
        depth += 1
    code = 0
    for k in range(1, n):
        code |= best[k] << (k * (k - 1) // 2)
    return code


@njit(cache=True)
def decode(n, code, out):
    for v in range(n):
        out[v] = 0
    b = 0
    for j in range(1, n):
        for i in range(j):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            b += 1


@njit(cache=True)
def extend(n, reps, res):
    """n = new order; reps = codes of order n - 1; res receives canonical codes of the min-score extensions."""
    m = n - 1
    out = np.zeros(n, np.int64)
    base = np.zeros(m, np.int64)
    cnt = 0
    for idx in range(reps.shape[0]):
        decode(m, reps[idx], base)
        for D in range(1 << m):
            d = popcount(D)
            ok = True
            for u in range(m):
                s = popcount(base[u]) + (0 if (D >> u) & 1 else 1)
                if d > s:
                    ok = False
                    break
            if not ok:
                continue
            for u in range(m):
                out[u] = base[u]
                if not ((D >> u) & 1):
                    out[u] |= 1 << m
            out[m] = D
            res[cnt] = canon(n, out)
            cnt += 1
    return cnt


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    known = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536, 10: 9733056}
    reps = np.array([0], np.int64)
    np.savetxt("cls_1.txt", reps, fmt="%d")
    print("order 1: 1 class")
    for n in range(2, nmax + 1):
        t0 = time.time()
        res = np.zeros(reps.shape[0] * (1 << (n - 1)), np.int64)
        cnt = extend(n, reps, res)
        reps = np.unique(res[:cnt])
        np.savetxt(f"cls_{n}.txt", reps, fmt="%d")
        print(f"order {n}: {cnt} min-score extensions, {reps.shape[0]} classes (known {known[n]}): "
              f"{'OK' if reps.shape[0] == known[n] else 'MISMATCH'}  [{time.time() - t0:.1f} s]", flush=True)


if __name__ == "__main__":
    main()
