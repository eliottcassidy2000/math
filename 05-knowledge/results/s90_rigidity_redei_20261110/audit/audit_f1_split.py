#!/usr/bin/env python3
"""HYP-9657 (A), the "reading": triples with s >= 1 solve N = x1 x2 + 2^s y1 y2 (class 0 favoured), triples with
s = 0 solve N + 1 = x1 x2 + y1 y2 (class -1 favoured), each "with about 1/p more solutions" in its favoured class.
Test: split the packet's r = r0 + r1 (s = 0: both u, v even in the canon coordinates; s >= 1: otherwise) and print,
for each prime, the class means of both parts at c = 0, c = -1 and over the other classes, the observed relative
excess (m_c / m_others - 1) and the heuristic prediction for a generic sum of two products,
p^2/(p^3 - p) = 1/(p - 1/p)  (ratio (p^3 + p^2 - p)/(p^3 - p) of the solution counts of ab + cd = 0 and = t != 0).
Usage: python3 audit_f1_split.py [K=23] [lo=21]
"""
import sys

import numpy as np


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 23
    LO = int(sys.argv[2]) if len(sys.argv) > 2 else 21
    X = 1 << K
    V = X + 2
    g = np.zeros(V, dtype=np.int64)
    g[0] = 1
    k = 0
    while (1 << k) < V:
        lo, hi = 1 << k, min(1 << (k + 1), V)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        g[lo:hi] = np.where(idx & 1 == 1, g[half], g[half] + half + (half & 1))
        k += 1
    ev = ((np.arange(V) & 1) == 0).astype(np.int64)
    r0 = np.zeros(X + 1, dtype=np.int64)
    r1 = np.zeros(X + 1, dtype=np.int64)
    bN, b0, b1 = [], [], []
    size = 0

    def flush():
        nonlocal bN, b0, b1, size
        if size:
            Nn = np.concatenate(bN)
            r0[:] += np.bincount(Nn, weights=np.concatenate(b0), minlength=X + 1).astype(np.int64)
            r1[:] += np.bincount(Nn, weights=np.concatenate(b1), minlength=X + 1).astype(np.int64)
        bN, b0, b1, size = [], [], [], 0

    u = 0
    while (1 + u) * (1 + u) <= X:
        vmax = X // (1 + u) - 1
        vs = np.arange(u, vmax + 1, dtype=np.int64)
        N = (1 + u) * (1 + vs) + int(g[u]) * g[u:vmax + 1] - int(ev[u]) * ev[u:vmax + 1]
        keep = N <= X
        N = N[keep]
        vs = vs[keep]
        ue = 1 - (u & 1)
        ve = ((vs & 1) == 0).astype(np.int64)
        ne = (vs != u).astype(np.int64)
        # ordered (u, v) with v even: s = 0 iff u even;  ordered (v, u) with u even (v != u): s = 0 iff v even
        w0 = ve * ue + ne * ue * ve
        w1 = ve * (1 - ue) + ne * ue * (1 - ve)
        bN.append(N)
        b0.append(w0)
        b1.append(w1)
        size += len(N)
        if size > 6_000_000:
            flush()
        u += 1
    flush()
    odd = np.arange(1, X + 1, 2)
    win = odd[odd >= (1 << LO)]
    a0 = r0[win].astype(float)
    a1 = r1[win].astype(float)
    print(f"odd N in [2^{LO}, 2^{K}): mean of the s = 0 part {a0.mean():.4f}, of the s >= 1 part {a1.mean():.4f}, total "
          f"{(a0 + a1).mean():.4f}")
    print(" p | s = 0 part: excess at c=0, at c=-1 | s >= 1 part: excess at c=0, at c=-1 | generic two-product prediction")
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        res = win % p
        line = f"{p:3d} |"
        for arr in (a0, a1):
            sums = np.bincount(res, weights=arr, minlength=p)
            ns = np.bincount(res, minlength=p).astype(float)
            means = sums / ns
            others = list(range(1, p - 1))
            mo = sums[others].sum() / ns[others].sum()
            line += f"   {means[0] / mo - 1:+.4f}  {means[p - 1] / mo - 1:+.4f}   |"
        line += f"   {p * p / (p ** 3 - p):.4f}"
        print(line)


if __name__ == "__main__":
    main()
