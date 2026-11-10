#!/usr/bin/env python3
"""HYP-9657 (A), second half: "the remaining p - 2 classes are equal up to sampling error; in particular the two
classes with Kronecker symbol (-N/p) = +1 and = -1 are not distinguished".  Test of that sentence.

r(N) = the packet's count (triples (h1, h2, s) = ordered pairs (u, v) with v even), computed from the canon form
F1(u, v) = (1+u)(1+v) + g(u)g(v) - e(u)e(v).  For each odd prime p <= 61 and for p = 73, 89, 127 on the window of odd
N in [2^lo, 2^K):
   * z0, z1 : (mean of class 0, resp. -1, minus mean of the other classes) in standard errors;
   * chi2   : sum over the classes c != 0, -1 of n_c (m_c - m_others)^2 / var(r), to be compared with a chi-square
              with p - 3 degrees of freedom (expected value p - 3, standard deviation sqrt(2(p - 3)));
   * zK     : (mean over classes with (-N/p) = +1) - (mean over classes with (-N/p) = -1), classes 0 and -1
              excluded, in standard errors;
   * z2     : (mean over the classes c in <2> mod p) - (mean over the other classes), 0 and -1 excluded, in standard
              errors (only when <2> is a proper subgroup that differs from the squares or p = 7).
Standard errors treat the r(N) as independent with the window variance; this is a heuristic yardstick, stated as such.
Usage: python3 audit_f1_classes.py [K=23] [lo=21]
"""
import sys
import time

import numpy as np


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 23
    LO = int(sys.argv[2]) if len(sys.argv) > 2 else 21
    X = 1 << K
    t0 = time.time()
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
    r = np.zeros(X + 1, dtype=np.int64)
    bufN = []
    bufW = []
    size = 0
    u = 0
    while (1 + u) * (1 + u) <= X:
        vmax = X // (1 + u) - 1
        vs = np.arange(u, vmax + 1, dtype=np.int64)
        N = (1 + u) * (1 + vs) + int(g[u]) * g[u:vmax + 1] - int(ev[u]) * ev[u:vmax + 1]
        keep = N <= X
        N = N[keep]
        vs = vs[keep]
        w = ((vs & 1) == 0).astype(np.int64) + np.where(vs != u, 1 - (u & 1), 0)
        bufN.append(N)
        bufW.append(w)
        size += len(N)
        if size > 6_000_000:
            r += np.bincount(np.concatenate(bufN), weights=np.concatenate(bufW), minlength=X + 1).astype(np.int64)
            bufN, bufW, size = [], [], 0
        u += 1
    if size:
        r += np.bincount(np.concatenate(bufN), weights=np.concatenate(bufW), minlength=X + 1).astype(np.int64)
    odd = np.arange(1, X + 1, 2)
    win = odd[odd >= (1 << LO)]
    rw = r[win].astype(np.float64)
    mean = rw.mean()
    var = rw.var()
    print(f"odd N in [2^{LO}, 2^{K}): {len(win)} integers, mean r = {mean:.4f}, sd = {var ** 0.5:.3f}, misses {(rw == 0).sum()}",
          f"[{time.time()-t0:.0f}s]")
    print(" p   class0  class-1  others(min..max)        z0      z1    chi2/dof(dof)    zK      z2   [<2> index]")
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 73, 89, 127):
        res = win % p
        sums = np.bincount(res, weights=rw, minlength=p)
        ns = np.bincount(res, minlength=p).astype(np.float64)
        means = sums / ns
        others = list(range(1, p - 1))
        mo = sums[others].sum() / ns[others].sum() if len(others) else float("nan")
        se1 = (var / ns.mean()) ** 0.5
        z0 = (means[0] - mo) / (se1 * (1 + 1 / max(len(others), 1)) ** 0.5)
        z1 = (means[p - 1] - mo) / (se1 * (1 + 1 / max(len(others), 1)) ** 0.5)
        chi2 = float(((means[others] - mo) ** 2 * ns[others]).sum() / var)
        dof = p - 3
        plus = [c for c in others if pow((-c) % p, (p - 1) // 2, p) == 1]
        minus = [c for c in others if c not in plus]
        if plus and minus:
            mp = sums[plus].sum() / ns[plus].sum()
            mm = sums[minus].sum() / ns[minus].sum()
            zK = (mp - mm) / (var * (1 / ns[plus].sum() + 1 / ns[minus].sum())) ** 0.5
            zKs = f"{zK:+7.1f}"
        else:
            zKs = "    n/a"
        sub2 = set()
        x = 1
        while x not in sub2:
            sub2.add(x)
            x = 2 * x % p
        idx2 = (p - 1) // len(sub2)
        in2 = [c for c in others if c in sub2]
        out2 = [c for c in others if c not in sub2]
        if in2 and out2:
            m_in = sums[in2].sum() / ns[in2].sum()
            m_out = sums[out2].sum() / ns[out2].sum()
            z2 = (m_in - m_out) / (var * (1 / ns[in2].sum() + 1 / ns[out2].sum())) ** 0.5
            z2s = f"{z2:+7.1f}"
        else:
            z2s = "    n/a"
        rng = f"{means[others].min():.4f}..{means[others].max():.4f}" if len(others) else "n/a"
        chis = f"{chi2/dof:8.2f} ({dof:3d})" if dof > 0 else "     n/a      "
        print(f"{p:3d}  {means[0]:7.4f} {means[p-1]:7.4f}  {rng:22s} {z0:+7.1f} {z1:+7.1f}  {chis}  {zKs} {z2s}   [{idx2}]")
        if p in (7, 11, 13, 31):
            print("       class means c = 0..p-1:", " ".join(f"{m:.3f}" for m in means))
    print(f"done [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
