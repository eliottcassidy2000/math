#!/usr/bin/env python3
"""Verifier (lane digits), check V3: one-apex words.

(a) closed form of the stay block (TT_k, W) against the product of the one-vertex matrices a, c: all words k <= 14.
(b) F(W) = V U - (2^k - 1)(Wd + sigma eps) from the DEFINITIONS for every nonempty W, k <= 24:
    zeros, and the minimum over proper nonempty W.
(c) the identity used in the verifier's proof of uniqueness for all k:
    F(W) = 2^(k-1) |W| + sum_{p<q in W} (2^(k-1-(q-p)) + 2^(q-p-1)) - (2^k - 1)(c1 + sigma eps),
    and the lower bound  F(W) >= sum over runs (2^(k-r) + 2^r - 2) - (2^k - 1) sigma eps, checked for k <= 16.
"""
import sys, time
import numpy as np
from numba import njit, prange


@njit(cache=True)
def F_def(k, mask):
    V = mask
    U = 0
    for p in range(1, k + 1):
        if (mask >> (p - 1)) & 1:
            U += 1 << (k - p)
    Wd = 0
    for p in range(1, k):
        if (mask >> (p - 1)) & 1 and (mask >> p) & 1:
            Wd += 1
    for e in range(2, k + 1):
        if not (mask >> (e - 1)) & 1:
            continue
        for f in range(1, e):
            if (mask >> (f - 1)) & 1:
                Wd += 1 << (e - f - 1)
    sig = mask & 1
    eps = (mask >> (k - 1)) & 1
    return V * U - ((1 << k) - 1) * (Wd + sig * eps)


@njit(parallel=True, cache=True)
def scan(k):
    full = (1 << k) - 1
    nz = 0
    mn = np.int64(1) << 62
    nth = 4
    zs = np.zeros(nth, np.int64)
    mns = np.full(nth, np.int64(1) << 62)
    neg = np.zeros(nth, np.int64)
    for t in prange(nth):
        for mask in range(1 + t, full + 1, nth):
            F = F_def(k, mask)
            if F == 0:
                zs[t] += 1
                if mask != full:
                    neg[t] += 1000000          # a proper zero would be a counterexample
            elif mask != full:
                if F < mns[t]:
                    mns[t] = F
                if F < 0:
                    neg[t] += 1
    return zs.sum(), mns.min(), neg.sum()


def closed_form(k, W):
    sig = 1 if 1 in W else 0
    eps = 1 if k in W else 0
    V = sum(1 << (p - 1) for p in W)
    U = sum(1 << (k - p) for p in W)
    Wd = sum(1 for p in W if p + 1 in W) + sum(1 << (e - f - 1) for e in W for f in W if e > f)
    return [[1, -sig, V, -Wd], [0, 0, 0, eps], [0, 0, 1 << k, -U], [0, 0, 0, 1]], (sig, eps, V, U, Wd)


def mm(A, B):
    return [[sum(A[i][t] * B[t][j] for t in range(4)) for j in range(4)] for i in range(4)]


def main():
    t0 = time.time()
    a = [[1, -1, 1, 0], [0, 0, 0, 1], [0, 0, 2, -1], [0, 0, 0, 1]]
    c = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]
    bad = 0
    n = 0
    badid = 0
    badlb = 0
    for k in range(1, 15):
        for mask in range(1 << k):
            W = [p for p in range(1, k + 1) if (mask >> (p - 1)) & 1]
            M = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
            for p in range(1, k + 1):
                M = mm(a if p in W else c, M)        # letter p applied after letters 1..p-1
            N, (sig, eps, V, U, Wd) = closed_form(k, W)
            n += 1
            if M != N:
                bad += 1
            if k <= 16 and W:
                F = V * U - ((1 << k) - 1) * (Wd + sig * eps)
                c1 = sum(1 for p in W if p + 1 in W)
                F2 = (1 << (k - 1)) * len(W) * 2
                # doubled to stay in integers: 2 F = 2^k |W| + sum (2^(k-d) + 2^d) - 2 (2^k - 1)(c1 + sig eps)
                for i, p in enumerate(W):
                    for q in W[i + 1:]:
                        d = q - p
                        F2 += (1 << (k - d)) + (1 << d)
                F2 -= 2 * ((1 << k) - 1) * (c1 + sig * eps)
                if F2 != 2 * F:
                    badid += 1
                runs = []
                r = 0
                for p in range(1, k + 2):
                    if p in W:
                        r += 1
                    elif r:
                        runs.append(r)
                        r = 0
                lb = sum((1 << (k - r)) + (1 << r) - 2 for r in runs) - ((1 << k) - 1) * sig * eps
                if F < lb:
                    badlb += 1
    print(f"(a) closed form vs product of one-vertex matrices: {n} words (k <= 14, incl. empty W), mismatches {bad}")
    print(f"(c) identity for F(W): violations {badid}; lower bound by runs: violations {badlb}  (k <= 14)")
    print("(b) k: zeros of F among nonempty W / min F over proper nonempty W / proper W with F <= 0")
    for k in range(1, 25):
        z, mn, neg = scan(k)
        print(f"    k = {k:2d}: zeros {z}, min proper F = {mn if k > 1 else '-'}, proper W with F <= 0: {neg}   [{time.time() - t0:.1f} s]", flush=True)


if __name__ == "__main__":
    main()
