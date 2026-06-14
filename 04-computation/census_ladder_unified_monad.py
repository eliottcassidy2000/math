#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The unified closed-walk census ladder  (monad-explorer-2026-06-14).

ONE principle generates every trace decomposition tr(A^k) for k <= 8:

  tr(A^k) = SUM over connected cycle-configurations C of total length k of
            (k / period(C)) * (#embeddings of C in T).

Configurations (parts >= 3, no 1/2-cycles in a tournament):
  - simple k-cycle:                 period k   -> coeff k        -> k * c_k
  - (k/2)-cycle traversed twice:    period k/2 -> coeff k/2      -> (k/2) * c_{k/2}   [k even]
  - two DISTINCT a,b-cycles (a+b=k) glued at >=1 shared vertex (asymmetric):
                                    period k   -> coeff k        -> k * O(a,b)

So, verified here exactly:
  tr A^6 = 6 c6 + 3 c3 + 6 p33            (p33 = overlapping triangle pairs)
  tr A^7 = 7 c7 + 7 TQ                    (TQ  = overlapping (tri,4cyc) pairs)
  tr A^8 = 8 c8 + 4 c4 + 8 Q44 + 8 TF     (Q44 = overlapping 4cyc pairs; TF = (tri,5cyc))

CONSERVATION COROLLARY: within a cospectral class (tr A^k fixed for all k, and all
spectral lower counts fixed) the SIMPLE cycle count trades 1-for-1 against the total
OVERLAP count:
  c6 + p33          = (tr A^6 - 3 c3)/6          (spectral) -> const on cospectral class
  c7 + TQ           =  tr A^7 / 7                (spectral) -> const
  c8 + Q44 + TF     = (tr A^8 - 4 c4)/8          (spectral) -> const
This is the precise mechanism of c_k non-spectrality (k>=6): the spectrum fixes the
SUM (closed walks) but not the SPLIT (simple vs reducible/overlapping).
"""

import sys, itertools, random
import numpy as np
from collections import defaultdict

random.seed(14)


def random_tournament(n):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if random.getrandbits(1):
                A[i, j] = 1
            else:
                A[j, i] = 1
    return A


def simple_cycles_by_len(A, n, maxL):
    adj = [[j for j in range(n) if A[i, j]] for i in range(n)]
    out = defaultdict(list)
    for start in range(n):
        stack = [(start, 1 << start, (start,))]
        while stack:
            v, vis, path = stack.pop()
            for w in adj[v]:
                if w == start and len(path) >= 3:
                    out[len(path)].append(path)
                elif w > start and not (vis >> w) & 1 and len(path) < maxL:
                    stack.append((w, vis | (1 << w), path + (w,)))
    return out


def opairs(la, lb, same):
    sa = [frozenset(x) for x in la]
    cnt = 0
    if same:
        for i in range(len(sa)):
            for j in range(i + 1, len(sa)):
                if sa[i] & sa[j]:
                    cnt += 1
    else:
        sb = [frozenset(x) for x in lb]
        for x in sa:
            for y in sb:
                if x & y:
                    cnt += 1
    return cnt


def trace(A, k):
    return int(np.trace(np.linalg.matrix_power(A.astype(object), k)))


def main():
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    NS = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    bad6 = bad7 = bad8 = 0
    cons6 = defaultdict(set); cons7 = defaultdict(set); cons8 = defaultdict(set)
    for _ in range(NS):
        A = random_tournament(n)
        cyc = simple_cycles_by_len(A, n, 8)
        c3 = len(cyc[3]); c4 = len(cyc[4]); c5 = len(cyc[5])
        c6 = len(cyc[6]); c7 = len(cyc[7]); c8 = len(cyc[8])
        p33 = opairs(cyc[3], None, True)
        TQ = opairs(cyc[3], cyc[4], False)
        Q44 = opairs(cyc[4], None, True)
        TF = opairs(cyc[3], cyc[5], False)
        tr6, tr7, tr8 = trace(A, 6), trace(A, 7), trace(A, 8)
        if tr6 != 6 * c6 + 3 * c3 + 6 * p33: bad6 += 1
        if tr7 != 7 * c7 + 7 * TQ: bad7 += 1
        if tr8 != 8 * c8 + 4 * c4 + 8 * Q44 + 8 * TF: bad8 += 1
        # spectral signature for conservation test
        sig = tuple(trace(A, k) for k in range(3, n + 1))
        cons6[sig].add(c6 + p33 - (tr6 - 3 * c3) // 6)
        cons7[sig].add(c7 + TQ - tr7 // 7)
        cons8[sig].add(c8 + Q44 + TF - (tr8 - 4 * c4) // 8)

    print(f"=== UNIFIED CENSUS LADDER  (n={n}, {NS} random tournaments) ===", flush=True)
    print(f"  tr A^6 = 6 c6 + 3 c3 + 6 p33        : exact {NS-bad6}/{NS}", flush=True)
    print(f"  tr A^7 = 7 c7 + 7 TQ                : exact {NS-bad7}/{NS}", flush=True)
    print(f"  tr A^8 = 8 c8 + 4 c4 + 8 Q44 + 8 TF : exact {NS-bad8}/{NS}", flush=True)
    # the conservation-offset must be identically 0 (algebraic), confirming RHS is spectral
    z6 = all(v == {0} for v in cons6.values())
    z7 = all(v == {0} for v in cons7.values())
    z8 = all(v == {0} for v in cons8.values())
    print(f"  conservation: c6+p33      = (trA^6-3c3)/6   identically: {z6}", flush=True)
    print(f"  conservation: c7+TQ       =  trA^7/7        identically: {z7}", flush=True)
    print(f"  conservation: c8+Q44+TF   = (trA^8-4c4)/8   identically: {z8}", flush=True)


if __name__ == "__main__":
    main()
