#!/usr/bin/env python3
"""Lane "digits": which words in the letters a (ray step) and c (sink), read as ONE stay-operation of the one-mark
system, have an exact affine class h -> 2^k h + b containing a moving actual state?

A word of length k in {a, c} is the stay block (TT_k, W), W = positions of the letters a (position p = 1 is the
first new vertex, the source of the block).  Its matrix on (h, e, B, g) is (PROOFS.txt, P1)

    [ 1  -sigma   V   -Wd ]        sigma = [1 in W],  eps = [k in W],
    [ 0    0      0   eps ]        V  = sum_{p in W} 2^(p-1),   U = sum_{p in W} 2^(k-p),
    [ 0    0     2^k  -U  ]        Wd = #{p, p+1 in W} + sum_{e > f in W} 2^(e-f-1).
    [ 0    0      0    1  ]

Criterion (PROOFS.txt, P4): a state with g != 0 lies in the class space of the multiplier 2^k iff
    V U = (2^k - 1) (Wd + sigma eps)      and      (sigma = 0  or  e = eps g).
This script (1) checks the matrix formula against the interpolated block matrices of onemark_matrices_5.json,
(2) lists, for every k <= KMAX, all nonempty W satisfying the criterion.

Usage: oneapex_exact_words.py KMAX
"""
import json
import sys
import time

import numpy as np
from numba import njit


@njit(cache=True)
def scan(k, out, maxout):
    full = (1 << k) - 1
    cnt = 0
    for W in range(1, 1 << k):
        V = W
        U = 0
        for p in range(k):
            if (W >> p) & 1:
                U |= 1 << (k - 1 - p)
        sigma = W & 1
        eps = (W >> (k - 1)) & 1
        Wd = 0
        for e in range(k):
            if (W >> e) & 1:
                for f in range(e):
                    if (W >> f) & 1:
                        Wd += 1 << (e - f - 1)
                        if e == f + 1:
                            Wd += 1
        if V * U == full * (Wd + sigma * eps):
            if cnt < maxout:
                out[cnt] = W
            cnt += 1
    return cnt


def matrix(k, W):
    V = W
    U = sum(1 << (k - 1 - p) for p in range(k) if (W >> p) & 1)
    sigma = W & 1
    eps = (W >> (k - 1)) & 1
    Wd = 0
    for e in range(k):
        if (W >> e) & 1:
            for f in range(e):
                if (W >> f) & 1:
                    Wd += 1 << (e - f - 1)
                    if e == f + 1:
                        Wd += 1
    return [[1, -sigma, V, -Wd], [0, 0, 0, eps], [0, 0, 1 << k, -U], [0, 0, 0, 1]]


def main():
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    t0 = time.time()
    # (1) check against the interpolated matrices: blocks X = TT_k in the labelling of onemark.py
    ops = json.load(open("onemark_matrices_5.json"))
    have = {tuple(tuple(r) for r in o["N"]) for o in ops}
    ok = True
    n = 0
    for k in range(1, 6):
        for W in range(1 << k):
            n += 1
            if tuple(tuple(r) for r in matrix(k, W)) not in have:
                ok = False
    print(f"(1) closed form of the {n} stay blocks (TT_k, W), k <= 5, found among the interpolated block matrices: {ok}")
    # (2)
    print("(2) nonempty W with V U = (2^k - 1)(Wd + sigma eps):")
    out = np.zeros(64, np.int64)
    for k in range(1, kmax + 1):
        c = scan(k, out, 64)
        sols = [format(int(out[i]), f"0{k}b")[::-1] for i in range(min(c, 64))]
        print(f"   k = {k:2d}: {c} solution(s): {sols}   (letters in order of creation; 1 = a, 0 = c)   [{time.time() - t0:.1f} s]", flush=True)


if __name__ == "__main__":
    main()
