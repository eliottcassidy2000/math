#!/usr/bin/env python3
"""Verifier (lane digits), check V8: actual two-marked tournaments and exact affine pairs of the two-mark rule system.

Own search logic.  States of all two-marked tournaments of order <= NMAX by own direct counts (17 state words of
twomark_rep_5_3.json); words act through the 17 x 17 matrices of that file (exact for every base: check V7).
Suffix words p: all 12 letters and all 144 two-letter words.  The multiplier is read off the data,
    a_p = (H(spp) - H(sp)) / (H(sp) - H(s))  (integer, H(sp) != H(s)),
single test: H(s p^(j+1)) = a_p H(s p^j) + b_p, j = 0, 1, 2;  pair test: the rules hold at every node of
s {p, q}^j, j <= 2, and the two maps do not commute.

Usage: v8_twomark_pairs.py NMAX
"""
import json
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit

sys.path.insert(0, sys.path[0] + "/..")
from tour import iso_classes  # noqa: E402

from v7_twomark_trace import OPS, H_word  # noqa: E402


@njit(cache=True)
def scan(S, M, hist, movers):
    ns, d = S.shape
    nc = M.shape[0]
    chain = np.zeros((4, d), np.int64)
    tree = np.zeros((15, d), np.int64)
    okc = np.zeros(nc, np.int64)
    ac = np.zeros(nc, np.int64)
    bc = np.zeros(nc, np.int64)
    pre = 0
    npass = 0
    for si in range(ns):
        nok = 0
        for c in range(nc):
            for k in range(d):
                chain[0, k] = S[si, k]
            for j in range(3):
                for i in range(d):
                    t = 0
                    for k in range(d):
                        t += M[c, i, k] * chain[j, k]
                    chain[j + 1, i] = t
            okc[c] = 0
            d0 = chain[1, 0] - chain[0, 0]
            if d0 == 0:
                continue
            d1 = chain[2, 0] - chain[1, 0]
            if d1 % d0 != 0:
                continue
            a = d1 // d0
            if chain[3, 0] - chain[2, 0] != a * d1:
                continue
            okc[c] = 1
            ac[c] = a
            bc[c] = chain[1, 0] - a * chain[0, 0]
            nok += 1
            movers[c] += 1
            if 0 <= a < 12:
                hist[c, a] += 1
        if nok < 2:
            continue
        for c1 in range(nc):
            if okc[c1] == 0:
                continue
            for c2 in range(c1 + 1, nc):
                if okc[c2] == 0:
                    continue
                if bc[c1] * (ac[c2] - 1) == bc[c2] * (ac[c1] - 1):
                    continue
                pre += 1
                for k in range(d):
                    tree[0, k] = S[si, k]
                good = True
                for i in range(7):
                    for r in range(d):
                        t1 = 0
                        t2 = 0
                        for k in range(d):
                            t1 += M[c1, r, k] * tree[i, k]
                            t2 += M[c2, r, k] * tree[i, k]
                        tree[2 * i + 1, r] = t1
                        tree[2 * i + 2, r] = t2
                    if tree[2 * i + 1, 0] != ac[c1] * tree[i, 0] + bc[c1] or tree[2 * i + 2, 0] != ac[c2] * tree[i, 0] + bc[c2]:
                        good = False
                        break
                if good:
                    npass += 1
    return pre, npass


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    t0 = time.time()
    import os
    rep = json.load(open(os.path.dirname(os.path.abspath(__file__)) + "/../twomark_rep_5_3.json"))
    A = [np.array([[int(Fraction(x)) for x in row] for row in M], dtype=np.int64) for M in rep["mats"]]
    sw = [tuple(w) for w in rep["state_words"]]
    classes = iso_classes(nmax)
    states = set()
    nb = 0
    for n in range(2, nmax + 1):
        for key, adj in classes[n].items():
            adj = np.array(adj, np.int64)
            for z1 in range(n):
                for z2 in range(n):
                    if z1 != z2:
                        nb += 1
                        states.add(tuple(H_word(n, adj, z1, z2, w) for w in sw))
        print(f"   order {n}: marked bases so far {nb}, distinct states {len(states)}   [{time.time() - t0:.1f} s]", flush=True)
    S = np.array(sorted(states), dtype=np.int64)
    names = [(o,) for o in range(12)] + [(o1, o2) for o1 in range(12) for o2 in range(12)]
    M = np.array([A[w[0]] if len(w) == 1 else A[w[1]] @ A[w[0]] for w in names], dtype=np.int64)
    hist = np.zeros((len(names), 12), np.int64)
    movers = np.zeros(len(names), np.int64)
    pre, npass = scan(S, M, hist, movers)
    print(f"{len(S)} distinct states of the {nb} two-marked tournaments of orders 2..{nmax}; suffix words {len(names)}")
    print(f"non-commuting pairs with both single rays exact to depth 2: {pre}; EXACT NON-COMMUTING PAIRS to depth 2: {npass}")

    def wn(w):
        return ".".join(f"({OPS[o][0]:02b},{OPS[o][1]})" for o in w)

    print("one-letter maps exact to depth 2 with H moving: op: states; multipliers (a: count)")
    for c, w in enumerate(names):
        if len(w) == 1 and movers[c]:
            print(f"   {wn(w)}: {int(movers[c])}; { {a: int(hist[c, a]) for a in range(12) if hist[c, a]} }")
    print(f"two-letter words with at least one moving exact state: {sum(1 for c, w in enumerate(names) if len(w) == 2 and movers[c])}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
