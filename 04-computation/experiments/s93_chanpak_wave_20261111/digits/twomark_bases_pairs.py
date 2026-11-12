#!/usr/bin/env python3
"""Lane "digits", finishing pass: do ACTUAL two-marked tournaments of orders <= NMAX carry an exact non-commuting
affine pair of one-vertex rules?  (The lane found 177 "virtual" class spaces for words of length <= 2 and no census
state inside; its census had the 302 marked bases of orders 2..5 and their images under short words.  Here: every
two-marked tournament (T; z1, z2) of order <= NMAX, one representative per isomorphism class of T, all ordered pairs.)

State of a base = (H(w_k(T; z1, z2)))_k for the 17 state words of twomark_rep_5_3.json (exact direct counts);
words act by the 17 x 17 integer matrices of that file (checked here against direct counts on a random sample).
Candidates: the 12 one-vertex operations with multiplier a in {1, 2, 3}; the 144 two-letter words with a in
{1, 3, 4, 9} (all integer eigenvalues >= 1 of words of length <= 2, twomark_classify_2.out).
single test: b_p = H(s p) - a H(s) constant on s, s p, s p p (depth 2);  pair test: b_p, b_q constant on every node
of s {p, q}^j, j <= 2, and b_p (a_q - 1) != b_q (a_p - 1).

Usage: twomark_bases_pairs.py NMAX
"""
import json
import random
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit

from tour import iso_classes
from twomark import OPS, H_of


@njit(cache=False)
def scan(S, M, cw, ca, depth, single_pass, single_moving, res, maxres):
    ns, d = S.shape
    nc = ca.shape[0]
    chain = np.zeros((depth + 2, d), np.int64)
    nn = (1 << (depth + 2)) - 1
    tree = np.zeros((nn, d), np.int64)
    okc = np.zeros(nc, np.int64)
    bc = np.zeros(nc, np.int64)
    npass = 0
    for si in range(ns):
        h = S[si, 0]
        nok = 0
        for c in range(nc):
            for k in range(d):
                chain[0, k] = S[si, k]
            good = True
            b = 0
            for j in range(depth + 1):
                for i in range(d):
                    t = 0
                    for k in range(d):
                        t += M[cw[c], i, k] * chain[j, k]
                    chain[j + 1, i] = t
                bj = chain[j + 1, 0] - ca[c] * chain[j, 0]
                if j == 0:
                    b = bj
                elif bj != b:
                    good = False
                    break
            okc[c] = 1 if good else 0
            bc[c] = b
            if good:
                nok += 1
                single_pass[c] += 1
                if b != (1 - ca[c]) * h:
                    single_moving[c] += 1
        if nok < 2:
            continue
        for c1 in range(nc):
            if okc[c1] == 0:
                continue
            for c2 in range(c1 + 1, nc):
                if okc[c2] == 0 or cw[c1] == cw[c2]:
                    continue
                if bc[c1] * (ca[c2] - 1) == bc[c2] * (ca[c1] - 1):
                    continue
                for k in range(d):
                    tree[0, k] = S[si, k]
                good = True
                for i in range((1 << (depth + 1)) - 1):
                    for r in range(d):
                        t1 = 0
                        t2 = 0
                        for k in range(d):
                            t1 += M[cw[c1], r, k] * tree[i, k]
                            t2 += M[cw[c2], r, k] * tree[i, k]
                        tree[2 * i + 1, r] = t1
                        tree[2 * i + 2, r] = t2
                    if tree[2 * i + 1, 0] - ca[c1] * tree[i, 0] != bc[c1] or tree[2 * i + 2, 0] - ca[c2] * tree[i, 0] != bc[c2]:
                        good = False
                        break
                if good:
                    if npass < maxres:
                        res[npass, 0] = si
                        res[npass, 1] = c1
                        res[npass, 2] = c2
                        res[npass, 3] = bc[c1]
                        res[npass, 4] = bc[c2]
                    npass += 1
    return npass


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    t0 = time.time()
    rep = json.load(open("twomark_rep_5_3.json"))
    d = rep["d"]
    mats = [np.array([[int(Fraction(x)) for x in row] for row in M], dtype=np.int64) for M in rep["mats"]]
    swords = [tuple(w) for w in rep["state_words"]]
    classes = iso_classes(nmax)
    states = {}
    nb = 0
    rng = random.Random(931)
    bad = 0
    ntest = 0
    for n in range(2, nmax + 1):
        for key, adj in classes[n].items():
            for z1 in range(n):
                for z2 in range(n):
                    if z1 == z2:
                        continue
                    nb += 1
                    s = tuple(H_of(n, adj, z1, z2, w) for w in swords)
                    states.setdefault(s, (n, [int(x) for x in adj], z1, z2))
                    if rng.random() < 0.01:
                        w = [rng.randrange(12) for _ in range(5)]
                        x = np.array(s, dtype=np.int64)
                        for o in w:
                            x = mats[o] @ x
                        ntest += 1
                        if int(x[0]) != H_of(n, adj, z1, z2, w):
                            bad += 1
        print(f"   order {n}: marked bases so far {nb}, distinct states {len(states)}  [{time.time() - t0:.1f} s]", flush=True)
    print(f"representation against direct counts: {ntest} random (base, word of length 5) checks, failures {bad}")
    S = np.array(sorted(states), dtype=np.int64)
    # words of length <= 2 as matrices; word (o1, o2) = apply o1 first
    W = []
    names = []
    for o in range(12):
        W.append(mats[o])
        names.append((o,))
    for o1 in range(12):
        for o2 in range(12):
            W.append(mats[o2] @ mats[o1])
            names.append((o1, o2))
    M = np.array(W, dtype=np.int64)
    cw = []
    ca = []
    for i, w in enumerate(names):
        for a in ((1, 2, 3) if len(w) == 1 else (1, 3, 4, 9)):
            cw.append(i)
            ca.append(a)
    cw = np.array(cw, dtype=np.int64)
    ca = np.array(ca, dtype=np.int64)
    sp = np.zeros(len(ca), np.int64)
    sm = np.zeros(len(ca), np.int64)
    res = np.zeros((200, 5), np.int64)
    npass = scan(S, M, cw, ca, 2, sp, sm, res, 200)

    def wn(w):
        return ".".join(f"({OPS[o][0]:02b},{OPS[o][1]})" for o in w)

    print(f"{len(S)} distinct states of the {nb} two-marked tournaments of orders 2..{nmax}; candidates {len(ca)}")
    print(f"NON-COMMUTING EXACT PAIRS (state, pair) passing to depth 2: {npass}")
    for i in range(min(npass, 20)):
        si, c1, c2, b1, b2 = [int(x) for x in res[i]]
        print(f"   PASS: base {states[tuple(int(x) for x in S[si])]} H = {int(S[si, 0])}: {wn(names[cw[c1]])} (a = {int(ca[c1])}, b = {b1}) & {wn(names[cw[c2]])} (a = {int(ca[c2])}, b = {b2})")
    print("single one-letter maps exact to depth 2: op (S, mu), a: states passing / with H moving")
    for c in range(len(ca)):
        if len(names[cw[c]]) == 1 and sp[c]:
            print(f"   {wn(names[cw[c]])}, a = {int(ca[c])}: {int(sp[c])} / {int(sm[c])}")
    mv2 = sum(1 for c in range(len(ca)) if len(names[cw[c]]) == 2 and sm[c] > 0)
    print(f"two-letter candidates with at least one moving exact state: {mv2}")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
