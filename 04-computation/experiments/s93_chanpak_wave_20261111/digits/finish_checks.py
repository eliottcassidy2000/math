#!/usr/bin/env python3
"""Lane "digits", finishing pass (independent checks of the interrupted lane; written by the finishing agent).

(V1) Integer eigenvalues of every stored operation matrix, recomputed EXACTLY (no floating point): the
     characteristic polynomial is monic with integer coefficients, so every rational eigenvalue is an integer
     dividing the lowest nonzero coefficient; all divisors are tested.  Compared with affine.int_eigs (which
     uses the window [-8, 64] plus float roots).  Files: gadget_matrices_{3,4,5,6}.json, onemark_matrices_{3,4,5}.json,
     twomark_rep_5_3.json (single letters and all words of length 2).
(V2) Closed form of the stay block (TT_k, W) of oneapex_exact_words.py against the product of the one-vertex
     matrices a (ray step) and c (sink, stay) for every word of length k <= KMAX2.
(V3) The one-vertex binary tree H -> 2H+1, 2H-1 recomputed with different code (hcount on the explicit
     extension, no ext_values, no canonical forms): for every isomorphism class T of order <= 6, is T in K_2?
     Also: which classes of order <= 5 have a one-vertex / a two-vertex step to 2H+1 and to 2H-1.

Usage: finish_checks.py [KMAX2]
"""
import itertools
import json
import sys
import time
from fractions import Fraction

import numpy as np

from affine import charpoly, int_eigs, matmul, polyval
from tour import extend, hcount, iso_classes


def divisors(n):
    n = abs(n)
    out = set()
    i = 1
    while i * i <= n:
        if n % i == 0:
            out.add(i)
            out.add(n // i)
        i += 1
    return out


def exact_int_eigs(M):
    c = charpoly(M)                 # c[0..d], monic
    k = 0
    out = set()
    while c[k] == 0:
        k += 1
    if k > 0:
        out.add(0)
    low = c[k]
    for dv in divisors(low):
        for t in (dv, -dv):
            if polyval(c, t) == 0:
                out.add(t)
    return sorted(out)


def v1():
    print("(V1) exact integer eigenvalues (divisor test) against affine.int_eigs")
    for name in ["gadget_matrices_3.json", "gadget_matrices_4.json", "gadget_matrices_5.json", "gadget_matrices_6.json"]:
        gad = json.load(open(name))
        diff = 0
        cnt = {}
        for g in gad:
            e1 = [a for a in exact_int_eigs(g["G"]) if a >= 1]
            e2 = [a for a in int_eigs(g["G"]) if a >= 1]
            if e1 != e2:
                diff += 1
            for a in e1:
                cnt[a] = cnt.get(a, 0) + 1
        print(f"   {name}: {len(gad)} matrices, disagreements {diff}, exact integer eigenvalues >= 1: {dict(sorted(cnt.items()))}")
    for name in ["onemark_matrices_3.json", "onemark_matrices_4.json", "onemark_matrices_5.json"]:
        ops = json.load(open(name))
        diff = 0
        cnt = {}
        for o in ops:
            e1 = [a for a in exact_int_eigs(o["N"]) if a >= 1]
            e2 = [a for a in int_eigs(o["N"]) if a >= 1]
            if e1 != e2:
                diff += 1
            for a in e1:
                cnt[a] = cnt.get(a, 0) + 1
        print(f"   {name}: {len(ops)} matrices, disagreements {diff}, exact integer eigenvalues >= 1: {dict(sorted(cnt.items()))}")
    rep = json.load(open("twomark_rep_5_3.json"))
    mats = [[[int(Fraction(x)) for x in row] for row in M] for M in rep["mats"]]
    diff = 0
    cnt = {}
    n = 0
    for w in list(itertools.product(range(12), repeat=1)) + list(itertools.product(range(12), repeat=2)):
        M = mats[w[0]]
        for o in w[1:]:
            M = matmul(mats[o], M)
        n += 1
        e1 = [a for a in exact_int_eigs(M) if a >= 1]
        e2 = [a for a in int_eigs(M) if a >= 1]
        if e1 != e2:
            diff += 1
        for a in e1:
            cnt[(len(w), a)] = cnt.get((len(w), a), 0) + 1
    print(f"   twomark_rep_5_3.json: {n} words of length <= 2, disagreements {diff}, (length, eigenvalue): count = {dict(sorted(cnt.items()))}")


A = [[1, -1, 1, 0], [0, 0, 0, 1], [0, 0, 2, -1], [0, 0, 0, 1]]
C = [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 2, 0], [0, 0, 0, 1]]


def closed(k, W):
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


def v2(kmax):
    print(f"(V2) closed form of the stay block (TT_k, W) against the product of the one-vertex matrices a, c, k <= {kmax}")
    bad = 0
    n = 0
    prev = {0: [[1 if i == j else 0 for j in range(4)] for i in range(4)]}
    for k in range(1, kmax + 1):
        cur = {}
        for W in range(1 << k):
            # bit p of W = letter number p + 1 (first new vertex = bit 0); the last letter is bit k - 1
            last = A if (W >> (k - 1)) & 1 else C
            M = matmul(last, prev[W & ((1 << (k - 1)) - 1)])
            cur[W] = M
            n += 1
            if M != closed(k, W):
                bad += 1
                if bad < 4:
                    print("      MISMATCH", k, format(W, f"0{k}b"), M, closed(k, W))
        prev = cur
    print(f"   {n} words checked, mismatches: {bad}")


def ext_h(n, adj):
    """H of every one-vertex extension, by a direct count on the explicit tournament (hcount)."""
    return [int(hcount(n + 1, extend(n, adj, O))) for O in range(1 << n)]


def v3():
    print("(V3) one-vertex tree by direct counts (hcount on explicit extensions), roots = all classes of order <= 6")
    classes = iso_classes(6)
    t0 = time.time()
    memo_both = {}

    def both(n, adj):
        """does (n, adj) have a one-vertex extension with 2H+1 and one with 2H-1?  (h, plus list, minus list)"""
        ev = ext_h(n, adj)
        h = ev[0]
        return h, [O for O in range(1 << n) if ev[O] == 2 * h + 1], [O for O in range(1 << n) if ev[O] == 2 * h - 1]

    k1 = {}
    k2 = {}
    for n in range(1, 7):
        c1 = 0
        c2 = 0
        for key, adj in classes[n].items():
            h, plus, minus = both(n, adj)
            if plus and minus:
                c1 += 1
                okp = False
                for O in plus:
                    a2 = extend(n, adj, O)
                    h2, p2, m2 = both(n + 1, a2)
                    assert h2 == 2 * h + 1
                    if p2 and m2:
                        okp = True
                        break
                okm = False
                for O in minus:
                    a2 = extend(n, adj, O)
                    h2, p2, m2 = both(n + 1, a2)
                    assert h2 == 2 * h - 1
                    if p2 and m2:
                        okm = True
                        break
                if okp and okm:
                    c2 += 1
        k1[n] = c1
        k2[n] = c2
    print(f"   classes in K_1 by order: {k1}   (lane: 0, 1, 1, 2, 5, 18)")
    print(f"   classes in K_2 by order: {k2}   (lane: none)   [{time.time() - t0:.1f} s]")
    # two-vertex steps at order <= 5
    t0 = time.time()
    res = {}
    for n in range(1, 6):
        for key, adj in classes[n].items():
            h = int(hcount(n, adj))
            one = set(ext_h(n, adj))
            two = set()
            for O in range(1 << n):
                a1 = extend(n, adj, O)
                two.update(ext_h(n + 1, a1))
            res.setdefault(n, []).append((h, (2 * h + 1 in one, 2 * h - 1 in one), (2 * h + 1 in two, 2 * h - 1 in two)))
    for n in range(1, 6):
        r = sorted(res[n])
        c1 = sum(1 for (h, o, t) in r if o[0] and o[1])
        c12 = sum(1 for (h, o, t) in r if (o[0] or t[0]) and (o[1] or t[1]))
        gain = [(h, o, t) for (h, o, t) in r if (t[0] and not o[0]) or (t[1] and not o[1])]
        print(f"   order {n}: classes with both targets by one-vertex steps: {c1}; by one- or two-vertex steps: {c12}; "
              f"classes where a two-vertex step reaches a target that no one-vertex step reaches: {gain}")
    print(f"   [{time.time() - t0:.1f} s]")


if __name__ == "__main__":
    kmax2 = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    v1()
    v2(kmax2)
    v3()
