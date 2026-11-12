#!/usr/bin/env python3
"""Verifier (lane digits), check V9: class spaces of the two-mark system (and of the two-apex family inside it).

Words of length <= 2 in the 12 operations, 17 x 17 integer matrices of twomark_rep_5_3.json (exact for every base:
check V7).  Integer eigenvalues >= 1: EXACT (Berkowitz characteristic polynomial, exact integer root search).
Class spaces and the categories (i) C = 0, (ii) kappa = 0 on C, (iii) kappa != 0 somewhere on C: linear algebra
modulo two 31-bit primes (ranks mod p; the two primes must agree) - a check of the lane's exact rational run.

The two-apex family G(2,0) is the orbit of (TT_2; source, sink) under the four operations with mu = 0, so its pair
statistics are those of the words in these four letters (letter (a, b) of G(2,0) = operation S = (1-a) + 2(1-b)).
"""
import json
import os
import sys
import time
from collections import Counter
from fractions import Fraction

import numpy as np
from numba import njit

from vlin import charpoly_int
from vroots import int_roots

OPS = [(S, mu) for S in range(4) for mu in range(3)]


@njit(cache=True)
def add_row(B, piv, nb, v, p):
    """reduce v by the echelon basis B[:nb] (pivot columns piv), add if nonzero.  Returns new nb (unchanged if dependent)."""
    d = v.shape[0]
    for r in range(nb):
        c = v[piv[r]]
        if c != 0:
            for k in range(d):
                v[k] = (v[k] - c * B[r, k]) % p
    pc = -1
    for k in range(d):
        if v[k] != 0:
            pc = k
            break
    if pc < 0:
        return nb
    # normalise
    inv = 1
    base = v[pc]
    e = p - 2
    while e > 0:
        if e & 1:
            inv = (inv * base) % p
        base = (base * base) % p
        e >>= 1
    for k in range(d):
        v[k] = (v[k] * inv) % p
    for r in range(nb):
        c = B[r, pc]
        if c != 0:
            for k in range(d):
                B[r, k] = (B[r, k] - c * v[k]) % p
    for k in range(d):
        B[nb, k] = v[k]
    piv[nb] = pc
    return nb + 1


@njit(cache=True)
def classify(M1, M2, lam1, lam2, kappa, p):
    """returns (dim V, 1 if kappa in V else 0)."""
    d = M1.shape[0]
    B = np.zeros((d, d), np.int64)
    piv = np.zeros(d, np.int64)
    queue = np.zeros((d + 4, d), np.int64)
    nq = 0
    nb = 0
    v = np.zeros(d, np.int64)
    w = np.zeros(d, np.int64)
    # seeds lambda_i (M_j - I)
    for i in range(2):
        for j in range(2):
            for k in range(d):
                t = 0
                for q in range(d):
                    if i == 0:
                        t = (t + lam1[q] * (M1[q, k] if j == 0 else M2[q, k])) % p
                    else:
                        t = (t + lam2[q] * (M1[q, k] if j == 0 else M2[q, k])) % p
                v[k] = (t - (lam1[k] if i == 0 else lam2[k])) % p
            for k in range(d):
                w[k] = v[k]
            nb2 = add_row(B, piv, nb, w, p)
            if nb2 > nb:
                nb = nb2
                for k in range(d):
                    queue[nq, k] = v[k]
                nq += 1
    while nq > 0 and nb < d:
        nq -= 1
        r = queue[nq].copy()
        for j in range(2):
            for k in range(d):
                t = 0
                for q in range(d):
                    t = (t + r[q] * (M1[q, k] if j == 0 else M2[q, k])) % p
                v[k] = t
            for k in range(d):
                w[k] = v[k]
            nb2 = add_row(B, piv, nb, w, p)
            if nb2 > nb:
                nb = nb2
                for k in range(d):
                    queue[nq, k] = v[k]
                nq += 1
            if nb == d:
                break
    for k in range(d):
        w[k] = kappa[k] % p
    for r in range(nb):
        c = w[piv[r]]
        if c != 0:
            for k in range(d):
                w[k] = (w[k] - c * B[r, k]) % p
    inV = 1
    for k in range(d):
        if w[k] != 0:
            inV = 0
    return nb, inV


def main():
    t0 = time.time()
    here = os.path.dirname(os.path.abspath(__file__))
    rep = json.load(open(here + "/../twomark_rep_5_3.json"))
    d = rep["d"]
    A = [[[int(Fraction(x)) for x in row] for row in M] for M in rep["mats"]]

    def mm(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(d)) for j in range(d)] for i in range(d)]

    names = [(o,) for o in range(12)] + [(o1, o2) for o1 in range(12) for o2 in range(12)]
    mats = [A[w[0]] if len(w) == 1 else mm(A[w[1]], A[w[0]]) for w in names]
    eigs = [[a for a in int_roots(charpoly_int(M)) if a >= 1] for M in mats]
    ec = Counter()
    for w, e in zip(names, eigs):
        for a in e:
            ec[(len(w), a)] += 1
    print(f"words of length <= 2 with an integer eigenvalue >= 1 (exact): {sum(1 for e in eigs if e)}; (length, eigenvalue): count = {dict(sorted(ec.items()))}   [{time.time() - t0:.1f} s]")
    primes = [2147483647, 2147483629]
    L = [1] + [0] * (d - 1)
    inc = [(i, a) for i in range(len(names)) for a in eigs[i]]
    lam = {}
    for (i, a) in inc:
        row = mats[i][0]                       # L M
        lam[(i, a)] = [row[k] - a * L[k] for k in range(d)]
    results = []
    for p in primes:
        Mp = [np.array([[x % p for x in r] for r in M], dtype=np.int64) for M in mats]
        lp = {k2: np.array([x % p for x in v], dtype=np.int64) for k2, v in lam.items()}
        zero = np.zeros(d, np.int64)
        # single class spaces: use the pair routine with M2 = M1, lam2 = lam1 (same seeds twice)
        single = {}
        for (i, a) in inc:
            nb, _ = classify(Mp[i], Mp[i], lp[(i, a)], lp[(i, a)], zero, p)
            single[(i, a)] = d - nb
        keep = [k2 for k2 in inc if single[k2] > 0]
        cat = {"all": Counter(), "len1": Counter(), "stay": Counter()}
        bytype = Counter()
        for x in range(len(keep)):
            i, a1 = keep[x]
            for y in range(x + 1, len(keep)):
                j, a2 = keep[y]
                if i == j:
                    continue
                kap = np.array([((a2 - 1) * u - (a1 - 1) * v) % p for u, v in zip(lam[(i, a1)], lam[(j, a2)])], dtype=np.int64)
                nb, inV = classify(Mp[i], Mp[j], lp[(i, a1)], lp[(j, a2)], kap, p)
                key = "(i) C = 0" if nb == d else ("(ii) commute on all of C" if inV else "(iii) kappa != 0 somewhere on C")
                cat["all"][key] += 1
                if len(names[i]) == 1 and len(names[j]) == 1:
                    cat["len1"][key] += 1
                if all(OPS[o][1] == 0 for o in names[i] + names[j]):
                    cat["stay"][key] += 1
                if key.startswith("(iii)"):
                    bytype[(a1, a2, d - nb)] += 1
        results.append((len(keep), {k2: dict(sorted(v.items())) for k2, v in cat.items()}, dict(sorted(bytype.items()))))
        print(f"prime {p}: (word, a) with C != 0: {len(keep)}   [{time.time() - t0:.1f} s]")
        print(f"   all words of length <= 2: {results[-1][1]['all']}")
        print(f"   one-letter words only:    {results[-1][1]['len1']}")
        print(f"   words in the four mu = 0 letters (two-apex family G(2,0)): {results[-1][1]['stay']}")
        print(f"   (iii) by (a1, a2, dim C): {results[-1][2]}")
    print("the two primes agree:", results[0] == results[1])
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
