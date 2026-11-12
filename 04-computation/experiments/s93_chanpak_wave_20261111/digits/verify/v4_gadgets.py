#!/usr/bin/env python3
"""Verifier (lane digits), check V4: the gadget system of THM-5130 (A) (claims D5, D6), own code.

(1) gadget law s+ = G s by direct counts (own DP) for the 8 gadget matrices with an integer eigenvalue >= 1 and a
    random sample of the others, on random bases; the sink operation sigma by direct counts.
(2) exact integer eigenvalues >= 1 of all 3,207 gadget matrices (no floating point), class spaces of the 8, all pairs.
(3) words: alphabet sigma + gadgets (X of order <= NMAX), ALL words of length <= LEN - nothing omitted (the lane's
    run for NMAX = 4 dropped the letter sigma itself and sigma^2, its run for NMAX = 3 dropped sigma^2, sigma^3).
    Exact integer eigenvalues >= 1, class spaces for the readouts L_h and L_1, all pairs of (word, eigenvalue).

Usage: v4_gadgets.py
"""
import json
import random
import sys
import time
from collections import Counter

import numpy as np
from numba import njit

from vlin import charpoly_int, class_space, dot, matmul, rowmat
from vroots import int_roots

SIGMA = [[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]]


def eig_ge1(M):
    return [a for a in int_roots(charpoly_int(M)) if a >= 1]


@njit(cache=True)
def end_table(n, adj):
    full = 1 << n
    E = np.zeros((full, n), np.int64)
    for v in range(n):
        E[1 << v, v] = 1
    for m in range(1, full):
        for v in range(n):
            c = E[m, v]
            if c == 0:
                continue
            for u in range(n):
                if (adj[v] >> u) & 1 and not (m >> u) & 1:
                    E[m | (1 << u), u] += c
    return E


@njit(cache=True)
def pair_state_direct(n, adj, O, s):
    """(r, M, h, E) of the pair (B, O) from the definitions of THM-5130."""
    full = (1 << n) - 1
    E = end_table(n, adj)
    radj = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(n):
            if (adj[i] >> j) & 1:
                radj[j] |= 1 << i
    R = end_table(n, radj)             # R[mask, v]: paths of B[mask] starting at v
    h = 0
    EO = 0
    for v in range(n):
        h += E[full, v]
        if (O >> v) & 1:
            EO += E[full, v]
    M = 0
    W = 0
    for S in range(0, full):          # S = vertex set of P1 (possibly empty), P2 on the complement, starts in O
        comp = full ^ S
        st = 0
        for f in range(n):
            if (O >> f) & 1 and (comp >> f) & 1:
                st += R[comp, f]
        if st == 0:
            continue
        if S == 0:
            M += st
        else:
            hs = 0
            eo = 0
            for e in range(n):
                hs += E[S, e]
                if (O >> e) & 1:
                    eo += E[S, e]
            M += hs * st
            W += eo * st
    s[0] = h - W
    s[1] = M
    s[2] = h
    s[3] = EO


def closing(nB, adjB, O, nX, adjX):
    n = nB + nX + 1
    a = np.zeros(n, np.int64)
    xmask = ((1 << nX) - 1) << nB
    v = n - 1
    for i in range(nB):
        a[i] = int(adjB[i]) | xmask
        if not (O >> i) & 1:
            a[i] |= 1 << v
    for i in range(nX):
        a[nB + i] = (int(adjX[i]) << nB) | (1 << v)
    a[v] = O
    return n, a


def rand_tournament(n, rng):
    a = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                a[i] |= 1 << j
            else:
                a[j] |= 1 << i
    return a


def words_run(gad, maxlen, t0):
    letters = [("s", SIGMA)] + [(f"g{i}", g["G"]) for i, g in enumerate(gad)]
    words = {}
    cur = {(): [[1 if i == j else 0 for j in range(4)] for i in range(4)]}
    ec = Counter()
    nw = 0
    for ln in range(1, maxlen + 1):
        nxt = {}
        for w, M in cur.items():
            for name, G in letters:
                nxt[w + (name,)] = matmul(G, M)
        cur = nxt
        for w, M in cur.items():
            nw += 1
            e = eig_ge1(M)
            if e:
                words[w] = (M, e)
                for a in e:
                    ec[(ln, a)] += 1
        print(f"      length {ln}: {len(cur)} words; words with an integer eigenvalue >= 1 so far {len(words)}  [{time.time() - t0:.1f} s]", flush=True)
    print(f"      (length, eigenvalue): count = {dict(sorted(ec.items()))}")
    for Lname, L in (("L_h", [0, 0, 1, 0]), ("L_1", [1, 2, 0, 0])):
        keep = []
        for w, (M, e) in words.items():
            for a in e:
                V, lam = class_space(L, [M], [a])
                if V.dim() < 4:
                    keep.append((w, a))
        cat = Counter()
        catsig = Counter()
        free = []
        for x in range(len(keep)):
            w1, a1 = keep[x]
            for y in range(x + 1, len(keep)):
                w2, a2 = keep[y]
                if w1 == w2:
                    continue
                V, lam = class_space(L, [words[w1][0], words[w2][0]], [a1, a2])
                sigpow = all(c == "s" for c in w1) or all(c == "s" for c in w2)
                if V.dim() == 4:
                    key = "(i) C = 0"
                else:
                    kappa = [(a2 - 1) * p - (a1 - 1) * q for p, q in zip(lam[0], lam[1])]
                    if V.contains(kappa):
                        key = "(ii) commute on all of C"
                    else:
                        key = "(iii) kappa != 0 somewhere on C"
                        free.append((w1, a1, w2, a2, 4 - V.dim()))
                cat[key] += 1
                if sigpow:
                    catsig[key] += 1
        print(f"      readout {Lname}: (word, a) with C != 0: {len(keep)}; pairs: {dict(sorted(cat.items()))}")
        print(f"          of these, pairs in which one word is a power of sigma (not all examined by the lane): {dict(sorted(catsig.items()))}")
        for f in free[:12]:
            print("          (iii):", f)


def main():
    t0 = time.time()
    here = sys.path[0]
    gad6 = json.load(open(here + "/../gadget_matrices_6.json"))
    mats = [g["G"] for g in gad6]
    eigs = [eig_ge1(G) for G in mats]
    special = [i for i in range(len(mats)) if eigs[i]]
    print(f"(2) {len(mats)} gadget matrices (X of order <= 6); with an integer eigenvalue >= 1 (exact): {len(special)}: "
          f"{[(i, gad6[i]['n'], gad6[i]['adj'], format(gad6[i]['P'], 'b'), eigs[i]) for i in special]}   [{time.time() - t0:.1f} s]")
    L = [0, 0, 1, 0]
    for i in special:
        for a in eigs[i]:
            V, lam = class_space(L, [mats[i]], [a])
            mu = [x - y for x, y in zip(rowmat(L, mats[i]), L)]
            print(f"      #{i}: a = {a}, dim C = {4 - V.dim()}, offset identically 0 on C: {V.contains(lam[0])}, H fixed on C: {V.contains(mu)}")
    cat = Counter()
    for x in range(len(special)):
        for y in range(x + 1, len(special)):
            i, j = special[x], special[y]
            for a1 in eigs[i]:
                for a2 in eigs[j]:
                    V, lam = class_space(L, [mats[i], mats[j]], [a1, a2])
                    if V.dim() == 4:
                        cat["(i) C = 0"] += 1
                        continue
                    kappa = [(a2 - 1) * p - (a1 - 1) * q for p, q in zip(lam[0], lam[1])]
                    cat["(ii) commute on all of C" if V.contains(kappa) else "(iii) kappa != 0 somewhere on C"] += 1
    print(f"      pairs among them: {dict(cat)}")
    # (1) direct checks of the law
    rng = random.Random(77)
    s = np.zeros(4, np.int64)
    s2 = np.zeros(4, np.int64)
    fails = 0
    nt = 0
    todo = [(i, 12) for i in special] + [(rng.randrange(len(mats)), 1) for _ in range(250)]
    for i, reps in todo:
        g = gad6[i]
        for _ in range(reps):
            nB = rng.randint(1, 5)
            adjB = rand_tournament(nB, rng)
            O = rng.randint(1, (1 << nB) - 1)
            pair_state_direct(nB, adjB, O, s)
            nn, a = closing(nB, adjB, O, g["n"], np.array(g["adj"], np.int64))
            pair_state_direct(nn, a, g["P"] << nB, s2)
            pred = [sum(g["G"][r][q] * int(s[q]) for q in range(4)) for r in range(4)]
            nt += 1
            if pred != [int(x) for x in s2]:
                fails += 1
    sf = 0
    for _ in range(200):
        nB = rng.randint(1, 7)
        adjB = rand_tournament(nB, rng)
        O = rng.randint(1, (1 << nB) - 1)
        pair_state_direct(nB, adjB, O, s)
        a = np.zeros(nB + 1, np.int64)
        for i in range(nB):
            a[i] = int(adjB[i]) | (1 << nB)
        pair_state_direct(nB + 1, a, O, s2)
        if [int(s[0]), 2 * int(s[1]), int(s[2]), 0] != [int(x) for x in s2]:
            sf += 1
    print(f"(1) gadget law by direct counts (own DP): {nt} (base, gadget) tests incl. 12 per special gadget, failures {fails}; "
          f"sink operation sigma: 200 tests, failures {sf}   [{time.time() - t0:.1f} s]")
    # (3) words
    for nmax, maxlen in ((4, 2), (3, 3)):
        gad = json.load(open(here + f"/../gadget_matrices_{nmax}.json"))
        print(f"(3) alphabet sigma + {len(gad)} gadget matrices (X of order <= {nmax}), all words of length <= {maxlen}:")
        words_run(gad, maxlen, t0)
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
