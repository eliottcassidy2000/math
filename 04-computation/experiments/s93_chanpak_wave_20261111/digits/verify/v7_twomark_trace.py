#!/usr/bin/env python3
"""Verifier (lane digits), check V7: is the 17-dimensional two-mark representation (twomark_rep_5_3.json) EXACT for
every two-marked tournament, and not only on the bases of order <= 5 on which it was interpolated?

Trace lemma (path surgery, as in the one-mark law).  (T; z1, z2), |T| >= 2, z1 != z2; a word w adds vertices X; the
only arcs from X to T end at z1 or z2.  A Hamiltonian path of T u X meets T in at most three maximal segments: at
most one initial segment (free start), the others entered from X, hence starting at z1 or at z2.  The admissible
completions on X depend only on the kind of each segment (initial / entered at z1 / entered at z2), their order, and
the end type of each segment (z1, z2, other).  So  H(T.w) = sum_k c_w(k) t_k(T)  with the 25 trace counts t_k below.

Consequence: a linear identity between functionals T -> H(T.w) that holds on a set of bases whose trace vectors span
Q^25 holds for EVERY base.  This script computes
 (a) the rank of the trace vectors of all two-marked tournaments of order <= 5 (and <= 6);
 (b) a numerical test of the trace lemma: H(T.w) for all 157 words of length <= 2 is a linear function of t(T)
     (solved on 25 independent bases, verified on all the others);
 (c) the identities x(T.o) = A_o x(T) (17 state words, 12 operations) by own direct counts on all bases of order <= 5.
"""
import itertools
import json
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit

sys.path.insert(0, sys.path[0] + "/..")
from tour import iso_classes  # noqa: E402

from vlin import Space  # noqa: E402

OPS = [(S, mu) for S in range(4) for mu in range(3)]


@njit(cache=True)
def path_table(n, adj):
    full = 1 << n
    P = np.zeros((full, n, n), np.int64)
    for v in range(n):
        P[1 << v, v, v] = 1
    for m in range(1, full):
        for s in range(n):
            if not (m >> s) & 1:
                continue
            for e in range(n):
                c = P[m, s, e]
                if c == 0:
                    continue
                for u in range(n):
                    if (adj[e] >> u) & 1 and not (m >> u) & 1:
                        P[m | (1 << u), s, u] += c
    return P


@njit(cache=True)
def trace(n, adj, z1, z2, t):
    P = path_table(n, adj)
    full = (1 << n) - 1
    for k in range(25):
        t[k] = 0
    # ---- one segment
    h = 0
    for s in range(n):
        for e in range(n):
            h += P[full, s, e]
        t[0] += P[full, s, z1]
        t[1] += P[full, s, z2]
    t[2] = h - t[0] - t[1]
    t[3] = P[full, z1, z2]
    t[5] = P[full, z2, z1]
    for e in range(n):
        if e != z2:
            t[4] += P[full, z1, e]
        if e != z1:
            t[6] += P[full, z2, e]
    # ---- (initial I, entered Z_i)
    for i in range(2):
        zi = z1 if i == 0 else z2
        zj = z2 if i == 0 else z1
        base = 7 + 5 * i
        for A in range(1, full + 1):
            if (A >> zi) & 1:
                continue
            B = full ^ A
            Izj = 0
            Io = 0
            for s in range(n):
                for e in range(n):
                    c = P[A, s, e]
                    if c:
                        if e == zj:
                            Izj += c
                        else:
                            Io += c
            Zzi = 1 if B == (1 << zi) else 0
            Zzj = 0
            Zo = 0
            for e in range(n):
                c = P[B, zi, e]
                if c and e != zi:
                    if e == zj:
                        Zzj += c
                    else:
                        Zo += c
            t[base + 0] += Izj * Zzi
            t[base + 1] += Izj * Zo
            t[base + 2] += Io * Zzi
            t[base + 3] += Io * Zzj
            t[base + 4] += Io * Zo
    # ---- (entered Z_1, entered Z_2) and (I, Z_1, Z_2)
    for A in range(1, full + 1):
        if not (A >> z1) & 1 or (A >> z2) & 1:
            continue
        a_s = 1 if A == (1 << z1) else 0
        a_o = 0
        for e in range(n):
            if e != z1:
                a_o += P[A, z1, e]
        rest = full ^ A
        B = rest
        while B > 0:
            if (B >> z2) & 1:
                b_s = 1 if B == (1 << z2) else 0
                b_o = 0
                for e in range(n):
                    if e != z2:
                        b_o += P[B, z2, e]
                C = rest ^ B
                if C == 0:
                    t[17] += a_s * b_s
                    t[18] += a_s * b_o
                    t[19] += a_o * b_s
                    t[20] += a_o * b_o
                else:
                    ic = 0
                    for s in range(n):
                        for e in range(n):
                            ic += P[C, s, e]
                    t[21] += ic * a_s * b_s
                    t[22] += ic * a_s * b_o
                    t[23] += ic * a_o * b_s
                    t[24] += ic * a_o * b_o
            B = (B - 1) & rest


@njit(cache=True)
def hcount(n, adj):
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
    tot = 0
    for v in range(n):
        tot += E[full - 1, v]
    return tot


def H_word(n, adj, z1, z2, word):
    N = n + len(word)
    a = np.zeros(N, np.int64)
    for i in range(n):
        a[i] = int(adj[i])
    cur = n
    for o in word:
        S, mu = OPS[o]
        v = cur
        out = 0
        if S & 1:
            out |= 1 << z1
        if S & 2:
            out |= 1 << z2
        a[v] = out
        for i in range(cur):
            if not (out >> i) & 1:
                a[i] |= 1 << v
        if mu == 1:
            z1 = v
        elif mu == 2:
            z2 = v
        cur += 1
    return int(hcount(N, a))


def solve_rows(rows, rhs_list):
    """rows: 25 independent integer vectors (as equations c . row = rhs); returns list of solutions c (Fractions)."""
    n = len(rows)
    M = [[Fraction(x) for x in r] + [Fraction(rh[i]) for rh in rhs_list] for i, r in enumerate(rows)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [[M[i][n + k] for i in range(n)] for k in range(len(rhs_list))]


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    t0 = time.time()
    classes = iso_classes(nmax)
    bases = []
    for n in range(2, nmax + 1):
        for key, adj in classes[n].items():
            for z1 in range(n):
                for z2 in range(n):
                    if z1 != z2:
                        bases.append((n, np.array(adj, np.int64), z1, z2))
    tv = np.zeros(25, np.int64)
    T = []
    for (n, adj, z1, z2) in bases:
        trace(n, adj, z1, z2, tv)
        T.append([int(x) for x in tv])
    # sanity: H = t0 + t1 + t2
    assert all(T[i][0] + T[i][1] + T[i][2] == H_word(*bases[i], ()) for i in range(0, len(bases), 7))
    V = Space(25)
    pivrows = []
    rank_by_order = {}
    for i, (n, adj, z1, z2) in enumerate(bases):
        if V.add(T[i]):
            pivrows.append(i)
        rank_by_order[n] = V.dim()
    print(f"(a) two-marked bases of orders 2..{nmax}: {len(bases)}; rank of the 25 trace counts by order: {rank_by_order}   [{time.time() - t0:.1f} s]")
    d = V.dim()
    small = [i for i in range(len(bases)) if bases[i][0] <= 5]
    # (b) trace lemma: H(T.w) linear in t(T) for all words of length <= 2 (bases of order <= 5)
    words = [()] + [(o,) for o in range(12)] + list(itertools.product(range(12), repeat=2))
    if d == 25:
        piv = [i for i in pivrows if bases[i][0] <= 5]
        if len(piv) == 25:
            rhs = [[H_word(*bases[i], w) for i in piv] for w in words]
            sols = solve_rows([T[i] for i in piv], rhs)
            bad = 0
            for w, c in zip(words, sols):
                for i in small:
                    if sum(c[k] * T[i][k] for k in range(25)) != H_word(*bases[i], w):
                        bad += 1
            print(f"(b) H(T.w) = c_w . t(T): {len(words)} words of length <= 2, coefficients solved on 25 bases, verified on all "
                  f"{len(small)} bases of order <= 5: failures {bad}   [{time.time() - t0:.1f} s]")
        else:
            print(f"(b) skipped: only {len(piv)} independent trace vectors among the bases of order <= 5")
    # (c) the identities of the representation by own direct counts
    rep = json.load(open(sys.path[0] + "/twomark_rep_5_3.json"))
    A = [[[int(Fraction(x)) for x in row] for row in M] for M in rep["mats"]]
    sw = [tuple(w) for w in rep["state_words"]]
    bad = 0
    ntest = 0
    for i in small:
        b = bases[i]
        x = [H_word(*b, w) for w in sw]
        for o in range(12):
            for k in range(17):
                ntest += 1
                if H_word(*b, (o,) + sw[k]) != sum(A[o][k][j] * x[j] for j in range(17)):
                    bad += 1
    print(f"(c) x(T.o) = A_o x(T), own direct counts, all {len(small)} two-marked bases of order <= 5: {ntest} identities, failures {bad}   [{time.time() - t0:.1f} s]")
    print(f"time {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
