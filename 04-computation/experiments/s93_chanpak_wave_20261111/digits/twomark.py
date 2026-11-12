#!/usr/bin/env python3
"""Lane "digits", part (1), two marks: the one-vertex rule system with two marks (z1, z2).

Marked tournament (T; z1, z2), z1 != z2.  Operation (S, mu): a new vertex v beats z_i for i in S (S a subset of
{1, 2}) and loses to every other vertex; then mu = 0: marks stay; mu = 1: z1 := v; mu = 2: z2 := v.   12 operations.
(All new vertices are "below" everything except the current marks; this is the class of rule sets with a
finite-dimensional state space, PROOFS.txt P2.)

This script finds the minimal linear representation numerically-exactly (Hankel method):
  bases  = every (T; z1, z2) with T of order <= NB (one representative per isomorphism class of T, all ordered pairs),
  suffix = every word of length <= LW in the 12 operations,
  Hankel matrix  Hk[base, word] = H(word(base)); its rank d (mod a 61-bit prime and over Q on the chosen minor)
  is the dimension of the smallest closed state space; the state of a base is (H(w_k(base)))_k for d suffix words
  w_k (the first one empty, so H is the first coordinate); every operation acts by a rational d x d matrix, solved
  on d bases and verified on all the others.
Output: twomark_rep.json (words w_k, matrices as fractions).

Usage: twomark.py NB LW
"""
import itertools
import json
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit

from tour import hcount, iso_classes

OPS = [(S, mu) for S in range(4) for mu in range(3)]     # S bitmask: bit0 -> beats z1, bit1 -> beats z2


def apply_word(n, adj, z1, z2, word):
    """returns (n', adj') after the word (list of op indices); marks tracked internally."""
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
    return N, a


def H_of(n, adj, z1, z2, word):
    N, a = apply_word(n, adj, z1, z2, word)
    return int(hcount(N, a))


P = (1 << 61) - 1


def rank_mod(M):
    """rank of an integer matrix mod P and the pivot columns (python ints, small sizes)."""
    A = [[x % P for x in row] for row in M]
    nr, nc = len(A), len(A[0])
    piv = []
    r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if A[i][c]), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = pow(A[r][c], P - 2, P)
        A[r] = [(x * inv) % P for x in A[r]]
        for i in range(nr):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(x - f * y) % P for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == nr:
            break
    return r, piv


def solve_exact(A, b):
    """solve A x = b (A square invertible, integers) over Q."""
    n = len(A)
    M = [[Fraction(x) for x in row] + [Fraction(y)] for row, y in zip(A, b)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]


def main():
    nb = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    lw = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    t0 = time.time()
    classes = iso_classes(nb)
    bases = []
    for n in range(2, nb + 1):
        for key, adj in classes[n].items():
            for z1 in range(n):
                for z2 in range(n):
                    if z1 != z2:
                        bases.append((n, adj, z1, z2))
    words = [()]
    for ln in range(1, lw + 1):
        words.extend(itertools.product(range(12), repeat=ln))
    print(f"{len(bases)} bases (orders 2..{nb}), {len(words)} suffix words (length <= {lw})", flush=True)
    Hk = [[H_of(n, adj, z1, z2, w) for w in words] for (n, adj, z1, z2) in bases]
    print(f"Hankel matrix computed  [{time.time() - t0:.1f} s]", flush=True)
    # rank growth with the word length
    for ln in range(0, lw + 1):
        ncol = sum(12 ** k for k in range(ln + 1))
        r, piv = rank_mod([row[:ncol] for row in Hk])
        print(f"   rank with suffix words of length <= {ln} ({ncol} columns): {r}")
    d, pivcols = rank_mod(Hk)
    # rank growth with the base order
    for nn in range(2, nb + 1):
        rows = [Hk[i] for i, b in enumerate(bases) if b[0] <= nn]
        r, _ = rank_mod(rows)
        print(f"   rank with bases of order <= {nn} ({len(rows)} rows): {r}")
    print(f"Hankel rank d = {d}; state words: {[''.join('%x' % o for o in words[c]) or 'empty' for c in pivcols]}")
    assert pivcols[0] == 0
    # d independent bases (rows) for the chosen columns
    sub = [[Hk[i][c] for c in pivcols] for i in range(len(bases))]
    _, pivrows = rank_mod([list(col) for col in zip(*sub)])
    assert len(pivrows) == d
    Ssq = [sub[i] for i in pivrows]
    # matrices of the operations: state(o(B))_k = H(w_k(o(B))) = H((o w_k)(B))
    mats = []
    bad = 0
    for o in range(12):
        rows = []
        for k, c in enumerate(pivcols):
            w = (o,) + tuple(words[c])
            tgt_all = [H_of(n, adj, z1, z2, w) for (n, adj, z1, z2) in bases]
            coef = solve_exact(Ssq, [tgt_all[i] for i in pivrows])
            for i in range(len(bases)):
                if sum(coef[j] * sub[i][j] for j in range(d)) != tgt_all[i]:
                    bad += 1
            rows.append(coef)
        mats.append(rows)
        print(f"   op {o} = (S = {OPS[o][0]:02b}, mu = {OPS[o][1]}): matrix solved, verification failures so far {bad}  [{time.time() - t0:.1f} s]", flush=True)
    dens = sorted(set(x.denominator for M in mats for row in M for x in row))
    print(f"denominators: {dens}; verification failures over all {len(bases)} bases: {bad}")
    json.dump({"d": d, "state_words": [list(words[c]) for c in pivcols],
               "mats": [[[str(x) for x in row] for row in M] for M in mats],
               "base_states": [[sub[i][j] for j in range(d)] for i in range(len(bases))]},
              open(f"twomark_rep_{nb}_{lw}.json", "w"))
    print(f"-> twomark_rep_{nb}_{lw}.json   [{time.time() - t0:.1f} s]")


if __name__ == "__main__":
    main()
