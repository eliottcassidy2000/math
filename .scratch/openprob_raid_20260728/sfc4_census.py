#!/usr/bin/env python3
"""
sfc4_census.py — certified census for the Strong Factorial Conjecture SFC(4):
f = a z^p + b z^q + c z^r + e z^s (0 <= p<q<r<s), window k:
L(f^{k+1}) = ... = L(f^{k+4}) = 0 => f = 0 ?

Same Macaulay full-graded-rank certificate as THM-2836, now in P^3:
degrees (k+1..k+4), D = sum(deg) - 3 = 4k+7; full rank of the degree-D
multiplication matrix mod a 30-bit prime => ideal contains m^{>=D} =>
projective variety empty over C.  Degenerate zeros (some coordinate zero)
are 3-term solutions of a 4-window, hence of a 3-window, excluded by the
THM-2836 SFC(3) census provided supports <= 12 and k+1 <= 8 — so stay in
that envelope: supports <= 8, k <= 2 here (also the matrix-size sweet spot).

Usage: python3 sfc4_census.py [SMAX] [KMAX]
"""
import sys
from math import comb, factorial
from itertools import combinations
import numpy as np

PRIMES = [1000000007, 998244353, 754974721]

def monomials4(deg):
    out = []
    for i in range(deg + 1):
        for j in range(deg + 1 - i):
            for l in range(deg + 1 - i - j):
                out.append((i, j, l, deg - i - j - l))
    return out

def moment_form4(p, q, r, s, m):
    out = {}
    fm = factorial(m)
    for (i, j, l, t) in monomials4(m):
        coef = fm // (factorial(i) * factorial(j) * factorial(l) * factorial(t))
        out[(i, j, l, t)] = coef * factorial(i * p + j * q + l * r + t * s)
    return out

def rank_mod(rows, ncols, prime):
    P = prime
    mat = np.array([[v % P for v in row] for row in rows], dtype=np.int64)
    nrows = mat.shape[0]
    rank = 0
    col = 0
    while rank < nrows and col < ncols:
        pivs = np.nonzero(mat[rank:, col])[0]
        if pivs.size == 0:
            col += 1
            continue
        piv = rank + int(pivs[0])
        if piv != rank:
            mat[[rank, piv]] = mat[[piv, rank]]
        inv = pow(int(mat[rank, col]), P - 2, P)
        mat[rank] = (mat[rank] * inv) % P
        colvals = mat[:, col].copy()
        colvals[rank] = 0
        nz = np.nonzero(colvals)[0]
        if nz.size:
            mat[nz] = (mat[nz] - colvals[nz, None] * mat[rank][None, :]) % P
        rank += 1
        col += 1
    return rank

def cell(p, q, r, s, k):
    degs = [k + 1, k + 2, k + 3, k + 4]
    D = sum(degs) - 3
    monsD = {m: i for i, m in enumerate(monomials4(D))}
    target = len(monsD)
    rows = []
    for m in degs:
        F = moment_form4(p, q, r, s, m)
        for alpha in monomials4(D - m):
            row = [0] * target
            for (i, j, l, t), cv in F.items():
                mono = (i + alpha[0], j + alpha[1], l + alpha[2], t + alpha[3])
                row[monsD[mono]] += cv
            rows.append(row)
    for prime in PRIMES:
        if rank_mod(rows, target, prime) == target:
            return True
    return False

def main():
    SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    bad = []
    tot = 0
    for (p, q, r, s) in combinations(range(SMAX + 1), 4):
        for k in range(KMAX + 1):
            tot += 1
            if not cell(p, q, r, s, k):
                bad.append((p, q, r, s, k))
                print("FLAGGED", (p, q, r, s, k), flush=True)
        print(f"done support ({p},{q},{r},{s})", flush=True)
    print(f"SFC4 census: {tot} cells, certified-empty {tot - len(bad)}, "
          f"flagged {len(bad)}: {bad}", flush=True)

if __name__ == "__main__":
    main()
