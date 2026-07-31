#!/usr/bin/env python3
"""
sfc3_census.py — bounded-exact census for the Strong Factorial Conjecture SFC(3)
(Edo–van den Essen), complementing THM-2812's consecutive-support k=0 theorem.

SFC(3): if f = a z^p + b z^q + c z^r (0 <= p < q < r) and
L(f^{k+1}) = L(f^{k+2}) = L(f^{k+3}) = 0 for some k >= 0 (L(z^n) = n!),
then f = 0.

For fixed (p,q,r,k) the three moments F_{k+1}, F_{k+2}, F_{k+3} are homogeneous
forms in (a,b,c) of degrees k+1,k+2,k+3 with integer coefficients:
    F_m(a,b,c) = sum_{i+j+l=m} m!/(i!j!l!) * (ip+jq+lr)! * a^i b^j c^l.

CERTIFICATE (rigorous, one-way): let D = (k+1)+(k+2)+(k+3) - 2 = 3k+4.
If the multiplication maps {x^alpha F_m : |alpha| = D - deg F_m} span the FULL
space of degree-D forms over F_prime (rank = C(D+2,2)), then they span over Q,
so the ideal contains m^{>=D} and the projective variety V(F1,F2,F3) in P^2 is
EMPTY over C.  Since any common zero with a coordinate = 0 would violate the
PROVED SFC(2) (or hit the never-zero pure-power axis moments), emptiness is
exactly "SFC(3) holds for this (p,q,r,k)".  Rank deficiency mod prime is
inconclusive: escalate to other primes, and report cells that stay deficient.

Output: verdict per cell over the box p<q<r <= RMAX, k <= KMAX.
"""
import sys
from math import comb, factorial
from itertools import combinations

PRIMES = [1000000007, 998244353, 754974721]  # < 2^31 so int64 products are safe

def moment_form(p, q, r, m):
    """coefficients dict {(i,j,l): int} of F_m."""
    out = {}
    for i in range(m + 1):
        for j in range(m + 1 - i):
            l = m - i - j
            coef = factorial(m) // (factorial(i) * factorial(j) * factorial(l))
            out[(i, j, l)] = coef * factorial(i * p + j * q + l * r)
    return out

def monomials(deg):
    return [(i, j, deg - i - j) for i in range(deg + 1) for j in range(deg + 1 - i)]

def rank_mod(rows, ncols, prime):
    """Vectorized Gaussian elimination rank mod prime (numpy object-free):
    use int64 with prime < 2^31 so products fit; reduce eagerly."""
    import numpy as np
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

def cell_verdict(p, q, r, k):
    degs = [k + 1, k + 2, k + 3]
    D = sum(degs) - 2
    target = comb(D + 2, 2)
    mons_D = {m: idx for idx, m in enumerate(monomials(D))}
    rows = []
    for m in degs:
        F = moment_form(p, q, r, m)
        for alpha in monomials(D - m):
            row = [0] * target
            for (i, j, l), cval in F.items():
                mono = (i + alpha[0], j + alpha[1], l + alpha[2])
                row[mons_D[mono]] += cval
            rows.append(row)
    for prime in PRIMES:
        if rank_mod(rows, target, prime) == target:
            return "EMPTY(certified mod %d)" % prime, True
    return "DEFICIENT(all primes) -- needs exact analysis", False

def main():
    RMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    bad = []
    total = 0
    for (p, q, r) in combinations(range(RMAX + 1), 3):
        for k in range(KMAX + 1):
            total += 1
            verdict, ok = cell_verdict(p, q, r, k)
            if not ok:
                bad.append((p, q, r, k))
                print(f"(p,q,r,k)=({p},{q},{r},{k}): {verdict}", flush=True)
    print(f"census done: {total} cells, certified-empty {total - len(bad)}, "
          f"flagged {len(bad)}: {bad}", flush=True)

if __name__ == "__main__":
    main()
