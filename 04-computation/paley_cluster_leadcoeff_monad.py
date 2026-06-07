#!/usr/bin/env python3
"""
paley_cluster_leadcoeff_monad.py
monad-explorer-2026-06-07 (deep-research, 4th session)

PIN DOWN the TRUE leading coefficient of A_{2k}/p^{k+1}, exactly, to settle
whether it is the Catalan number C_k (prior claim, THM-438) or the value
S'_k = (-1)^k * sum_{even cacti} mu  =  6, 25, 132 ... (this session's census).

Method: A_{2k} is computed EXACTLY at each prime via the flow-Moebius identity
    A_{2k} = sum_{set partitions rho of {0..2k}} mu(0,rho) * M_rho,
    M_rho = (-1)^k p^{V-k} F(rho),   F(rho) = sum_{F_p flows} prod chi(t_e),
which is fast (Bell(2k+1) partitions, each a small flow sum) -- NO p^{2k}
brute force.  We then fit  A_{2k}/p^{k+1} = c0 + c1/p + c2/p^2 + ...  by
solving the Vandermonde system on many Paley primes, getting c0 to high
precision.  Compare c0 to C_k and to S'_k.
"""
import math
from itertools import product as iproduct
from fractions import Fraction
import numpy as np


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def set_partitions(c):
    c = list(c)
    if len(c) == 1:
        yield [c]
        return
    f = c[0]
    for sm in set_partitions(c[1:]):
        for i, s in enumerate(sm):
            yield sm[:i] + [[f] + s] + sm[i + 1:]
        yield [[f]] + sm


def mu_partition(blocks):
    m = 1
    for B in blocks:
        b = len(B)
        m *= ((-1) ** (b - 1)) * math.factorial(b - 1)
    return m


def build_graph(blocks, L):
    pos2blk = {}
    for bi, B in enumerate(blocks):
        for pos in B:
            pos2blk[pos] = bi
    return [(pos2blk[i], pos2blk[i + 1]) for i in range(L)], len(blocks)


def nullspace_basis(edges, nb, p):
    E = len(edges)
    A = np.zeros((nb, E), dtype=np.int64)
    for ei, (u, v) in enumerate(edges):
        A[v, ei] += 1
        A[u, ei] -= 1
    A %= p
    rows, cols = A.shape
    pivots = []
    r = 0
    for cc in range(cols):
        piv = next((rr for rr in range(r, rows) if A[rr, cc] % p != 0), None)
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, cc]), p - 2, p)) % p
        for rr in range(rows):
            if rr != r and A[rr, cc] % p != 0:
                A[rr] = (A[rr] - A[rr, cc] * A[r]) % p
        pivots.append(cc)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        vec = np.zeros(cols, dtype=np.int64)
        vec[fc] = 1
        for ri, pc in enumerate(pivots):
            vec[pc] = (-A[ri, fc]) % p
        basis.append(vec)
    return np.array(basis, dtype=np.int64) if basis else np.zeros((0, cols), dtype=np.int64)


def flow_sum(edges, nb, p):
    basis = nullspace_basis(edges, nb, p)
    m = basis.shape[0]
    E = len(edges)
    chi = np.array([legendre(z, p) for z in range(p)], dtype=np.int64)
    if m == 0:
        return (1 if E == 0 else 0)
    # iterate flows; for m up to ~ (2k - V + 1) this is p^m -- keep p moderate
    grids = np.array(np.meshgrid(*[range(p)] * m, indexing='ij')).reshape(m, -1).T
    T = (grids @ basis) % p
    V = chi[T]
    return int(V.prod(axis=1).sum())


def A_via_moebius(k, p, partitions):
    L = 2 * k
    tot = 0
    sign = (-1) ** k
    for blocks in partitions:
        mu = mu_partition(blocks)
        if mu == 0:
            continue
        edges, nb = build_graph(blocks, L)
        if any(u == v for (u, v) in edges):
            continue
        fs = flow_sum(edges, nb, p)
        if fs == 0:
            continue
        tot += mu * sign * (p ** (nb - k)) * fs
    return tot


def catalan(k):
    return math.comb(2 * k, k) // (k + 1)


# Paley primes p = 3 mod 4
def paley_primes(lo, hi):
    out = []
    for p in range(lo, hi + 1):
        if p % 4 != 3:
            continue
        if all(p % d for d in range(2, int(p ** 0.5) + 1)):
            out.append(p)
    return out


def fit_leading(ks_vals, ncoef):
    """ks_vals: list of (p, value=A/p^{k+1}); fit value = sum_{j=0}^{ncoef-1} c_j / p^j."""
    P = [p for p, _ in ks_vals]
    Y = [v for _, v in ks_vals]
    # use the ncoef largest primes
    idx = sorted(range(len(P)), key=lambda i: -P[i])[:ncoef]
    Pn = [P[i] for i in idx]
    Yn = [Y[i] for i in idx]
    M = [[Fraction(1, Pn[r] ** j) for j in range(ncoef)] for r in range(ncoef)]
    # solve M c = Yn over rationals (Yn are exact rationals)
    Yn = [Fraction(y).limit_denominator(10 ** 12) for y in Yn]
    # gaussian elimination
    A = [row[:] + [Yn[r]] for r, row in enumerate(M)]
    n = ncoef
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        inv = A[col][col]
        A[col] = [x / inv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][c] - f * A[col][c] for c in range(n + 1)]
    return [A[r][n] for r in range(n)]


for k in [2, 3, 4]:
    L = 2 * k
    print("=" * 70)
    print(f"k={k}: A_{2*k}/p^{k+1} leading coefficient   (C_k={catalan(k)})")
    parts = list(set_partitions(range(L + 1)))
    # choose primes: keep p^m small. m <= k+1. Use moderate primes.
    if k == 2:
        primes = paley_primes(7, 120)
    elif k == 3:
        primes = paley_primes(11, 59)        # max cycle rank m=3, p^3 manageable
    else:
        primes = paley_primes(11, 31)        # k=4: max m=4, p^4<=10^6
    data = []
    for p in primes:
        A = int(A_via_moebius(k, p, parts))
        val = Fraction(A, p ** (k + 1))
        data.append((p, val))
        print(f"   p={p:4d}   A/p^{k+1} = {float(val):.6f}")
    # fit with increasing #coeffs to see c0 stabilize
    print(f"   --- extrapolated leading coeff c0 (Vandermonde on largest primes) ---")
    for nc in range(3, min(7, len(data)) + 1):
        c = fit_leading(data, nc)
        print(f"     using {nc} primes:  c0 = {float(c[0]):.6f}"
              + (f"  ~ {c[0].limit_denominator(50)}" if hasattr(c[0], 'limit_denominator') else ""))
    print(f"   Catalan C_{k} = {catalan(k)}    S'_k=(-1)^k*sum_ec mu (this-session census) for ref")
