#!/usr/bin/env python3
"""Second, independent implementation of the Koszul homology of (d^n f)_{n>=0}, f = x_0^N,
on C[x_0, x_1, ...] (wt x_j = j+2), by weight.  Different basis enumeration (itertools) and
mod-P ranks (two primes), plus an explicit test of whether the product s_{2N+3} * s_{2N+4} of the
two lowest H_1 classes is a boundary.   Usage: python koszul_check2.py N Wmax
"""
import sys
from itertools import combinations
import flint

N = int(sys.argv[1])
Wmax = int(sys.argv[2])
PR = [2147483647, 2147483629]


def parts(D, mx=None):
    if mx is None:
        mx = D
    if D == 0:
        return [()]
    out = []
    for k in range(min(D, mx), 1, -1):
        for r in parts(D - k, k):
            out.append((k,) + r)
    return out


def mono_mul(a, b):
    return tuple(sorted(a + b, reverse=True))


def dpoly(p):
    out = {}
    for m, c in p.items():
        for i, k in enumerate(m):
            l = list(m)
            l[i] = k + 1
            nm = tuple(sorted(l, reverse=True))
            out[nm] = out.get(nm, 0) + c * (k - 1)
    return {m: c for m, c in out.items() if c}


df = [{(2,) * N: 1}]
for _ in range(Wmax):
    df.append(dpoly(df[-1]))
PB = {w: parts(w) for w in range(Wmax + 1)}


def Kbasis(i, w):
    out = []
    maxn = w
    gens = [n for n in range(0, w + 1) if 2 * N + n <= w]
    for S in combinations(gens, i):
        ws = sum(2 * N + n for n in S)
        if ws <= w:
            for m in PB[w - ws]:
                out.append((S, m))
    return out


def bd(i, w):
    tgt = Kbasis(i - 1, w)
    idx = {b: j for j, b in enumerate(tgt)}
    rows = []
    for S, m in Kbasis(i, w):
        row = {}
        for pos, n in enumerate(S):
            sg = -1 if pos % 2 else 1
            rest = S[:pos] + S[pos + 1:]
            for mm, c in df[n].items():
                key = idx[(rest, mono_mul(m, mm))]
                row[key] = row.get(key, 0) + sg * c
        rows.append(row)
    return rows, len(tgt)


def rk(rows, ncols, P):
    if not rows or ncols == 0:
        return 0
    M = flint.nmod_mat(len(rows), ncols, P)
    for i, r in enumerate(rows):
        for j, c in r.items():
            M[i, j] = c % P
    return M.rank()


for P in PR:
    H1, H2 = [], []
    for w in range(Wmax + 1):
        d1, n0 = bd(1, w) if Kbasis(1, w) else ([], 0)
        d2, n1 = bd(2, w) if Kbasis(2, w) else ([], 0)
        d3, n2 = bd(3, w) if Kbasis(3, w) else ([], 0)
        k1, k2 = len(Kbasis(1, w)), len(Kbasis(2, w))
        r1 = rk(d1, len(PB[w]), P)
        r2 = rk(d2, k1, P)
        r3 = rk(d3, k2, P)
        H1.append(k1 - r1 - r2)
        H2.append(k2 - r2 - r3)
    print(f"N={N} P={P}: H_1 by weight 0..{Wmax}: {H1}")
    print(f"N={N} P={P}: H_2 by weight 0..{Wmax}: {H2}")
    fH2 = next((w for w in range(Wmax + 1) if H2[w]), None)
    print(f"   first nonzero H_2 weight: {fH2}   (4N+7 = {4 * N + 7})")
