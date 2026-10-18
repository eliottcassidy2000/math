#!/usr/bin/env python3
"""Derived arc space of the fat point X = Spec C[x]/(x^N) (wt x = 2): the Koszul complex of f_n = d^n(x_0^N), n >= 0,
over R = C[x_0, x_1, ...] (wt x_j = j+2).  Weight by weight, compute dim H_0 (= arc space), dim H_1, and the Euler
characteristic, and compare with 1/prod_{m=2}^{2N-1}(1-q^m).  H_1 = ker(K_1 -> K_0)/im(K_2 -> K_1).
Ranks mod P (flint)."""
import sys, functools
import flint
from vir_ss import parts_ge2, deriv
P = 2147483647
N, Wmax = int(sys.argv[1]), int(sys.argv[2])

@functools.lru_cache(maxsize=None)
def monos(w):  # monomials of weight w in R (partitions of w into parts >= 2)
    return tuple(sorted(parts_ge2(w, w))) if w >= 0 else ()

def mul(a, b):  # monomial product (partitions)
    return tuple(sorted(a + b, reverse=True))

# f_n as dict partition -> coeff
f = [{(2,) * N: 1}]
while 2 * N + len(f) - 1 < Wmax:
    f.append({k: v % P for k, v in deriv(f[-1]).items()})
wt = lambda n: 2 * N + n

def K1_basis(w):   # pairs (n, monomial m) with wt(m) + wt(n) = w
    return [(n, m) for n in range(len(f)) if wt(n) <= w for m in monos(w - wt(n))]

def K2_basis(w):
    return [(n1, n2, m) for n1 in range(len(f)) for n2 in range(n1 + 1, len(f)) if wt(n1) + wt(n2) <= w
            for m in monos(w - wt(n1) - wt(n2))]

def rank(rows, ncols):
    if not rows or not ncols: return 0
    M = flint.nmod_mat(len(rows), ncols, [x for r in rows for x in r], P)
    return M.rank()

out = []
for w in range(0, Wmax + 1):
    R0 = monos(w); idx0 = {m: i for i, m in enumerate(R0)}
    B1 = K1_basis(w); idx1 = {b: i for i, b in enumerate(B1)}
    # d1: e_n * m -> f_n * m
    rows1 = []
    for (n, m) in B1:
        row = [0] * len(R0)
        for mono, co in f[n].items():
            row[idx0[mul(mono, m)]] = (row[idx0[mul(mono, m)]] + co) % P
        rows1.append(row)
    rk1 = rank(rows1, len(R0))
    # d2: e_n1 ^ e_n2 * m -> f_n1 * e_n2 m - f_n2 * e_n1 m
    B2 = K2_basis(w)
    rows2 = []
    for (n1, n2, m) in B2:
        row = [0] * len(B1)
        for mono, co in f[n1].items():
            j = idx1[(n2, mul(mono, m))]; row[j] = (row[j] + co) % P
        for mono, co in f[n2].items():
            j = idx1[(n1, mul(mono, m))]; row[j] = (row[j] - co) % P
        rows2.append(row)
    rk2 = rank(rows2, len(B1))
    H0 = len(R0) - rk1
    H1 = (len(B1) - rk1) - rk2
    out.append((w, len(R0), H0, H1))
# Euler characteristic target
targ = [0] * (Wmax + 1); targ[0] = 1
for m in range(2, 2 * N):
    for i in range(m, Wmax + 1): targ[i] += targ[i - m]
print(f"N={N}: weight, dim R, dim H0 (arc space), dim H1, [1/prod_(m=2)^(2N-1)(1-q^m)] coefficient")
for (w, r, h0, h1), t in zip(out, targ):
    print(f"  w={w:2d}  R={r:4d}  H0={h0:4d}  H1={h1:4d}  H0-H1={h0-h1:4d}  target={t}")
