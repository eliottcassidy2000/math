# The census tower of Paley tournaments P_p (p = 3 mod 4 prime): which subtournament counts are polynomials in p and
# which carry arithmetic?  3-census: fixed by regularity; 4-census: fixed by double regularity (= g^2 = -p);
# 5-census: TT5 = p * TT4(T[QR]) by vertex-transitivity (a TT5 has a unique source), computed exactly here for all
# p < Pmax and tested against a polynomial in p.  mac-mini-2026-10-25-pyth.
import sys
from math import comb
import numpy as np
from sympy import primerange

def paley(p):
    sq = {(x * x) % p for x in range(1, p)}
    return np.array([[1 if (j - i) % p in sq else 0 for j in range(p)] for i in range(p)], dtype=np.int64), sq

def t3(A, X):
    """transitive triples in T[X] = sum over u of C(outdeg in X, 2)"""
    sub = A[np.ix_(X, X)]; d = sub.sum(axis=1)
    return int(sum(comb(int(k), 2) for k in d))

def tt4(A, X):
    """transitive 4-sets in T[X]: unique source w, the rest a transitive triple in N+(w) within X"""
    sub = A[np.ix_(X, X)]; tot = 0
    for a in range(len(X)):
        Y = [X[b] for b in range(len(X)) if sub[a, b]]
        tot += t3(A, Y)
    return tot

Pmax = int(sys.argv[1]) if len(sys.argv) > 1 else 160
rows = []
for p in primerange(7, Pmax):
    if p % 4 != 3: continue
    A, sq = paley(p)
    QR = sorted(sq)
    TT5 = p * tt4(A, QR)
    TT4 = p * t3(A, QR)
    rows.append((p, TT4, TT5))
    n = p; m = (n - 1) // 2
    print(f"p={p:4d}: TT4 = {TT4} (DRT formula n(n-1)(n-3)(n-7)/64 = {n * (n - 1) * (n - 3) * (n - 7) // 64}), TT5 = {TT5}", flush=True)
# polynomial test for TT5(p)/p: degree-4 polynomial in p? fit on all but the last 3, test on the last 3
ps = np.array([r[0] for r in rows], dtype=float); y = np.array([r[2] / r[0] for r in rows], dtype=float)
for deg in (4,):
    c = np.polyfit(ps[:-3], y[:-3], deg)
    pred = np.polyval(c, ps)
    print(f"degree-{deg} fit of TT5/p on p <= {int(ps[-4])}: max |residual| in-sample {np.max(np.abs(y[:-3] - pred[:-3])):.3e}, "
          f"out-of-sample {np.abs(y[-3:] - pred[-3:])}")
# exact candidate: the random-like expectation for a doubly regular tournament is not known a priori; print residuals of
# the best rational quartic (least squares) times 1024 to look for an arithmetic error term
c = np.polyfit(ps, y, 4); res = (y - np.polyval(c, ps))
print("TT5/p minus best quartic (residual):", [f"{r:+.3f}" for r in res])
