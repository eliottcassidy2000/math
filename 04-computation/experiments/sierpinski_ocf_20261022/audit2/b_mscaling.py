# Audit2 (B, C), second independent method: the excess grading by blow-up scaling.
# For M = Shat (n x n, zero diagonal) and the "replica" matrix M' = M (x) J_m / m on n*m vertices, every connected
# cluster of excess e satisfies Z_pi(M') = m^(#blocks - #arcs) Z_pi(M) = m^(-e) Z_pi(M), hence EXACTLY
#     log Psi_(M')(z) = sum_e m^(-e) X_e(M; z),      Psi_M(z) = sum over directed linear forests F of z^|F| prod M_e.
# Psi_(M') is computed without any cluster theory: a forest on [n] x [m] projects to a tuple of walks on [n]; the
# number of injective lifts is prod_v (m)_(c_v) (c_v = visits to v).  So, with W(t) = sum over walks (>= 1 arc) of
# z^len prod S prod t_(visited) and Lambda: t^c -> prod_v (m)_(c_v) (linear),
#     L! n^L m^L [z^L] Psi_(M') = Lambda[ L! [z^L] exp(W) ],
# an integer polynomial in m.  Then log in z with Laurent-polynomial coefficients in m; X_e^(L) = coefficient of m^(-e).
# Multivariate polynomial arithmetic: python-flint fmpz_mpoly (exact integers).
import sys
from fractions import Fraction
from math import factorial
from collections import defaultdict
import numpy as np
import flint

def falling_polys(K):
    """F[k] = coefficient list (index = power of m) of m(m-1)...(m-k+1)"""
    F = [[1]]
    for k in range(1, K + 1):
        prev = F[-1]; new = [0] * (len(prev) + 1)
        for i, c in enumerate(prev):
            new[i + 1] += c          # * m
            new[i] -= (k - 1) * c    # * (-(k-1))
        F.append(new)
    return F

def polymul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] += x * y
    return r

def mscaling(S, Lmax, verbose=True):
    S = [[int(x) for x in row] for row in np.asarray(S)]
    n = len(S)
    ctx = flint.fmpz_mpoly_ctx.get(('t', n), 'lex')
    t = ctx.gens()
    # walk polynomials by endpoint
    P = [t[v] for v in range(n)]
    W = {}
    for l in range(1, Lmax + 1):
        P = [t[w] * sum((S[v][w] * P[v] for v in range(n) if S[v][w] != 0), ctx.from_dict({})) for w in range(n)]
        W[l] = sum(P, ctx.from_dict({}))
        if verbose: print(f"   W_{l}: {len(W[l])} monomials{' (zero, as it must be for odd l)' if l % 2 and W[l] == 0 else ''}", flush=True)
    E = {0: ctx.from_dict({(0,) * n: 1})}
    F = falling_polys(2 * Lmax + 2)
    cache = {}
    def lam_of(exps):
        key = tuple(sorted(e for e in exps if e))
        if key not in cache:
            p = [1]
            for e in key: p = polymul(p, F[e])
            cache[key] = p
        return cache[key]
    psi = {}
    for L in range(1, Lmax + 1):
        acc = ctx.from_dict({})
        for l in range(1, L + 1):
            if W[l] == 0 or E[L - l] == 0: continue
            acc += (factorial(L - 1) // factorial(L - l)) * l * W[l] * E[L - l]
        E[L] = acc                                   # = L! [z^L] exp(W)
        # Lambda
        grouped = defaultdict(int)
        for exps, c in acc.terms():                     # streamed (no dict of all monomials)
            grouped[tuple(sorted(e for e in exps if e))] += int(c)
        poly = [0]
        for key, c in grouped.items():
            if c == 0: continue
            p = lam_of(key)
            if len(p) > len(poly): poly += [0] * (len(p) - len(poly))
            for i, x in enumerate(p): poly[i] += c * x
        den = factorial(L) * n ** L
        # psi_L as Laurent polynomial in m: dict power -> Fraction ; power = i - L
        psi[L] = {i - L: Fraction(c, den) for i, c in enumerate(poly) if c}
        if verbose: print(f"   L={L}: E_L has {len(acc)} monomials; psi_L(m) has powers {sorted(psi[L])}", flush=True)
    # log in z
    lam = {}
    for L in range(1, Lmax + 1):
        s = defaultdict(Fraction)
        for p, c in psi[L].items(): s[p] += c
        for k in range(1, L):
            for p1, c1 in lam[k].items():
                for p2, c2 in psi[L - k].items():
                    s[p1 + p2] -= Fraction(k, L) * c1 * c2
        lam[L] = {p: c for p, c in s.items() if c != 0}
    # X_e^(L) = coefficient of m^(-e)
    return {L: {-p: c for p, c in lam[L].items()} for L in lam}

if __name__ == '__main__':
    sys.path.insert(0, '.')
    from common import rot_A, skew, X0_series, X1_series
    S = skew(rot_A(7))
    X = mscaling(S, 6)
    x0 = X0_series(S, 6); x1 = X1_series(S, 6)
    for L in range(1, 7):
        print(L, {e: float(v) for e, v in sorted(X[L].items())}, X[L].get(0, 0) == x0[L], X[L].get(1, 0) == x1[L])
