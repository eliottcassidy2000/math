# mac-mini-2026-11-06-microcosm. THM-4900 (6): sw_1 = 2 c3 for all n; exact tower OCF H(Gamma^r T) = sum_j sw_j(T) H(Gamma^r TT_(n-2j)).
import itertools, math, random
from fractions import Fraction
from collections import Counter
import os; exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cover_tower_values.py')).read().split('cores = {')[0])   # forward_counts, P_value, tower_H
def P_coeffs(n, arc):
    """coefficients of P_T(q) via interpolation of exact values at q = 0..n."""
    cnt = forward_counts(n, arc)
    import sympy as sp
    q = sp.symbols('q')
    pts = [(k, P_value(cnt, n, k)) for k in range(n + 1)]
    poly = sp.interpolate(pts, q)
    return sp.Poly(poly, q).all_coeffs()[::-1], cnt   # ascending
def c3(n, arc):
    s = [sum(arc[i]) for i in range(n)]
    return math.comb(n, 3) - sum(math.comb(x, 2) for x in s)
random.seed(11); ok = True; okc3 = True
TT = lambda m: [[i < j for j in range(m)] for i in range(m)]
for trial in range(40):
    n = random.randint(2, 7)
    arc = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5: arc[i][j] = True
            else: arc[j][i] = True
    co, cnt = P_coeffs(n, arc)
    okc3 &= (co[n - 2] == 2 * c3(n, arc))                       # sw_1 = 2 c3, all n
    for r in range(0, 5):
        lhs = tower_H(cnt, n, r)
        rhs = 0
        for jj in range(0, n // 2 + 1):
            d = n - 2 * jj
            if d < 0 or co[d] == 0: continue
            if d == 0: val = tower_H(Counter({0: 1}), 0, r) if False else None
            if d == 0:   # P = 1 (empty core): L^r(1)(1)
                lo, hi = 1 - r, 1 + r; V = {k: 1 for k in range(lo, hi + 1)}
                for i in range(r):
                    lo += 1; hi -= 1; V = {k: k * (V[k + 1] + k * V[k] - V[k - 1]) for k in range(lo, hi + 1)}
                val = V[1]
            else:
                val = tower_H(forward_counts(d, TT(d)), d, r)
            rhs += co[d] * val
        ok &= (lhs == rhs)
print('sw_1 = coefficient of q^(n-2) = 2 c3(T) for 40 random cores of orders 2..7:', okc3)
print('tower OCF  H(Gamma^r T) = sum_j sw_j(T) H(Gamma^r TT_(n-2j))  for r = 0..4:', ok)
