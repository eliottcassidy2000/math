#!/usr/bin/env python3
"""Checks for THM-4735 Theorem C: for every first-exit word (sources 3..LIMIT): lambda_w < 1; every split value
V(u,v) = (D_u/(lambda_u-1))(1-lambda_v) + D_v (lambda_u > 1, lambda_v < 1) satisfies V >= D_w; and pays <=> P_w < (1-lambda_w) z."""
from fractions import Fraction as F
import sys
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 100001
def word(n):
    x = n; w = []
    while True:
        y = 3 * x + 1; a = 0
        while y % 2 == 0: y //= 2; a += 1
        w.append(a); x = y
        if x < n: return w, x
def affine(w):
    lam = F(1); P = F(0); D = F(0)
    for a in w:
        mu = F(3, 2 ** a); lam *= mu; P *= mu; D *= mu
        if a == 1: P += 1
        elif a >= 3: D += 1 - F(4, 2 ** a)
    return lam, P, D
bad_lambda = []; bad_V = 0; n_splits = 0; pays_mismatch = 0; cnt = 0; min_gap = None
for n in range(3, LIMIT + 1, 2):
    w, h = word(n); m = len(w); lam, P, D = affine(w); z = n - 1
    cnt += 1
    if lam >= 1: bad_lambda.append((n, w))
    pays = D < n - h; alt = P < (1 - lam) * z
    if pays != alt: pays_mismatch += 1
    pref = [(F(1), F(0))]; lu, du = F(1), F(0)
    for a in w:
        mu = F(3, 2 ** a); lu *= mu; du *= mu
        if a >= 3: du += 1 - F(4, 2 ** a)
        pref.append((lu, du))
    for i in range(1, m):
        lu, du = pref[i]; lv = lam / lu
        if lu <= 1 or lv >= 1: continue
        _, _, dv = affine(w[i:])
        V = du / (lu - 1) * (1 - lv) + dv; n_splits += 1
        if V < D: bad_V += 1
        gap = V - D
        if min_gap is None or gap < min_gap: min_gap = gap
print(f"sources {cnt}: lambda_w >= 1 cases {len(bad_lambda)} {bad_lambda[:5]}; splits checked {n_splits}, V < D_w violations {bad_V}, min (V - D_w) = {float(min_gap) if min_gap is not None else None:.6f}; pays <=> P < (1-lambda) z mismatches {pays_mismatch}")
