#!/usr/bin/env python3
"""Debt cones with terminal packets (charge_palettes_boundary_20261020g/charge.md conventions).
Odd source n > 1: odd word w = (a_1..a_m), a_i = v_2(3x_i+1), until the first odd iterate h < n (first exit).
z = n-1 coordinates: z_end = lambda z + P - D, lambda = 3^m/2^K, P = sum_{a_i=1} lambda_(i+1:m), D = sum_{a_i>=3} (1-4/2^a_i) lambda_(i+1:m).
Payment ("positive charge pays"): D_w < n - h. Last-step cone (c <= 1): prefix u = w[:-1] has D_u <= lambda_u - 1.
Sharpened last-step cone: D_u <= c_max(a_m)(lambda_u - 1), c_max(a) = (2^a+4)/(2^a-3) (a >= 3), 8 (a = 2), from c(1-mu)+D_a <= 2.
Optimal split: min over splits (u = w[:i], v = w[i:], lambda_v < 1, lambda_u >= 1) of ratio_u (1-lambda_v) + D_v where
ratio_u = D_u/(lambda_u-1) (0 if u empty); covered iff <= 2. Reports counts, 63's data, worst cases, and the
distribution of the minimal split value."""
from fractions import Fraction as F
import sys, json
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 100001
def word(n):
    x = n; w = []
    while True:
        y = 3 * x + 1; a = 0
        while y % 2 == 0: y //= 2; a += 1
        w.append(a); x = y
        if x < n: return w, x
def affine(w):
    # returns lambda, P, D as Fractions for the word w
    lam = F(1); P = F(0); D = F(0)
    for a in w:  # prepend-free: process left to right, constants multiplied by later multipliers
        mu = F(3, 2 ** a)
        lam *= mu; P *= mu; D *= mu
        if a == 1: P += 1
        elif a >= 3: D += 1 - F(4, 2 ** a)
    return lam, P, D
def cmax(a):
    return F(8) if a == 2 else F(2 ** a + 4, 2 ** a - 3)
stats = dict(sources=0, pays=0, last_step_cone=0, sharp_last_step=0, split_covered=0, no_reset=0, three_plus_resets=0)
worst = []; hist = {}; unresolved_sharp = []
for n in range(3, LIMIT + 1, 2):
    w, h = word(n); m = len(w)
    lam, P, D = affine(w)
    stats['sources'] += 1
    resets = sum(1 for a in w if a >= 3)
    if resets == 0: stats['no_reset'] += 1
    if resets >= 3: stats['three_plus_resets'] += 1
    if D < n - h: stats['pays'] += 1
    # prefix data for all splits
    lam_u, P_u, D_u = F(1), F(0), F(0)
    best = None; last_ok = False; sharp_ok = False
    pref = [(F(1), F(0))]
    for i, a in enumerate(w):
        mu = F(3, 2 ** a); lam_u *= mu; D_u *= mu
        if a >= 3: D_u += 1 - F(4, 2 ** a)
        pref.append((lam_u, D_u))
    lam_w = pref[-1][0]
    for i in range(0, m):  # split: u = w[:i], v = w[i:], v nonempty
        lu, du = pref[i]
        lv = lam_w / lu  # lambda_v
        if lv >= 1: continue
        # D_v: debt of the suffix alone
        _, _, dv = affine(w[i:])
        if i == 0:
            val = dv
        else:
            if lu <= 1: continue
            val = du / (lu - 1) * (1 - lv) + dv
        if best is None or val < best[0]: best = (val, i)
        if i == m - 1:
            ratio = (du / (lu - 1)) if lu > 1 else None
            if ratio is not None and ratio <= 1: last_ok = True
            if ratio is not None and ratio <= cmax(w[-1]): sharp_ok = True
    if m == 1: last_ok = sharp_ok = True; stats['single_step'] = stats.get('single_step', 0) + 1
    if last_ok: stats['last_step_cone'] += 1
    if sharp_ok: stats['sharp_last_step'] += 1
    else: unresolved_sharp.append(n)
    if D >= 2: stats['debt_at_least_2'] = stats.get('debt_at_least_2', 0) + 1; print(f'  D_w >= 2: source {n}, D = {float(D):.4f}, n-h = {n-h}, m = {m}, resets = {resets}, word = {w}')
    if best is not None and best[0] <= 2: stats['split_covered'] += 1
    bv = float(best[0]) if best else None
    b = round(bv, 1) if bv is not None else None
    hist[b] = hist.get(b, 0) + 1
    if bv is not None and bv > 1.5: worst.append((bv, n, best[1], m, h))
    if n == 63 or n == 55:
        print(f"source {n}: word {w} (m={m}, resets={resets}), h={h}, lambda={float(lam):.4f}, D={float(D):.4f}, n-h={n-h}, pays={D < n-h}")
        print(f"  last-step ratio D_u/(lam_u-1) = {float(pref[-2][1]/(pref[-2][0]-1)) if pref[-2][0]>1 else None}  (cone c<=1: {last_ok}; sharpened c_max(a_m={w[-1]})={float(cmax(w[-1])):.3f}: {sharp_ok})")
        vals = []
        for i in range(0, m):
            lu, du = pref[i]; lv = lam_w / lu
            if lv >= 1: vals.append((i, None)); continue
            _, _, dv = affine(w[i:])
            val = dv if i == 0 else (du / (lu - 1) * (1 - lv) + dv if lu > 1 else None)
            vals.append((i, None if val is None else round(float(val), 4)))
        print(f"  split values (i, ratio_u(1-lam_v)+D_v): {vals}")
        print(f"  best split {best[1]} value {float(best[0]):.4f}")
print(json.dumps(stats))
worst.sort(reverse=True)
print('worst 15 minimal split values (value, n, split, m, h):', [(round(a, 3), b, c, d, e) for a, b, c, d, e in worst[:15]])
print('histogram of minimal split value (rounded 0.1):', dict(sorted((k, v) for k, v in hist.items() if k is not None)))
print('first 30 sources unresolved by the sharpened last-step cone:', unresolved_sharp[:30], 'count', len(unresolved_sharp))
json.dump({'stats': stats, 'worst': worst[:200], 'unresolved_sharp': unresolved_sharp}, open('cone_splits.json', 'w'))
