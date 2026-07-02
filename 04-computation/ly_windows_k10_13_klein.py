#!/usr/bin/env python3
"""klein-2026-07-02-S96 -- HYP-4004: N1 leg (a) completed at k=10..13 (full windows,
float scan + exact top-8 confirm; spread caps per THM-534 conventions)."""
from fractions import Fraction as F
import itertools
exec(open("04-computation/ly_far_element_rate_lemma_klein.py").read().split('print("="*92)')[0])
def avoid_float(E, A):
    bad = []
    for e in E:
        if e == 0: continue
        for j in A:
            lo, hi = j/7.0, (j+1)/7.0
            for a in range(e):
                bad.append(((a+lo)/e, (a+hi)/e))
    bad.sort()
    tot, cl, ch = 0.0, None, None
    for l, h in bad:
        if ch is None: cl, ch = l, h
        elif l <= ch: ch = max(ch, h)
        else: tot += ch - cl; cl, ch = l, h
    if ch is not None: tot += ch - cl
    return 1.0 - tot
def Ly_f(E, k):
    v = 1.0
    for r, y in DUALS[k]:
        v += float(y) * sum(avoid_float(E, A) for A in itertools.combinations(range(1,7), r))
    return v
for k, spread in [(10, 14), (11, 13), (12, 13), (13, 13)]:
    rows = []
    for tail in itertools.combinations(range(1, spread+1), k-1):
        E = [0]+list(tail)
        rows.append((Ly_f(E, k), E))
    rows.sort(reverse=True)
    exact = sorted(((Ly(E, k), E) for _, E in rows[:8]), reverse=True)
    mx, arg = exact[0]
    overf = sum(1 for v, _ in rows if v > float(CAPS[k]) + 1e-9)
    print(f"k={k} spread<={spread}: {len(rows)} shapes; max L_y = {mx} = {float(mx):.6f} at {arg}; "
          f"cap = {float(CAPS[k]):.6f}; float-over-cap = {overf}; exact max <= cap: {mx <= CAPS[k]}")
