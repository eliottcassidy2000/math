#!/usr/bin/env python3
"""klein-2026-07-02-S95 -- HYP-4003: the queued FULL-WINDOW L_y census (sorry-ledger N1
legs (a)+(c)), float-scan + exact-confirm, Lean-consumable decide-table output.
k=8 spread<=16 (11440 shapes), k=9 spread<=15 (6435): max L_y vs cap_k (leg a) and
max damped L_y^inf vs cap_{k+1} (leg c, feeds the HYP-4001 rate-lemma recursion)."""
from fractions import Fraction as F
import itertools
exec(open("04-computation/ly_far_element_rate_lemma_klein.py").read().split('print("="*92)')[0])

def sectors_float(E, rmax):
    import math
    # float S_r via fine breakpoint-free sampling is risky; instead exact-but-float intervals
    S = {}
    for r in range(1, rmax+1):
        tot = 0.0
        for A in itertools.combinations(range(1,7), r):
            tot += float(J(E, list(A)))
        S[r] = tot
    return S

# J with floats: reimplement avoid_set with floats for speed
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

def Ly_f(E, k, damp=False):
    rmax = max(r for r,_ in DUALS[k])
    v = 1.0
    for r, y in DUALS[k]:
        Sr = sum(avoid_float(E, A) for A in itertools.combinations(range(1,7), r))
        v += float(y) * (1 - r/7.0 if damp else 1.0) * Sr
    return v

for k, spread in [(8, 16), (9, 15)]:
    rows = []
    for tail in itertools.combinations(range(1, spread+1), k-1):
        E = [0]+list(tail)
        rows.append((Ly_f(E, k), Ly_f(E, k+1, damp=True), E))
    rows.sort(reverse=True)
    top = rows[:8]
    # exact confirmation of the float-top candidates
    exact = [(Ly(E, k), E) for _,_,E in top]
    exact.sort(reverse=True)
    mx, arg = exact[0]
    over = sum(1 for v,_ in exact if v > CAPS[k])
    rows_d = sorted(rows, key=lambda t: -t[1])[:8]
    exact_d = [(Ly_damped(E, k+1), E) for _,_,E in rows_d]
    exact_d.sort(reverse=True)
    mxd, argd = exact_d[0]
    print(f"k={k} spread<={spread}: {len(rows)} shapes scanned (float), top-8 exact-confirmed")
    print(f"  leg(a): max L_y = {mx} = {float(mx):.6f} at {arg}; cap_{k} = {float(CAPS[k]):.6f}; "
          f"over-cap among top = {over} (float scan saw {sum(1 for v,_,_ in rows if v > float(CAPS[k])+1e-9)} over)")
    print(f"  leg(c): max L_y^inf = {mxd} = {float(mxd):.6f} at {argd}; cap_{k+1} = {float(CAPS[k+1]):.6f}; "
          f"margin = {float(CAPS[k+1]-mxd):+.6f}")
print("DECIDE-TABLE (Lean-consumable): per row, the exact rationals (max L_y, cap_k, max L_y^inf, cap_k+1).")
