#!/usr/bin/env python3
"""
Signed-LRC pairwise gap: VERIFY the r_min lower bounds and REFINE inf_S Gstar.
monad-explorer-2026-06-06-S3.

THEOREM candidate (THM-427):
  Let r_min(S) = min over cuts (A,B) of #distinct nonzero relative speeds in
  W(A,B) = {|v_i-v_j| same side} u {v_i+v_j across}.  Then
    (i)  Gstar(S) >= 1/(2 r_min)            [UNCONDITIONAL, measure/union bound]
    (ii) Gstar(S) >= 1/(r_min + 1)          [conditional on LRC for r_min runners;
                                             unconditional when r_min <= 7]
  COROLLARY: r_min(S) <= n-1  =>  Gstar(S) >= 1/n   (the observer floor survives).
             So Gstar < 1/n  =>  r_min >= n  (every cut has >= n distinct rel speeds).

This script:
  A) verifies (i) for ALL gcd-1 sets (must hold; sanity on the measure bound),
  B) verifies (ii) for all sets with r_min <= 7 (LRC-proven range),
  C) verifies the COROLLARY: every set with r_min <= n-1 has Gstar >= 1/n,
     and reports whether EVERY below-floor set has r_min >= n,
  D) refines inf_S Gstar at higher speed bound B (min-tracking only, fast).
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def candidate_times(Ws):
    T = set()
    for w in Ws:
        for a in range(0, w):
            T.add(F(2 * a + 1, 2 * w))
    for x, y in combinations(Ws, 2):
        for d in (x + y, x - y):
            if d == 0: continue
            d = abs(d)
            for a in range(1, d):
                T.add(F(a, d))
    return T

def maximin_gap(W):
    Ws = sorted(set(abs(w) for w in W if w != 0))
    if not Ws:
        return F(1, 2)
    best = F(0)
    for t in candidate_times(Ws):
        ok = True; m = F(1, 2)
        for w in Ws:
            nv = norm(w * t)
            if nv <= best:
                ok = False; break
            if nv < m: m = nv
        if ok and m > best:
            best = m
    return best

def cut_W(V, side):
    r = len(V); W = []
    for i, j in combinations(range(r), 2):
        W.append(abs(V[i]-V[j]) if side[i]==side[j] else V[i]+V[j])
    return W

def gstar_and_rmin(V):
    r = len(V)
    best = F(-1); rmin = None
    for tail in product([0,1], repeat=r-1):
        side = (0,)+tail
        W = cut_W(V, side)
        k = len(set(w for w in W if w != 0))
        if rmin is None or k < rmin: rmin = k
        g = maximin_gap(W)
        if g > best: best = g
    return best, rmin

def enum_sets(r, B):
    return [c for c in combinations(range(1,B+1), r) if reduce(gcd,c)==1]

def main():
    print("="*78)
    print("SIGNED-LRC r_min LOWER BOUNDS + refined inf   (monad-S3)")
    print("="*78)
    configs = {2:18, 3:16, 4:13, 5:11, 6:10, 7:9}   # r:B
    for r in range(2,8):
        B = configs[r]; n = r+1; floor = F(1,n)
        sets = enum_sets(r,B)
        bound_i_fail = 0; bound_ii_fail = 0; cor_fail = 0
        below = []; inf_g = F(2); inf_set = None
        rmin_le_nm1_count = 0
        belowfloor_rmin_ge_n = 0; belowfloor_total = 0
        max_rmin = 0
        for V in sets:
            g, rmin = gstar_and_rmin(V)
            max_rmin = max(max_rmin, rmin)
            if g < F(1, 2*rmin): bound_i_fail += 1
            if rmin <= 7 and g < F(1, rmin+1): bound_ii_fail += 1
            if rmin <= n-1:
                rmin_le_nm1_count += 1
                if g < floor: cor_fail += 1
            if g < floor:
                belowfloor_total += 1
                if rmin >= n: belowfloor_rmin_ge_n += 1
                below.append((V,g,rmin))
            if g < inf_g: inf_g = g; inf_set = (V,rmin)
        print(f"\n--- n={n} (r={r}), B<={B}, {len(sets)} sets, floor=1/{n} ---")
        print(f"  inf_S Gstar = {inf_g} = {float(inf_g):.5f}  (n*inf={float(n*inf_g):.4f}); "
              f"minimizer {inf_set[0]} r_min={inf_set[1]}")
        print(f"  max_S r_min = {max_rmin}  (C(r,2)={r*(r-1)//2})")
        print(f"  (i)  Gstar>=1/(2 r_min): {len(sets)-bound_i_fail}/{len(sets)} ok "
              f"({'PASS' if bound_i_fail==0 else f'FAIL {bound_i_fail}'})")
        print(f"  (ii) Gstar>=1/(r_min+1) for r_min<=7: "
              f"{'PASS' if bound_ii_fail==0 else f'FAIL {bound_ii_fail}'}")
        print(f"  COROLLARY r_min<=n-1 => Gstar>=1/n: {rmin_le_nm1_count} sets have r_min<=n-1, "
              f"{'PASS' if cor_fail==0 else f'FAIL {cor_fail}'}")
        print(f"  below-floor sets: {belowfloor_total}; with r_min>=n: {belowfloor_rmin_ge_n} "
              f"({'ALL have r_min>=n' if belowfloor_total==belowfloor_rmin_ge_n else 'SOME r_min<n!'})")
        if below:
            below.sort(key=lambda x:x[1])
            for V,g,rmin in below[:5]:
                print(f"      {V}: Gstar={g}={float(g):.4f}, r_min={rmin}")

if __name__ == "__main__":
    main()
