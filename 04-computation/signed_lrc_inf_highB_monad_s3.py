#!/usr/bin/env python3
"""
Signed-LRC pairwise gap: FOCUSED high-B search for inf_S Gstar(S).
monad-explorer-2026-06-06-S3.

Speedup: SET-LEVEL early-out.  We track the running minimum Gstar over sets.
For a new set, as soon as ANY cut gives M >= running_inf, this set cannot lower
the inf, so we abandon it (skip remaining cuts).  Only genuine inf-candidates
pay the full 2^{r-1}-cut cost.  This lets us reach larger B.

Goal: is n * inf_S Gstar(S) bounded below (Theta(1/n)) or decaying (toward 1/n^2)?
Reports inf and the minimizer at each (n, B); compare across B to see if the inf
has stabilized (true inf) or is still dropping with B.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce
import sys

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def candidate_times(Ws):
    T = set()
    for w in Ws:
        for a in range(0, w):
            T.add(F(2*a+1, 2*w))
    for x, y in combinations(Ws, 2):
        for d in (x+y, x-y):
            if d == 0: continue
            d = abs(d)
            for a in range(1, d):
                T.add(F(a, d))
    return T

def maximin_at_least(Ws, thresh):
    """Return (M, reached) where reached=True once M > thresh (strict): then this
    cut already exceeds the running inf, enough to abandon the set."""
    best = F(0)
    for t in candidate_times(Ws):
        ok = True; m = F(1, 2)
        for w in Ws:
            nv = norm(w*t)
            if nv <= best:
                ok = False; break
            if nv < m: m = nv
        if ok and m > best:
            best = m
            if best > thresh:
                return best, True
    return best, False

def cut_Ws(V, side):
    r = len(V); W = []
    for i, j in combinations(range(r), 2):
        W.append(abs(V[i]-V[j]) if side[i]==side[j] else V[i]+V[j])
    return sorted(set(w for w in W if w != 0))

def gstar_capped(V, inf_so_far):
    """Gstar(V), but abandon early once a cut reaches inf_so_far."""
    r = len(V); best = F(0)
    for tail in product([0,1], repeat=r-1):
        side = (0,)+tail
        Ws = cut_Ws(V, side)
        m, reached = maximin_at_least(Ws, inf_so_far)
        if m > best: best = m
        if reached:
            return best, False   # Gstar > inf_so_far: not a candidate
    return best, True            # fully evaluated; best is the true Gstar (<= inf)

def enum_sets(r, B):
    return [c for c in combinations(range(1,B+1), r) if reduce(gcd,c)==1]

def search(n, B):
    r = n-1
    sets = enum_sets(r, B)
    inf_g = F(1, 2); inf_sets = []
    for V in sets:
        g, is_cand = gstar_capped(V, inf_g)
        if is_cand and g < inf_g:
            inf_g = g; inf_sets = [V]
        elif is_cand and g == inf_g:
            inf_sets.append(V)
    return inf_g, inf_sets, len(sets)

def main():
    # (n, [B values to compare convergence])
    plan = [(6, [10, 13, 16]), (7, [9, 11, 13]), (8, [9, 11])]
    if len(sys.argv) > 1:
        # custom: n B
        n = int(sys.argv[1]); B = int(sys.argv[2])
        g, sets, ns = search(n, B)
        print(f"n={n} B<={B} ({ns} sets): inf Gstar={g}={float(g):.6f}  n*inf={float(n*g):.4f}")
        for V in sets[:8]: print(f"   minimizer {V}")
        return
    print("="*78)
    print("SIGNED-LRC FOCUSED inf SEARCH (high B, set-level pruning)  monad-S3")
    print("="*78)
    for n, Bs in plan:
        print(f"\n=== n={n} (floor 1/n={float(F(1,n)):.5f}) ===")
        for B in Bs:
            g, sets, ns = search(n, B)
            flag = "FLOOR-BREAK" if g < F(1,n) else ""
            print(f"  B<={B:2d} ({ns:5d} sets): inf={g}={float(g):.6f}  n*inf={float(n*g):.4f}  {flag}")
            print(f"       minimizers ({len(sets)}): {sets[:6]}")
            sys.stdout.flush()

if __name__ == "__main__":
    main()
