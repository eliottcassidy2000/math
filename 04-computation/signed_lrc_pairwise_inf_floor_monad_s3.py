#!/usr/bin/env python3
"""
Signed LRC pairwise gap: the INF floor inf_S Gstar(S) vs n, tight-set census,
and the LRC-reduction lower bound.   monad-explorer-2026-06-06-S3.

Builds on THM-426 (sign patterns = cuts of K_{n-1}) and HYP-2293 (REFUTED:
Gstar can be < 1/n; V=(2,3,4,6,8) -> 3/19 at n=6).  The live handoff (T764):

   Gstar(S) = max over cuts (A,B) of  M( W(A,B) ),
   W(A,B) = { |v_i-v_j| : same side } u { v_i+v_j : across },
   M(W) = max_t min_{w in W} ||w t||   (the LRC loneliness of the relative-
          speed SET; only distinct values matter).

Questions:
  (Q1) inf_S Gstar(S) as a function of n (over gcd-1 sets, speeds <= B).
       Track the MINIMIZER(s).  Does it scale ~ c/n, ~ c/n^2, or ->0?
  (Q2) Characterize the TIGHT sets Gstar(S) = 1/n exactly.  Are they
       consecutive-block / pairwise worry-sets?
  (Q3) LOWER BOUND: Gstar(S) = max_cut M(W) >= 1/(2 r_min) where r_min =
       min over cuts of #distinct(W) (union/measure bound, UNCONDITIONAL).
       Also report 1/(r_min+1) (the LRC-conjecture bound).  How do these
       compare to the actual inf?

Also computes the "FULL" mutual variant including the OBSERVER (speed 0):
   W_full(A,B) = {v_i : all i}  u  W(A,B)
(observer-mover pairs give the bare speed |v_i|, sign-invariant).
Gstar_full(S) <= min(M_obs(S), Gstar(S)).  This is "all n runners mutually
lonely with the best signs".
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd
from functools import reduce

def norm(x):
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def candidate_times(W):
    """Local maxima of t->min_w ||w t|| occur only at:
       single-tent peaks  t=(2a+1)/(2w)   and
       tent crossings     t=a/(w_i +- w_j).
    This is the exact (minimal-ish) candidate set."""
    W = sorted(set(abs(w) for w in W if w != 0))
    T = set()
    for w in W:                       # single-tent peaks (odd / 2w)
        for a in range(0, w):
            T.add(F(2 * a + 1, 2 * w))
    for x, y in combinations(W, 2):   # tent crossings integer/(x+-y)
        for d in (x + y, x - y):
            if d == 0:
                continue
            for a in range(1, abs(d)):
                T.add(F(a, abs(d)))
    return T

def maximin_gap(W):
    """max_t min_{w in W} ||w t|| exactly.  Only distinct |w| matter.
    Early-termination prune: speeds sorted ascending; break once a norm
    falls at/below the running best (cannot improve)."""
    Ws = sorted(set(abs(w) for w in W if w != 0))
    if not Ws:
        return F(1, 2), None
    best = F(0)
    bestt = None
    for t in candidate_times(Ws):
        m = F(1, 2)
        ok = True
        for w in Ws:
            nv = norm(w * t)
            if nv <= best:
                ok = False
                break
            if nv < m:
                m = nv
        if ok and m > best:
            best = m; bestt = t
    return best, bestt

def cut_speeds(V, side):
    """side: tuple in {0,1}^r giving the side of each mover. Return W multiset."""
    r = len(V)
    W = []
    for i, j in combinations(range(r), 2):
        if side[i] == side[j]:
            W.append(abs(V[i] - V[j]))
        else:
            W.append(V[i] + V[j])
    return W

def gstar(V, include_obs=False):
    """max over cuts of M(W). Returns (gstar, best_side, best_W_distinct)."""
    r = len(V)
    best = F(-1); best_side = None; best_t = None
    # fix side[0]=0 to kill global flip symmetry
    for tail in product([0, 1], repeat=r - 1):
        side = (0,) + tail
        W = cut_speeds(V, side)
        if include_obs:
            W = list(V) + W
        g, t = maximin_gap(W)
        if g > best:
            best = g; best_side = side; best_t = t
    return best, best_side, best_t

def r_min_cut(V):
    """min over cuts of #distinct relative speeds (for the LRC lower bound)."""
    r = len(V)
    best = None
    for tail in product([0, 1], repeat=r - 1):
        side = (0,) + tail
        W = set(w for w in cut_speeds(V, side) if w != 0)
        if best is None or len(W) < best:
            best = len(W)
    return best

def enum_sets(r, B):
    out = []
    for combo in combinations(range(1, B + 1), r):
        if reduce(gcd, combo) == 1:
            out.append(combo)
    return out

def is_consecutive_block(V):
    return all(V[i+1] - V[i] == 1 for i in range(len(V)-1))

def main():
    print("=" * 78)
    print("SIGNED-LRC PAIRWISE INF FLOOR  inf_S Gstar(S)  vs n   (monad-S3, exact)")
    print("=" * 78)
    configs = {2: 16, 3: 14, 4: 11, 5: 9, 6: 8}   # r:B  (n=r+1)
    for r in range(2, 7):
        B = configs[r]
        n = r + 1
        floor = F(1, n)
        sets = enum_sets(r, B)
        inf_g = F(2); inf_sets = []
        tight = []          # Gstar == 1/n
        below_floor = []    # Gstar < 1/n
        # also track the lower-bound diagnostics on the minimizer
        for V in sets:
            g, side, t = gstar(V)
            if g < inf_g:
                inf_g = g; inf_sets = [(V, side, t)]
            elif g == inf_g:
                inf_sets.append((V, side, t))
            if g == floor:
                tight.append(V)
            if g < floor:
                below_floor.append((V, g))
        print(f"\n--- n={n} (r={r} movers), B<={B}, {len(sets)} gcd-1 sets, floor 1/n={floor}={float(floor):.4f} ---")
        print(f"  inf_S Gstar = {inf_g} = {float(inf_g):.5f}   (n*inf = {float(n*inf_g):.4f})")
        for V, side, t in inf_sets[:6]:
            A = [V[i] for i in range(r) if side[i] == 0]
            Bs = [V[i] for i in range(r) if side[i] == 1]
            rmin = r_min_cut(V)
            print(f"     minimizer V={V}: best cut A={A}|B={Bs}, t*={t}; "
                  f"r_min={rmin}, LB 1/(2 r_min)={F(1,2*rmin)}={float(F(1,2*rmin)):.4f}, "
                  f"1/(r_min+1)={float(F(1,rmin+1)):.4f}")
        print(f"  #below floor (Gstar<1/n): {len(below_floor)};  #tight (Gstar=1/n): {len(tight)}")
        if tight:
            ncb = sum(1 for V in tight if is_consecutive_block(V))
            print(f"     tight sets: {len(tight)} total, {ncb} are consecutive blocks")
            for V in tight[:10]:
                tag = "CONSEC" if is_consecutive_block(V) else ""
                print(f"        {V} {tag}")
        if below_floor:
            below_floor.sort(key=lambda x: x[1])
            print(f"     lowest below-floor sets:")
            for V, g in below_floor[:8]:
                print(f"        {V}: Gstar={g}={float(g):.4f}")

if __name__ == "__main__":
    main()
