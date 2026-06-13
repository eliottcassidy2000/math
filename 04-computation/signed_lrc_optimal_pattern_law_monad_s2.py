#!/usr/bin/env python3
"""
Signed LRC pairwise gap: structural law of the OPTIMAL sign pattern.
monad-explorer-2026-06-06-S2b.

Tests, exhaustively at small n, candidate explanations for which sign pattern
eps maximizes the pairwise gap  G_pair(eps) = max_t min_{i<j} ||(eps_i v_i - eps_j v_j) t||:

  C1  argmax G_pair  ==  argmax min_{i<j} |eps_i v_i - eps_j v_j|
        (maximize the SMALLEST relative-speed magnitude -- avoid small clocks)
  C2  the all-+ pattern is optimal  <=>  no sign flip raises the min |rel speed|
  C3  relation of optimal G_pair to the observer floor 1/n.

Also: distribution of the IMPROVEMENT ratio gmax/allpos, and a head-to-head of
G_pair against the trivial single-speed bound and the per-set max relative speed.
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
    W = [w for w in W if w != 0]
    ms = set(abs(w) for w in W)
    for a, b in combinations(W, 2):
        ms.add(abs(a + b)); ms.add(abs(a - b))
    ms.discard(0)
    T = set()
    for m in ms:
        for a in range(1, 2 * m):
            T.add(F(a, 2 * m))
    return T

def maximin_gap(W):
    W = [w for w in W if w != 0]
    if not W: return F(1, 2)
    best = F(0)
    for t in candidate_times(W):
        m = min(norm(w * t) for w in W)
        if m > best: best = m
    return best

def signed_diff_set(V, eps):
    return [eps[i]*V[i] - eps[j]*V[j] for i, j in combinations(range(len(V)), 2)]

def enum_sets(r, B):
    return [c for c in combinations(range(1, B+1), r) if reduce(gcd, c) == 1]

def main():
    print("="*78)
    print("OPTIMAL SIGN-PATTERN LAW for the signed pairwise gap  monad-S2b")
    print("="*78)
    configs = {3:9, 4:8, 5:7}
    for r in range(3, 6):
        B = configs[r]; n = r+1
        sets = enum_sets(r, B)
        c1_ok = c1_bad = 0       # C1: argmax G == argmax min|rel|
        c2_ok = c2_bad = 0       # C2: all-+ optimal <=> no flip raises min|rel|
        above_floor = at_floor = below_floor = 0
        n_opt_dist = {}
        ratios = []
        floor = F(1, n)
        for V in sets:
            gvals = {}; minrel = {}
            for tail in product([1,-1], repeat=r-1):
                eps = (1,)+tail
                D = signed_diff_set(V, eps)
                gvals[eps] = maximin_gap(D)
                minrel[eps] = min(abs(d) for d in D)
            gmax = max(gvals.values())
            allpos = gvals[tuple([1]*r)]
            ratios.append((gmax, allpos, V))
            opt_set = {e for e,g in gvals.items() if g == gmax}
            n_opt_dist[len(opt_set)] = n_opt_dist.get(len(opt_set),0)+1
            # C1: every min|rel| maximizer is a G_pair maximizer?  (and vv weaker)
            mr_max = max(minrel.values())
            mr_argmax = {e for e,m in minrel.items() if m == mr_max}
            # is at least one min|rel|-maximizer also a G_pair-maximizer?
            if mr_argmax & opt_set:
                c1_ok += 1
            else:
                c1_bad += 1
            # C2: all-+ optimal in G  <=>  all-+ maximizes min|rel|
            allpos_eps = tuple([1]*r)
            g_allpos_opt = (allpos == gmax)
            mr_allpos_opt = (minrel[allpos_eps] == mr_max)
            if g_allpos_opt == mr_allpos_opt:
                c2_ok += 1
            else:
                c2_bad += 1
            # C3
            if gmax > floor: above_floor += 1
            elif gmax == floor: at_floor += 1
            else: below_floor += 1
        print(f"\n--- r={r} (n={n}), {len(sets)} sets, floor 1/n={floor} ---")
        print(f"  C1 (a min|rel|-maximizer is also a G_pair-maximizer): "
              f"{c1_ok}/{len(sets)} hold, {c1_bad} fail")
        print(f"  C2 (all-+ G-optimal  <=>  all-+ min|rel|-optimal):    "
              f"{c2_ok}/{len(sets)} hold, {c2_bad} fail")
        print(f"  C3 max_eps G_pair  vs observer floor 1/n: "
              f"above={above_floor}, equal={at_floor}, below={below_floor}")
        print(f"  #optimal-pattern distribution: {dict(sorted(n_opt_dist.items()))}")
        # extremes of improvement
        ratios.sort(key=lambda x: x[0]/x[1], reverse=True)
        print(f"  biggest improvements gmax/allpos:")
        for gm, ap, V in ratios[:4]:
            print(f"     V={V}: {ap} -> {gm}   (x{float(gm/ap):.3f})")
    # ---- C1 sharper: is min|rel| a sufficient PROXY? show failures explicitly
    print("\n" + "="*78)
    print("C1 FAILURE INSPECTION (r=4, where min|rel| proxy fails)")
    print("="*78)
    r, B = 4, 8
    fails = 0
    for V in enum_sets(r,B):
        gvals={}; minrel={}
        for tail in product([1,-1],repeat=r-1):
            eps=(1,)+tail; D=signed_diff_set(V,eps)
            gvals[eps]=maximin_gap(D); minrel[eps]=min(abs(d) for d in D)
        gmax=max(gvals.values()); mr_max=max(minrel.values())
        opt={e for e,g in gvals.items() if g==gmax}
        mr_arg={e for e,m in minrel.items() if m==mr_max}
        if not (mr_arg & opt) and fails < 8:
            fails+=1
            e0=next(iter(opt)); e1=next(iter(mr_arg))
            print(f"  V={V}: G-opt {e0} G={gvals[e0]} minrel={minrel[e0]} | "
                  f"minrel-opt {e1} G={gvals[e1]} minrel={minrel[e1]}")

if __name__ == "__main__":
    main()
