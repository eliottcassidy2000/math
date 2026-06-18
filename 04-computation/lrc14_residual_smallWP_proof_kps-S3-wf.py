#!/usr/bin/env python3
"""
LRC(14) S3 -- a FULLY RIGOROUS lemma for the small-W_P regime, plus the bounded
reduction structure.

The worst global-witness margins live where the small-part widest safe arc I_P is SHORT
relative to the cluster period 1/Vmax (i.e. Vmax*W_P small). In that regime the cluster
plants FEW teeth in I_P, and we can bound the largest free gap RIGOROUSLY.

RIGOROUS LEMMA (short-arc regime).  Let I_P=(alpha,beta) be a small-part safe arc of
width W_P. Let L={u_1<...<u_c} be the cluster (all in [Vmin,Vmax]).
  - Each runner u plants at most floor(W_P * u) + 1 tooth-centers in I_P (centers j/u are
    spaced 1/u apart, so an interval of length W_P holds <= floor(W_P*u)+1 of them).
  - Hence the cluster plants at most  N := sum_{u in L} (floor(W_P*u)+1) teeth in I_P.
    These merge into <= N danger components, so I_P has <= N+1 free gaps.
  - Each tooth contributes danger measure <= 2/(14u) = 1/(7u) <= 1/(7 Vmin).
    Total danger in I_P <= sum_u (floor(W_P*u)+1)/(7u) <= sum_u (W_P*u+1)/(7u)
      = c*W_P/7 + (1/7) sum_u 1/u.
  - Free measure in I_P >= W_P - c*W_P/7 - (1/7) H_L,  H_L := sum_{u in L} 1/u.
  - Largest free gap >= free_measure / (#gaps) >= [W_P(1 - c/7) - H_L/7] / (N+1).

For this to beat 1/(7 Vmax) we need the numerator positive and N small. When
Vmax*W_P < 1 (short arc), floor(W_P*u)=0 for all u (since W_P*u <= W_P*Vmax < 1), so:
  - N = sum_u 1 = c  (each runner at most ONE tooth in I_P!)
  - #gaps <= c+1
  - danger in I_P <= sum_u (1)/(7u) = H_L/7 <= c/(7 Vmin)
  - free >= W_P - H_L/7 >= W_P - c/(7 Vmin)
  - largest gap >= (W_P - H_L/7)/(c+1).

CLAIM (short-arc lemma): if  W_P*Vmax < 1  and
        (W_P - H_L/7)/(c+1) > 1/(7 Vmax)
   <=>  7 Vmax (W_P - H_L/7) > c+1
   <=>  7 Vmax W_P - Vmax H_L > c+1,
   then a global witness exists (M(S) >= 1/14). RIGOROUS.

We TEST: (a) is the short-arc regime (Vmax*W_P<1) actually where the worst margins are?
         (b) does the short-arc lemma inequality hold there?
         (c) what fraction of ALL S3 does the short-arc lemma + single-gap lemma close?
         (d) for the rest, are they bounded (finite check)?

NOTE W_P here = width of a small-part safe arc; we use the WIDEST. But the lemma applies
to ANY small arc, so we may pick the arc maximizing the bound.
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
from itertools import combinations
import random
from collections import Counter

def merge(ivs):
    ivs = sorted(ivs); m = []
    for a, b in ivs:
        if m and a <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    return m

def small_safe_arcs(P, h=F(1, 14)):
    iv = []
    for u in P:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    m = merge(iv)
    safe = []; prev = F(0)
    for a, b in m:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def teeth_in(u, lo, hi, h=F(1, 14)):
    out = []
    jmin = int(lo * u) - 1; jmax = int(hi * u) + 1
    for j in range(jmin, jmax + 1):
        c = F(j, u); a = c - h/u; b = c + h/u
        aa = max(a, lo); bb = min(b, hi)
        if aa < bb: out.append((aa, bb))
    return out

def actual_global_gap(S):
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    safe = small_safe_arcs(P); best = F(0)
    for (a, b) in safe:
        dgr = merge([t for u in L for t in teeth_in(u, a, b)])
        prev = a
        for (x, y) in dgr:
            if x > prev and x - prev > best: best = x - prev
            prev = max(prev, y)
        if b > prev and b - prev > best: best = b - prev
    return best

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'
def missing_q(P): return [q for q in range(2, 15) if not any(v % q == 0 for v in P)]

def single_gap_closable(S):
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    Vmin, Vmax = min(L), max(L); s = Vmax - Vmin
    safe = small_safe_arcs(P); K = 0
    while (13 * Vmin - Vmax > 14 * K * s) if s > 0 else (K == 0):
        lo = F(14*K+1, 14)/Vmin; hi = F(14*K+13, 14)/Vmax
        if lo < hi:
            for (a, b) in safe:
                if max(a, lo) < min(b, hi): return True
        if hi >= 1 or (s == 0 and K > 0): break
        K += 1
        if K > 14*Vmax: break
    return False

def short_arc_lemma_closes(S):
    """Apply the rigorous short-arc lemma over ALL small arcs; return True if the
    RIGOROUS inequality fires for some arc (=> proof). Uses provable bound only."""
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    Vmax = max(L); HL = sum(F(1, u) for u in L); c = len(L)
    for (a, b) in small_safe_arcs(P):
        WP = b - a
        # general provable bound: N=sum_u(floor(WP*u)+1); danger<=sum_u (floor(WP*u)+1)/(7u)
        N = sum(floor(WP * u) + 1 for u in L)
        danger_ub = sum(F(floor(WP * u) + 1, 7 * u) for u in L)
        free_lb = WP - danger_ub
        if free_lb <= 0: continue
        gap_lb = free_lb / (N + 1)
        if gap_lb > F(1, 7 * Vmax): return True
    return False

def gen_constructive(seed=0, target=2000, Vrange=(50, 2000)):
    rng = random.Random(seed); out = []; tries = 0
    smalls = []
    base = list(range(1, 14))
    for sz in (11, 10, 9, 8, 7):
        for P in combinations(base, sz): smalls.append(list(P))
    while len(out) < target and tries < target * 300:
        tries += 1
        P = rng.choice(smalls); c = 13 - len(P)
        if c < 2: continue
        miss = missing_q(P); V = rng.randint(*Vrange)
        spread = rng.choice([2, 3, 5, 7, 10, 14, 20, 28, 40, 56, 70, 90])
        window = list(range(V, V + spread + 1))
        if len(window) < c: continue
        cluster = rng.sample(window, c)
        if not all(any(v % q == 0 for v in cluster) for q in miss): continue
        S = sorted(set(P) | set(cluster))
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    allS = []
    for (lo, hi) in [(50, 200), (200, 800), (800, 3000), (3000, 10000)]:
        allS += gen_constructive(seed=lo+3, target=1200, Vrange=(lo, hi))
    n = len(allS)
    print(f"n S3 = {n}")
    sg = 0; sa = 0; both = 0; neither = []
    for S in allS:
        a = single_gap_closable(S); b = short_arc_lemma_closes(S)
        if a: sg += 1
        if b: sa += 1
        if a or b: both += 1
        else: neither.append(S)
    print(f"single-gap lemma closes:        {sg}/{n} ({100*sg/n:.1f}%)")
    print(f"short-arc rigorous lemma closes:{sa}/{n} ({100*sa/n:.1f}%)")
    print(f"EITHER (rigorous proof):        {both}/{n} ({100*both/n:.1f}%)")
    print(f"neither (residual):             {len(neither)}/{n} ({100*len(neither)/n:.1f}%)")
    # residual: bounded Vmax? distribution of Vmax, c, s among neither
    if neither:
        vmaxs = sorted(max(u for u in S if u > 13) for S in neither)
        cs = Counter(sum(1 for u in S if u > 13) for S in neither)
        ss = Counter((max(u for u in S if u>13)-min(u for u in S if u>13)) for S in neither)
        print(f"\nresidual Vmax range: [{vmaxs[0]}, {vmaxs[-1]}]  median~{vmaxs[len(vmaxs)//2]}")
        print(f"residual cluster sizes c: {dict(sorted(cs.items()))}")
        # check all residual still have a global witness (sanity)
        gw_ok = sum(1 for S in neither if actual_global_gap(S) > 0)
        print(f"residual sets with an actual global witness (meas>0): {gw_ok}/{len(neither)}")
        print("sample residual sets:")
        for S in neither[:10]:
            L = [u for u in S if u>13]
            print(f"   {S}  c={len(L)} s={max(L)-min(L)} Vmax={max(L)}")
