#!/usr/bin/env python3
"""
lrc14_wsb_convexorder_stage2_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

Stage 2 of the moment-convex-order angle.  Stage 1 established:
  * exact law of N engine works (cross-checked vs Sr).
  * P1: THM-536 B2 is a POINTWISE coupling N_E(x) >= N_{AP_{N+1}}(x) when E subset {0..N}.
  * P2: usual stochastic order does NOT make consec_k extremal among same-k sets (1/792).

NEW QUESTIONS (stage 2):
 (Q1) THE SPAN-CAP TABLE.  For each k in {8,9,10}: tabulate meas(S7(AP_{N+1})) for N=k-1,k,...
      and find N*(k) = largest N with meas(S7(AP_{N+1})) <= cap_k.  THM-536 says N*=7,8,10.
      This certifies (via pointwise dom) ALL E of span<=N*.  The residual = span>N* shapes.
      => quantify the residual: how many primitive k-sets have span in (N*, B] for various B,
         and is meas(S7) of EACH actually <= cap (the finite-check claim)?
 (Q2) CONVEX ORDER.  g(t) is NOT convex.  But L_y = E[g(N)].  Decompose g into a part handled
      by stochastic order + a residual.  Check: is N_consec <=_cx or >=_cx N_E for same k?
      (convex order ignores mean, compares spread).  Compute means E[N] and variances; does
      consec minimize/maximize E[N] (=6 - avg#hit = S_1-ish)?  Where does g's non-monotonicity bite?
 (Q3) THE KEY STRUCTURAL TEST: is there a UNIVARIATE function phi and order such that
      L_y(E) = E[g(N_E)] is controlled by a MONOTONE functional of the law that consec extremizes?
      Test candidate orders: (a) increasing-convex (icx), (b) the "g-cone" = functions h with
      same sign pattern as g.  For each, check whether consec is the extremal law in the bank.
 (Q4) PARTIAL SUMS / Sturmian for AP: the law of N_AP comes from a mechanical word cover-time.
      Inspect the law p_i(AP_k) as a function of k -- is it a nice closed shape (e.g. p_6 = 1/(7*(k-1))?)
"""
import sys, itertools, functools
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def law_N(E):
    E = sorted(set(E)); bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for j in range(7):
            c = F(j, 7)
            for m in range(e): bps.add((c + m) / e)
    bps = sorted(z for z in bps if F(0) <= z < F(1))
    p = [F(0)] * 7
    for i in range(len(bps)):
        x0 = bps[i]; x1 = bps[i + 1] if i + 1 < len(bps) else F(1)
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hit = set()
        for e in E:
            fr = (e * xm) % 1
            hit.add((fr.numerator * 7) // fr.denominator)
        nempty = len([j for j in range(1, 7) if j not in hit])
        p[nempty] += x1 - x0
    return p

def Ly_from_law(p, k):
    if k == 8:   return p[0] + F(1, 10) * p[3] + p[6]
    if k in (9, 10):
        g = lambda t: F(-(t - 2) * (t - 3) * (t - 6), 36)
        return sum(g(t) * p[t] for t in range(7))
    g = lambda t: F((t - 3) * (t - 4), 12)
    return sum(g(t) * p[t] for t in range(7))
def g_of(k):
    if k==8: return lambda t: F((t-1)*(t-2)*(t-4)*(t-5),40)
    if k in (9,10): return lambda t: F(-(t-2)*(t-3)*(t-6),36)
    return lambda t: F((t-3)*(t-4),12)

H = F(1, 14)
def danger(u):
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - H / u) % 1; b = (c + H / u) % 1
        if a < b: iv.append((a, b))
        else: iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def measGP(P):
    dz = sorted([iv for u in P for iv in danger(u)]); s = F(0); prev = F(0)
    for a, b in dz:
        if a > prev: s += a - prev
        prev = max(prev, b)
    if prev < 1: s += 1 - prev
    return s
@functools.lru_cache(None)
def cap(k): return min(measGP(P) for P in itertools.combinations(range(1, 14), 13 - k))

def mean_var(p):
    m = sum(F(t)*p[t] for t in range(7))
    v = sum(F(t)**2*p[t] for t in range(7)) - m*m
    return m, v
def tail(p): return [sum(p[t:],F(0)) for t in range(7)]

print("="*78)
print("(Q1) SPAN-CAP table: meas(S7(AP_{N+1})) vs cap_k; residual span shapes far under cap?")
print("="*78)
for k in (8,9,10):
    capk = cap(k)
    print(f"--- k={k}  cap={float(capk):.5f} ---")
    Nstar = None
    for N in range(k-1, k+10):
        AP = list(range(N+1))
        m = law_N(AP)[0]
        flag = "<=cap" if m <= capk else "OVER"
        if m <= capk: Nstar = N
        print(f"    AP_{N+1} (span {N}): meas(S7)={float(m):.5f}  {flag}")
        if m > capk and N > Nstar + 1: break
    print(f"    => N*({k}) = {Nstar}  (subset-dom certifies all span<= {Nstar})")
    # residual: among primitive k-sets with span in (N*, 13], is each meas(S7) <= cap?
    over = 0; tot = 0; maxm = F(0)
    for r in itertools.combinations(range(1,14), k-1):
        E = [0]+list(r)
        if max(E) <= Nstar: continue
        tot += 1
        m = law_N(E)[0]
        if m > maxm: maxm = m
        if m > capk: over += 1
    print(f"    residual span in ({Nstar},13]: {tot} sets, max meas(S7)={float(maxm):.5f}, OVER cap: {over}")

print()
print("="*78)
print("(Q2) means/variances; does consec extremize E[N]?  where g non-monotonicity bites")
print("="*78)
k=8; pc=law_N(list(range(k)))
mc,vc = mean_var(pc)
print(f"  consec_8: E[N]={float(mc):.4f} ({mc})  Var={float(vc):.4f}  L_y={float(Ly_from_law(pc,k)):.5f}")
bank=[[0]+list(r) for r in itertools.combinations(range(1,13),k-1)]
mn_lt=0; mn_gt=0   # E[N_E] < / > consec
both_cx=0
for E in bank:
    p=law_N(E); m,v=mean_var(p)
    if m < mc - F(1,10**12): mn_lt+=1
    if m > mc + F(1,10**12): mn_gt+=1
print(f"  bank {len(bank)}: E[N_E] < consec for {mn_lt}; > consec for {mn_gt}")
print(f"  (E[N]=6-avg#hit; consec MIN E[N] would mean consec hits most sectors on avg)")

print()
print("="*78)
print("(Q3) g-cone order: is consec the extremal law for EVERY h with g's sign/curvature class?")
print("="*78)
# The dual g for k=8: g=(t-1)(t-2)(t-4)(t-5)/40, values:
gk = g_of(8)
gvals = [gk(t) for t in range(7)]
print(f"  g(0..6) (k=8) = {[str(x) for x in gvals]}")
# decompose g in the basis of 'convex order generators': 1, t, and (t-c)_+ ... we instead test
# whether L_y(E) <= L_y(consec) follows from a FINITE set of linear functionals consec extremizes.
# Specifically: find all extreme rays of the cone of laws q (prob on {0..6}) and see if consec's law
# is on the boundary maximizing <g, q>. We instead do: among the bank, compute the support of
# laws and check if consec is a vertex of the achievable moment region for (S1,S2,S3,S4).
# Simpler diagnostic: the 'pointwise-dom certified' shapes have N_E >= N_consec? NO (P2). So what
# weaker order holds for the WIDE-SPAN dangerous shapes specifically?
# Print the law of the few shapes whose L_y is CLOSEST to consec (the tightest competitors).
scored = []
for E in bank:
    p=law_N(E); ly=Ly_from_law(p,k)
    scored.append((ly, E, p))
scored.sort(reverse=True)
print("  top-6 L_y competitors among same-k bank (consec should be #1):")
for ly,E,p in scored[:6]:
    m,v=mean_var(p)
    print(f"    L_y={float(ly):.6f}  E={E}  E[N]={float(m):.4f} Var={float(v):.4f}  p0={float(p[0]):.4f} p3={float(p[3]):.4f} p6={float(p[6]):.4f}")

print()
print("="*78)
print("(Q4) shape of AP law as function of k: is p_6(AP_k) = 1/(7(k-1))?  p_5? etc.")
print("="*78)
for k in range(7,15):
    p=law_N(list(range(k)))
    # test p_6 = 1/(7(k-1))
    guess6 = F(1,7*(k-1))
    print(f"  k={k}: p6={p[6]}  guess 1/(7(k-1))={guess6}  match={p[6]==guess6}  "
          f"p5={p[5]} p0={p[0]}")
print("\nDONE stage 2.")
