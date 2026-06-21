#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A — THE CLEAN RUNG: prove consec maximizes the even band B_2.

B_2 = S_0 - S_1 + S_2 = 1 - S_1 + S_2 = E[g_2(N)],  g_2(N)=1 - N + N(N-1)/2,
N = #empty INNER sectors among the 6 inner sectors {1..6} (sector 0 always hit).
WAIT: the IE is over ALL 7 sectors but sector 0 is always hit (0 in E), so the
miss events live on the 6 inner sectors.  Let M = N = #missed inner sectors in
{0..6} (sector 0 never missed).  S_r = E[C(N,r)] where N = #missed sectors
(0..7, but realized in 0..6).  measS7 = P(N=0).

g_2(N) = 1 - N + C(N,2) = (N^2 - 3N + 2)/2 = (N-1)(N-2)/2.  CONVEX in N (2nd diff=1>0).
g_2 >= 0, and g_2(0)=1, g_2(1)=0, g_2(2)=0, g_2(N)>0 for N>=3.

So B_2 = E[(N-1)(N-2)/2].  Maximizing B_2 over E is a CONVEX-ORDER extremal problem.

THE RELATION-CODE / MACWILLIAMS READ of B_2:
  S_1 = E[N] = sum over inner sectors s of P(s missed).
  S_2 = E[C(N,2)] = sum over pairs {s,t} of inner sectors of P(s,t both missed).
  Each P(...) is a Weyl integral = (uniform/iid value) + (relation-code correction).
  iid: under k iid uniform residues, E[N_iid] and E[C(N,2)_iid] are FIXED constants
       (depend only on k), so the iid part of B_2 is shape-independent.
  The SHAPE-DEPENDENT part of B_2 = (relation-code corrections to S_1) and
  (relation-code corrections to S_2).  consec-max of B_2 <=> consec extremizes
  this relation-code combination.

CLAIM TO TEST (the clean rung):
  (i)  B_2 = E[g_2(N)], g_2 convex deg-2.  [exact]
  (ii) consec maximizes B_2 over all primitive shapes (exhaustive box).  [VERIFY]
  (iii) decompose B_2's shape-dependence: -E[N] favors LOW mean empty (good
        equidistribution = consec via 3-distance theorem), +E[C(N,2)] favors HIGH
        pair-co-missing = positive association (FKG).  Show consec wins BOTH or
        wins the combination.
  (iv) MacWilliams/relation reading: the corrections to S_1 are SUPPORT->? and
        to S_2 are support->? in Lambda(E); B_2 = 1 - corr_S1 + corr_S2 packaged
        as a relation weight-enumerator with the OUTER Krawtchouk g_2 signs.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

def occupancy_law(E):
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    occ=defaultdict(lambda:F(0)); pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2
        hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator)
            hit.add((v.numerator*7)//v.denominator)
        occ[frozenset(hit)]+=hi-lo; pi[len(hit)]+=hi-lo
    return dict(occ),pi

def S_r_list(E):
    occ,pi=occupancy_law(E)
    Sr=[sum(pi[h]*comb(7-h,r) for h in range(8)) for r in range(8)]
    return Sr,pi

def B_J(E,J):
    Sr,pi=S_r_list(E)
    return sum((-1)**r*Sr[r] for r in range(J+1)),Sr,pi

def g2(N): return F((N-1)*(N-2),2)

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

if __name__=="__main__":
    print("="*78); print("(i) B_2 = E[g_2(N)], g_2(N)=(N-1)(N-2)/2 convex deg-2"); print("="*78)
    for N in range(8):
        print(f"   N={N}: g_2={g2(N)}  (>=0, convex: 2nd diff const =1)")
    print()
    for E in [list(range(8)),[0,2,3,4,5,6,7,8]]:
        B2,Sr,pi=B_J(E,2)
        Eg2=sum(pi[7-N]*g2(N) for N in range(8))  # pi indexed by h=7-N
        # actually pi[h]=P(h hit)=P(N=7-h missed); E[g2(N)]=sum_h pi[h] g2(7-h)
        Eg2b=sum(pi[h]*g2(7-h) for h in range(8))
        print(f" E={E}: B_2={float(B2):.6f}  E[g2(N)]=sum_h pi[h]g2(7-h)={float(Eg2b):.6f}  match={B2==Eg2b}")
    print()
    print("="*78); print("(ii) EXHAUSTIVE: does consec maximize B_2 over primitive shapes?"); print("="*78)
    for k,W in [(8,11),(8,12),(8,13),(9,11),(9,12)]:
        C=consec(k); B2c,_,_=B_J(C,2)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        bank=[E for E in bank if primitive(E)]
        beaters=[]; best=B2c; bestE=C
        for E in bank:
            B2,_,_=B_J(list(E),2)
            if B2>B2c+F(1,10**15): beaters.append((E,float(B2)))
            if B2>best: best=B2; bestE=list(E)
        tag="CONSEC MAX" if not beaters else f"{len(beaters)} BEATERS!"
        print(f" k={k} span<= {W}: {len(bank)} shapes, B_2(consec)={float(B2c):.6f}  best={float(best):.6f} by {bestE}  -> {tag}")
        for E,b in beaters[:5]: print(f"      beater {list(E)}: B_2={b:.6f}")

# ============================================================================
# (iii) DECOMPOSE B_2 = 1 - S_1 + S_2 into equidistribution + association.
# ============================================================================
# S_1 = sum_{s=0}^{6} P(sector s missed) = E[N].  Sector 0 never missed (0 in E),
#       so S_1 = sum_{s=1}^{6} P(sector s missed) = E[N_inner], N_inner in 0..6.
# S_2 = sum_{s<t} P(sectors s,t both missed) = E[C(N,2)].
# B_2 = 1 - E[N] + E[C(N,2)].
#
# EXACT FORM of P(sector s missed):  P(s missed) = meas{x : no e_i x lands in
#   [s/7,(s+1)/7) mod 1}.  By symmetry of the torus, for a FIXED shape E the
#   single-sector miss probabilities P(s missed), s=1..6, are NOT all equal
#   (sector 0 is special), but S_1 = their sum.
# KEY: S_1 is a SUM of measures of "gap" events.  For consec E={0..k-1}, the
#   residues frac(e_i x)=frac(i x) are an ARITHMETIC PROGRESSION on the torus
#   -> three-distance theorem -> most EVEN spacing -> MINIMAL clumping ->
#   minimal sum of empty-sector measure?  Test: does consec MINIMIZE S_1 (=E[N])?
# And does consec MAXIMIZE S_2 (=E[C(N,2)], positive association)?
# B_2 consec-max could come from EITHER or the COMBINATION.

def S1_S2(E):
    Sr,pi=S_r_list(E)
    return Sr[1],Sr[2],pi

if __name__=="__main__":
    print()
    print("="*78); print("(iii) decompose B_2 = 1 - S_1 + S_2  (S_1=E[N], S_2=E[C(N,2)])"); print("="*78)
    for k,W in [(8,12),(9,12)]:
        C=consec(k); s1c,s2c,_=S1_S2(C)
        bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
        bank=[E for E in bank if primitive(E)]
        s1_lower=0; s2_higher=0   # consec MIN S1? consec MAX S2?
        s1_min=s1c; s2_max=s2c
        for E in bank:
            s1,s2,_=S1_S2(list(E))
            if s1<s1c-F(1,10**15): s1_lower+=1
            if s2>s2c+F(1,10**15): s2_higher+=1
            if s1<s1_min: s1_min=s1
            if s2>s2_max: s2_max=s2
        print(f" k={k} span<= {W}: {len(bank)} shapes")
        print(f"   S_1(consec)=E[N]={float(s1c):.6f}; #shapes with LOWER S_1 = {s1_lower}  (consec min? {s1_lower==0}) global-min={float(s1_min):.6f}")
        print(f"   S_2(consec)=E[C(N,2)]={float(s2c):.6f}; #shapes with HIGHER S_2 = {s2_higher}  (consec max? {s2_higher==0}) global-max={float(s2_max):.6f}")
        print(f"   -> B_2 consec-max source: {'S_1-min only' if s1_lower==0 and s2_higher>0 else ('S_2-max only' if s2_higher==0 and s1_lower>0 else ('BOTH' if s1_lower==0 and s2_higher==0 else 'COMBINATION (neither alone)'))}")
