#!/usr/bin/env python3
"""closed_form_hunt_kps_S128c81.py -- kind-pasteur S128 cont.81.
CROSS-AGENT CLOSED-FORM HUNT.

boxeph-S120's located-maximizer theorem says the optimum is attained at a STRADDLING ACTIVE
PAIR and equals
        M = |v_i a_j - v_j a_i| / (v_i + v_j)
for integers a_i, a_j.  My thread has produced several EXACT M values -- 1/7 for {2..12},
the k/(7k+1) family, and the clustered worst cases -- so this is a sharp cross-check, and if
it holds it hands my families closed forms.
Also test my own guesses: |B| = 1/297 and concentration ratio 198/7 = 28.2857."""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def Mexact(V,qhi=4000):
    best=(F(0),None)
    for q in range(2,qhi+1):
        for a in range(1,q//2+1):
            if gcd(a,q)!=1: continue
            t=F(a,q); mn=min(nd(v*t) for v in V)
            if mn>best[0]: best=(mn,t)
    return best
def boxeph_pair(V,t,M):
    """find a straddling active pair realising M = |v_i a_j - v_j a_i|/(v_i+v_j)"""
    act=[v for v in V if nd(v*t)==M]
    out=[]
    for vi,vj in itertools.combinations(sorted(set(V)),2):
        s=vi+vj
        for ai in range(-3,4):
            for aj in range(-3,4):
                if F(abs(vi*aj-vj*ai),s)==M and (vi in act or vj in act):
                    out.append((vi,vj,ai,aj,s))
    return act,out[:3]
print("### boxeph-S120's located-maximizer formula vs my exact values ###")
FAMS=[("{2..12}",list(range(2,13))),
      ("{2..12} u {157,314}",list(range(2,13))+[157,314]),
      ("{2..12} u {169,182}",list(range(2,13))+[169,182]),
      ("{1..11} u {312,364}",list(range(1,12))+[312,364]),
      ("{1..12}",list(range(1,13))),
      ]
for name,V in FAMS:
    M,t=Mexact(V)
    act,pairs=boxeph_pair(V,t,M)
    ok = len(pairs)>0
    print("  %-22s M = %-9s at t = %-9s  active speeds %s"%(name,M,t,act))
    if pairs:
        vi,vj,ai,aj,s = pairs[0]
        print("        boxeph pair (v_i,v_j) = (%d,%d), (a_i,a_j) = (%d,%d): |%d*%d - %d*%d|/%d = %s  MATCH: %s"%(
            vi,vj,ai,aj,vi,aj,vj,ai,s,F(abs(vi*aj-vj*ai),s),F(abs(vi*aj-vj*ai),s)==M))
    else:
        print("        no straddling pair found within |a|<=3")
print()
print("### the closed forms this gives my clustered families ###")
print("  family                       M         v_i+v_j     M*(v_i+v_j)   note")
for name,V in FAMS:
    M,t=Mexact(V)
    act,pairs=boxeph_pair(V,t,M)
    if pairs:
        vi,vj,ai,aj,s=pairs[0]
        print("  %-28s %-9s %-11d %-13s %s"%(name,M,s,M*s,"= |v_i a_j - v_j a_i|"))
print()
print("### my own guesses: |B| and the concentration ratio ###")
print("  measured |B| (140^3 grid) = 33/9800 = %.7f"%(33/9800))
for cand,lab in [(F(1,297),"1/297"),(F(1,343),"1/343 (interior only)"),(F(33,9800),"33/9800 (grid)"),(F(1,294),"1/294"),(F(2,588),"2/588")]:
    print("    %-22s = %.7f   ratio (2/21)/cand = %s = %.4f"%(lab,float(cand),F(2,21)/cand,float(F(2,21)/cand)))
print()
print("  measured concentration 28.28 ; candidates:")
for cand,lab in [(F(198,7),"198/7"),(F(2,21)/F(1,297),"(2/21)/(1/297)"),(F(2,21)/F(33,9800),"(2/21)/(33/9800)")]:
    print("    %-18s = %.5f"%(lab,float(cand)))
print("DONE")
