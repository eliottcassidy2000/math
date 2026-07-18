#!/usr/bin/env python3
"""small_modulus_avoid_kps_S128c58.py -- kind-pasteur S128 cont.58.
THE SMALL-MODULUS AVOIDANCE CRITERION.  For q <= 28, la(s,q) >= 2 gives ||.|| >= 2/q >= 1/14.
So it SUFFICES to find q <= 28 and a with  v*a mod q  outside {0, 1, q-1}  for every speed v.
This tolerates WRAPAROUND (unlike the no-wrap band lemma), which is exactly what the 21
cont.58 failures needed.  Test the criterion on the full covering wide-cluster census.
PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def avoid_witness(V,qlo=15,qhi=28):
    """find (q,a) with v*a mod q not in {0,1,q-1} for all v; returns all hits"""
    hits=[]
    for q in range(qlo,qhi+1):
        for a in range(1,q):
            ok=True
            for v in V:
                r=(v*a)%q
                if r==0 or r==1 or r==q-1: ok=False; break
            if ok: hits.append((q,a))
    return hits
print("### criterion check on ALL covering wide clusters (core in {1..12}, killers 13a & 14b) ###")
tot=0; ok=0; fails=[]
byq={}
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]
    mu,M=min(P),max(P); span=M-mu
    for a13 in range(1,40):
        k1=13*a13
        if k1<=13*M: continue
        for b14 in range(1,40):
            k2=14*b14
            if k2<=13*M or k2==k1: continue
            if abs(k2-k1)<=span: continue
            K=sorted([k1,k2]); V=P+K
            if len(set(V))!=13 or K[0]<=13*M: continue
            tot+=1
            h=avoid_witness(V)
            if h:
                ok+=1; byq[h[0][0]]=byq.get(h[0][0],0)+1
            else: fails.append((drop,P,K))
print("  covering wide clusters: %d ; small-modulus-avoidance certified: %d ; failures: %d"%(tot,ok,len(fails)))
print("  first successful modulus distribution:", dict(sorted(byq.items())))
if fails:
    print("  failures (first 6):")
    for drop,P,K in fails[:6]: print("    drop=%d core=%s K=%s"%(drop,P,K))
print()
print("### the 21 former band-failures, under the new criterion ###")
FAILS=[(2,[308,338]),(2,[294,351]),(2,[351,378]),(3,[308,338]),(3,[294,351]),(3,[351,378]),
       (4,[308,338]),(4,[294,351]),(4,[351,378]),(5,[308,338]),(5,[294,351]),(5,[351,378])]
for drop,K in FAILS:
    P=[x for x in range(1,13) if x!=drop]; V=P+K
    h=avoid_witness(V)
    if h:
        q,a=h[0]; t=F(a,q); mn=min(nd(v*t) for v in V)
        print("  drop=%-2d K=%-12s -> q=%-3d a=%-3d t=%-7s min||vt||=%-7s >=1/14: %s  (#valid (q,a): %d)"%(
            drop,str(K),q,a,t,mn,mn>=F(1,14),len(h)))
    else:
        print("  drop=%-2d K=%-12s -> NO (q,a) in [15,28]"%(drop,str(K)))
print()
print("### how much slack?  number of valid (q,a) pairs across the census ###")
import statistics
cnts=[]
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]; M=max(P); span=M-min(P)
    for a13 in range(13,30):
        k1=13*a13
        if k1<=13*M: continue
        for b14 in range(12,30):
            k2=14*b14
            if k2<=13*M or k2==k1 or abs(k2-k1)<=span: continue
            K=sorted([k1,k2]); V=P+K
            if len(set(V))!=13 or K[0]<=13*M: continue
            cnts.append(len(avoid_witness(V)))
cnts.sort()
print("  families: %d ; valid (q,a) pairs: min=%d  p10=%d  median=%d  max=%d"%(
    len(cnts),cnts[0],cnts[len(cnts)//10],cnts[len(cnts)//2],cnts[-1]))
print("  families with ZERO valid pairs: %d"%sum(1 for c in cnts if c==0))
print("DONE")
