#!/usr/bin/env python3
"""shallow_witness_lemma_kps_S128c56.py -- kind-pasteur S128 cont.56.
THE SHALLOW-WITNESS COUNTING LEMMA for the d_min<=5 covering stratum.

CONSTRUCTION. V = P u {v1, v2}, P subset {1..12}, v2 = v1 + d (d <= 5), v1 >= 156.
Pick m and set q = ceil(v1/m).  Then v1 = m*q - s with 0 <= s < m, so v1 = -s (mod q)
and v2 = d - s (mod q): BOTH KILLERS COLLAPSE TO SMALL SPEEDS mod q.
Let e1 = |least-abs residue of v1|, e2 = same for v2, E = P u {e1,e2} (all <= 12ish).
BAND: if e_max <= 13 e_min, then every integer a with
      q/(14 e_min)  <=  a  <=  13q/(14 e_max)
gives t = a/q with p*a in [q/14, 13q/14] for EVERY effective speed, hence ||v t|| >= 1/14
for all v in V.  The interval has length q(13/e_max - 1/e_min)/14, so an integer exists
once that exceeds 1.
VERIFY exactly on the 8 critical covering families (and a wider sweep). Data only."""
import sys
from fractions import Fraction as F
from math import gcd, ceil, floor
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def least_abs(v,q):
    r=v%q
    return min(r,q-r)
def build_witness(P,v1,v2,mlist=(2,3,4,5)):
    """return (q,a,margin) or None"""
    for m in mlist:
        q=ceil(v1/m)
        if q<2: continue
        e1=least_abs(v1,q); e2=least_abs(v2,q)
        if e1==0 or e2==0: continue
        E=sorted(set(list(P)+[e1,e2]))
        emin,emax=min(E),max(E)
        if emax>13*emin: continue
        lo=F(q,14*emin); hi=F(13*q,14*emax)
        if hi<lo: continue
        a=ceil(lo)
        if a>hi: continue
        if a<=0: continue
        return q,a,e1,e2,E,float(hi-lo)
    return None
print("=== the 8 critical covering d_min<=5 families ===")
crit=[(168,169),(195,196),(208,210),(221,224),(234,238),(247,252),(294,299),(308,312)]
P=list(range(2,13))
allok=True
print("  %-11s %-5s %-5s %-8s %-8s %-13s %-9s"%("killers","q","a","e1","e2","min||v t||",">=1/14"))
for v1,v2 in crit:
    r=build_witness(P,v1,v2)
    if r is None:
        print("  %-11s NO WITNESS FOUND"%("%d,%d"%(v1,v2))); allok=False; continue
    q,a,e1,e2,E,ilen=r
    V=sorted(P+[v1,v2])
    t=F(a,q)
    mn=min(nd(v*t) for v in V)
    ok = mn>=F(1,14)
    if not ok: allok=False
    print("  %-11s %-5d %-5d %-8d %-8d %-13s %-9s"%("%d,%d"%(v1,v2),q,a,e1,e2,mn,ok))
print("  all 8 witnessed at level 1/14:",allok)
print()
print("=== wider sweep: P = 11-subsets of {1..12}, v1 in [156,900], d in 1..5 ===")
import random
random.seed(56)
tot=0; okc=0; fails=[]
for _ in range(400):
    drop=random.randint(1,12)
    P2=[x for x in range(1,13) if x!=drop]
    v1=random.randint(157,900); d=random.randint(1,5); v2=v1+d
    if v1<=13*max(P2): continue
    V=sorted(P2+[v1,v2])
    if len(set(V))!=13: continue
    tot+=1
    r=build_witness(P2,v1,v2)
    if r is None:
        fails.append((tuple(P2),v1,v2,'no-q')); continue
    q,a,e1,e2,E,ilen=r
    t=F(a,q); mn=min(nd(v*t) for v in V)
    if mn>=F(1,14): okc+=1
    else: fails.append((tuple(P2),v1,v2,str(mn)))
print("  witnessed %d / %d"%(okc,tot))
if fails:
    print("  failures (%d), first 5: %s"%(len(fails),fails[:5]))
print("DONE")
