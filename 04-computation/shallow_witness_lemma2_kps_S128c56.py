#!/usr/bin/env python3
"""shallow_witness_lemma2_kps_S128c56.py -- the band construction with a WIDE q-scan.
Fix: the interval [q/(14 emin), 13q/(14 emax)] has length q(13/emax - 1/emin)/14, which is
q/168 when emin=1,emax=12 (too short for small q) but q/24 when emin=2,emax=12.
So SCAN q and require the reduced killer residues to keep emin >= 2 (when 1 is not in P).
Report: does a valid q always exist, and in what range?  Data only."""
import sys, random
from fractions import Fraction as F
from math import ceil
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(v,q):
    r=v%q; return min(r,q-r)
def find_q(P,v1,v2,qlo=20,qhi=None):
    """scan q; return first (q,a,e1,e2) whose band interval contains an integer."""
    if qhi is None: qhi=v1
    for q in range(qlo,qhi+1):
        e1=la(v1,q); e2=la(v2,q)
        if e1==0 or e2==0: continue
        E=sorted(set(list(P)+[e1,e2]))
        emin,emax=min(E),max(E)
        if emax>13*emin: continue
        lo=F(q,14*emin); hi=F(13*q,14*emax)
        if hi<lo: continue
        a=ceil(lo)
        if a>hi or a<=0: continue
        return q,a,e1,e2,emin,emax
    return None
print("=== the 8 critical families, WIDE q-scan ===")
crit=[(168,169),(195,196),(208,210),(221,224),(234,238),(247,252),(294,299),(308,312)]
P=list(range(2,13))
allok=True
print("  %-11s %-5s %-5s %-4s %-4s %-6s %-12s %-6s"%("killers","q","a","e1","e2","emin","min||v t||",">=1/14"))
for v1,v2 in crit:
    r=find_q(P,v1,v2)
    if r is None:
        print("  %-11s NO q FOUND"%("%d,%d"%(v1,v2))); allok=False; continue
    q,a,e1,e2,emin,emax=r
    V=sorted(P+[v1,v2]); t=F(a,q)
    mn=min(nd(v*t) for v in V); ok=mn>=F(1,14)
    if not ok: allok=False
    print("  %-11s %-5d %-5d %-4d %-4d %-6d %-12s %-6s"%("%d,%d"%(v1,v2),q,a,e1,e2,emin,mn,ok))
print("  all 8 witnessed:",allok)
print()
print("=== wide sweep: P = 11-subsets of {1..12}, v1 in [157,900], d in 1..5 ===")
random.seed(561)
tot=0; okc=0; qs=[]; fails=[]
for _ in range(300):
    drop=random.randint(1,12)
    P2=[x for x in range(1,13) if x!=drop]
    v1=random.randint(157,900); d=random.randint(1,5); v2=v1+d
    if v1<=13*max(P2): continue
    V=sorted(P2+[v1,v2])
    if len(set(V))!=13: continue
    tot+=1
    r=find_q(P2,v1,v2)
    if r is None: fails.append((v1,v2,drop,'no-q')); continue
    q,a,e1,e2,emin,emax=r
    t=F(a,q); mn=min(nd(v*t) for v in V)
    if mn>=F(1,14): okc+=1; qs.append(q)
    else: fails.append((v1,v2,drop,str(mn)))
print("  witnessed %d / %d"%(okc,tot))
if qs: print("  q used: min %d, median %d, max %d  (v1 range 157..900)"%(min(qs),sorted(qs)[len(qs)//2],max(qs)))
if fails: print("  failures %d, first 5: %s"%(len(fails),fails[:5]))
print("DONE")
