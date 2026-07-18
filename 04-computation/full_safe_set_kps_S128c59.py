#!/usr/bin/env python3
"""full_safe_set_kps_S128c59.py -- kind-pasteur S128 cont.59.
CLOSING THE MIXED CASE with the FULL core-safe set instead of one interval round 1/13.
S(P) = {t : ||p t|| >= 1/14 for all p in P} is an exact finite union of rational intervals.
For a small killer k1 remove Bad_k1 EXACTLY; let L be the largest surviving interval.
A second killer k2 then only needs  2/(7 k2) < L*(6/7), i.e.  k2 > 1/(3L).
If max over cores and small k1 of 1/(3L) is below 874, the finite check already done
(both killers < 874, 41986/41986) closes everything.  PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
    """exact union of intervals in [0,1] where every p in P has ||p t|| >= 1/14"""
    bps={F(0),F(1)}
    for p in P:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                v=F(j,p)+s
                if 0<=v<=1: bps.add(v)
    B=sorted(bps); out=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        mid=(a+b)/2
        if all(nd(p*mid)>=F(1,14) for p in P):
            if out and out[-1][1]==a: out[-1]=(out[-1][0],b)
            else: out.append((a,b))
    return out
def remove_bad(iv,k):
    """subtract {t : ||k t|| < 1/14} from a list of intervals"""
    out=[]
    for a,b in iv:
        jlo=int(a*k)-1; jhi=int(b*k)+1
        bad=[]
        for j in range(jlo,jhi+1):
            x=F(j,k)-F(1,14*k); y=F(j,k)+F(1,14*k)
            if y<=a or x>=b: continue
            bad.append((max(x,a),min(y,b)))
        bad.sort(); cur=a
        for x,y in bad:
            if x>cur: out.append((cur,x))
            cur=max(cur,y)
        if cur<b: out.append((cur,b))
    return out
print("### the full core-safe set S(P) ###")
print("  drop  #intervals  total measure   largest interval   1/(3L)")
SS={}
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]
    iv=safe_set(P); SS[drop]=(P,iv)
    tot=sum(b-a for a,b in iv); L=max(b-a for a,b in iv)
    print("  %-5d %-11d %-15.7f %-18.7f %.1f"%(drop,len(iv),float(tot),float(L),float(F(1,3)/L)))
print()
print("### after removing ONE small killer k1: worst-case largest surviving interval ###")
print("  drop   worst k1   largest surviving L    implied k2 threshold 1/(3L)")
globalmax=0; argmax=None
for drop in range(1,13):
    P,iv=SS[drop]; M=max(P)
    worst=None
    for k1 in range(13*M+1,874):
        rem=remove_bad(iv,k1)
        L=max((b-a for a,b in rem),default=F(0))
        if L==0:
            worst=(F(0),k1); break
        if worst is None or L<worst[0]: worst=(L,k1)
    L,k1=worst
    thr=float(F(1,3)/L) if L>0 else float('inf')
    if L>0 and thr>globalmax: globalmax=thr; argmax=(drop,k1)
    print("  %-6d %-10d %-22.8f %s"%(drop,k1,float(L),"INFINITE (k1 swallows S)" if L==0 else "%.1f"%thr))
print()
print("  worst k2 threshold over all cores and all small k1: %.1f at %s"%(globalmax,argmax))
print("  finite check already covers both killers < 874")
print("  => gap remains iff  threshold > 874  for some (core, k1)")
print()
print("### direct verification on CONSTRUCTED covering mixed families ###")
print("  (k1 small, k2 large, built so 13| and 14| are supplied)")
tested=0; ok=0; bad=[]
for drop in range(1,13):
    P,iv=SS[drop]; M=max(P)
    for k1 in [13*M+1,169,182,200,260,364,500,700,860]:
        if k1<=13*M: continue
        for mult in [874,1000,5000,50000,500000]:
            for k2 in range(mult,mult+40):
                V=P+[k1,k2]
                if len(set(V))!=13: continue
                if not all(any(v%q==0 for v in V) for q in range(2,15)): continue
                tested+=1
                rem=remove_bad(iv,k1); got=None
                for a,b in rem:
                    if b-a<=0: continue
                    N=400
                    for j in range(N+1):
                        t=a+(b-a)*F(j,N)
                        if min(nd(v*t) for v in V)>=F(1,14): got=t; break
                    if got: break
                if got: ok+=1
                else: bad.append((drop,k1,k2))
                break
print("  constructed covering mixed families: %d ; certified by a t in S(P)\Bad_k1: %d"%(tested,ok))
if bad: print("  grid missed:",bad[:6])
print("DONE")
