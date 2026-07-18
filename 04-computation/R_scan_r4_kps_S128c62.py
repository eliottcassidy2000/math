#!/usr/bin/env python3
"""R_scan_r4_kps_S128c62.py -- kind-pasteur S128 cont.62.  r=4 R scan, small end only
(where r=2 and r=3 both put the worst case).  R = T/k_max-removed; R<1 => measure horn alone.
Trend so far: 0.519 (r=2), 0.734 (r=3) -- rising, so whether it stays below 1 matters."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
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
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((float(a),float(b)))
    mg=[]
    for a,b in out:
        if mg and abs(mg[-1][1]-a)<1e-15: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return mg
def rm(iv,k):
    out=[]; w=1.0/(14.0*k)
    for a,b in iv:
        jlo=int(a*k)-1; jhi=int(b*k)+1; cur=a
        for j in range(jlo,jhi+1):
            x=j/k-w; y=j/k+w
            if y<=a or x>=b: continue
            x=max(x,a); y=min(y,b)
            if x>cur: out.append((cur,x))
            cur=max(cur,y)
        if cur<b: out.append((cur,b))
    return out
def thresh(iv):
    if not iv: return None
    mu=sum(b-a for a,b in iv); N=len(iv); L=max(b-a for a,b in iv)
    if mu<=0 or L<=0: return None
    return min(N/(6*mu), 1.0/(3*L))
C9=[sorted(c) for c in itertools.combinations(range(1,13),9)]
print("### r=4 : 220 nine-speed cores, killers in [lo, lo+55) exhaustive ###")
worst=None; zero=0
for P in C9:
    M=max(P); lo=13*M+1; iv0=safe_set(P)
    ks=list(range(lo,lo+55))
    for i in range(len(ks)):
        r1=rm(iv0,ks[i])
        if not r1: zero+=1; continue
        for j in range(i+1,len(ks)):
            r2=rm(r1,ks[j])
            if not r2: zero+=1; continue
            for h in range(j+1,len(ks)):
                T=thresh(rm(r2,ks[h]))
                if T is None: zero+=1; continue
                R=T/ks[h]
                if worst is None or R>worst[0]: worst=(R,tuple(P),ks[i],ks[j],ks[h],T)
R,Pw,a,b,c,T=worst
print("  MAX R = %.5f at core %s, killers (%d,%d,%d), T = %.1f"%(R,list(Pw),a,b,c,T))
print("  swallow cases: %d ; R < 1: %s"%(zero,R<1))
print()
print("### the trend in r ###")
print("  r=2  max R = 0.51852  (exhaustive to 4000)")
print("  r=3  max R = 0.73375")
print("  r=4  max R = %.5f"%R)
print("  margin to 1 at r=4: %.5f"%(1-R))
print("DONE")
