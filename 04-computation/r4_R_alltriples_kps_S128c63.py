#!/usr/bin/env python3
"""r4_R_alltriples_kps_S128c63.py -- kind-pasteur S128 cont.63 (background).
SETTLE r=4's R.  cont.62 scanned only killers in [lo, lo+55).  Here: ALL triples with
k3 <= lo+100 exhaustively (the worst case provably sits at the bottom for r=2/3/4), plus a
decay check showing R falls as k3 grows, so the bottom window is the whole story.
R = T/k3 with T = min(N/(6 mu), 1/(3L)).  PRINT DATA ONLY."""
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
W=100
print("### r=4 : ALL triples with k3 <= lo+%d, over all %d cores ###"%(W,len(C9)))
worst=None; cnt=0; zero=0
for ci,P in enumerate(C9):
    M=max(P); lo=13*M+1; iv0=safe_set(P)
    ks=list(range(lo,lo+W))
    for i in range(len(ks)):
        r1=rm(iv0,ks[i])
        if not r1: zero+=1; continue
        for j in range(i+1,len(ks)):
            r2=rm(r1,ks[j])
            if not r2: zero+=1; continue
            for h in range(j+1,len(ks)):
                T=thresh(rm(r2,ks[h]))
                cnt+=1
                if T is None: zero+=1; continue
                R=T/ks[h]
                if worst is None or R>worst[0]: worst=(R,tuple(P),ks[i],ks[j],ks[h],T)
    if ci%50==0: print("  ... core %d/%d  triples=%d  maxR=%.5f"%(ci,len(C9),cnt,worst[0] if worst else 0))
R,Pw,a,b,c,T=worst
print()
print("  triples tested: %d ; swallow cases: %d"%(cnt,zero))
print("  MAX R = %.5f at core %s, killers (%d,%d,%d), T = %.1f"%(R,list(Pw),a,b,c,T))
print("  R < 1: %s   margin %.5f"%(R<1,1-R))
print()
print("### decay check: R as k3 grows with k1,k2 pinned at the worst pair ###")
P=list(Pw); iv0=safe_set(P); r2=rm(rm(iv0,a),b)
for k3 in [c,c+50,c+150,c+400,c+1000,c+3000,c+9000]:
    T3=thresh(rm(r2,k3))
    if T3: print("  k3=%-8d T=%-10.1f R=%.5f"%(k3,T3,T3/k3))
print("DONE")
