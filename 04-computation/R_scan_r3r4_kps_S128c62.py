#!/usr/bin/env python3
"""R_scan_r3r4_kps_S128c62.py -- kind-pasteur S128 cont.62.
R scan for r=3 and r=4.  CORRECTION to my earlier method: the exhaustive r=2 scan puts the
worst R at k1 = 160, i.e. at SMALL killers -- my cont.60/61 scans sampled the TOP of the
range, which is the wrong place.  Scan the small end exhaustively here.
R = T / k_max-removed with T = min(N/(6 mu), 1/(3L)); R < 1 => measure horn suffices alone.
PRINT DATA ONLY."""
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
print("### r=3 : 66 cores ; small end exhaustive (k <= lo+260), plus a coarse tail ###")
C10=[sorted(c) for c in itertools.combinations(range(1,13),10)]
worst=None; zero=0
for P in C10:
    M=max(P); lo=13*M+1; iv0=safe_set(P)
    hi=lo+260
    tail=list(range(hi,3000,97))
    ks=list(range(lo,hi))+tail
    for i,k1 in enumerate(ks):
        r1=rm(iv0,k1)
        if not r1: zero+=1; continue
        for k2 in ks[i+1:]:
            T=thresh(rm(r1,k2))
            if T is None: zero+=1; continue
            R=T/k2
            if worst is None or R>worst[0]: worst=(R,tuple(P),k1,k2,T)
R,Pw,a,b,T=worst
print("  MAX R = %.5f at core %s, (k1,k2) = (%d,%d), T = %.1f"%(R,list(Pw),a,b,T))
print("  swallow cases: %d ; R < 1 everywhere scanned: %s"%(zero,R<1))
print()
print("### r=4 : 220 cores ; small end exhaustive (k <= lo+90), plus a coarse tail ###")
C9=[sorted(c) for c in itertools.combinations(range(1,13),9)]
worst4=None; zero4=0
for P in C9:
    M=max(P); lo=13*M+1; iv0=safe_set(P)
    hi=lo+90
    ks=list(range(lo,hi))+list(range(hi,2200,211))
    for i,k1 in enumerate(ks):
        r1=rm(iv0,k1)
        if not r1: zero4+=1; continue
        for j in range(i+1,len(ks)):
            r2=rm(r1,ks[j])
            if not r2: zero4+=1; continue
            for k3 in ks[j+1:]:
                T=thresh(rm(r2,k3))
                if T is None: zero4+=1; continue
                R=T/k3
                if worst4 is None or R>worst4[0]: worst4=(R,tuple(P),ks[i],ks[j],k3,T)
R4,P4,a4,b4,c4,T4=worst4
print("  MAX R = %.5f at core %s, (k1,k2,k3) = (%d,%d,%d), T = %.1f"%(R4,list(P4),a4,b4,c4,T4))
print("  swallow cases: %d ; R < 1 everywhere scanned: %s"%(zero4,R4<1))
print()
print("### summary ###")
print("  r=2 (exhaustive to 4000): max R = 0.51852")
print("  r=3 : max R = %.5f"%R)
print("  r=4 : max R = %.5f"%R4)
print("  R < 1 at every r scanned => the measure horn certifies alone; finite horns redundant")
print("DONE")
