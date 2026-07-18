#!/usr/bin/env python3
"""R_exhaustive_kps_S128c62.py -- kind-pasteur S128 cont.62.
SCANNING THE RIGHT QUANTITY.  The measure horn removes the r-1 smaller killers and needs the
largest, k_r, to exceed T = min(N/(6 mu), 1/(3L)) of the surviving set.  Since k_r exceeds
the largest REMOVED killer, it suffices that
        R := T / k_max-removed  <  1 .
If R < 1 holds for every core and every choice of removed killers, then the measure horn
certifies EVERY family on its own and the finite horns of THM-1051/1061/1071 are REDUNDANT.
r=2 is cheap enough to scan EXHAUSTIVELY.  PRINT DATA ONLY."""
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
print("### r=2 : EXHAUSTIVE over all 12 cores and every removed killer k1 ###")
KHI=4000
worst=None; zeros=[]
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]; M=max(P)
    iv0=safe_set(P)
    for k1 in range(13*M+1,KHI):
        T=thresh(rm(iv0,k1))
        if T is None: zeros.append((drop,k1)); continue
        R=T/k1
        if worst is None or R>worst[0]: worst=(R,drop,k1,T)
print("  killers scanned: (13*maxP, %d)  -- %d cores"%(KHI,12))
print("  L = 0 cases (removal swallows S(P)): %d %s"%(len(zeros),zeros[:6]))
R,drop,k1,T=worst
print("  MAX R = %.5f at core-drop %d, k1 = %d (T = %.1f)"%(R,drop,k1,T))
print("  R < 1 everywhere: %s"%(R<1))
print()
print("### r=3 : all 66 cores, exhaustive pairs in (13maxP, %d) ###"%1200)
KHI3=1200
worst3=None; z3=0
C10=[sorted(c) for c in itertools.combinations(range(1,13),10)]
for P in C10:
    M=max(P); iv0=safe_set(P); lo=13*M+1
    for k1 in range(lo,KHI3):
        r1=rm(iv0,k1)
        if not r1: z3+=1; continue
        for k2 in range(k1+1,KHI3):
            T=thresh(rm(r1,k2))
            if T is None: z3+=1; continue
            R=T/k2
            if worst3 is None or R>worst3[0]: worst3=(R,tuple(P),k1,k2,T)
R3,P3,a3,b3,T3=worst3
print("  MAX R = %.5f at core %s, (k1,k2) = (%d,%d), T = %.1f"%(R3,list(P3),a3,b3,T3))
print("  L = 0 cases: %d ; R < 1 everywhere: %s"%(z3,R3<1))
print("DONE-PART1")
