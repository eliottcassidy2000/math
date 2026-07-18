#!/usr/bin/env python3
"""r6_maxT_chunk_kps_S128c64.py -- kind-pasteur S128 cont.64 (background).
STEP 1 FIRST, per the rule I got wrong at r=5: the finite-horn bound is set by max T over
the region where the measure horn FAILS (R >= 1), NOT by the largest removed killer.
r=6: cores are 7-subsets of {1..12} (792).  Remove five killers exactly; the sixth needs
k6 > T = min(N/(6 mu), 1/(3L)).  R = T/k5.  Report max R, max T over R>=1, and whether the
worst case sits strictly inside the scan window (if it sits at the boundary the window is
too small and the answer is not trustworthy).  PRINT DATA ONLY."""
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
import os
START=int(os.environ.get('C0','0')); END=int(os.environ.get('C1','792'))
C7=[sorted(c) for c in itertools.combinations(range(1,13),7)][START:END]
W=16
print("### r=6 : %d seven-speed cores ; window [lo, lo+%d) ###"%(len(C7),W))
st=[]
for P in C7[:40]:
    iv=safe_set(P); st.append(sum(b-a for a,b in iv))
print("  S(P) measure (sample): min %.5f max %.5f"%(min(st),max(st)))
print()
worstR=None; maxT=0.0; maxk=0; cases=0; tested=0
for ci,P in enumerate(C7):
    M=max(P); lo=13*M+1; iv0=safe_set(P)
    ks=list(range(lo,lo+W))
    for i in range(len(ks)):
        r1=rm(iv0,ks[i])
        if not r1: continue
        for j in range(i+1,len(ks)):
            r2=rm(r1,ks[j])
            if not r2: continue
            for h in range(j+1,len(ks)):
                r3=rm(r2,ks[h])
                if not r3: continue
                for m in range(h+1,len(ks)):
                    r4=rm(r3,ks[m])
                    if not r4: continue
                    for n in range(m+1,len(ks)):
                        T=thresh(rm(r4,ks[n]))
                        tested+=1
                        if T is None: continue
                        R=T/ks[n]
                        if worstR is None or R>worstR[0]: worstR=(R,tuple(P),ks[i],ks[j],ks[h],ks[m],ks[n],T)
                        if R>=1.0:
                            cases+=1
                            if T>maxT: maxT=T
                            if ks[n]>maxk: maxk=ks[n]
    if ci%60==0: print("  ... core %d/%d tested=%d maxR=%.4f maxT=%.1f"%(ci,len(C7),tested,worstR[0] if worstR else 0,maxT))
R,Pw,a,b,c,d,e,T=worstR
print()
print("  quintuples tested: %d"%tested)
print("  MAX R = %.5f at core %s, killers (%d,%d,%d,%d,%d), T = %.1f"%(R,list(Pw),a,b,c,d,e,T))
print("  quintuples with R >= 1: %d ; largest killer among them: %d ; MAX T = %.1f"%(cases,maxk,maxT))
lo_w=13*max(Pw)+1
print("  worst case at offset %d within window of width %d -- %s"%(
    e-lo_w,W,"INSIDE (trustworthy)" if e-lo_w<W-3 else "AT BOUNDARY (window too small)"))
print()
print("  => r=6 finite-horn bound KB = max T + 25 = %d"%(int(maxT)+25))
print("DONE")
