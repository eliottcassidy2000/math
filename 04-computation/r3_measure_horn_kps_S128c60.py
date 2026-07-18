#!/usr/bin/env python3
"""r3_measure_horn_kps_S128c60.py -- kind-pasteur S128 cont.60.
THE r=3 MEASURE HORN.  Sharper than THM-1051's crude 2r/(L(7-r)) formula: remove the
(r-1) SMALL killers from S(P) EXACTLY, take the largest surviving interval L, and the LAST
killer only needs k > 1/(3L) -- the threshold never depends on r, only on L.
Core now has 10 speeds (13-3), so S(P) is LARGER than in the r=2 case.
Find the worst L over cores and small-killer pairs -> the split point B.  PRINT DATA ONLY."""
import sys, itertools, random
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
        mid=(a+b)/2
        if all(nd(p*mid)>=F(1,14) for p in P):
            if out and out[-1][1]==a: out[-1]=(out[-1][0],b)
            else: out.append((a,b))
    return out
def remove_bad(iv,k):
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
CORES=[sorted(c) for c in itertools.combinations(range(1,13),10)]
print("### S(P) for 10-speed cores (r=3 setting): %d cores ###"%len(CORES))
stats=[]
for P in CORES:
    iv=safe_set(P)
    stats.append((max(b-a for a,b in iv),sum(b-a for a,b in iv),len(iv),tuple(P)))
stats.sort()
print("  largest-component: min %.6f  median %.6f  max %.6f"%(
    float(stats[0][0]),float(stats[len(stats)//2][0]),float(stats[-1][0])))
print("  total measure:     min %.6f  max %.6f"%(
    float(min(s[1] for s in stats)),float(max(s[1] for s in stats))))
print("  worst core (smallest largest-component): %s  L=%.6f -> 1/(3L)=%.1f"%(
    list(stats[0][3]),float(stats[0][0]),float(F(1,3)/stats[0][0])))
print()
print("### worst L after removing TWO small killers exactly ###")
print("  scanning k1<k2 in (13*maxP, KB); reporting the worst over cores")
KB=900
random.seed(60)
worst=None
for P in CORES:
    iv=safe_set(P); M=max(P); lo=13*M+1
    ks=list(range(lo,KB))
    # full scan is 66 x ~280k pairs; sample densely plus the known-bad tail
    cand=ks[-60:]+random.sample(ks,min(120,len(ks)))
    for i in range(len(cand)):
        for j in range(i+1,len(cand)):
            k1,k2=cand[i],cand[j]
            if k1==k2: continue
            rem=remove_bad(remove_bad(iv,k1),k2)
            L=max((b-a for a,b in rem),default=F(0))
            if worst is None or L<worst[0]: worst=(L,tuple(P),k1,k2)
L,Pw,k1w,k2w=worst
print("  worst L = %.8f at core %s, killers (%d,%d)"%(float(L),list(Pw),k1w,k2w))
print("  => measure-horn threshold for the third killer: 1/(3L) = %.1f"%(float(F(1,3)/L) if L>0 else -1))
print()
print("### exhaustive worst-case on the WORST core only (all pairs) ###")
P=list(Pw); iv=safe_set(P); M=max(P)
w2=None
for k1 in range(13*M+1,KB):
    r1=remove_bad(iv,k1)
    if not r1: 
        w2=(F(0),k1,None); break
    for k2 in range(k1+1,KB):
        rem=remove_bad(r1,k2)
        L2=max((b-a for a,b in rem),default=F(0))
        if w2 is None or L2<w2[0]: w2=(L2,k1,k2)
print("  core %s : worst L over ALL pairs = %.8f at (%s,%s) -> threshold %.1f"%(
    P,float(w2[0]),w2[1],w2[2],float(F(1,3)/w2[0]) if w2[0]>0 else -1))
print("DONE")
