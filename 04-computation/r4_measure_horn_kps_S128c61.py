#!/usr/bin/env python3
"""r4_measure_horn_kps_S128c61.py -- kind-pasteur S128 cont.61 (background).
r=4 MEASURE HORN: remove the three smaller killers from S(P) exactly; the fourth needs
k4 > 1/(3L).  Core is now a 9-subset of {1..12} (220 cores) -- LARGER safe set than r=3.
Also size the finite horn that must follow.  PRINT DATA ONLY."""
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
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((a,b))
    mg=[]
    for a,b in out:
        if mg and mg[-1][1]==a: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return [(float(a),float(b)) for a,b in mg]
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
CORES=[sorted(c) for c in itertools.combinations(range(1,13),9)]
print("### r=4 : %d nine-speed cores ###"%len(CORES))
st=[]
for P in CORES:
    iv=safe_set(P)
    st.append((max(b-a for a,b in iv),sum(b-a for a,b in iv),tuple(P)))
st.sort()
print("  S(P) largest component: min %.6f  median %.6f  max %.6f"%(st[0][0],st[len(st)//2][0],st[-1][0]))
print("  S(P) total measure:     min %.6f  max %.6f"%(min(x[1] for x in st),max(x[1] for x in st)))
print()
print("### worst L after removing THREE killers exactly ###")
random.seed(61); KB=900; TOP=45
worst=None
for P in CORES:
    ivF=safe_set(P); M=max(P); lo=13*M+1
    ks=list(range(lo,KB))
    pool=sorted(set(ks[-TOP:]+random.sample(ks,min(35,len(ks)))))
    w=None
    for i in range(len(pool)):
        r1=rm(ivF,pool[i])
        if not r1: w=(0.0,pool[i],None,None); break
        for j in range(i+1,len(pool)):
            r2=rm(r1,pool[j])
            if not r2: continue
            for h in range(j+1,len(pool)):
                r3=rm(r2,pool[h])
                L=max((b-a for a,b in r3),default=0.0)
                if w is None or L<w[0]: w=(L,pool[i],pool[j],pool[h])
    if w and (worst is None or w[0]<worst[0]): worst=(w[0],tuple(P),w[1],w[2],w[3])
L,Pw,a,b,cc=worst
print("  GLOBAL WORST L = %.8f at core %s killers (%s,%s,%s)"%(L,list(Pw),a,b,cc))
thr = 1.0/(3*L) if L>0 else -1
print("  => r=4 measure-horn threshold for the fourth killer: 1/(3L) = %.1f"%thr)
print()
print("### sizing the r=4 finite horn at KB = %d ###"%(int(thr)+1))
KB4=int(thr)+1
tot=0
for P in CORES:
    M=max(P); lo=13*M+1
    if lo>=KB4: continue
    rng=list(range(lo,KB4))
    A=[k for k in rng if k%13==0]; B2=[k for k in rng if k%14==0]
    tot+=len(A)*len(B2)*len(rng)*len(rng)//4
print("  crude upper bound on quadruples to test: ~%.2e  (covering forces 13| and 14| among killers)"%tot)
print("DONE")
