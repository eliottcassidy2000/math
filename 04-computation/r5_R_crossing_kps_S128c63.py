#!/usr/bin/env python3
"""r5_R_crossing_kps_S128c63.py -- kind-pasteur S128 cont.63.
r=5 : find where R crosses BACK below 1.  The ladder (0.519, 0.734, 0.985) predicts R > 1
at r=5 for small killers, but R decays as killers grow (asymptotically ~1/6), so there is a
crossing.  The finite horn only has to cover killers BELOW the crossing.
Cores are 8-subsets of {1..12} (495).  R = T/k4 with T = min(N/(6mu), 1/(3L)).
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
C8=[sorted(c) for c in itertools.combinations(range(1,13),8)]
print("### r=5 : %d eight-speed cores ; S(P) profile ###"%len(C8))
st=[]
for P in C8:
    iv=safe_set(P); st.append((sum(b-a for a,b in iv),len(iv),max(b-a for a,b in iv),tuple(P)))
st.sort()
print("  measure: min %.5f  max %.5f ; components: %d..%d"%(
    st[0][0],st[-1][0],min(x[1] for x in st),max(x[1] for x in st)))
print()
print("### max R at the bottom window (all quadruples with k4 <= lo+W) ###")
W=26
worst=None
for P in C8:
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
                    T=thresh(rm(r3,ks[m]))
                    if T is None: continue
                    R=T/ks[m]
                    if worst is None or R>worst[0]: worst=(R,tuple(P),ks[i],ks[j],ks[h],ks[m],T)
R,Pw,a,b,c,d,T=worst
print("  MAX R = %.5f at core %s, killers (%d,%d,%d,%d), T = %.1f"%(R,list(Pw),a,b,c,d,T))
print("  R < 1 at the bottom: %s"%(R<1))
print()
print("### where does R cross below 1?  (worst core, killers scaled up together) ###")
P=list(Pw); iv0=safe_set(P)
print("  scale s   killers (a+s,b+s,c+s,d+s)     T          R")
cross=None
for s in [0,10,25,50,100,200,400,800,1600,3200]:
    ks=[a+s,b+s,c+s,d+s]
    iv=iv0
    ok=True
    for k in ks[:-1]:
        iv=rm(iv,k)
        if not iv: ok=False; break
    if not ok: continue
    T2=thresh(rm(iv,ks[-1]))
    if T2 is None: continue
    R2=T2/ks[-1]
    if cross is None and R2<1: cross=ks[-1]
    print("  %-9d %-30s %-10.1f %.5f"%(s,str(tuple(ks)),T2,R2))
print()
print("  first scanned scale with R < 1: k4 = %s"%cross)
print("  => the r=5 finite horn must cover killers below roughly that point")
print("DONE")
