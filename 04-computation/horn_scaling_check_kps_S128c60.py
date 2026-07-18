#!/usr/bin/env python3
"""horn_scaling_check_kps_S128c60.py -- kind-pasteur S128 cont.60.
COMPLETENESS CHECK for THM-1051 and the r=3 result.  Both measure horns were scanned with
the REMOVED killers bounded (874 / 900).  Beyond that bound the surviving interval L shrinks
like ~1/k, so the threshold 1/(3L) GROWS -- but so does the next killer.  The horns overlap
for ALL scales iff  1/(3L) < k_max_removed  (since the next killer exceeds it).
Measure the ratio  R = (1/(3L)) / k_max_removed  across scales.  R < 1 everywhere is what
the argument needs.  PRINT DATA ONLY."""
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
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P):
            if out and out[-1][1]==a: out[-1]=(out[-1][0],b)
            else: out.append((out[-1][0] if False else a,b))
    # merge adjacent
    mg=[]
    for a,b in out:
        if mg and mg[-1][1]==a: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return [(float(a),float(b)) for a,b in mg]
def remove_bad_f(iv,k):
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
random.seed(6001)
C11=[sorted(c) for c in itertools.combinations(range(1,13),11)]
C10=[sorted(c) for c in itertools.combinations(range(1,13),10)]
print("### r=2 : remove ONE killer k1, ratio R = (1/(3L)) / k1 ###")
print("  scale        samples   max R     (R < 1 => threshold below k1 < k2, horns overlap)")
for lo,hi in [(157,900),(900,3000),(3000,10000),(10000,60000),(60000,400000)]:
    mx=0; arg=None
    for _ in range(140):
        P=random.choice(C11); k1=random.randint(max(lo,13*max(P)+1),hi)
        iv=safe_set(P); rem=remove_bad_f(iv,k1)
        L=max((b-a for a,b in rem),default=0.0)
        if L<=0: mx=float('inf'); arg=(P,k1); break
        R=(1.0/(3*L))/k1
        if R>mx: mx=R; arg=(tuple(P),k1)
    print("  [%-7d,%-7d] %-9d %-9.4f %s"%(lo,hi,140,mx,"OK" if mx<1 else "*** R>=1"))
print()
print("### r=3 : remove TWO killers k1<k2, ratio R = (1/(3L)) / k2 ###")
print("  scale        samples   max R")
for lo,hi in [(131,900),(900,3000),(3000,10000),(10000,60000),(60000,400000)]:
    mx=0; arg=None
    for _ in range(120):
        P=random.choice(C10)
        base=max(lo,13*max(P)+1)
        k1=random.randint(base,hi); k2=random.randint(base,hi)
        if k1==k2: continue
        k1,k2=min(k1,k2),max(k1,k2)
        iv=safe_set(P); rem=remove_bad_f(remove_bad_f(iv,k1),k2)
        L=max((b-a for a,b in rem),default=0.0)
        if L<=0: mx=float('inf'); arg=(P,k1,k2); break
        R=(1.0/(3*L))/k2
        if R>mx: mx=R; arg=(tuple(P),k1,k2)
    print("  [%-7d,%-7d] %-9d %-9.4f %s   worst=%s"%(lo,hi,120,mx,"OK" if mx<1 else "*** R>=1",
        (list(arg[0]),arg[1],arg[2]) if arg and mx>=1 else ""))
print()
print("### adversarial: k1,k2 sharing structure (k2 = 2k1, k2 = k1+1, both highly composite) ###")
for tag,mk in [("k2=2k1",lambda k:2*k),("k2=k1+1",lambda k:k+1),("k2=k1+2",lambda k:k+2)]:
    mx=0; arg=None
    for _ in range(120):
        P=random.choice(C10); k1=random.randint(max(131,13*max(P)+1),200000); k2=mk(k1)
        iv=safe_set(P); rem=remove_bad_f(remove_bad_f(iv,k1),k2)
        L=max((b-a for a,b in rem),default=0.0)
        if L<=0: mx=float('inf'); arg=(P,k1,k2); break
        R=(1.0/(3*L))/max(k1,k2)
        if R>mx: mx=R; arg=(tuple(P),k1,k2)
    print("  %-10s max R = %-9.4f %s"%(tag,mx,"OK" if mx<1 else "*** R>=1 at %s"%(str(arg))))
print("DONE")
