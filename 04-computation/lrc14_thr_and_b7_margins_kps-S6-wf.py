#!/usr/bin/env python3
"""
Compute thr_k EXACTLY for k=8..13 as thr_k = 1 - min_{|P|=13-k} meas(G_P),
where P ranges over size-(13-k) subsets of the SMALL primes/integers that can
appear (P = small part of a primitive covering 13-set; p in {1..13}, but p=1
forces ||x||>=1/14 nontrivial; the small co-prime structure). Then check
whether the TRUE (deterministic, exact) B7 minorant at consec_k is >= thr_k for
ALL k=8..12 (the route's actual requirement), and report the real margins.

We DON'T know the exact admissible P-universe from the prompt, but thr_k is the
'1 - min meas(G_P)' and the prompt PROVIDES thr_8=3637/5880 and thr_12=1/7.
We reverse-engineer: min_{|P|=5} meas(G_P) = 1-thr_8 = 2243/5880; and for k=12,
|P|=1: min_p meas(G_{p}) = 1-1/7 = 6/7, achieved at p=7 (||7x||>=1/14 keeps 6/7).
We compute meas(G_P) over P subsets of {1,...,13} (the natural small universe for
a covering 13-set, since elements are 1..13 = "small" part) of each size 13-k and
take the min, to get thr_k exactly, then compare to B7(consec_k).

KEY QUESTION for the route's viability: is B7(consec_k) >= thr_k for all k=8..12?
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

def meas_GP(P):
    P=sorted(set(p for p in P if p!=0))
    if not P: return F(1)
    bp=set([F(0),F(1)])
    for p in P:
        for m in range(0,p+1):
            bp.add(F(14*m+1,14*p)); bp.add(F(14*m+13,14*p)); bp.add(F(m,p))
    bp=sorted(b for b in bp if 0<=b<=1); tot=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; ok=True
        for p in P:
            t=(p*mid)%1; nd=min(t,1-t)
            if nd<F(1,14): ok=False;break
        if ok: tot+=(b-a)
    return tot

# B7 exact
def some_fixed_arc_empty(E,x):
    arcs=[(F(2*i+1,14),F(2*i+3,14)) for i in range(7)]
    pts=[(e*x)%1 for e in E]
    for (lo,hi) in arcs:
        has=False
        for p in pts:
            if hi<=1:
                if lo<=p<hi: has=True;break
            else:
                if p>=lo or p<(hi-1): has=True;break
        if not has: return True
    return False
def b7_bps(E):
    E=sorted(set(E)); bp=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,e+1):
            for i in range(7):
                bp.add((F(m)+F(2*i+1,14))/e); bp.add((F(m)+F(2*i+3,14))/e)
    return sorted(b for b in bp if 0<=b<=1)
def B7(E):
    bp=b7_bps(E); t=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        if some_fixed_arc_empty(E,(a+b)/2): t+=(b-a)
    return t

# mu engine
def mu_theta(E, theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[(E[order[t]]*mid).__floor__() for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t];o2=order[(t+1)%n];k1=ks[t];k2=ks[(t+1)%n];wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1];c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s);  subs.append((lo,b)) if lo<b else None
            else:
                hi=min(b,(theta-c)/s);  subs.append((a,hi)) if a<hi else None
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

print("=== thr_k = 1 - min_{|P|=13-k, P subset {1..13}} meas(G_P) ===")
thr={}
for k in range(8,14):
    sz=13-k
    if sz==0:
        thr[k]=F(0)
        print(f"k={k}: |P|=0 -> thr=0")
        continue
    best=None; argbest=None
    for P in combinations(range(1,14), sz):
        m=meas_GP(P)
        if best is None or m<best:
            best=m; argbest=P
    thr[k]=1-best
    print(f"k={k}: |P|={sz}  min meas(G_P)={best}={float(best):.5f} at P={argbest}  "
          f"thr_k=1-min={thr[k]}={float(thr[k]):.5f}")

print()
print("=== ROUTE CHECK: B7(consec_k) >= thr_k for k=8..12 (the route's requirement)? ===")
print("   (B7 = deterministic exact minorant; thr from above)")
allok=True
for k in range(8,13):
    E=list(range(k))
    b=B7(E); t=thr[k]; m=mu_theta(E,F(1,7))
    ok = b>=t
    allok &= ok
    print(f"k={k}: B7(consec)={float(b):.5f}  thr_k={float(t):.5f}  "
          f"B7-thr={float(b-t):+.5f}  {'OK B7>=thr' if ok else '**B7 < thr: B7-route FAILS even at consec**'}")
    print(f"      (mu(consec)={float(m):.5f} >= thr? {m>=t}  margin {float(m-t):+.5f})")
print()
print(f"B7-minorant route survives at consec for all k=8..12? {allok}")
print("NOTE: even if B7(consec)>=thr, the route also needs B7(E)>=thr for ALL E,")
print("AND consec to be (near) the B7-minimizer; and a uniform analytic bound.")
