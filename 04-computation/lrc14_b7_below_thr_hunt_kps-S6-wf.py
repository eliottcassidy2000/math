#!/usr/bin/env python3
"""
FINAL adversarial check on the B7-minorant route. The route requires the TRUE
deterministic minorant B7(E) >= thr_k for ALL primitive E (0 in E, |E|=k).
Since B7(consec_8)=0.736 leaves only 0.118 margin to thr_8=0.6185 (NOT the 0.30+
the prompt claimed), hunt for ANY E with B7(E) < thr_k -- which would BREAK the
B7 lower-bound route directly (B7 too weak to reach thr). Also find the actual
B7-minimizer at k=8,9,10 over small spread, and stress large spread.

Note: B7(E) < thr_k does NOT refute LRC(14) (mu could still be >= thr; B7 is only
a lower bound). It would only show the B7 route is INSUFFICIENT. We ALSO track
whether mu itself ever drops below thr (that WOULD be a true counterexample).
"""
from fractions import Fraction as F
from itertools import combinations
import random

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

thr={8:F(3637,5880),9:F(2025,4004),10:F(36,91),11:F(25,91),12:F(1,7)}

print("=== B7-minimizer hunt + B7<thr / mu<thr detection ===")
b7_below=[]; mu_below=[]
# Exhaustive small spread for k=8 (the binding case)
print("-- k=8 exhaustive spread<=13 --")
best=(None,F(2)); cnt=0
for MAX in range(7,14):
    for combo in combinations(range(1,MAX),6):
        E=(0,)+combo+(MAX,); cnt+=1
        b=B7(list(E))
        if b<best[1]: best=(E,b)
        if b<thr[8]: b7_below.append((E,b))
print(f"  tested {cnt}; min B7={float(best[1]):.5f} at {best[0]}; thr_8={float(thr[8]):.5f}")
print(f"  any B7<thr_8? {len(b7_below)>0}  (these would only show B7-route too weak)")
# verify mu on the B7-minimizer
if best[0]:
    mm=mu_theta(list(best[0]),F(1,7))
    print(f"  mu at B7-minimizer {best[0]} = {float(mm):.5f}  >= thr_8? {mm>=thr[8]}")

# random large-spread descent minimizing B7 (k=8..12)
print("-- large-spread random descent minimizing B7, k=8..12 --")
random.seed(2026)
for k in range(8,13):
    bestk=(tuple(range(k)),B7(list(range(k))))
    for _ in range(25):
        W=random.choice([k+3, 2*k, 4*k, 60, 120])
        E=sorted(set([0]+random.sample(range(1,W+1),k-1)))
        if len(E)!=k: continue
        cur=B7(E); improved=True; steps=0
        while improved and steps<40:
            improved=False; steps+=1
            for idx in range(1,len(E)):
                for d in (-2,-1,1,2,-5,5):
                    cand=E[idx]+d
                    if cand<=0 or cand in E: continue
                    nE=sorted(set(E[:idx]+[cand]+E[idx+1:]))
                    if len(nE)!=k or nE[0]!=0: continue
                    bv=B7(nE)
                    if bv<bestk[1]: bestk=(tuple(nE),bv)
                    if bv<thr[k]: b7_below.append((tuple(nE),bv,k))
                    if bv<cur:
                        cur=bv; E=nE; improved=True; break
                if improved: break
    m=mu_theta(list(bestk[0]),F(1,7))
    print(f"  k={k}: min B7 found={float(bestk[1]):.5f} thr={float(thr[k]):.5f} "
          f"B7>=thr? {bestk[1]>=thr[k]}  | mu there={float(m):.5f} mu>=thr? {m>=thr[k]} "
          f"at E={bestk[0] if len(bestk[0])<=10 else str(bestk[0][:8])+'...'}")
    if m<thr[k]: mu_below.append((bestk[0],m,k))

print()
print("=== SUMMARY ===")
print(f"E with B7 < thr_k (B7-route insufficiency witnesses): {len(b7_below)}")
for x in b7_below[:8]: print("   ", x)
print(f"E with mu < thr_k (TRUE LRC(14) counterexamples): {len(mu_below)}")
for x in mu_below[:8]: print("   ", x)
if not mu_below:
    print("  >>> No true counterexample (mu>=thr held everywhere tested). LRC(14) target survives.")
