#!/usr/bin/env python3
"""LRC(14) covering MARGIN, lean+fast (kps 2026-06-27).
Q: does ANY covering 13-set have M < 1/13?  (conjecture C1: covering => M>=1/13)
Plus the aliasing table explaining why {1..12,14} is the minimizer (14==1 mod 13)."""
from fractions import Fraction as F
import random, math

def nrm(x):
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)
def nrm_f(x):
    x = x - math.floor(x); return min(x, 1.0-x)
def M_float(S):
    S = sorted(set(S)); n=len(S); c=set()
    for v in S:
        for k in range(v): c.add((2*k+1,2*v))
    for i in range(n):
        for j in range(i+1,n):
            a,b=S[i],S[j]
            for D in (a+b,b-a):
                if D>0:
                    for k in range(1,D//2+1): c.add((k,D))
    best=0.0
    for num,den in c:
        t=num/den
        if 0<t<=0.5:
            d=min(nrm_f(v*t) for v in S)
            if d>best: best=d
    return best
def M_exact(S):
    S=sorted(set(S)); n=len(S); c=set()
    for v in S:
        for k in range(v): c.add(F(2*k+1,2*v))
    for i in range(n):
        for j in range(i+1,n):
            a,b=S[i],S[j]
            for D in (a+b,b-a):
                if D>0:
                    for k in range(1,D//2+1): c.add(F(k,D))
    best=F(0)
    for t in c:
        if 0<t<=F(1,2):
            d=min(nrm(v*t) for v in S)
            if d>best: best=d
    return best
def is_cov(S): return any(v%14==0 for v in S)
T13=1.0/13

print("="*66); print("LRC(14) COVERING MARGIN (lean)"); print("="*66)
print(f"\n[S1] Any covering 13-set with M < 1/13={T13:.5f} ?")
lo=1.0; lo_set=None; below=[]; N=0
def consider(S):
    global lo,lo_set,N
    S=tuple(sorted(set(S)))
    if len(S)!=13 or not is_cov(S): return
    N+=1
    m=M_float(list(S))
    if m<lo-1e-12: lo=m; lo_set=S
    if m<T13-1e-9: below.append((m,S))

# structured families most likely to be near-tight:
# (a) GW-type {1..11,13,X}, X possibly a mult of 14
for X in range(12,70):
    consider(list(range(1,12))+[13,X])
# (b) dilated AP with one entry -> a multiple of 14
for d in range(1,5):
    base=[d*i for i in range(1,14)]
    for j in range(13):
        for mlt in range(14,14*25,14):
            consider(base[:j]+base[j+1:]+[mlt])
# (c) {1..12} (LRC13-tight) + any single extra speed that is a multiple of 14
for mlt in range(14,14*60,14):
    consider(list(range(1,13))+[mlt])
# (d) {1..k} consecutive core + fillers + one mult of 14, moderate range
rng=random.Random(4)
for _ in range(15000):
    mlt=rng.choice([14,28,42,56])
    rest=rng.sample([v for v in range(1,45) if v%14!=0],12)
    consider(rest+[mlt])

print(f"   examined {N} covering sets; global min M (float)={lo:.6f} at {lo_set}")
if below:
    print(f"   *** {len(below)} sets with M<1/13 (C1 REFUTED):")
    for m,S in sorted(below)[:8]:
        print(f"       M={M_exact(list(S))}={m:.5f}  {S}")
else:
    me=M_exact(list(lo_set))
    print(f"   NONE below 1/13. min exact = {me} = {float(me):.6f} at {lo_set}")
    print(f"   => C1 holds on search: covering => M >= 1/13, attained at {lo_set}")

print("\n[S3] Aliasing: R={1..12}*d is LRC13-tight (M=1/13). Add 14*j; stays>=1/13?")
for d in (1,2,3):
    R=[d*i for i in range(1,13)]
    print(f"   d={d}: M(R)={M_exact(R)} ; adding 14*j (j=1..13):")
    drop=[]
    for j in range(1,14):
        mlt=14*j
        if mlt in R: continue
        m=M_exact(R+[mlt])
        if m<F(1,13): drop.append((j,mlt,float(m)))
    if drop:
        print(f"      drops below 1/13 at: {drop}")
        print(f"      (14j == j mod 13; unsafe iff optimum kills it, i.e. 13|j)")
    else:
        print(f"      stays >= 1/13 for all j=1..13.")
print("="*66); print("DONE")
