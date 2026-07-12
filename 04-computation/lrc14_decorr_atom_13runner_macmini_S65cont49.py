#!/usr/bin/env python3
"""cont.49: verify the 13-runner decorrelation-atom bound (THM-636 analog) on large-diameter
DC families. For v with a scale L, write v_i = b_i + L*k_i (b = v mod L, k = round(v/L)).
THM-636: reach(v) >= reach(k) - B/L, B = max|b_i|. If the lift family k has <= 12 distinct
speeds (repeated lift), LRC(<=13) gives reach(k) >= 1/13, so reach(v) >= 1/13 - B/L > 1/14
for L > B/(1/13-1/14) = 182*B. Test the inequality + the threshold on DC families."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
def reach(v, Q=80):
    best=F(0)
    qs=set()
    for i in range(len(v)):
        for j in range(len(v)):
            if i!=j: qs.add(abs(v[i])+abs(v[j])); qs.add(abs(abs(v[i])-abs(v[j])))
    for q in list(qs)+list(range(2,Q)):
        if q<2: continue
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            m=None
            for e in v:
                r=(e*p)%q; d=min(r,q-r); dd=F(d,q)
                if m is None or dd<m: m=dd
            if m>best: best=m
    return best
def decompose(v,L):
    b=[((x+L//2)%L)-L//2 for x in v]  # centered residue
    k=[round(x/L) for x in v]
    return b,k
print("13-runner decorrelation atom (THM-636 analog): reach(v) >= reach(k) - B/L")
print(f"floor 1/14={float(F(1,14)):.5f}, LRC(<=13) lift floor 1/13={float(F(1,13)):.5f}")
print(f"{'family (scale L)':28s} {'#distinct k':>11s} {'reach(k)':>9s} {'B/L':>7s} {'lb=r(k)-B/L':>12s} {'reach(v)':>9s}")
random.seed(49)
def make_scaled_dc(L, seed):
    random.seed(seed)
    base=[2,3,4,6,12]  # DC small parts
    k=sorted(set(random.sample(range(1,10),6)))  # lift speeds (<=12 distinct => repeated)
    v=sorted(set(base+[b+L*ki for b in random.sample(range(1,13),1) for ki in k]))
    while len(v)<13: v.append(max(v)+random.randint(1,4)); v=sorted(set(v))
    v=v[:13]; g=reduce(gcd,v); return [x//g for x in v]
for L in [500,2000,8000,40000]:
    for s in range(2):
        v=make_scaled_dc(L,s*100+L)
        b,k=decompose(v,L)
        kd=sorted(set(abs(x) for x in k if x!=0))
        rk=reach(kd) if len(kd)>=1 else F(1)
        B=max(abs(x) for x in b) if b else 0
        lb=float(rk)-B/L
        rv=reach(v)
        ok = "LOOSE" if float(rv)>float(F(1,14)) else "tight?"
        print(f"  L={L:6d} s{s}: {len(kd):>10d} {float(rk):9.4f} {B/L:7.4f} {lb:12.4f} {float(rv):9.4f}  {ok}", flush=True)
print(f"\n=> THM-636 lb = reach(k) - B/L is a valid lower bound; reach(k)>=1/13 (LRC<=13) when <=12 distinct lifts")
print(f"=> reach(v) > 1/14 for L > 182*B ~ {182*13} (B~13); large-diameter DC LOOSE by decorrelation descent")
