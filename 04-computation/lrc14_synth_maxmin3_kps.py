from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import random

# Use a guaranteed-correct but faster max-min: max_t min_v ||v t||.
# Equivalent: 1/14 - (lonely-set existence). Actually we only need to know if max-min < 1/14.
# max-min >= 1/14  <=>  there exists t with all ||v t|| >= 1/14 <=> the OPEN danger sets
# (||v t||<1/14) do NOT cover [0,1) ... but max-min is about closed >=1/14.
# CLEAN: max-min = max over breakpoints. Breakpoints of f(t)=min_v||vt|| are where some ||vt||
# hits 0 (t=k/v) or two sawtooths cross. Restrict candidates to crossings t=(i+j)/(a+b) and
# t=(i-j)/(a-b) AND t=k/v peaks. We just need the MAX so test all; but prune.

def maxmin(S):
    cset=set()
    Sl=sorted(S)
    for v in Sl:
        for k in range(v):
            cset.add(Fr(2*k+1,2*v))
    for i in range(len(Sl)):
        a=Sl[i]
        for j in range(i+1,len(Sl)):
            b=Sl[j]
            s=a+b; d=b-a
            for ka in range(a):
                # crossing up: a t - ka = -(b t - kb) => t=(ka+kb)/s ; and a t-ka=b t-kb => t=(kb-ka)/d
                for kb in range(b):
                    t=Fr(ka+kb,s)
                    if 0<t<1: cset.add(t)
                    if d!=0:
                        t2=Fr(kb-ka,d)
                        if 0<t2<1: cset.add(t2)
    best=Fr(0)
    for t in cset:
        m=Fr(1)
        ok=True
        for v in Sl:
            x=v*t; lo=x.numerator//x.denominator
            d1=x-lo; r=(lo+1)-x
            dv=d1 if d1<r else r
            if dv<m: m=dv
            if m<=best: ok=False; break
        if ok and m>best: best=m
    return best

print("AP       max-min =", maxmin(list(range(1,14))), "== 1/14?", maxmin(list(range(1,14)))==Fr(1,14))
print("sporadic max-min =", maxmin([1,2,3,4,5,6,7,8,9,10,11,13,24]), "== 1/14?", maxmin([1,2,3,4,5,6,7,8,9,10,11,13,24])==Fr(1,14))
print("12->36   max-min =", maxmin([1,2,3,4,5,6,7,8,9,10,11,13,36]))

random.seed(99)
target=Fr(1,14); mn=None; below=0; N=1500
for _ in range(N):
    S=sorted(random.sample(range(1,101),13))
    g=reduce(gcd,S); S=[x//g for x in S]
    mm=maxmin(S)
    if mn is None or mm<mn[0]: mn=(mm,tuple(S))
    if mm<target: below+=1
print(f"Random {N} 13-sets [1..100]: min max-min={mn[0]}={float(mn[0]):.6f}, #<1/14: {below}, min==1/14? {mn[0]==Fr(1,14)}")
