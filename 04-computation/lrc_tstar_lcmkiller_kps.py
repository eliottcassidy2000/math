import sys
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
def fn(x):
    r=x-(x.numerator//x.denominator); return min(r,1-r)
def maxmin(S):
    S=sorted(set(S)); cand=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for d in (S[i]+S[j],S[i]-S[j]):
                if d==0: continue
                d=abs(d)
                for m in range(1,d): cand.add(Fraction(m,d))
    best=Fraction(0); bt=Fraction(0)
    for t in cand:
        mv=min(fn(v*t) for v in S)
        if mv>best: best,bt=mv,t
    return best,bt
print("=== optimal time t* and binding pair for F(k,a)={1..k-2,k,a(k-1)} at dip-k ===")
for (a,k) in [(2,5),(3,7),(3,13),(4,31)]:
    S=sorted(set(list(range(1,k-1))+[k,a*(k-1)]))
    M,t=maxmin(tuple(S))
    binding=[v for v in S if fn(v*t)==M]
    print(f"  a={a} k={k:3d}: M={M} t*={t}  binding speeds={binding}  (W=a(k-1)={a*(k-1)}, q={M.denominator})")
print()
print("=== E1 own-check: can an lcm-KILLER config beat the primorial family at smaller k? ===")
print("    {1,...,k-1} + W where W=lcm(1..m) kills clocks t=j/d for ALL d<=m")
for k in [11,13,16,21]:
    floor=Fraction(1,k+1); med=Fraction(2,2*k+1)
    best=(med, None)
    # try W = lcm(1..m)*c for various m,c, replacing nothing (use {1..k-1}+W) -> k speeds
    for m in range(2, k+2):
        L=lcm(*range(1,m+1))
        for c in range(1, 6):
            W=L*c
            if W in range(1,k): continue
            S=tuple(sorted(set(list(range(1,k))+[W])))
            if len(S)!=k: continue
            M,t=maxmin(S)
            if floor<M<best[0]: best=(M,S,W,m)
    if best[1]:
        M=best[0]; g=M-floor
        a = M.numerator if (M.numerator*(k+1)-M.denominator)==1 else None
        print(f"  k={k}: best lcm-killer M={M} g*k^2={float(g)*k*k:.4f} a={a} W={best[2]}(=lcm(1..{best[3]})) S={best[1]}")
    else:
        print(f"  k={k}: no lcm-killer below mediant (best stays {med})")
