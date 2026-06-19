import sys
from fractions import Fraction
from functools import reduce
from math import gcd
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
print("=== A'_k = {1,...,k-2, k, 3(k-1)} : does it give 3/(3k+2)?  k mod 6 noted ===")
for k in range(5,26):
    S=tuple(sorted(set(list(range(1,k-1))+[k,3*(k-1)])))
    if len(S)!=k: 
        print(f"  k={k}: size {len(S)}!=k (collision), skip"); continue
    M,t=maxmin(S)
    floor=Fraction(1,k+1); med=Fraction(2,2*k+1); a3=Fraction(3,3*k+2)
    tag = "=3/(3k+2) DIP a=3" if M==a3 else ("=mediant a=2" if M==med else ("BELOW med" if M<med else "above med"))
    print(f"  k={k:2d} (k%6={k%6}): M(A'_k)={M} ({float(M):.5f})  3/(3k+2)={a3} med={med} [{tag}]  S={S if k<=14 else '...'}")
