import sys
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
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
    best=Fraction(0)
    for t in cand:
        mv=min(fn(v*t) for v in S)
        if mv>best: best=mv
    return best
# F(k,a) = {1,...,k-2, k, a(k-1)}.  When does M = a/(a(k+1)-1) (the level-a mediant)?
print("=== F(k,a)={1..k-2,k,a(k-1)} : k values where M = a/(a(k+1)-1) (level-a dip) ===")
for a in range(2,8):
    hits=[]
    for k in range(a+2, 160):
        S=sorted(set(list(range(1,k-1))+[k,a*(k-1)]))
        if len(S)!=k: continue
        M=maxmin(tuple(S))
        target=Fraction(a, a*(k+1)-1)
        if M==target:
            hits.append(k)
    # detect modular pattern
    diffs=[hits[i+1]-hits[i] for i in range(len(hits)-1)] if len(hits)>1 else []
    print(f"  a={a}: k with M=a/(a(k+1)-1): {hits[:20]}  (gaps {sorted(set(diffs))[:5]})")
print()
print("=== so the SMALLEST k achieving level a (=> a_max can reach a) ===")
for a in range(2,8):
    first=None
    for k in range(a+2,200):
        S=sorted(set(list(range(1,k-1))+[k,a*(k-1)]))
        if len(S)!=k: continue
        if maxmin(tuple(S))==Fraction(a,a*(k+1)-1):
            first=k; break
    if first:
        gk2=float(Fraction(a,a*(first+1)-1)-Fraction(1,first+1))*first*first
        print(f"  a={a}: first k={first}  (k%6={first%6}, k%2={first%2})  g*k^2={gk2:.4f} -> ~1/a={1/a:.4f}")
