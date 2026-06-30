"""Thorough covering-min search at n=7: enumerate covering 6-sets (elements <=21), find true min M."""
import math
from fractions import Fraction
from itertools import combinations
def M_exact(S,Qmax):
    best=Fraction(0)
    for q in range(2,Qmax+1):
        bb=0
        for a in range(1,q):
            if math.gcd(a,q)!=1: continue
            m=q
            for s in S:
                r=(s*a)%q; d=r if r<=q-r else q-r
                if d<m:m=d
                if m<=bb:break
            if m>bb:bb=m
        v=Fraction(bb,q)
        if v>best:best=v
    return best
n=7; Phi6=n*n-n+1; constr=Fraction(n,Phi6); Q=2*Phi6
pool=list(range(1,22))
best=Fraction(1); bestS=None; below=[]
cnt=0
for combo in combinations(pool,n-1):
    S=list(combo)
    if not all(any(x%q==0 for x in S) for q in range(2,n+1)): continue
    cnt+=1
    M=M_exact(S,Q)
    if M<best: best=M; bestS=S
    if M<constr: below.append((S,M))
print(f"n=7: tested {cnt} covering 6-sets (elements<=21). construction=7/43={float(constr):.5f}")
print(f"   TRUE covering-min in this range = {best}={float(best):.5f} at {bestS}")
print(f"   construction beaten: {best<constr}  | #sets below construction: {len(below)}")
from collections import Counter
if below:
    dist=Counter(str(M) for _,M in below)
    print(f"   M-values below construction: {dict(dist)}")
    for S,M in sorted(below,key=lambda t:t[1])[:6]:
        print(f"      M={str(M):>7}={float(M):.5f}  S={S}")
print(f"   claimed beater 2/13={float(Fraction(2,13)):.5f}: achieved? {any(M==Fraction(2,13) for _,M in below)}")
