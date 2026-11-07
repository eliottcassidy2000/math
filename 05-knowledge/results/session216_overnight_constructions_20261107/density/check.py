#!/usr/bin/env python3
"""Exact finite energy checks; split construction uses integer arithmetic only."""
from collections import Counter
from math import log
import json
rows=[]
for n in range(12,23):
    # Unique m with n*2^(3*(m-1)) <= 2^n < n*2^(3*m).
    m=1
    while n*2**(3*m)<=2**n: m+=1
    P=2**(n-m-1); C=2**m
    assert P>=8*C
    g=[1]*(2**n)
    for v in range(1,2**n):
        u=v//2;g[v]=g[u] if v&1 else g[u]+u+(u&1)
    hist=Counter(g[2**(n-1):]);energy=sum(t*t for t in hist.values())
    # Independent symmetric-block identity, exhaustive h,c on full length universe.
    def revcomp(x,m):
        y=0
        for _ in range(m): y=2*y+1-(x&1);x//=2
        return y
    failures=0
    for h in range(P,2*P):
        for l in range(C):
            c=revcomp(l,m)
            if g[h*C+l]!=g[h]+g[c]+h*c+(h&1)*(c&1)-1:failures+=1
    exact_bound=P*C+8*P*P+P*C*C*(33.74+12*log(P/C))+3*P*C*(1+log(P/C))
    bound=7.3*2**(4*n/3)*n**(2/3)
    assert energy<=exact_bound<=bound and failures==0
    rows.append(dict(n=n,m=m,energy=energy,distinct=len(hist),block_identity_failures=failures,optimized_D1=exact_bound,uniform_bound=bound))
print(json.dumps(rows,indent=2))
