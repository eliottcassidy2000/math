"""Exact structural controls; not a Chan-Pak census."""
from itertools import combinations
from math import comb
import json

def chain(r,pad=0):
    n=3*r+pad
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for b in range(r):
        i=3*b;j=i+2
        rows[i]^=1<<j;rows[j]^=1<<i
    return rows

def h(rows,subset=None):
    if subset is None:subset=list(range(len(rows)))
    m=len(subset)
    if not m:return 1
    dp=[[0]*m for _ in range(1<<m)]
    for i in range(m):dp[1<<i][i]=1
    for mask in range(1,1<<m):
        for i in range(m):
            if dp[mask][i]:
                for j in range(m):
                    if not(mask>>j&1) and rows[subset[i]]>>subset[j]&1:
                        dp[mask|1<<j][j]+=dp[mask][i]
    return sum(dp[-1])

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

out=[]
for r,pad in [(1,0),(1,7),(2,0),(2,4),(3,0),(3,3)]:
    rows=chain(r,pad);n=len(rows)
    assert h(rows)==3**r
    poly=[1]
    for _ in range(r):poly=mul(poly,[1,3,3])
    for _ in range(pad):poly=mul(poly,[1,1])
    for k in range(3,min(6,n)+1):
        actual=sum(h(rows,S)==1 for S in combinations(range(n),k))
        predicted=poly[k] if k<len(poly) else 0
        assert actual==predicted
        failures=comb(n,k)-actual
        assert failures<=r*comb(n-3,k-3)
        out.append(dict(r=r,padding=pad,n=n,k=k,H=3**r,nontransitive_samples=failures,total=comb(n,k)))
print(json.dumps(dict(status='PASS',controls=out),indent=2))
