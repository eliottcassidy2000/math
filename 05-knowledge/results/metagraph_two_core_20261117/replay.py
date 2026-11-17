from math import isqrt
from functools import cache

def need(x):
    if not x: raise AssertionError

def extreme(j):
    n=j+2
    A=[((1<<n)-1)^((2<<i)-1) for i in range(n)]
    A[0]^=1<<(n-1);A[-1]^=1
    return A

def read(A,B,v,u): # B is converse of displayed right input
    n=len(A)
    X=list(A)
    for i,r in enumerate(B):
        rev=sum(1<<j for j,t in enumerate(B) if t>>i&1)
        X.append((rev<<n)|((1<<n)-1))
    X[v]|=1<<(n+u);X[n+u]^=1<<v
    return X

def H(rows):
    @cache
    def d(S,v):
        R=S^(1<<v)
        return 1 if not R else sum(d(R,u) for u,r in enumerate(rows) if R>>u&1 and r>>v&1)
    return sum(d((1<<len(rows))-1,v) for v in range(len(rows)))

for j,k,p in [(3,2,13),(3,5,89),(6,5,613)]:
    A=extreme(j);W=read(A,extreme(k),1,0)
    n=len(A);sc=sorted(x.bit_count() for x in W)
    need(sum(sc[:n])-n*(n-1)//2==1 and sc[n-1]<sc[n])
    need(H(W)==7*p==(2**(j+1)+3)*(2**k+1)-2**(j-1))
    need(all(p%i for i in range(2,isqrt(p)+1)))
    print(j,k,p,W)
need(H([234,240,3,6,109,12,172,60])==343)
print('343 verified')
