#!/usr/bin/env python3
"""Independent finite replay of endpoint phases and additive long-guard H."""
from itertools import product
from random import Random

def literal(w,m):
    n=len(w); h=1; q=pow(2,n+1,m)
    for i,a in enumerate(w):
        if a!='O':continue
        for j in range(i+1,n):
            if w[j]=='I':
                d=j-i-1
                h+=pow(2,d,m)+(d==0)
                q+=pow(2,n-j+i,m)*(pow(3,d,m)+1)
    return h%m,q%m
def boundary(w,k):
    q=1<<k;n=len(w)
    def E(side,x):
        return sum(c==z and pow(x,i,q) or 0 for i,(c,z) in enumerate(side))%q
    A=list(zip(w[:k],'O'*k));B=list(zip(w[::-1][:k],'I'*k))
    X=2*pow(3,-1,q)%q
    return 2*(pow(3,n-2,q)*E(A,X)*E(B,X)+E(A,2)*E(B,2))%q

rnd=Random(11);cases=0
for k in range(1,10):
    mod=1<<k;g=max(k,4)
    for n in range(2*k+4,2*k+36):
        for _ in range(8):
            w='O'+''.join(rnd.choice('OI') for j in range(n-2))+'I'
            H,Q=literal(w,mod)
            assert Q==boundary(w,k)
            assert (Q-4)%min(16,mod)==0 if n%2==0 else Q%min(8,mod)==0
            cases+=1
    C='O'*g+'I'*g; tiles=[C,'O'+C,C+'I']
    for m in range(1,5):
        for word in product(tiles,repeat=m):
            assert literal(''.join(word),mod)[0]==(1+2*m)%mod
    for parity in range(2):
        got=set()
        for n in range(2*k+4+parity,2*k+4+parity+2*mod,2):
            first=['O'+'I'*(k-1)]
            if parity and k>=2:first.append('OO'+'I'*(k-2))
            last='O'*(k-1)+'I'
            for left in first:
                w=left+'O'*(n-2*k)+last
                got.add(boundary(w,k))
        expect={q for q in range(mod) if ((q-4)%min(16,mod)==0 if parity==0 else q%min(8,mod)==0)}
        assert got==expect
    n=2*k+8
    hi='O'*k+'O'*(n-2*k)+'I'*k
    lo='O'+'I'*(k-1)+'O'*(n+2-2*k)+'O'*(k-1)+'I'
    assert boundary(hi,k)==boundary(lo,k)==2*(pow(3,n,mod)+1)%mod
print(cases,'literal endpoint phases; all image cosets, long-guard phases, and near-power boundary coincidences through 2^9')
