from itertools import product
from math import ceil,log2
from cmath import exp,pi

def K(d):
    return (1<<(d-1))+(d==1)

def positions(L,d):
    U=[];V=[]
    for k in range(d):
        a=list(range(k,L,d))
        for i,j in zip(a[::2],a[1::2]):U.append(i);V.append(j)
    return U,V

# Exact column cancellation certificate, with the earliest changed V bit.
for j in range(1,5):
    d=j+1;mod=1<<j
    for L in range(d+1,2*d+4):
        U,V=positions(L,d)
        assert 2*len(U)>=L-d
        for w in product((-1,0,1),repeat=len(V)):
            if not any(w):continue
            v=min(i for i,b in zip(V,w) if b)
            u=v-d
            assert u in U
            assert sum((K(abs(i-u))//2)*b for i,b in zip(V,w))%mod==mod//2

def H(b,n):
    S=[i for i in range(n) if b>>i&1]
    rev=int(f'{b:0{n}b}'[::-1],2)
    J=sum(K(j-i) for l,i in enumerate(S) for j in S[l+1:])
    return -11+82*((1<<(n+1))-rev)+9*(2*b-J)

# Small-scale check of full CRT support under fixed exterior bits.
for k,q,L,n,lo,hi in [(2,3,12,16,1,3),(3,5,14,18,3,1),(4,3,15,19,0,2)]:
    S={((H((hi<<(L+2))+(b<<2)+lo,n)-1)//2%(1<<k),
        ((hi<<(L+2))+(b<<2)+lo)%q) for b in range(1<<L)}
    assert len(S)==(1<<k)*q
print('Exact triangular column cancellations and pinned-exterior CRT controls pass.')
