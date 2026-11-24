from collections import Counter
from itertools import product
from math import prod

def runs(bits):
    t=0
    while t<len(bits):
        if not bits[t]:t+=1;continue
        u=t
        while t<len(bits) and bits[t]:t+=1
        yield u,t

def poly(gaps):
    ans=Counter()
    for bits in product(range(2), repeat=len(gaps)):
        dim=0
        for u,v in runs(bits):
            g=sum(gaps[u:v])
            if not g and (v-u)%2:break
            dim+=g or 1
        else:ans[dim]+=1
    return ans

def rows(gaps):
    n=1+len(gaps)+sum(gaps)
    A=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    i=0
    for g in gaps:
        j=i+1+g;A[i]^=1<<j;A[j]^=1<<i;i=j
    return A

def H(A):
    n=len(A);N=1<<n
    D=[[0]*n for _ in range(N)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,N):
        for i,c in enumerate(D[s]):
            if not c:continue
            w=A[i]&~s
            while w:
                b=w&-w;w-=b;D[s|b][b.bit_length()-1]+=c
    return sum(D[-1])

for t in range(7):
    for g in product(range(3),repeat=t):
        if 1+t+sum(g)>10:continue
        P=poly(g)
        assert H(rows(g))==sum(v*2**k for k,v in P.items())
        if all(g):
            for q in (1,2,3,4,5):
                assert sum(v*q**k for k,v in P.items())==prod(1+q**h for h in g)

def K(t,q):
    def one(a,b):
        ans=0
        for bits in product(range(2),repeat=t):
            p=1
            for u,v in runs((a,)+bits+(b,)):
                if u and v<t+2:p*=q if (v-u)%2==0 else 0
            ans+=p
        return ans
    return [[one(a,b) for b in (0,1)] for a in (0,1)]

def det(a):return a[0][0]*a[1][1]-a[0][1]*a[1][0]
F=[0,1]
for i in range(15):F.append(sum(F[-2:]))
for t in range(12):
    assert K(t,1)==[[F[t+1],F[t+2]],[F[t+2],F[t+3]-1]]
    assert det(K(t,1))==(-1)**t-F[t+1]
for q in range(1,10):
    assert K(3,q)==[[1+2*q,4+q],[4+q,7]]
    assert det(K(3,q))==-(q-3)**2

assert H(rows((1,3)))==27 and H(rows((2,2)))==25
assert sum(poly((1,3)).values())==sum(poly((2,2)).values())==4
for g in range(1,8):
    pieces=[g];counts=[]
    while True:
        assert sum(2**t for t in pieces)==2**g
        counts.append(len(pieces)+1) # one extra empty-cycle stratum
        if max(pieces)==0:break
        d=max(pieces);pieces.remove(d);pieces.extend([d-1,d-1])
    assert counts==list(range(2,2**g+2))
print('checked point polynomials, Fibonacci determinants, and both nonrecoverability families')
