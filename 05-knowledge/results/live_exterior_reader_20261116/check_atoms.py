from itertools import combinations, permutations, product
from functools import cache

@cache
def hp(a, mask):
    ids = [i for i in range(len(a)) if mask>>i&1]
    return sum(all(a[i]>>j&1 for i,j in zip(p,p[1:])) for p in permutations(ids))

def slots(a,t):
    F=(1<<len(a))-1
    total=0
    for word in product(range(t),repeat=len(a)):
        weight=1
        for j in range(t):
            M=sum(1<<i for i,k in enumerate(word) if k==j)
            weight *= hp(a, M)
        total += weight
    return total

def stats(a,r):
    F=(1<<len(a))-1
    R=S=0
    for M in range(1,F+1):
        if not M>>r&1: continue
        ids=[i for i in range(len(a)) if (M^(1<<r))>>i&1]
        count=sum(all(a[i]>>j&1 for i,j in zip((r,)+p,p)) for p in permutations(ids))
        count*=hp(a,F^M)
        if M.bit_count()%2: S+=count
        else: R+=count
    return hp(a,F),hp(a,F^(1<<r)),R,S

A=((14,12,8,0),(14,4,8,2),(10,12,9,0),(10,4,9,2))
for a,c in zip(A,(0,1,1,2)):
    for t in range(1,8):
        assert slots(a,t)==t**4+2*c*t**2
assert [stats(a,1) for a in A[1:3]]==[(3,1,1,2),(3,1,3,3)]

def rev(b,n):
    return int(f'{b:0{n}b}'[::-1],2)
def H(b,n):
    sites=[i for i in range(n) if b>>i&1]
    J=sum((1<<(j-i-1))+(j-i==1) for i,j in combinations(sites,2))
    return -11+82*((2<<n)-rev(b,n))+9*(2*b-J)
for n in range(2,15):
    m=n//2
    B=1<<m
    k=n-m
    C=1<<k
    for q in range(m):
        L=1<<q
        for U in range(0,C,L):
            for X in range(0,B,L):
                b=B*U+rev(X,m)
                assert rev(b,n)==C*X+rev(U,k)
                H0=-11+(82*C+9*U)*(2*B-X)
                assert abs(H0-H(b,n))<=220*B//L+9*n
print('The atom polynomials, root profiles and paired-address approximation replayed.')
