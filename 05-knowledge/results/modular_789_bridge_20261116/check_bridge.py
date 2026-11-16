from itertools import combinations, permutations, product
from fractions import Fraction as F
from math import isqrt,prod

def prime(n):return n>=2 and all(n%i for i in range(2,isqrt(n)+1))
def genus(n):
    mu=F(n**3,2)
    for p in range(2,n+1):
        if n%p==0 and prime(p):mu*=1-F(1,p*p)
    g=1+mu*F(n-6,12*n)
    assert g.denominator==1
    return int(g)

assert [genus(n) for n in range(3,14)]==[0,0,0,1,3,5,10,13,26,25,50]
a=(102,252,248,241,225,192,128,1)
assert all(not(a[i]>>i&1) for i in range(8))
assert all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i,j in combinations(range(8),2))
assert sum(all(a[i]>>j&1 for i,j in zip(p,p[1:])) for p in permutations(range(8)))==63
p=(2,3,7,43)
assert {prod(S)+1 for k in range(1,5) for S in combinations(p,k) if prime(prod(S)+1)}=={3,7,43}

f=lambda x: 2*x*x*(x-2)+3
m=7;r=3
for k,A in enumerate((2,26,236,236),1):
    assert pow(3,A-1,m)==r and f(r)%m==0
    r=next(r+m*j for j in range(7) if f(r+m*j)%(7*m)==0)
    m*=7
for a in range(2,61):
    v=f(3**(a-1))
    assert v%72==(21 if a%2==0 else 57)
    assert (v%7==0)==(a%6==2)
    assert (v%49==0)==(a%42==26)
assert [j for j in range(7) if ((j*j-1)*(j-6)+36)%7==0]==[0]

def word(x,k,s):
    w=[]
    for _ in range(k):
        w.append(x%2)
        x=(3*x+s)//2 if x%2 else x//2
    return tuple(w),x
for s,k,q in product((-1,1),range(1,8),(1,3,7,9,11)):
    for r in range(q):
        T=[word(r+j*q,k,s) for j in range(1<<k)]
        assert len({w for w,_ in T})==1<<k
        for j,(w,y) in enumerate(T):
            w2,z=word(-(r+j*q),k,-s)
            assert w2==w and z==-y
print('verified genus table, H63, 7-adic lifts, prime support closure and +/- parity lifts')
