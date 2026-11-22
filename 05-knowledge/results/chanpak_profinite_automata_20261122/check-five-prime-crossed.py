"""Exact transfer and simultaneous-phase replay for the 110011100 1^k crossed rail."""
from collections import defaultdict
from itertools import product

def stream(tags, swapped, forward_arc):
    d={(True,0,False):1}
    for p,q in tags:
        if swapped:p,q=q,p
        e=defaultdict(int)
        for (left,mid,right),n in d.items():
            e[not p,mid,right]+=n
            if mid or p:e[left,1 if q else 2,right]+=n
            if right or q:e[left,mid,True]+=n
        d=e
    return sum(n for (left,mid,right),n in d.items()
               if left and (mid==2 or (mid==0 and forward_arc)))

def counts(k,arc):
    W=tuple(map(int,"10110110"+"1"*k+"0011"))
    tags=tuple((v,v) for v in W)+((0,0),(1,0),(0,1),(0,0))+tuple((1-v,1-v) for v in reversed(W))
    return stream(tags,False,arc)+stream(tags,True,not arc)

def midpoint(X,Z):
    return 3412542-982640804*X+2328668912*Z+71433464766*X*X-1453657649592*X*Z+8410552316658*Z*Z


def crt(a,b,m,n):
    return a+m*((b-a)*pow(m,-1,n)%n) if n>1 else a

def phase(s,t,u,r,j,R):
    assert R%64==0 and R%(9*25*7*11)==1
    W=0
    for i in range(s-6):
        if ((10644932798157+37848067592189*W+33642209266632*W*W-R//64)>>i)&1:W+=1<<i
    want=9+16*W
    v2=0
    for i in range(4,s):
        m=1<<(i+1)
        if ((pow(3,22+60*v2,m)-want)>>i)&1:v2+=1<<(i-4)
    def lift(p,start,precision,onlyx=False):
        v=0
        for i in range(start,precision):
            m=p**(i+1)
            v=next(a for d in range(p)
                   if midpoint(pow(2,22+60*(a:=v+d*p**(i-start)),m),
                       0 if onlyx else pow(3,22+60*a,m))%m==R%m)
        return v
    v3,v5,v7,v11=lift(3,2,t,True),lift(5,2,u),lift(7,1,r),lift(11,1,j)
    V,M=0,1
    for a,m in ((v2,2**(s-4)),(v3,3**(t-2)),(v5,5**(u-2)),
                (v7,7**(r-1)),(v11,11**(j-1))):
        V,M=crt(V,a,M,m),M*m
    return 22+60*(V+M*max(s,t,u,r,j))

E=2**60*(2**61-1)
for k in range(1,10):
    M=midpoint(2**k,3**k)
    assert (counts(k,True),counts(k,False))==(M-1,M+1)
for p,start,a,b in ((5,2,1,51),(7,1,8,36),(11,1,56,100)):
    m=p**(start+1)
    assert [midpoint(pow(2,k,m),pow(3,k,m))%m for k in (22,82)]==[a,b]
vals={midpoint(pow(2,k,19),pow(3,k,19))%19 for k in range(18)}
assert E%19==15 and E%19 not in vals
for s,t,u,r,j in product((6,13,31),(2,8),(2,9),(1,6),(1,7)):
    for R in (E,E+64*9*25*7*11):
        k=phase(s,t,u,r,j,R)
        for m in (2**s,3**t,5**u,7**r,11**j):
            assert midpoint(pow(2,k,m),pow(3,k,m))%m==R%m
print("simultaneous five-prime phase replay PASS; E61 excluded at 19")
