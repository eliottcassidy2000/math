from functools import cache
from fractions import Fraction
from math import isqrt

witness = {3:('1',1), 5:('10',1), 11:('110',1), 13:('120',0),
    17:('1000',1),19:('1010',1),31:('1201',1),37:('1200',1),
    41:('3100',1),61:('100010',1),73:('11200',0),151:('33100',1),241:('112210',1)}

def count(word, ab):
    n=len(word); edges=[0]*(n+2)
    for i,z in enumerate(map(int,word)):
        for j in range(i+1,n):edges[i]|=1<<j
        for j,c in ((n,z//2),(n+1,z%2)):
            if c:edges[j]|=1<<i
            else:edges[i]|=1<<j
    a,b=(n,n+1) if ab else (n+1,n)
    edges[a]|=1<<b
    @cache
    def paths(m,i):
        rem=m^(1<<i)
        return sum(paths(rem,j) for j in range(n+2) if rem & edges[i] & (1<<j)) if rem else 1
    return sum(paths((1<<(n+2))-1,i) for i in range(n+2))

def prime(p):
    return p>=2 and all(p%i for i in range(2,isqrt(p)+1))

def order(p):
    return next(d for d in range(1,p) if pow(2,d,p)==1)

def odd(d):
    return d//(d&-d)

covers = [
    (72,66837283,[(3,2,0),(5,4,1),(17,8,3),(13,12,7),(241,24,15),
                 (73,9,5),(37,36,35),(19,18,11)]),
    (120,57013502549,[(3,2,1),(5,4,2),(17,8,0),(13,12,8),(241,24,12),
                     (31,5,1),(11,10,2),(41,20,8),(61,60,40),(151,15,4)])]

for p,rep in witness.items(): assert prime(p) and count(*rep)==p
for N,K,rows in covers:
    assert K%2 and K-1>21*max(p for p,_,_ in rows)
    assert all(p in witness and order(p)==d and N%d==0 and K*pow(2,r,p)%p==1
               for p,d,r in rows)
    assert all(any(e%d==r for p,d,r in rows) for e in range(N))
    print('Exact covered ray and physical labels',N,K)
orders={p:order(p) for p in range(3,242,2) if p!=7 and prime(p)}
eligible=[]
for N in range(1,72):
    ds=[d for d in orders.values() if N%d==0]
    assert sum(Fraction(1,d) for d in ds if odd(d)==1)<1
    if sum(Fraction(1,d) for d in ds)>=1 and sum(Fraction(1,odd(d)) for d in ds if odd(d)>1)>=1:
        eligible.append(N)
assert eligible==[48,60]
assert sorted(d for d in orders.values() if 48%d==0)==[2,4,8,12,24,48]
for N,cap in ((60,241),(120,240)):
    ds=[odd(d) for p,d in orders.items() if p<=cap and N%d==0]
    assert [ds.count(d) for d in (3,5,15)]==([1,3,2])
assert (3-1)*(5-3)>2
print('Minimal capped 7-free period 72; least maximum label for period 120 is 241.')
