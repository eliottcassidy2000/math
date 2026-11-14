from replay_forest import orders
from itertools import product
from functools import cmp_to_key

def price(q,p,i):
    return ((2*q,p+q),(p-q,2*p),(p+q,2*p))[i]
def inverse_price(q,p):
    digit=0 if q%2==0 else (1 if 2*q<p else 2)
    e,o=(q,p) if q%2==0 else (p,q)
    return tuple(sorted((e//2,abs(o-e//2)))),digit
def berggren(m,n,i):
    return ((n,2*n-m),(n,2*n+m),(m,2*m+n))[i]
def fold(fun,w):
    u=(1,2)
    for i in w:u=fun(*u,i)
    return u
def cmp(w,v):
    sign=1
    for i,j in zip(reversed(w),reversed(v)):
        if i!=j:return sign*(i-j)
        if i==1:sign=-sign
    return 0

from fractions import Fraction
from math import gcd
for k in range(1,8):
    W=list(product(range(3),repeat=k))
    P=[fold(price,w) for w in W]
    assert all(gcd(q,p)==1 and (p-q)%2 for q,p in P)
    assert len(set(P))==3**k
    assert all(inverse_price(*p)==(fold(price,w[:-1]),w[-1]) for p,w in zip(P,W))
    assert all((p-q==1)==(1 not in w) for (q,p),w in zip(P,W))
    assert sorted(W,key=cmp_to_key(cmp))==sorted(W,key=lambda w:Fraction(*reversed(fold(berggren,w))))
print('Price inverse, binary subtree, and Berggren reflected rank checked through depth seven')

W=(1,0,2,0,0)
anchors={(3,2,4,0):set(),(1,2,4,0):set()}
for p in orders(W):
    if p.index(5)>p.index(6):continue
    a,b=p.index(5),p.index(6)
    runs=(p[:a],p[a+1:b],p[b+1:])
    if not all(runs):continue
    R,S,T=runs
    u=(R[-1],S[0],S[-1],T[0])
    if u not in anchors:continue
    (free,)=set(range(5))-set(u)
    anchors[u].update(i for i,r in enumerate(runs) if free in r)
assert list(anchors.values())==[{0,2},{1,2}]
print('Verified literal same-tournament frames:',anchors)
