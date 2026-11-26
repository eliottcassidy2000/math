from itertools import combinations_with_replacement, product
from math import isqrt

P=(5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61)
absent={
5:{(0,b) for b in range(1,5)},
7:set(product((0,3,5,6),(2,5))),
11:{(q,b%11) for q,n in ((1,5),(3,2),(4,4),(5,1),(9,3)) for b in (n,-n)},
13:{(0,b%13) for b in (2,3,-2,-3)}
}
bad0={7,19,23,47}

def F(x): return (x**3-x)//6
def moments(y):
    assert sum(y)==0
    return sum(v*v for v in y)//2,sum(v**3 for v in y)//6

def table(p):
    R={}
    sq=[i*i for i in range(p)]; cb=[i**3 for i in range(p)]
    i2,i6=pow(2,-1,p),pow(6,-1,p)
    for a,b,c,d in combinations_with_replacement(range(p),4):
        e=(-a-b-c-d)%p
        if d>e: continue
        y=(a,b,c,d,e)
        ij=(sum(sq[v] for v in y)*i2%p,sum(cb[v] for v in y)*i6%p)
        R[ij]=R.get(ij,False) or len(set(y))>=3
    return R

for p in P:
    tab=table(p)
    assert set(product(range(p),repeat=2))-tab.keys()==absent.get(p,set())
    assert {k for k,v in tab.items() if not v}==({(0,0)} if p in bad0 else set())
    assert all(any(tab.get((q,(r-c*q)%p),False) for q in range(p)) for r,c in product(range(p),repeat=2))
    print(p, 'absent',len(absent.get(p,set())), 'singular-only',p in bad0)

def bring_count(p):
    n=0
    i2,i6=pow(2,-1,p),pow(6,-1,p)
    for i in range(4):
        for tail in product(range(p),repeat=3-i):
            x=(0,)*i+(1,)+tail
            y=(*x,-sum(x))
            q,b=moments(y)
            n+=not(q%p or b%p)
    return n
for p in P[1:]:
    n=bring_count(p)
    ne=1+sum((y*y+x*y+y-x**3+126*x+552)%p==0 for x,y in product(range(p),repeat=2))
    assert n==p+1-4*(p+1-ne)
    assert (n==0)==(p in bad0)
print('Bring point counts independently matched conductor-50 trace for all 15 tested primes > 5')

gold=tuple(range(5))
g=tuple((x-2) for x in gold)
from collections import Counter
mod5=Counter()
for ds in product(range(5),repeat=4):
    y=(*ds,-sum(ds))
    ij=tuple(v%5 for v in moments(y))
    if len(set(v%5 for v in y))>=3: mod5[ij]+=1
counts=Counter()
goldcount=Counter()
for ds in product(range(25),repeat=4):
    y=(*ds,-sum(ds))
    ij=tuple(v%25 for v in moments(y))
    rs=tuple(v%5 for v in y)
    if len(set(rs))>=3: counts[ij]+=1
    if rs==tuple(v%5 for v in g): goldcount[ij]+=1
assert len(goldcount)==25 and set(goldcount.values())=={25}
assert all(counts[q,b]==25*mod5[q%5,b%5] for q,b in product(range(25),repeat=2))
assert all(mod5[0,b]==120*(b==0) for b in range(5))
assert all(mod5[q,b]==(20 if q in (1,4) else 30) for q,b in product(range(1,5),range(5)))
print('390625 mod-25 tuples; designated golden tuple lifts to all 25 moment pairs with 25 lifts each')

def scalar(y,c): q,b=moments(y); return c*q+b
witness=(
 ((1,0,0,0,-1),(2,2,0,0,-4),(1,1,0,0,-2)),
 ((1,1,1,1,-4),(1,0,0,0,-1),(2,1,1,0,-4)),
 ((2,2,2,1,-7),(2,1,1,0,-4),(1,0,0,0,-1)))
for c,r in product(range(3),repeat=2):
    y=witness[c][r]
    assert scalar(y,c)%3==r
    assert len({(c*a+2*a*a)%3 for a in y})>1
for c in range(12):
    for k in range(1,8):
        a={sum(F(c+x) for x in (2*t,1-2*t,-1,0,0))%(2**k) for t in range(2**k)}
        assert len(a)==2**k
print('scalar Hensel seeds checked at 2 and 3')

def selector(c,Q,B):
    # Return five offsets. All restrictions apply before evaluating F.
    if c<1 or Q<0: return None
    if Q==0: return (0,)*5 if B==0 else None
    W=isqrt(2*Q)
    lo=max(1-c,-W)
    hi=min(4*(c-1),W)
    for u in range(lo,hi+1):
        for s in range(max(2*lo,-u-2*hi),min(2*hi,-u-2*lo)+1):
            t=-u-s
            if s<=t: continue
            v=24*B-4*u**3-s**3-t**3
            if v%3: continue
            A=4*Q-2*u*u-s*s-t*t
            d=s-t
            a2,rem=divmod(v//3-t*A,d)
            b2=A-a2
            if rem or min(a2,b2)<0: continue
            a,b=isqrt(a2),isqrt(b2)
            if a*a!=a2 or b*b!=b2 or (s-a)%2 or (t-b)%2: continue
            y=(u,(s+a)//2,(s-a)//2,(t+b)//2,(t-b)//2)
            if min(y)>=1-c:
                assert moments(y)==(Q,B)
                return y

def compositions_sum(n,k,m=1):
    if k==1:
        if n>=m: yield (n,)
        return
    for i in range(m,n//k+1):
        for xs in compositions_sum(n-i,k-1,i): yield (i,*xs)
allpairs=0
for c in range(1,14):
    ys=[tuple(x-c for x in xs) for xs in compositions_sum(5*c,5)]
    oracle={moments(y) for y in ys}
    assert all(selector(c,Q,B) is not None for Q,B in oracle)
    allpairs+=len(oracle)
    for Q,B in list(oracle)[::13]:
        for eps in (1,5,19):
            assert (selector(c,Q,B+eps) is not None)==((Q,B+eps) in oracle)
print('selector checked on every realized pair for centers 1..13:',allpairs)
for c in range(2,14):
    a=c-1
    for xs in compositions_sum(5*c,5):
        Q,B=moments(tuple(x-c for x in xs))
        assert Q<=10*a*a and 10*B*B<=Q**3
        assert 15*a*B>=-5*a*a*Q+2*Q*Q
print('sharp real cusp and positive-window constraints checked through center 13')
