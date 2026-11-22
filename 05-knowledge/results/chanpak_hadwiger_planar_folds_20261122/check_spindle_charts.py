from fractions import Fraction as F
from itertools import combinations, permutations
from math import prod

BASE=(80,73,35,5,14,27,60)
MARK=((0,1),(2,6),(3,4))
U=(F(1),F(0));V=(F(0),F(1));W=(F(3,5),F(4,5));WT=(F(5,13),F(12,13))
M=tuple(4*((b>>0)&1)+7*((b>>1)&1)+8*((b>>2)&1) for b in range(8))

def rows(b):
    a=list(BASE)
    for j,(u,v) in enumerate(MARK):
        if b>>j&1:a[u]^=1<<v;a[v]^=1<<u
    return tuple(a)

def reads(a):
    n=len(a); full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,full):
        for i,v in enumerate(dp[S]):
            if not v:continue
            r=a[i]&~S
            while r:
                j=(r&-r).bit_length()-1;r&=r-1
                dp[S|1<<j][j]+=v
    p=[sum(r) for r in dp];p[0]=1
    return p[-1],sum(p[S]*p[full^S] for S in range(full+1))

def strong(a):
    n=len(a);full=(1<<n)-1
    for i in range(n):
        S=1<<i
        while True:
            T=S
            for j in range(n):
                if S>>j&1:T|=a[j]
            if S==T:break
            S=T
        if S!=full:return False
    return True

def point(b,unit):
    return tuple(sum(t[i] for j,t in enumerate(unit) if b>>j&1) for i in (0,1))

def sq(P,Q):return sum((x-y)**2 for x,y in zip(P,Q))

def project(P):x,y=P;return 4*x+7*y,F(y,4)
def unproject(P):m,q=P;return F(m-28*q,4),4*q

def poly(m,values):
    return sum(values[i]*prod(F(m-M[j],M[i]-M[j]) for j in range(8) if j!=i) for i in range(8))

P=[point(b,(U,V,W)) for b in range(8)]
PT=[point(b,(U,V,WT)) for b in range(8)]
vals=[F(b>>2) for b in range(8)]
codes=[F(project(z)[1])-(b.bit_count()&1) for b,z in enumerate(P)]
assert len(set(P))==len(set(PT))==len(set(M))==8

for b,c in combinations(range(8),2):
    edge=(b^c).bit_count()==1
    assert (sq(P[b],P[c])==1)==edge==(sq(PT[b],PT[c])==1)
    changed=sum(((a^d).bit_count()) for a,d in zip(rows(b),rows(c)))//2
    assert changed==(b^c).bit_count()

for b in range(8):
    p=P[b]; pt=PT[b];m,q=project(p);mm,qq=project(pt)
    parity=b.bit_count()&1
    assert m==mm==M[b] and q+(F(2,65))*poly(m,vals)==qq
    assert unproject((m,qq-F(2,65)*poly(m,vals)))==p
    assert unproject((m,parity+poly(m,codes)))==p
    assert q-poly(m,codes)==parity
    assert strong(rows(b))
    assert reads(rows(b))==(141+2*parity,1656+4*parity)
assert sq(P[1],P[4])==F(4,5) and sq(PT[1],PT[4])==F(16,13)
assert sum(vals[i]/prod(M[i]-M[j] for j in range(8) if j!=i) for i in range(8))==F(17,300960)

def iso(a,b):
    return next(pi for pi in permutations(range(7)) if all((a[u]>>v&1)==(b[pi[u]]>>pi[v]&1) for u,v in combinations(range(7),2)))
for b in range(8):
    iso(rows(b),rows(b.bit_count()&1))
print("Eight strong physical tournaments, two isomorphism classes, 12 exact unit-distance seam edges.")
print("Both reciprocal shears, bit/height recovery and nonedge deformation verified over Q.")

# Exact isomorphism-quotient Moser spindle, distinct from the unit-gap cube.
MOSER=((0,1),(0,2),(0,4),(0,5),(1,2),(1,3),(2,3),(3,6),(4,5),(4,6),(5,6))
SM=(20,16,22,26,81,156,89)
SM_TARGET=(16,22,84,276,8208,2064,4118,90,337,89,4252)
ED=list(combinations(range(6),2))
P6=list(permutations(range(6)))

def rows6(m):
    a=[0]*6
    for j,(u,v) in enumerate(ED):
        if m>>j&1:a[u]|=1<<v
        else:a[v]|=1<<u
    return tuple(a)
def orbit(m):
    a=rows6(m)
    return {sum(1<<j for j,(u,v) in enumerate(ED) if a[pi[u]]>>pi[v]&1) for pi in P6}

OR=[orbit(m) for m in SM]
assert all(map(strong,map(rows6,SM)))
assert all(not OR[i]&OR[j] for i,j in combinations(range(7),2))
assert [reads(rows6(m))[0] for m in SM]==[19,17,15,19,23,33,25]
for i,j in combinations(range(7),2):
    assert any(SM[i]^(1<<b) in OR[j] for b in range(15))==((i,j) in MOSER)
for (i,j),v in zip(MOSER,SM_TARGET):
    assert (SM[i]^v).bit_count()==1 and v in OR[j]
assert sum(abs(reads(rows6(SM[i]))[0]-reads(rows6(SM[j]))[0])==2 for i,j in MOSER)==4
a,b=rows6(8208),rows6(22)
pi=(1,2,3,4,5,0)
assert all((a[i]>>j&1)==(b[pi[i]]>>pi[j]&1) for i,j in ED)
assert sum(pi[i]>pi[j] for i,j in ED)%2==1

# Q(sqrt(3),sqrt(11)), basis 1,sqrt(3),sqrt(11),sqrt(33).
def k(v=0,index=0):return tuple(F(v) if i==index else F(0) for i in range(4))
def plus(a,b):return tuple(x+y for x,y in zip(a,b))
def minus(a,b):return tuple(x-y for x,y in zip(a,b))
def times(a,b):
    c=[F(0)]*4
    for i in range(4):
        for j in range(4):
            g=i&j
            c[i^j]+=a[i]*b[j]*(3 if g&1 else 1)*(11 if g&2 else 1)
    return tuple(c)
def qdist(a,b):return plus(*(times(minus(x,y),minus(x,y)) for x,y in zip(a,b)))
c,s=k(F(5,6)),k(F(1,6),2)
def rotate(z):
    x,y=z
    return minus(times(c,x),times(s,y)),plus(times(s,x),times(c,y))
A=(k(),k()); P=(k(F(1,2),1),k(F(1,2))); Q=(P[0],k(F(-1,2))); T=(k(1,1),k())
SP=[A,P,Q,T,rotate(P),rotate(Q),rotate(T)]
assert len(set(SP))==7
assert len({plus(*z) for z in SP})==7
for i,j in combinations(range(7),2):assert (qdist(SP[i],SP[j])==k(1))==((i,j) in MOSER)
print("Strong quotient Moser spindle: 7 different classes, exact 11 Euclidean unit edges, 4 unit-gap lifts.")
