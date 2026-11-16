from itertools import combinations,product
from math import isqrt,gcd
from functools import lru_cache

def read_profile(T):
    n=len(T); F=(1<<n)-1
    d=[[0]*n for _ in range(F+1)]
    for i in range(n):d[1<<i][i]=1
    for m in range(1,F+1):
        for i,c in enumerate(d[m]):
            if not c: continue
            s=T[i]&~m
            while s:
                b=s&-s;s-=b
                d[m|b][b.bit_length()-1]+=c
    h=[1]+[sum(a) for a in d[1:]]
    return h[F],sum(h[m]*h[F^m] for m in range(F+1))

def Htable(a):
    n=len(a);F=(1<<n)-1
    path=[[0]*n for _ in range(1<<n)]
    H=[1]+[0]*F
    for s in range(1,1<<n):
        for j in range(n):
            if s>>j&1:
                u=s^(1<<j)
                path[s][j]=sum(path[u][i] for i in range(n) if a[j]>>i&1) if u else 1
        H[s]=sum(path[s])
    return H,path

def profile(a,r):
    h,path=Htable(a);F=(1<<len(a))-1
    E=sum(path[s][r]*h[F^s] for s in range(1<<len(a)) if s.bit_count()%2==0)
    O=sum(path[s][r]*h[F^s] for s in range(1<<len(a)) if s.bit_count()%2==1)
    return h[F],E,O,h[F^(1<<r)]

from itertools import permutations

def join(a,b):
    n=len(a); mask=((1<<len(b))-1)<<n
    return tuple(x|mask for x in a)+tuple(x<<n for x in b)

def tr(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def all_blocks(cap):
    # Keep one representative per degree-ordered canonical mask.
    def canon(t):
        n=len(t)
        cells=[tuple(i for i in range(n) if t[i].bit_count()==d) for d in range(n)]
        return min(tuple(sum((t[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
            for z in product(*(permutations(x) for x in cells if x)) for p in [sum(z,())])
    layer={(0,)}
    while layer:
        nxt=set()
        for t in sorted(layer):
            u,v=read_profile(t)
            if v>cap:continue
            yield t,u,v
            if 2*v>cap:continue
            n=len(t)
            for m in range(1<<n):
                child=tuple(row|((1-(m>>i&1))<<n) for i,row in enumerate(t))+(m,)
                if read_profile(child)[1]<=cap:nxt.add(canon(child))
        layer=nxt

def prime(p):return p>1 and all(p%d for d in range(2,isqrt(p)+1))
W41=(30,60,120,112,0,81,19);W53=(62,92,120,16,32,74,25);W5=(2,4,9,3)
def pad(a,n):return join(a,tr(n))
BANK=[
 ((6,12,24,17,3),15,120),(tr(3),1,8),
 (pad(join(join(W41,W5),W5),10),1025,872415232),
 ((14,4,8,2),3,24),((30,92,56,48,96,67,13),121,1556),
 (W41,41,832),(pad((14,4,8,2),3),3,192),((126,28,88,48,96,70,10),43,848),
 (W53,53,844),((126,28,24,112,96,70,6),45,864),(tr(9),1,512),
 (pad(W41,6),41,53248),(tr(6),1,64),((22,12,24,17,2),13,104),
 (pad(join(join(W53,W5),W5),10),1325,884998144)]
for i,(b,u,v) in enumerate(BANK):
    if len(b)<=13:assert read_profile(b)==(u,v),(i,b)
    else:
        x=read_profile(W41 if i==2 else W53); y=read_profile(W5)
        assert x[0]*y[0]**2==u and x[1]*y[1]**2*(1<<10)==v
    assert v%7==u%7!=0 and (52*v+11*u)//7%15==i
    for a in range(1,4):
        m=7**a
        observed={((52*v*pow(4096,t,7*m)+11*u)//7)%m for t in range(m)}
        assert len(observed)==m
A=(2,16,3,23,165,79,31,111); root=1
def graft(a,r,b):
    n=len(a); mask=((1<<len(b))-1)<<n
    return tuple(v|(mask if j!=r else 0) for j,v in enumerate(a))+tuple((x<<n)|(1<<r) for x in b)
assert profile(A,root)==(63,26,37,11)
assert profile(graft(A,root,W53),root)==(44471,21944,22527,583)
assert read_profile(W41)==(41,832) and read_profile(W53)==(53,844)
assert prime(6353)
print('literal witnesses and all 15 profiles checked',flush=True)
min_any=10**9; min_prime=10**9; layer={}
for b,u,v in all_blocks(855):
    layer[len(b)]=layer.get(len(b),0)+1
    h=52*v+11*u
    if h%7 or (p:=h//7)%3!=2:continue
    min_any=min(min_any,p)
    if prime(p):min_prime=min(min_prime,p)
assert (min_any,min_prime)==(6245,6353),(min_any,min_prime)
print('complete bounded class counts:',layer,'least quotient/prime:',min_any,min_prime,flush=True)

def f(x):
    v=1
    while x:
        b,x=x&1,x>>1
        if not b:v+=x+(x&1)
    return v
lane=sorted((h//7,z) for z in range(1,34)
    if (h:=52*(z+1)+11*f(z))<=1799 and h%21==14)
assert lane==[(257,30)] and prime(257)
T=A
for bit in '10000':
    n=len(T)
    T=tuple(x|((int(bit) or j!=root)<<n) for j,x in enumerate(T))+(0 if bit=='1' else 1<<root,)
assert profile(T,root)==(1799,920,931,11)
for k in range(1,13):
    z=(1<<(k+1))-2
    assert 52*(z+1)+11*f(z)==115*(1<<k)-41
print('binary first hit: 257, word 10000, literal H=1799')
