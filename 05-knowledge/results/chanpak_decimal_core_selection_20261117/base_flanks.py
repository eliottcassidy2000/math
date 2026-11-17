from functools import lru_cache
from itertools import combinations,product
from random import Random

D=(56,33,35,6,14,24)
def require(x):
    if not x:raise ArithmeticError('control failed')
def htable(a):
    n=len(a);M=1<<n
    arr=[[0]*n for _ in range(M)]
    for i in range(n):arr[1<<i][i]=1
    for s in range(1,M):
        for i,z in enumerate(arr[s]):
            ll=a[i]&~s
            while z and ll:
                x=ll&-ll;ll-=x;arr[s|x][x.bit_length()-1]+=z
    return arr
def H(a):return sum(htable(a)[-1]) if a else 1
def TT(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def swap(a):return tuple(sum(((x>>i)&1)<<j for j,x in enumerate(a)) for i in range(len(a)))
def apex(d,tag):
    n=len(d);return tuple(r|((1-(tag>>i&1))<<n) for i,r in enumerate(d))+(tag,)
def join(A,B):return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)
def make(blocks):
    t=();mask=0
    for d,v in blocks:
        mask|=v<<len(t);t=join(t,d)
    return apex(t,mask)
def local(d,v):
    a,b=htable(d),htable(swap(d));M=len(a)-1;T=[1]+[sum(row) for row in a[1:]]
    L=[0,0];R=L.copy();P=[0,0]
    for s in range(M+1):
        j=s.bit_count()%2;h=T[M^s]
        P[j]+=T[s]*h
        L[j]+=sum(z for i,z in enumerate(b[s]) if v>>i&1)*h
        R[j]+=sum(z for i,z in enumerate(a[s]) if not v>>i&1)*h
    require(P[0]==P[1])
    return T[-1],sum(P),L,R,H(apex(d,v))-T[-1]
def law(blocks):
    from math import prod
    p=[local(*b) for b in blocks]; h=[z[0] for z in p]
    ans=prod(h)
    for i,z in enumerate(p):ans+=prod(h[:i]+h[i+1:])*z[4]
    for i,j in combinations(range(len(p)),2):
        L,R=p[i][2],p[j][3]
        if j==i+1:e=2*(L[0]*R[0]+L[1]*R[1])
        else:e=sum(L)*sum(R)*prod(z[1] for z in p[i+1:j])
        ans+=prod(h[:i]+h[j+1:])*e
    return ans
@lru_cache(None)
def F(t):
    if not t:return 1
    z=t//2
    return F(z) if t%2 else F(z)+z+(z%2)
def flanks(u,w,tag):
    p=[]
    if u:
        b=[int(x) for x in bin(u)[2:]]
        p.append((TT(len(b)),sum(x<<i for i,x in enumerate(b))))
    p.append((D,tag))
    if w:
        b=[1-int(x) for x in reversed(bin(w)[2:])]
        p.append((TT(len(b)),sum(x<<i for i,x in enumerate(b))))
    return p
for v,delta in ((6,92),(9,94)):
    h,p,L,R,g=local(D,v)
    require((h,p,sum(L),sum(R),g)==(45,432,129,258,delta))
    for u,w in product(range(0,13,2),repeat=2):
        q=flanks(u,w,v)
        r=45+delta+45*(F(u)+F(w)-2)+258*u+129*w+432*u*w
        require(law(q)==r==H(make(q)))
    for m,n in product(range(9),repeat=2):
        u,w=2*(4**m-1)//3,2*(4**n-1)//3
        require(45+delta+45*(F(u)+F(w)-2)+258*u+129*w+432*u*w==192*4**(m+n)-86*4**n+30*(m+n)+31+(delta-92))
rng=Random(503)
for _ in range(90):
    q=[]
    for z in range(rng.randrange(1,5)):
        n=rng.randrange(1,5);a=[0]*n
        for i,j in combinations(range(n),2):
            if rng.randrange(2):a[i]|=1<<j
            else:a[j]|=1<<i
        q.append((tuple(a),rng.randrange(1<<n)))
    require(law(q)==H(make(q)))
print('PASS: 98 literal variable-flank counts, 162 dyadic formula controls, 90 random interval laws')
