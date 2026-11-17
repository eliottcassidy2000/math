from math import isqrt
from functools import lru_cache
def TT(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def join(A,B):
    x=((1<<len(B))-1)<<len(A)
    return tuple(y|x for y in A)+tuple(y<<len(A) for y in B)
def bits(u):
    return bin(u)[2:] if u else ''
def make(A,r,u):
    w=bits(u);k=len(w);B=TT(k);S=join(B,A);i=len(S)
    tag=sum(int(b)<<j for j,b in enumerate(w))|(((1<<len(A))-1)^(1<<r))<<k
    return tuple(z|((1-(tag>>j&1))<<i) for j,z in enumerate(S))+(tag,),i
def contrast(a):
    return tuple(sum(((a[j]>>i)&1)<<j for j in range(len(a))) for i in range(len(a)))
def paths(a):
    n=len(a);U=1<<n;d=[[0]*n for _ in range(U)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,U):
        for i,z in enumerate(d[s]):
            if not z:continue
            b=a[i]&~s
            while b:
                x=b&-b;b-=x;d[s|x][x.bit_length()-1]+=z
    return [1]+[sum(b) for b in d[1:]],d
def profile(a,r):
    s,d=paths(a);_,e=paths(contrast(a));U=len(s)-1
    return (s[U],sum(d[j][r]*s[U^j] for j in range(U+1)),e[U][r],d[U][r])
def F(u):
    v=1
    while u:
        x,b=divmod(u,2)
        if not b:v+=x+x%2
        u=x
    return v
a=(20,1,2,7,14)
def ck(v):
    if not v:raise ArithmeticError('check failed')
def covers(a,q):
    H,_=paths(a);N=len(H)
    p=[1]+[0]*(N-1)
    for _ in range(q):
        R=[0]*N
        for s in range(N):
            t=s
            while True:
                R[s]+=H[t]*p[s^t]
                if not t:break
                t=(t-1)&s
        p=R
    return p
H,U=paths(a);_,V=paths(contrast(a));N=len(H)-1
for r,S,L in [(0,1,11),(3,3,23)]:
    ck(profile(a,r)==(9,11,S,1))
    ck(sum(V[j][r]*H[N^j] for j in range(N+1))==L)
ck(covers(a,2)[N]==84)
for q in range(1,5):
    def parity(r):
        P=covers(a,q);ans=[0,0]
        for j in range(1,N+1):ans[j.bit_count()%2]+=U[j][r]*P[N^j]
        return ans
    x,y=parity(0),parity(3)
    ck((y[0]-x[0],y[1]-x[1])==(-q**3,q**2))
for u in range(32):
    w=bits(u);k=len(w)
    x,_=make(a,0,u);y,_=make(a,3,u)
    ck(sum((v^z).bit_count() for v,z in zip(x,y))==4)
    X,Y=paths(x)[0][-1],paths(y)[0][-1]
    ck(X==9*F(u)+11*(u+1)-(1-u%2) and X-Y==2*(1-u%2))
    for q in (2,3) if k<=4 else (2,):
        A=q**k+sum(q**i*(q+1)**(k-i-1) for i,b in enumerate(w) if b=='1')
        B=-q**k+sum(q**i*(q-1)**(k-i-1) for i,b in enumerate(w) if b=='1')
        d=q**3*((1-q)*A+(1+q)*B)
        ck(covers(y,q)[-1]-covers(x,q)[-1]==d)
        if u%2:ck(d<0)
for u,z in [(12,239),(424,8191)]:
    ck(9*F(u)+11*(u+1)-3==z)
    if u==12:ck(all(all(x%d for d in range(2,isqrt(x)+1)) for x in (z,z+2)))
ck(all(9*F(u)+11*(u+1)-S!=65537 for u in range(0,65537//11+1,2) for S in (1,3)))
print('PASS rooted primitive pair, literal arc switches, higher-slot detection, and exact examples')
