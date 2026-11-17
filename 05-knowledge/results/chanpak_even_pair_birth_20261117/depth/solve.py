#!/usr/bin/env python3
"""Binary submodular cut decoder; rows and profiles are exact integers."""
from collections import deque
from math import gcd
from functools import cache
from itertools import combinations,product

@cache
def F(u):
    if not u:return 1
    p=u//2
    return F(p)+(0 if u%2 else p+p%2)

def count(S,u,w):
    h,P,L,R,ls,rs,k=S
    return k+h*(F(u)+F(w)-2)+R*u+L*w+P*u*w-rs*(u%2)-ls*(w%2)

def build(S,m,n):
    h,P,L,R,ls,rs,k=S
    lin=[R*2**i-rs*(i==0) for i in range(m)]+[-L*2**i+ls*(i==0) for i in range(n)]
    k+=L*((1<<n)-1)-ls*(n>0)
    arcs=[]
    for i,j in combinations(range(m),2):arcs.append((j,i,h*(2**(j-i-1)+(j==i+1))))
    for i,j in combinations(range(n),2):arcs.append((m+i,m+j,h*(2**(j-i-1)+(j==i+1))))
    for i in range(m):
        for j in range(n):arcs.append((i,m+j,P*2**(i+j)))
    return k,lin,arcs

def reduce(E,fixed):
    k,lin,arcs=E
    lin=lin.copy()
    keep=[]
    for i,j,c in arcs:
        a,b=fixed[i],fixed[j]
        if a==0 or b==1:continue
        if a==1 and b==0:k+=c
        elif a==1:k+=c;lin[j]-=c
        elif b==0:lin[i]+=c
        else:keep.append((i,j,c))
    for i,v in enumerate(fixed):
        if v!=-1:k+=v*lin[i];lin[i]=0
    return k,lin,keep

def minimum(E):
    k,lin,arcs=E
    n=len(lin);s,t=n,n+1;G=[[] for _ in range(n+2)]
    def arc(i,j,c):
        G[i].append([j,c,len(G[j])]);G[j].append([i,0,len(G[i])-1])
    for i,j,c in arcs:arc(i,j,c)
    for i,a in enumerate(lin):
        if a>0:arc(i,t,a)
        if a<0:k+=a;arc(s,i,-a)
    while True:
        dist=[-1]*(n+2);dist[s]=0;Q=deque([s])
        while Q:
            i=Q.popleft()
            for j,c,_ in G[i]:
                if c and dist[j]<0:dist[j]=dist[i]+1;Q.append(j)
        if dist[t]<0:return k
        p=[0]*(n+2)
        def send(i,f):
            if i==t:return f
            while p[i]<len(G[i]):
                e=G[i][p[i]];j,c,r=e
                if c and dist[j]==dist[i]+1:
                    g=send(j,min(f,c))
                    if g:e[1]-=g;G[j][r][1]+=g;return g
                p[i]+=1
            return 0
        bound=1+sum(e[1] for e in G[s])
        while (a:=send(s,bound)):k+=a

def dfs(S,N,m,n,full=False):
    """Search address lengths <=m,n; result (u,w) or None; no other cores."""
    E=build(S,m,n)
    fix=[-1]*(m+n); stats=[0,0,0,0]
    if full:
        if m:fix[m-1]=1
        if n:fix[m+n-1]=0
    def search():
        stats[0]+=1
        D=reduce(E,fix)
        k,Ls,A=D
        if N<k+sum(min(0,a) for a in Ls) or N>k+sum(max(0,a) for a in Ls)+sum(c for i,j,c in A):return
        if not A and not any(Ls):
            if N==k:
                z=[max(0,x) for x in fix]
                u=sum(z[i]*2**i for i in range(m));w=sum((1-z[m+j])*2**j for j in range(n))
                assert count(S,u,w)==N
                return u,w
            return
        linear=Ls.copy()
        g=0
        for i,j,c in A:
            linear[i]+=c;g=gcd(g,c)
        g=gcd(g,*linear)
        if g and (N-k)%g:
            stats[3]+=1;return
        stats[1]+=1
        if N<minimum(D):stats[2]+=1;return
        gain=[abs(a) for a in Ls]
        for i,j,c in A:gain[i]+=c;gain[j]+=c
        j=max((i for i in range(m+n) if fix[i]<0),key=lambda i:gain[i])
        for x in (0,1):
            fix[j]=x
            v=search()
            if v is not None:return v
        fix[j]=-1
    return search(),stats

def htable(rows):
    n=len(rows);A=[[0]*n for _ in range(1<<n)]
    for i in range(n):A[1<<i][i]=1
    for s in range(1,len(A)):
        for i,a in enumerate(A[s]):
            q=rows[i]&~s
            while a and q:
                b=q&-q;q-=b;A[s|b][b.bit_length()-1]+=a
    return A
def H(rows):return sum(htable(rows)[-1])
def converse(rows):return tuple(sum((r>>i&1)<<j for j,r in enumerate(rows)) for i in range(len(rows)))
def ordinal(A,B):return tuple(r|((1<<len(B))-1)<<len(A) for r in A)+tuple(r<<len(A) for r in B)
def TT(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def apex(D,q):
    n=len(D)
    return tuple(r|((1-(q>>i&1))<<n) for i,r in enumerate(D))+(q,)
def profile(rows,q):
    A,B=htable(rows),htable(converse(rows));M=len(A)-1
    h=[1]+[sum(z) for z in A[1:]];L=[0,0];R=[0,0];P=0
    for s in range(M+1):
        j=s.bit_count()%2
        P+=h[s]*h[M^s]
        L[j]+=sum(a for i,a in enumerate(B[s]) if q>>i&1)*h[M^s]
        R[j]+=sum(a for i,a in enumerate(A[s]) if not q>>i&1)*h[M^s]
    return h[-1],P,sum(L),sum(R),L[0]-L[1],R[0]-R[1],H(apex(rows,q))
def literal(D,q,u,w):
    a=tuple(map(int,bin(u)[2:])) if u else ()
    b=tuple(1-int(x) for x in reversed(bin(w)[2:])) if w else ()
    tag=sum(v<<i for i,v in enumerate(a))+(q<<len(a))+sum(v<<(len(a)+len(D)+i) for i,v in enumerate(b))
    return apex(ordinal(ordinal(TT(len(a)),D),TT(len(b))),tag)

def test():
    S=profile((4,1,2),1)
    for m,n in product(range(4),repeat=2):
        E=build(S,m,n)
        for fixed in product((-1,0,1),repeat=m+n):
            actual=[]
            for z in product((0,1),repeat=m+n):
                if any(v!=-1 and v!=z[i] for i,v in enumerate(fixed)):continue
                u=sum(z[i]*2**i for i in range(m));w=sum((1-z[m+j])*2**j for j in range(n))
                actual.append(count(S,u,w))
            assert actual and minimum(reduce(E,fixed))==min(actual)
        values={count(S,u,w) for u in range(1<<m) for w in range(1<<n)}
        for N in range(min(values),max(values)+1):
            ans,stats=dfs(S,N,m,n)
            assert (ans is not None)==(N in values)
    for u,w in product(range(8),repeat=2):assert count(S,u,w)==H(literal((4,1,2),1,u,w))
    print("PASS exact cut minima for all frozen cells of three+three bits; complete small search and literal transfer")
if __name__=='__main__':test()
