from math import prod
from collections import Counter
import argparse
M=(1<<60)-1
TARGET=M*((1<<61)+1)
B=(
(1321,775,5,(256,289,67,23,295,141,59,95,236)),
(331,671,3,(6,72,42,49,71,83,141,63)),
(451,1321,4,(768,257,3,7,815,399,831,863,524,46)),
(305,2317,2,(256,257,3,23,39,399,319,95,156)),
(3,3,2,(2,4,1)),
(3,3,2,(2,4,1)),
(5,1,5,(2,8,3,5)),
(1057,151,8,(256,257,515,263,783,31,831,127,164,427)),
(13,13,3,(8,5,17,6,11)),
(31,41,3,(4,33,18,7,11,29)),
)
U=tuple(x for *_,x in B)+((0,),)*24

def strong(A):
    n=len(A);f=(1<<n)-1
    for i in range(n):
        S=1<<i
        while True:
            z=S
            for j in range(n):
                if S>>j&1:z|=A[j]
            if z==S:break
            S=z
        if S!=f:return False
    return True

def ends(A):
    n=len(A);N=1<<n
    D=[[0]*n for _ in range(N)]
    for i in range(n):D[1<<i][i]=1
    for m in range(1,N):
        for i,h in enumerate(D[m]):
            if not h:continue
            E=A[i]&~m
            while E:
                e=E&-E;E-=e
                D[m|e][e.bit_length()-1]+=h
    return D

def reverse(A):
    n=len(A)
    return tuple(((1<<n)-1)^(1<<i)^A[i] for i in range(n))

def read(A,need=3,root=True):
    n=len(A);N=1<<n;F=N-1
    E=ends(A);H=[1]+list(map(sum,E[1:]));z=[[1]+[0]*F,H]
    for q in range(2,need+1):
        R=[0]*N
        for m in range(N):
            a=m;t=0
            while True:
                t+=H[a]*z[q-1][m^a]
                if not a:break
                a=(a-1)&m
            R[m]=t
        z.append(R)
    Z=tuple(a[F] for a in z)
    if not root:return Z
    S=ends(reverse(A))
    def rooted(D):
        return tuple(tuple(sum(D[m][i]*z[q-1][F^m] for m in range(N) if m>>i&1) for i in range(n)) for q in range(1,need+1))
    return (Z,)+tuple(rooted(D) for D in (S,E))

def flip(A,i,j):
    R=list(A);R[i]^=1<<j;R[j]^=1<<i
    return tuple(R)

def ordsum(L):
    p=0;total=sum(map(len,L));f=(1<<total)-1;ans=[]
    for A in L:
        n=len(A)
        ans.extend((a<<p)|(f-((1<<(p+n))-1)) for a in A)
        p+=n
    return tuple(ans)

def cyclic(A):
    n=len(A)
    return tuple(a|(1<<n) for a in A)+(1<<(n+1),(1<<n)-1)

def weighted_cross(P,i,j,u,v,q):
    fl=prod(Z[q] for Z,*_ in P[:i]+P[j+1:])
    zp=prod(Z[q+1] for Z,*_ in P[i+1:j])
    zm=prod(Z[q-1] for Z,*_ in P[i+1:j])
    Z,S,E=P[i]; Z1,S1,E1=P[j]
    return q*fl*(S[q][u]*E1[q][v]*zp-E[q-1][u]*S1[q-1][v]*zm)

def outercontact(P,i,u):
    Z,S,E=P[i]
    def ss(q):
        return S[q-1][u]*prod(Z[q-1] for Z,*_ in P[:i])*prod(Z[q] for Z,*_ in P[i+1:])
    def ee(q):
        return E[q-1][u]*prod(Z[q] for Z,*_ in P[:i])*prod(Z[q-1] for Z,*_ in P[i+1:])
    return ss(2)-ss(1)-ee(2)-ee(1),ee(2)-ee(1)-ss(2)-ss(1)

def scan(L):
    P=tuple(read(A) for A in L)
    H=prod(Z[1] for Z,*_ in P);Q=prod(Z[2] for Z,*_ in P)
    J=len(L);n=sum(map(len,L))
    def all_changes():
        yield ('outer',),-Q
        for i,A in enumerate(L):
            Z,*_=P[i]
            for u in range(len(A)):
                x,y=outercontact(P,i,u)
                yield ('A',i,u),x
                yield ('B',i,u),y
                for v in range(u+1,len(A)):
                    R=read(flip(A,u,v),2,False)
                    yield ('internal',i,u,v),H//Z[1]*(R[1]-Z[1])+Q//Z[2]*(R[2]-Z[2])
                for j,B in enumerate(L[i+1:],i+1):
                    for v in range(len(B)):
                        yield ('cross',i,j,u,v),weighted_cross(P,i,j,u,v,1)+weighted_cross(P,i,j,u,v,2)
    count=Counter(); nearest={}; se=[]
    for tag,y in all_changes():
        name=tag[0];count[name]+=1
        if (d:=nearest.get(name)) is None or (abs(y-2),y,tag)<d:nearest[name]=(abs(y-2),y,tag)
        if y==2:se.append(tag)
    assert sum(count.values())==(n+2)*(n+1)//2
    return H,Q,count,nearest,se

def controls():
    a=(2,4,1);b=(2,8,3,5);x=(0,)
    for L in ((a,x),(x,b),(a,x,x),(x,a,x),(a,x,b),(a,b),(x,x,x,x)):
        P=tuple(read(A) for A in L);D=ordsum(L);H=read(D,2,False)
        n=len(D);C=cyclic(D); base=read(C,1,False)[1]
        assert base==H[1]+H[2]
        R=[sum(map(len,L[:i])) for i in range(len(L))]
        for i,A in enumerate(L):
            for u in range(len(A)):
                p=R[i]+u
                x,y=outercontact(P,i,u)
                assert x==read(flip(C,p,n),1,False)[1]-base
                assert y==read(flip(C,p,n+1),1,False)[1]-base
                for j,B in enumerate(L[i+1:],i+1):
                    for v in range(len(B)):
                        q=R[j]+v
                        R0=read(flip(D,p,q),2,False)
                        for t in (1,2):
                            assert R0[t]-H[t]==weighted_cross(P,i,j,u,v,t)

if __name__=='__main__':
    controls()
    for h,odd,v,A in B:
        assert strong(A)
        _,H,Q=read(A,2,False)
        assert H==h and Q==odd<<v
    h,q,counts,nearest,se=scan(U)
    assert len(cyclic(ordsum(U)))==93 and h==M and q==(M<<61)
    print('strong on 93 vertices: H=',h+q)
    print('arc counts',dict(counts))
    print('nearest to +2',nearest)
    print('seams:',se)
