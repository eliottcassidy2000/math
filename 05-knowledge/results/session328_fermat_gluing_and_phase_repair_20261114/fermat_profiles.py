from profile_tools import *

def ends(A):
    n=len(A);N=1<<n;E=[[0]*n for _ in range(N)]
    for i in range(n):E[1<<i][i]=1
    for s in range(1,N):
        for i,e in enumerate(E[s]):
            if not e:continue
            r=A[i]&~s
            while r:
                b=r&-r;r-=b;E[s|b][b.bit_length()-1]+=e
    return E

def submasks(m):
    s=m
    while True:
        yield s
        if not s:return
        s=(s-1)&m

def rooted(A):
    m=(1<<len(A))-1
    E=ends(A);P=ends(tuple(m^r^(1<<i) for i,r in enumerate(A)))
    h=[1]+[sum(e) for e in E[1:]]
    q=[sum(h[s]*h[j^s] for s in submasks(j)) for j in range(m+1)]
    for r,o in enumerate(A):
        t=m^(1<<r);inside=t^o
        L=h[t];Ls=q[t];Es=h[t]
        for Q in submasks(t):
            if not Q:continue
            ps=sum(P[Q][i] for i in range(len(A)) if o>>i&1)
            ee=sum(E[Q][i] for i in range(len(A)) if inside>>i&1)
            L+=ps*h[t^Q]; Ls+=ps*q[t^Q]; Es+=ee*h[t^Q]
        yield (h[m],L,q[m],Ls,Es,r)

def glue(A,r,j):
    n=len(A); high=((1<<(j+1))-1)<<n
    return tuple(row|(high if i!=r else high^(1<<(n+j))) for i,row in enumerate(A))+tuple(
        (high&~((1<<(n+k+1))-1)) | ((1<<r) if k==j else 0) for k in range(j+1))
