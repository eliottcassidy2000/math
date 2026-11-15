from itertools import product
from verify_seam import *
from collections import defaultdict
import itertools
def prod(p,q):
    r=defaultdict(int)
    for (a,b),x in p.items():
        for(c,d),y in q.items():r[a+c,b+d]+=x*y
    return dict(r)
def packets(a,u,v):
    E=lambda i,j: (a[i]>>j)&1
    options=[(),(u,),(v,),((u,v) if E(u,v) else (v,u))]
    old=[x for x in range(len(a)) if x not in (u,v)]
    def gap(x,y):
        g={}
        for w in options:
            z=([x] if x is not None else [])+list(w)+([y] if y is not None else [])
            if all(E(a,b) for a,b in zip(z,z[1:])):g[(int(u in w),int(v in w))]=1
        return g
    G={(x,y):gap(x,y) for x in [None]+old for y in [None]+old if x is None or y is None or x!=y}
    def f(w,b):
        if not w:return G[None,None]
        p={(0,0):1} if b else G[None,w[0]]
        for x,y in zip(w,w[1:]):p=prod(p,G[x,y])
        return prod(p,G[w[-1],None])
    h=defaultdict(int);s={i:defaultdict(int) for i in old}
    for w in itertools.permutations(old):
        for k,v in f(w,False).items():h[k]+=v
        for t in range(1,len(w)+1):
            for k,v in prod(f(w[:t],True),f(w[t:],False)).items():s[w[0]][k]+=v
    return {k:x for k,x in h.items() if all(k)}, {i:{k:x for k,x in p.items() if all(k)} for i,p in s.items()}
from math import comb
H=((43,36,8,0),(54,78,40,6),(16,34,28,6))
S=((74,91,42,6),(66,106,66,12),(16,34,28,6))
C=((191,218,92,12),(186,290,172,30),(48,102,84,18))
a=tournament(6,7002); assert strong(a)
h,s=packets(a,4,5)
for i in range(1,7):
    for j in range(1,7):
        old=i<=3 and j<=4
        z=(i,j)
        hh=H[i-1][j-1] if old else 0;ss=S[i-1][j-1] if old else 0
        assert h.get(z,0)==hh and s[2].get(z,0)==ss
        assert s[1].get(z,0)==ss+(z==(1,1))
        if old:assert hh+2*ss==C[i-1][j-1]
def insert(a,f,D):
    n=len(a);k=len(D)
    chunks=[[(r,1) for r in D] if i==f else [(0,0)] for i in range(n)]
    masks=[];base=0
    for c in chunks:masks.append(((1<<len(c))-1)<<base);base+=len(c)
    rows=[];base=0
    for i,c in enumerate(chunks):
        cross=sum(x for j,x in enumerate(masks) if a[i]>>j&1)
        rows.extend((u<<base)|cross for u,_ in c)
        base+=len(c)
    return rows
def pows(a,k):
    h,_=stats(a);f=len(h)-1
    P=[0]*(k+1);term=[0]*(f+1);term[0]=1
    for j in range(1,k+1):
        nex=[0]*(f+1)
        for u in range(1,f+1):
            t=u
            while t:
                nex[u]+=h[t]*term[u^t];t=(t-1)&u
        term=nex;P[j]=term[f]
    return P
def f(m,n):
    a,b,c,d,e=2**m,2**n,3**m,3**n,4**n
    return 51-20*b-24*d+6*e-70*a+80*a*b+16*a*d-24*a*e+24*c-42*c*b+12*c*d+18*c*e
def bounds(N):
    m=1
    while 2*3**m*4<=3*N:
        lo,hi=1,2
        while f(m,hi)<N-2:hi*=2
        while lo<hi:
            j=(lo+hi)//2
            if f(m,j)<N-2:lo=j+1
            else:hi=j
        if f(m,lo) in (N,N-2):yield (m,lo,N-f(m,lo))
        m+=1
# Exact insertion replay. Insert V (vertex 5) first, leaving index 4.
operands=[tt(1),tt(2),tt(3),tournament(3,2),tournament(4,17)]
for U,V in product(operands,repeat=2):
    u,v=pows(U,3),pows(V,4)
    A=insert(insert(a,5,V),4,U)
    base=sum(C[i-1][j-1]*u[i]*v[j] for i in range(1,4) for j in range(1,5))
    for r,off in [(2,0),(1,2*u[1]*v[1])]:
        T=reverse(join(A,tt(2)),r,len(A)+1)
        assert strong(T) and stats(T)[0][-1]==base+off
    if U==tt(len(U)) and V==tt(len(V)):assert base==f(len(U),len(V))
for i in range(1,9):
    for j in range(1,9):
        for d in [0,2]:
            assert (i,j,d) in list(bounds(f(i,j)+d))
print('two-module coefficient table, 50 literal checks, 128 exact memberships')
