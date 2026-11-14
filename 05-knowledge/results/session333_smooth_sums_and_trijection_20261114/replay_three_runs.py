from collections import defaultdict,Counter
from itertools import permutations,product
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import strong, hp_all, ends

def require(z):
    if not z: raise RuntimeError("three-run replay failed")
def order_induced(A,V):
    W=sorted(V,key=lambda u:-sum((A[u]>>v)&1 for v in V))
    require(all(A[u]>>v&1 for i,u in enumerate(W) for v in W[i+1:]))
    return W
def join_edges(A, p):
    return all(A[a]>>b&1 for a,b in zip(p,p[1:]))
def inflate(A,w,m):
    p=[]
    for i in range(len(A)):p.extend([i]*(m if i==w else 1))
    return tuple(sum(1<<j for j,b in enumerate(p) if (a==b and i<j) or (a!=b and A[a]>>b&1)) for i,a in enumerate(p))

def linear(A,C,B):
    V=order_induced(A,C); k=len(B);r=k+1;total=0
    for P in permutations(B):
        ds={(0,)*r:1}
        for v in V:
            es=defaultdict(int)
            for u,s in ds.items():
                for i in range(r):
                    if not u[i] and i and not (A[P[i-1]]>>v&1):continue
                    q=list(u);q[i]=1 if i==k or A[v]>>P[i]&1 else 2
                    es[tuple(q)]+=s
            ds=es
        for u,s in ds.items():
            if 2 in u:continue
            if any(not u[i] and not (A[P[i-1]]>>P[i]&1) for i in range(1,k)):continue
            total+=s
    return total

def group(A,C,B):
    V=order_induced(A,C);R=len(B)+1
    index={v:i for i,v in enumerate(V)}
    gs=defaultdict(int)
    for P in permutations(B):
      for f in product(range(R),repeat=len(V)):
        runs=[[v for v,c in zip(V,f) if c==r] for r in range(R)]
        w=[v for i in range(R) for v in (runs[i]+([P[i]] if i<len(P) else []))]
        if join_edges(A,w):
            g=tuple((u[0],u[-1]) if u else () for u in runs)
            gs[(P,g)]+=1
    hist=Counter()
    for (P,g),size in gs.items():
        anchors={e for ep in g for e in ep}
        a=b=0
        for v in V:
            if v in anchors:continue
            n=sum(len(k)==2 and index[k[0]]<index[v]<index[k[1]] for k in g)
            require(n in [1,2,3,4])
            a+=(n==2)+2*(n==4);b+=n==3
        require(2**a*3**b==size);hist[a,b]+=1
    return hist

def check(A,B,decompose=False):
    C=[i for i in range(len(A)) if i not in B]
    H=hp_all(A)[-1]
    require(linear(A,C,B)==H)
    if decompose:
        e=group(A,C,B)
        require(sum((2**a*3**b)*n for (a,b),n in e.items())==H)
        return H,sorted(e.items())
    return H

examples=[
 ("H69",(2,8,3,21,103,79,143,63),(3,6)),
 ("H45",tuple(sum((1<<j) for j in range(6) if (i//2==j//2 and i<j) or (i//2+1)%3==j//2) for i in range(6)),(4,5)),
 ("frame135",(258,276,9,291,301,327,415,319,4),(0,3,5)),
 ("frame405",(262,28,424,337,293,331,279,379,2),(1,2,4)),
]
for name,A,ap in examples:
    print(name,check(A,ap,True))
for ext in range(1,5):
    _,A,ap=examples[0]
    print("inflated H69 frame",ext,check(inflate(A,7,ext),ap,True))
# A literal trijection fibre with a free vertex covered by all three runs.
n=9; x,y=7,8; A=[0]*n
for a in range(n):
  for b in range(a+1,n):
    rev=(b==x and a in (0,1)) or (b==y and a in (1,2))
    if rev:A[b]|=1<<a
    else:A[a]|=1<<b
A=tuple(A)
require(strong(A))
for s in range(3):
    runs=[[0,4],[1,5],[2,6]];runs[s].insert(1,3)
    R,S,T=runs
    u=R+[x]+S+[y]+T
    require(join_edges(A,u))
    e=set(zip(u,u[1:]));e.remove((x,S[0]));e.remove((y,T[0]));e.add((x,R[0]));e.add((y,S[0]))
    require(len(e)==8 and all(A[a]>>b&1 for a,b in e))
print("all three ternary fillings mapped to 3 components with exactly n−1 arcs")
# Arithmetic antichain expansion by literal Euclidean remainder test.
def smooth(N):
    r=[];a=0
    while N:
        while not N%2:a+=1;N//=2
        b=0;p=1
        while p*3<=N:p*=3;b+=1
        r.append((a,b));N-=p
    return r
for n in range(1,10000):
    z=smooth(n)
    require(sum(2**a*3**b for a,b in z)==n)
    require(all(a<c and b>d for (a,b),(c,d) in zip(z,z[1:])))
print({n:smooth(n) for n in (34,45,69)})

rails=[
 (135,(258,276,9,291,301,327,415,319,4),2,(115,41,20,176)),
 (405,(262,28,424,337,293,331,279,379,2),1,(379,83,56,488)),
 (1215,(760,1013,769,518,44,652,956,540,697,16),4,(1095,309,188,1524)),
]
for D,A,y,P in rails:
    n=len(A)-1;U=tuple(o&((1<<n)-1) for o in A[:-1]);F=hp_all(U);M=len(F)-1
    E=ends(U);R=ends(tuple(M^o^(1<<i) for i,o in enumerate(U)))
    X=(F[M],F[M^(1<<y)],E[M][y],sum(R[T][y]*F[M^T] for T in range(M+1)))
    require(X==P and A[-1]==1<<y and strong(A))
    H,h,e,L=X
    require((H+2*h-e-L,L-h)==(1,D))
    for s in (1,2,3,4):
        W=inflate(A,n,s)
        require(strong(W) and hp_all(W)[-1]==1+D*(1<<s))
require([1+D*(1<<s) for D,s in zip((15,45,135,405,1215),(6,5,4,3,2))]==[961,1441,2161,3241,4861])
print("135,405,1215 sole-exit rails and exact a=9 predecessor chain checked")

def decoder(r,bits,space=False,first=True):
    # chain y,c0,(z),c1,...,c_r,t, then apex x
    n=1+1+space+r+1;t=n-1;x=n
    firstbit=2+space
    losses={t}
    if not first:losses.add(1)
    losses.update(firstbit+i for i in range(r) if bits>>i&1)
    A=tuple(sum(1<<j for j in range(i+1,n))|((1<<x) if i in losses else 0) for i in range(n))+(sum(1<<j for j in range(n) if j not in losses),)
    return A,1,t
def start_end(A,s,t):
    n=len(A);N=1<<n;P=[[0]*n for _ in range(N)];P[1<<s][s]=1
    for M in range(N):
        for i,h in enumerate(P[M]):
            if not h or i==t:continue
            E=A[i]&~M
            while E:
                b=E&-E;E-=b;P[M|b][b.bit_length()-1]+=h
    return P[-1][t]
gaps=set()
for r in range(6):
  for N in range(1<<r):
    A,s,t=decoder(r,N)
    require(strong(A) and start_end(A,s,t)==N)
    U,*_=decoder(r,N,True,True);V,*_=decoder(r,N,True,False)
    require(strong(U) and strong(V) and sum((a^b).bit_count() for a,b in zip(U,V))==2)
    h=hp_all(V)[-1];delta=hp_all(U)[-1]-h
    require(delta==2*N+(2<<r)-2)
    if r:require(2*h>=5*delta+16)
    gaps.add(delta)
require(gaps==set(range(0,max(gaps)+1,2)))
print("all 63 marked ending inputs and all corresponding single-arc absolute gaps replayed")
