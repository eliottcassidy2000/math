# Replay guards are explicit even under Python -O.
def need(x):
    if not x: raise ArithmeticError("factor-seven certificate failed")

from itertools import combinations

def induced(T,order):
    return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(order)) for i in order)

def reachable(a,s,mask):
    r=1<<s
    while True:
        z=r
        for j in range(len(a)):
            if r>>j&1: z|=a[j]&mask
        if z==r: return r
        r=z

def scc(a,verts):
    verts=list(verts)
    mask=sum(1<<i for i in verts)
    grouped={}
    for i in verts:grouped.setdefault(reachable(a,i,mask),[]).append(i)
    return [grouped[z] for z in sorted(grouped,key=int.bit_count,reverse=True)]

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

def pad(a,r,m):
    for _ in range(m):
        a=tuple(v+((i!=r)<<len(a)) for i,v in enumerate(a))+(1<<r,)
    return a

def strong(a):
    F=(1<<len(a))-1
    return all(reachable(a,i,F)==F for i in range(len(a)))
from math import gcd,isqrt,prod
pieces=[
 [(11,800,(0,1,3,71,47,79,151,63)),
  (53,1640,(2,32,3,23,39,77,159,63)),
  (97,2584,(0,17,35,39,77,211,143,31))],
 [(77,1904,(2,4,9,19,103,79,143,63)),
  (53,1520,(2,16,3,7,45,79,159,63)),
  (97,4100,(256,1,11,275,327,287,431,319,6)),
  (1,2,(0,))]
]
def two_profile(T):
 need(all(not (r>>i&1) for i,r in enumerate(T)))
 need(all((T[i]>>j&1)^(T[j]>>i&1) for i in range(len(T)) for j in range(i)))
 H,_=Htable(T)
 F=(1<<len(T))-1
 return H[F],sum(H[s]*H[F^s] for s in range(F+1))
for P in pieces:
 for v,S,T in P:
  need(two_profile(T)==(v,S))
  print(T,'literal profile',(v,S))
profiles=[(prod(v for v,S,T in P),prod(S for v,S,T in P)) for P in pieces]
need(profiles==[(56551,3390208000),(395857,23731456000)])
need(profiles[1]==tuple(7*x for x in profiles[0]))
v,S=profiles[0]
need(gcd(v,S)==1)
def prime(x):return x>1 and all(x%d for d in range(2,isqrt(x)+1))
for n,p in [(1,3390264551),(4,50852384837),(6,213579654389)]:
 ans=[(S-v)*(1<<n)+3*v-S for v,S in profiles]
 need(ans==[p,7*p] and prime(p))
 print('transitive order',n,'prime count',p,'exact larger count',7*p)
need(two_profile((0,33,11,19,39,77,159,63))==(49,1624))
need(two_profile((26,60,25,32,8,21))==(23,272))
