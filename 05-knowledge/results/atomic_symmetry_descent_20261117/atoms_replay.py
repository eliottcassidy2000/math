def need(x):
    if not x: raise ArithmeticError("atomic factorization certificate failed")

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
from math import isqrt
for H,T in [(49,(2,8,3,21,71,31,175,63)),(63,(16,5,17,7,74,31,175,63)),(343,(4,37,112,23,99,201,139,31))]:
 need(strong(T))
 need(Htable(T)[0][-1]==H)
 print(H,T,'strong and literal Hamilton count verified')
def lengths(a,b,c):
 return sorted({b+c+x+y-z for x in range(a+1) for y in range(a+1) for z in range(min(a,b//2)+1) if 0<=a-2*x-3*y-z<=c})
need(lengths(2,0,2)==[2,3])
need(lengths(2,2,1)==[2,4])
need(lengths(2,4,0)==[2,5])
for n in range(2,23):
 need(lengths(n,2*n,0)[:2]==[n,n+3])
for a in range(15):
 for b in range(15):
  for c in range(15):
   L=lengths(a,b,c)
   need(all(1<=j-i<=3 for i,j in zip(L,L[1:])))
   if L and L[0]:need(2*L[-1]<=5*L[0])
   if a!=1:need(L[-1]==b+c+a//2)
print('factorization lengths, gaps and maximum receipts pass')
def two_profile(T):
 H,_=Htable(T)
 F=(1<<len(T))-1
 return H[F],sum(H[s]*H[F^s] for s in range(F+1))
need(two_profile((126,60,120,112,96,64,2))==(17,448))
need(two_profile((24,1,3,6,14))==(11,100))
v,S=56551,3390208000
u,R=17,448
x,T=11,100
p=292235682641783
D=p-2*v*x*u
a,b=R-u,2*S*T-v*x
need(a*b==D)
need(((a+u)//R,(b+v*x)//(S*T))==(1,2))
need(all(p%i for i in range(3,isqrt(p)+1,2)))
print('prescribed mixed-core inverse',p,'factor pair',a,b,'7p',7*p)
