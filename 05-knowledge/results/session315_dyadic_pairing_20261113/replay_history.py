from math import prod, gcd, isqrt
from fractions import Fraction
from itertools import combinations
from collections import Counter

pr=[p for p in range(2,98) if all(p%d for d in range(2,isqrt(p)+1))]
A=[]
for p in pr:
    if p==2: continue
    T=next(i for i in range(1,p) if pow(2,i,p)==1)
    A.extend((p,T,i) for i in range(T) if (99*pow(2,i,p)+23)%p==0)
assert A==[(5,4,3),(13,12,1),(17,8,1),(19,18,9),(29,28,27),(37,36,23),(53,52,25),(59,58,35),(61,60,0),(67,66,56),(83,82,15)]
live=lambda r:not any(r%t==a for _,t,a in A[:4])
ex=[(4,25,13),(2,35,29),(12,0,5),(6,56,11),(2,15,41)]
L=72*prod(p for _,_,p in ex)
N=[sum(prod(p-(r%m==a%m) for m,a,p in ex) for r in range(v,72,2) if live(r)) for v in (0,1)]
assert L==61209720 and N==[28657278,3696000]
assert [Fraction(x,L//2) for x in N]==[Fraction(103,110),Fraction(5600,46371)]
assert live(4) and live(5) and not live(17) and live((17+L//2)%72)
assert not live(45) and live((45+L//3)%72)
for m,a,p in ex:
    assert any(live(r) and r%m==a%m for r in range(72)) and gcd(L//p,p)==1
for Es,W in [((5,7,3,4,1,2,1,1,1,1),(92016,120,16,4,2,0,1)),((4,7,3,4,1,2,1,1,1,1),(76656,120,16,6,0,0,1))]:
    d=lambda j:prod(e//j+1 for e in Es)
    assert tuple(d(j)-d(j+1) for j in range(1,8))==W

def stats(A,cuts=None):
 n=len(A);N=(1<<n)-1
 s=[[0]*n for _ in range(N+1)];e=[[0]*n for _ in range(N+1)]
 for i in range(n):s[1<<i][i]=e[1<<i][i]=1
 for S in range(1,N+1):
  if S&(S-1)==0:continue
  for i in range(n):
   if S>>i&1:
    T=S^(1<<i)
    s[S][i]=sum(s[T][j] for j in range(n) if A[i]>>j&1)
    e[S][i]=sum(e[T][j] for j in range(n) if A[j]>>i&1)
 h=[1]+list(map(sum,s[1:]))
 ans={}
 for B in (range(N+1) if cuts is None else cuts):
  start=sum(s[N][i] for i in range(n) if B>>i&1)
  end=sum(e[N][i] for i in range(n) if not B>>i&1)
  sstar=sum(h[N^U]*sum(s[U][i] for i in range(n) if B>>i&1) for U in range(1,N+1))
  estar=sum(h[N^U]*sum(e[U][i] for i in range(n) if not B>>i&1) for U in range(1,N+1))
  split=sum(sum(s[U][i] for i in range(n) if B>>i&1)*sum(e[N^U][i] for i in range(n) if not B>>i&1) for U in range(1,N))
  ans[B]=((start,sstar),(end,estar,split))
 return h[N],sum(h[U]*h[N^U] for U in range(N+1)),ans

for A,p,H,Q,cuts in [((48,1,3,7,14,30),23,23,296,{1:((1,25),(22,244,24))}),((50,8,3,5,14,30),5,25,304,{1:((3,42),(23,252,39)),2:((2,28),(17,199,26))})]:
    h,q,D=stats(A,cuts)
    assert (h,q,D)==(H,Q,cuts)
    for i in range(6):
        K=1<<i
        for _ in range(6):
            for j in range(6):
                if K>>j&1: K|=A[j]
        assert K==63
    q2=q if p==23 else 2*q
    assert pow(q2,p-1,p*p)!=(1)
    assert all(pow(q2,i,p)!=1 for i in range(1,p-1))
    for _,(sS,eEb) in cuts.items():
        S=sS[1] if p==23 else 2*sS[1]
        E=eEb[1] if p==23 else eEb[1]+q
        assert S*E%p

def adj(n):
 return {i:tuple(j for j in range(1,n+1) if i!=j and isqrt(i+j)**2==i+j) for i in range(1,n+1)}
def paths(n):
 N=adj(n);A=(1<<n)-1
 def rec(P,U):
  if U==A:yield P;return
  for j in N[P[-1]]:
   if not U>>(j-1)&1:yield from rec(P+(j,),U|1<<(j-1))
 return list(rec((18,),1<<17))

_,_,padded_cuts=stats([r|64 for r in (50,8,3,5,14,30)]+[0], [1,2,5,12])
assert [(a[1],b[1]) for a,b in padded_cuts.values()]==[(84,556),(56,503),(152,506),(126,481)]
out=[]
for n in range(23,27):
    P=paths(n)
    D=Counter('end' if p[-1]==n else (p[p.index(n)-1],p[p.index(n)+1]) for p in P)
    out.append(D)
assert out==[{(2,13):2,(13,2):1},{},{(11,24):8,(24,11):2},{(10,23):6,(23,10):4,'end':2}]
print("S315 historical clock and transfer controls passed")

