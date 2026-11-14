from itertools import product
from collections import defaultdict
from functools import lru_cache
cores=[(14,28,24,16,1),(10,28,25,16,1)]
def profile(A,U):
 @lru_cache(None)
 def P(s,v,w):
  if s==(1<<v):return int(v==w)
  if v==w:return 0
  t=s^(1<<v)
  return sum(P(t,z,w) for z in range(len(A)) if (A[v]&t)>>z&1)
 @lru_cache(None)
 def H(s):
  return sum(P(s,v,w) for v in range(len(A)) for w in range(len(A)) if s>>v&1 and s>>w&1) if s else 1
 f=(1<<len(A))-1; es=os=ee=oe=c=0
 for s in range(1,f+1):
  for v,w in product(range(len(A)),repeat=2):
   if s>>v&1 and s>>w&1:
    wgt=P(s,v,w)*H(f^s)
    if s.bit_count()%2:
     if U>>v&1:os+=wgt
     if not U>>w&1:oe+=wgt
    else:
     if U>>v&1:es+=wgt
     if not U>>w&1:ee+=wgt
     if U>>v&1 and not U>>w&1:c+=wgt
 return H(f),sum(H(s)*H(f^s) for s in range(f+1)),es+os,os-es,ee+oe,oe-ee,H(f)+2*c
from math import isqrt
source=[
((254,220,248,240,224,194,128,0),182,9,89),
((254,252,184,240,224,192,132,0),202,9,91),
((238,252,248,240,225,192,128,0),169,7,93)]
def check(ok):
 if not ok: raise AssertionError("h9 bilateral replay")
for T,U,t,k in source:
 check(profile(T,U)==(9,672,245,9,209,t,k))
def g(z):
 if z<2:return 1
 x=z//2
 return g(x)+(x+x%2 if not(z%2) else 0)
def val(x,y,t,k):
 return k+9*(g(x)+g(y)-2)+209*x+245*y+672*x*y+t*(x%2)+9*(y%2)
def p2(x):
 return x.bit_length()-1 if x>0 and x&(x-1)==0 else None
def axis(N,a,b,c,d):
 M=N-d
 if M<=0:return set()
 t=(M&-M).bit_length()-1
 out=set()
 if t>=1:
  u=1<<t
  q,z=divmod(M-b*u,a*u+c)
  y=p2(q)
  if z==0 and y is not None and y>t:out.add((t,y))
  q,z=divmod(M-c*u,a*u+b)
  x=p2(q)
  if z==0 and x is not None and x>t:out.add((x,t))
 v=(b+c)&-(b+c)
 k=t-(v.bit_length()-1)
 if k>=1:
  u=1<<k
  if M==a*u*u+(b+c)*u:out.add((k,k))
 return out
for a,b,c,R,X,Y in [(13440,-11305,-2189,(1905,1907,1909),lambda i:10*(1<<i)-2,lambda j:2*(1<<j)-2),
 (4704,-463,-2989,(361,363),lambda i:(1<<i)-1,lambda j:7*(1<<j)-1)]:
 L=defaultdict(set)
 check((a&-a)>(abs(b+c)&-abs(b+c)))
 for k,m in product(range(1,32),repeat=2):
  for (T,U,t,H),d in zip(source,R):
   N=a*2**(k+m)+b*2**k+c*2**m+d
   check(N==val(X(k),Y(m),t,H))
   check((k,m) in axis(N,a,b,c,d))
   if N<100000:L[N,d].add((k,m))
 for N in range(1,100000,2):
  for d in R:check(axis(N,a,b,c,d)==L[N,d])
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
check(all(prime(p) for p in (374351,53479,8831,61819)))
check(374353==7*53479 and 61817==7*8831)
check(pow(2,42,147)==1)
print("Three full marked profiles, exact inverses and numerical prime guards verified.")
for m,ref in [(12,(9,84,34,4,41,5,29)),(13,(9,84,45,5,30,4,31))]:
 check(profile(cores[1],m)==ref)
for m,k in [(15,29),(23,31),(14,33)]:
 h,q,p,s,e,t,v=profile(cores[1],m)
 check((h,q,p+e,s+t,v)==(9,84,87,7,k))
for x in range(1,2048):
 for a,b,D in [(75,9,(11,13)),(87,7,(11,13,15))]:
  for d in D:
   N=84*x*x+a*x+18*g(x)+b*(x%2)+d
   check(isqrt(N//84)==x)
check(prime(73961) and prime(517729) and 7*73961==517727)
print("Diagonal unit flip, three-arm profiles, square-root inverses and extra prime pair verified.")
from itertools import permutations
T=cores[1]
anti=[r for r in permutations(range(5)) if all(
 ((T[i]>>j)&1)!=((T[r[i]]>>r[j])&1) for i in range(5) for j in range(i+1,5))]
check(anti==[(0,4,3,2,1)])
U=(12,13,15,23,14)
check([sum(1<<anti[0][i] for i in range(5) if not ((u>>i)&1)) for u in U]==[19,18,2,4,3])
print("Rigid-core anti-isomorphism and converse-disjoint masks verified.")
