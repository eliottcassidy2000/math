def require(x):
 if not x: raise AssertionError("lift check failed")
#!/usr/bin/env python3
from random import Random
from math import prod
Q=(11,4,18,7); B=(5,2,8,3); I=(1,0,0,1); z=(25,34); primes=(7,13,19,31)
def mul(A,B,m):
 a,b,c,d=A;e,f,g,h=B
 return ((a*e+b*g)%m,(a*f+b*h)%m,(c*e+d*g)%m,(c*f+d*h)%m)
def power(A,k,m):
 R=I
 while k:
  if k&1:R=mul(R,A,m)
  A=mul(A,A,m);k//=2
 return R
def inv(A,m):
 a,b,c,d=A;t=pow((a*d-b*c)%m,-1,m);return tuple(x*t%m for x in (d,-b,-c,a))
def act(A,v,m):
 a,b,c,d=A;x,y=v;return ((a*x+b*y)%m,(c*x+d*y)%m)
def loops(m):
 L=power(B,576,m); return L,mul(mul(Q,L,m),inv(Q,m),m)
C={}
for p in primes:
 pp=p*p; require(all(x%p==y for x,y in zip(power(B,576,pp),I)))
 C[p]=[tuple(((act(A,z,pp)[i]-z[i])%pp)//p for i in range(2)) for A in loops(pp)]
 (u,v),(w,t)=C[p];d=(u*t-v*w)%p;require(d)
 print('p, tangent columns, determinant:',p,C[p],d)
rng=Random(15531)
for depths in ((2,2,2,2),(3,2,2,2),(2,3,2,2),(2,2,3,2),(2,2,2,3)):
 n=prod(p**a for p,a in zip(primes,depths)); L,M=loops(n)
 for _ in range(30):
  target=tuple((z[i]+53599*rng.randrange(n//53599))%n for i in range(2));cur=z; lev={p:1 for p in primes}
  for p,a in zip(primes,depths):
   while lev[p]<a:
    r=lev[p]; k=prod(q**(lev[q]-1) for q in primes if q!=p); e=p**(r-1)*k
    (u,v),(w,t)=C[p]; rhs=[((target[i]-cur[i])//p**r)%p for i in range(2)]; det=(u*t-v*w)%p
    x=((t*rhs[0]-w*rhs[1])*pow(det*k,-1,p))%p
    y=((-v*rhs[0]+u*rhs[1])*pow(det*k,-1,p))%p
    cur=act(power(M,y*e,n),act(power(L,x*e,n),cur,n),n);lev[p]+=1
  if cur!=target: raise AssertionError((depths,cur,target))
 print('constructive seed-fibre tests',depths,': 30 passed')
