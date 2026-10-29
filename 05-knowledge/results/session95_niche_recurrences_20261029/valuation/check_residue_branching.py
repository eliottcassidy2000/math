#!/usr/bin/env python3
"""Exact finite checks; no empirical extrapolation and no external packages."""
from collections import defaultdict
from fractions import Fraction as Q
from math import comb,lcm
import json
CHECKS=0

def check(ok,context):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(context)

def selector(d):
 h=1
 while (2**h-1)%d: h+=1
 K=(2**h-1)//d
 return h,{r for r in range(1,h+1) if K & (1<<(h-r))}

def mul(a,b,T):
 c=[0]*(T+1)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b[:T+1-i]):
    if y:c[i+j]+=x*y
 return c

def power(a,k,T):
 b=[1]+[0]*T
 for _ in range(k):b=mul(b,a,T)
 return b

def formula(d,h,S,n,T):
 if (n-1)%d:return [0]*(T+1)
 m=(n-1)//d;L=n-m
 A=[int(1<=r<=h and r not in S) for r in range(T+1)]
 B=[int(r in S) for r in range(T+1)]
 P=mul(power(A,L,T),power(B,m,T),T)
 H=[comb(n+r//h-1,r//h) if r%h==0 else 0 for r in range(T+1)]
 C=mul(P,H,T); shapes=comb(n,m)//n
 check(comb(n,m)%n==0,('shape integrality',d,n))
 return [shapes*x for x in C]

def frontier(d,h,S,N,T):
 # Each exact valuation word of total D has probability 2^-D.
 states={(1,0):1};ans={}
 for n in range(1,N+1):
  nxt=defaultdict(int);out=[0]*(T+1)
  for (q,t),w in states.items():
   for a in range(1,T-t+1):
    r=(a-1)%h+1; qq=q-1+(d if r in S else 0)
    if qq==0:out[t+a]+=w
    elif qq<=N-n and t+a+qq<=T:nxt[qq,t+a]+=w
  ans[n]=out;states=nxt
 return ans

rows=[]
for d,N,T in [(3,13,25),(5,11,25),(7,15,28),(9,10,25),(11,12,30),(15,16,30)]:
 h,S=selector(d); bank=frontier(d,h,S,N,T)
 for n in range(1,N+1):
  f=formula(d,h,S,n,T)
  for t in range(T+1):check(bank[n][t]==f[t],('frontier/coefficient',d,n,t))
 rows.append({'arity':d,'period':h,'internal_residues':sorted(S),'max_vertices':N,'max_divisions':T})
for d in range(3,102,2):
 h,S=selector(d);M=2**h-1
 check(sum(2**(h-r) for r in S)*d==M,('critical',d))
 check(all((2**j-1)%d for j in range(1,h)),('minimal',d))
 for j in range(1,3*h+1):
  if j%h==0:
   K=(2**j-1)//d
   sj={r for r in range(1,j+1) if K & (1<<(j-r))}
   check(all((r in sj)==(((r-1)%h+1) in S) for r in range(1,j+1)),('period lift',d,j))
for d in range(3,42,2):
 h,S=selector(d);M=2**h-1
 probs={r:Q(2**(h-r),M) for r in range(1,h+1)}
 mi=sum(r*p for r,p in probs.items() if r in S)*d
 ml=sum(r*p for r,p in probs.items() if r not in S)*Q(d,d-1)
 vi=sum(r*r*p for r,p in probs.items() if r in S)*d-mi*mi
 vl=sum(r*r*p for r,p in probs.items() if r not in S)*Q(d,d-1)-ml*ml
 delta=mi-ml
 for m in range(5):
  n=d*m+1;L=n-m
  check(m*mi+L*ml+Q(n*h,M)==2*n-delta/d,('general mean',d,m))
  raw=m*vi+L*vl+Q(n*h*h*2**h,M*M)
  check(raw==n*(2-Q(d-1,d*d)*delta*delta)+(vl-vi)/d,('general variance',d,m))
for h in range(2,9):
 d=2**h-1;S={h};M=d
 weights={r:Q(2**(h-r),M) for r in range(1,h+1)}
 mu=sum(r*w for r,w in weights.items() if r<h)/(1-Q(1,d))
 var=sum(r*r*w for r,w in weights.items() if r<h)/(1-Q(1,d))-mu*mu
 check(mu==Q(2*(d-h),d-1),('leaf mean',h))
 check(var==Q(6*d-2*h*h-4*h,d-1)-mu*mu,('leaf variance',h))
 for m in range(9):
  n=d*m+1;L=n-m
  mean=h*m+L*mu+Q(n*h,d)
  check(mean==2*n+Q(d*(2-h)-h,d*(d-1)),('conditional mean',h,m))
check(formula(3,2,{2},1,8)[2]==0,'ternary parity')
check(formula(7,3,{3},1,8)[1:3]==[1,1],'seven parity hostile')
check(formula(15,4,{4},1,8)[1:3]==[1,1],'fifteen parity hostile')
for h in range(2,8):
 for k in range(2,8):
  H=lcm(h,k);M=2**H-1
  joint=sum(Q(2**(H-r),M) for r in range(1,H+1) if r%h==r%k==0)
  check(joint==Q(1,M),('shared residue',h,k))
check(Q(3*7,63)-1==Q(-2,3),'negative shared covariance')
check(Q(3*15,15)-1==2,'nested covariance')
# Literal native valuation cylinder for the inherited growth control.
a=[2,1,1,1];D=sum(a);matches=[]
for x in range(1,2**(D+1),2):
 y=x;observed=[]
 for _ in a:
  z=3*y+1;j=0
  while z%2==0:z//=2;j+=1
  observed.append(j);y=z
 if observed==a:matches.append(x)
check(matches==[41],('guard',matches))
B=t=0
for j in a:B=3*B+2**t;t+=j
check((B,t)==(103,5),'affine carry')
check((3**4*41+B)//2**t==107,'native growth')
print(json.dumps({'status':'PASS','checks':CHECKS,'frontier_comparisons':rows,'source':'iid valuations, probability 2^-a','no_native_convergence_claim':True},sort_keys=True))
