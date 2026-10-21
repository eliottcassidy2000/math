#!/usr/bin/env python3
"""Exact word-model controls. Synthetic alphabets are not new tournament censuses."""
import json, math, hashlib
from fractions import Fraction
from pathlib import Path
checks=0
def check(b,msg):
 global checks
 checks+=1
 if not b: raise RuntimeError(msg)
def dist(N,l):
 return [math.comb(l+N-2*m-1,l-1)*math.comb(N-2*m,m) for m in range(N//3+1)]
def conv(a,b,N):
 c=[0]*(N+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b[:N+1-i]): c[i+j]+=x*y
 return c
def coeffs(N,l):
 a=[0]*(N+1);a[0]=1
 for n in range(1,N+1): a[n]=a[n-1]+(a[n-3] if n>=3 else 0)
 out=[1]+[0]*N
 for _ in range(l):out=conv(out,a,N)
 return a,out
# Separate polynomial DP for gaps represented by coloured separator compositions.
for l in range(1,6):
 for N in range(0,35):
  words=[{0:1}]+[{} for _ in range(N)]
  for n in range(1,N+1):
   for m,v in words[n-1].items(): words[n][m]=words[n].get(m,0)+v
   if n>=3:
    for m,v in words[n-3].items(): words[n][m+1]=words[n].get(m+1,0)+v
  dp=[{0:1}]+[{} for _ in range(N)]
  for _ in range(l):
   nxt=[{} for _ in range(N+1)]
   for i in range(N+1):
    for j in range(N+1-i):
     for a,x in dp[i].items():
      for b,y in words[j].items():nxt[i+j][a+b]=nxt[i+j].get(a+b,0)+x*y
   dp=nxt
  d=dist(N,l)
  check(d==[dp[N].get(m,0) for m in range(N//3+1)],f'dist {N,l}')
  a,c=coeffs(N,l)
  check(sum(d)==c[N],f'count {N,l}')
  if l>1:
   _,rest=coeffs(N,l-1)
   moment=sum(g*a[g]*rest[N-g] for g in range(N+1))
   check(Fraction(moment,c[N])==Fraction(N,l),f'gap mean {N,l}')
# Large n are exact finite moment calculations, NOT proof of CLT or convergence.
rho=.6823278038280193;a=rho*rho
mu=a/(1+3*a);var=a*(1+a)/(1+3*a)**3
rows=[]
for l in [1,2,4]:
 for N in [60,180,540,1080]:
  d=dist(N,l);Z=sum(d)
  E=Fraction(sum(m*x for m,x in enumerate(d)),Z)
  E2=Fraction(sum(m*m*x for m,x in enumerate(d)),Z)
  V=E2-E*E
  row={'gaps':l,'neutral_vertices':N,'triangle_mean':float(E),'triangle_variance':float(V),'mean_minus_muN':float(E)-mu*N,'variance_minus_sigma2N':float(V)-var*N}
  check(abs(row['mean_minus_muN'])<10,'mean finite control')
  check(abs(row['variance_minus_sigma2N'])<10,'variance finite control')
  aa,cc=coeffs(N,l)
  if l==1: second=Fraction(1)
  else:
   _,rr=coeffs(N,l-1)
   second=Fraction(sum(g*g*aa[g]*rr[N-g] for g in range(N+1)),cc[N]*N*N)
  row['gap_second_moment_normalized']=float(second)
  row['dirichlet_target']=2/(l*(l+1))
  check(abs(float(second)-row['dirichlet_target'])<.03,'Dirichlet finite control')
  rows.append(row)
# Distinguish uniform isomorphism classes from labelled weighting: C3 has 2 labels,
# K1^3 has 6, so a 3-vertex neutral model has probabilities 1/2 vs1/4 for C3.
check(Fraction(1,2)!=Fraction(2,8),'labelling hostile')
# Explicit h=3^r dominant skeleton R5^r has total size5r; no invented atom census.
for r in range(4):
 N=18
 a,c=coeffs(N,r+1)
 check(c[N]==sum(dist(N,r+1)),f'R5 skeleton {r}')
# Exact mirror dependence controls on each fixed center parity and contribution.
for N in [0,1,3,12,60,180]:
 for l in [1,2,4]:
  d=dist(N,l);Z=sum(d)
  E=Fraction(sum(m*x for m,x in enumerate(d)),Z)
  V=Fraction(sum(m*m*x for m,x in enumerate(d)),Z)-E*E
  for c,e in [(0,0),(1,0),(3,1),(5,0)]:
   EM=Fraction(sum((2*m+e)*x for m,x in enumerate(d)),Z)
   VM=Fraction(sum((2*m+e)**2*x for m,x in enumerate(d)),Z)-EM*EM
   check(EM==2*E+e,'mirror mean')
   check(VM==4*V,'mirror variance')
   if V:check(VM!=2*V,'independent mirror hostile')
# Exact multiset sampler and converse rejection controls (actual K1/C3/R5 words).
from itertools import combinations
for r in range(3):
 for n in range(5*r,16):
  N=n-5*r;total=pal=0
  for m in range(N//3+1):
   a=N-3*m;L=a+m+r
   w=math.factorial(L)//(math.factorial(a)*math.factorial(m)*math.factorial(r))
   check(w==dist(N,r+1)[m],'sampler distribution')
   if sum(c%2 for c in (a,m,r))<=1:
    pc=math.factorial(L//2)//(math.factorial(a//2)*math.factorial(m//2)*math.factorial(r//2))
   else:pc=0
   actual=ap=0
   for rp in combinations(range(L),r):
    rem=[j for j in range(L) if j not in rp]
    for cp in combinations(rem,m):
     word=tuple(5 if j in rp else 3 if j in cp else 1 for j in range(L))
     actual+=1;ap+=word==word[::-1]
     # Sequential draws reproduce each individual word probability exactly.
     counts={1:a,3:m,5:r};prob=Fraction(1)
     for j,t in enumerate(word):prob*=Fraction(counts[t],L-j);counts[t]-=1
     check(prob==Fraction(1,w),'uniform sequential word')
   check(actual==w and ap==pc,'palindrome formula')
   total+=w;pal+=pc
  check((total-pal)%2==0,'converse pairs')
check(dist(0,2)==[1],'R5-only no nonpalindrome hostile')
root=Path(__file__).resolve().parent
proof=root/'RENEWAL_LIMIT_LAWS.md'
print(json.dumps({'checks':checks,'scope':'exact neutral-word coefficients and moments; no new strong-atom census; finite tests do not prove limits','proof_sha256':hashlib.sha256(proof.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rho':rho,'mu':mu,'sigma2':var,'moment_controls':rows},indent=2))
