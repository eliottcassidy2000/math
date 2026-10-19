#!/usr/bin/env python3
"""Exact controls for the common-center AP-petal capacity obstruction.
No contact search and no floating-point arithmetic enter this verifier.
"""
from math import gcd,isqrt,prod
from collections import Counter
import json
checks=Counter()
def check(c,n):
 checks[n]+=1
 if not c: raise RuntimeError(n)
def factor(n):
 out={};p=2
 while p*p<=n:
  while n%p==0: out[p]=out.get(p,0)+1;n//=p
  p+=1
 if n>1:out[n]=out.get(n,0)+1
 return out
def capacity(H):return (prod(2*e+1 for p,e in factor(H).items() if p%4==1)-1)//2
def petals(H):
 out=[]
 for n in range(1,isqrt(H)+1):
  for m in range(1,n):
   c=n*n+m*m
   if c>H:break
   if (n-m)%2==0 or gcd(n,m)>1 or H%c:continue
   a,b=n*n-m*m,2*n*m
   out.append((a,b,c))
 return out
def direct_ap(H):
 out=[];M=2*H*H
 for p in range(isqrt(M)+1,2*H):
  q=isqrt(8*H*H-p*p)
  if p*p+q*q==8*H*H:
   d=4*H*H-p*p
   out.append((M-d,M,M+d))
 return out
def from_ppt(H,T):
 a,b,c=T;g=H//c
 return tuple(2*g*g*z for z in (c*c-4*a*b,c*c,c*c+4*a*b))
samples=[]
for H in list(range(1,151))+[205,325,533,1105,1681,32045]:
 ts=petals(H);good=[t for t in ts if t[2]**2>4*t[0]*t[1]]
 check(len(ts)==capacity(H),'gaussian_product_vs_euclid')
 check(sorted(from_ppt(H,t) for t in good)==sorted(direct_ap(H)), 'direct_square_vs_ppt')
 for t in good:
  a,M,b=from_ppt(H,t)
  check(0<a<M<b<4*H*H and a+b==2*M,'positive_AP_labels')
  check(all(isqrt(z)**2==z for z in (a+M,a+b,M+b)),'all_three_square_edges')
 if H in [41,65,82,123,145,205,325,533,1105,1681,32045]:
  samples.append({'H':H,'unfiltered_capacity':len(ts),'positive_petals':len(good)})
# Exact capacity lower optimization. Rearranging exponents to descending order
# on the increasing primes 1 mod 4 cannot increase H and retains the capacity.
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
ps=[n for n in range(5,500) if n%4==1 and prime(n)]
def minimum_capacity_scale(k):
 target=2*k+1;initial=1;D=1;r=0
 while D<target:initial*=ps[r];D*=3;r+=1
 best=[initial,[1]*r]
 def visit(j,last,value,D,es):
  if D>=target:
   if value<best[0]:best[:]=[value,es]
   return
  if j>=len(ps):raise RuntimeError('prime budget exhausted')
  v=value
  for e in range(1,last+1):
   v*=ps[j]
   if v>=best[0]:break
   visit(j+1,e,v,D*(2*e+1),es+[e])
 bound=0;z=initial
 while z>=5:z//=5;bound+=1
 visit(0,bound,1,1,[])
 check(capacity(best[0])>=k,'optimization_positive_control')
 return {'k':k,'minimum_unfiltered_H':best[0],'exponents':best[1],
         'center_lower_bound':2*best[0]**2,
         'capacity_at_bound':capacity(best[0])}
opt=[minimum_capacity_scale(k) for k in [1,2,4,8,16,32,64,128,256,512,1024]]
# Independent exhaustive minimization for small targets.
for r in opt:
 if r['minimum_unfiltered_H']<=4000:
  check(all(capacity(H)<r['k'] for H in range(1,r['minimum_unfiltered_H'])),
        'optimization_vs_exhaustive')
# Explicit max count for a bounded common scale, via sieve (not pair contacts).
B=100000;counts=[0]*(B+1);unfiltered=[0]*(B+1)
for n in range(2,isqrt(B)+1):
 for m in range(1,n):
  c=n*n+m*m
  if c>B:break
  if (n-m)%2==0 or gcd(n,m)>1:continue
  good=c*c>4*(n*n-m*m)*(2*n*m)
  for H in range(c,B+1,c):
   unfiltered[H]+=1
   if good:counts[H]+=1
check(all(unfiltered[H]==capacity(H) for H in range(1,B+1)), 'capacity_sieve_complete_100000')
record=[];old=0
for H,n in enumerate(counts):
 if n>old:record.append({'H':H,'positive_petals':n,'unfiltered_capacity':unfiltered[H]});old=n
out={'status':'FINITE-EXACT; universal theorem has separate elementary proof',
     'controls':dict(checks),'sample_centers':samples,
     'minimum_capacity_scales':opt,
     'positive_capacity_records_through_100000':record,
     'scope':'Counts AP square-sum triangles with one common center. No inter-petal contacts tested.'}
print(json.dumps(out,indent=2))
