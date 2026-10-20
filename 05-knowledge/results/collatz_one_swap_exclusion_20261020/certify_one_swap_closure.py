#!/usr/bin/env python3
"""Exact rational/continued-fraction and direct small-word certificate.
External input: real Matveev bound as stated by Bugeaud--Mignotte--Siksek,
arXiv:math/0403046, Theorem 9.4. Internal input: proved sqrt(7) content gate.
No floating-point value is used for a decision.
"""
from fractions import Fraction as F
from math import factorial,gcd
from hashlib import sha256
from pathlib import Path
import json
CHECKS=0
def need(b,*args):
 global CHECKS
 CHECKS+=1
 if not b:raise RuntimeError(args)
def logs(t,N=220):
 s=2*sum((t**(2*j+1)/F(2*j+1) for j in range(N)),F(0))
 err=2*t**(2*N+1)/((2*N+1)*(1-t*t))
 return s,s+err
L2,U2=logs(F(1,3));L3,U3=logs(F(1,2));alpha_lo=L3/U2;alpha_hi=U3/L2
Cupper=F(7,5)*30**5*23*F(7,10)*F(11,10)
need(Cupper<10**9,'constant');need(U2<F(7,10) and U3<F(11,10),'heights')
L54,U54=logs(F(1,9))
need(U2+12*(3*U2+U54)<30,'cutoff logarithm')
need(F(10**12,9)>10**9*31,'global cutoff')
threshold=100
exp_lower=sum((F(100,9)**j/factorial(j) for j in range(60)),F(0))
need(F(2*threshold)/exp_lower<L2,'Legendre threshold')
# CF extraction is certified because both exact interval endpoints have same floor.
lo,hi=alpha_lo,alpha_hi;hprev,h=0,1;kprev,k=1,0
cf=[];rows=[]
for i in range(100):
 ai=lo.numerator//lo.denominator
 need(ai==hi.numerator//hi.denominator,'CF ambiguity',i)
 cf.append(ai);hprev,h=h,ai*h+hprev;kprev,k=k,ai*k+kprev
 q,p=h,k;low=q*L2-p*U3;high=q*U2-p*L3
 if threshold<=p<10**12:
  if high<0:status='negative gap excludes positive cycle'
  else:
   need(low>0,'ambiguous sign',q,p)
   cap=min(p//9,256);upper=F(1,2**cap)
   need(low>upper,'unexcluded convergent',q,p)
   status='positive gap exceeds exponential upper bound'
  rows.append({'q':q,'p':p,'status':status,'gap_lower':str(low),'gap_upper':str(high)})
 if p>=10**12:break
 lo,hi=1/(hi-ai),1/(lo-ai)
else:raise RuntimeError('no terminal convergent')
# Every actual unequal adjacent swap at the unique positive critical length,
# for every 2<=p<100 and gcd(p,q)=1. No omitted rotations or boundary cut.
small=[];swaps=0
for p0 in range(2,threshold):
 qlo=(p0*alpha_lo).__ceil__();qhi=(p0*alpha_hi).__ceil__()
 need(qlo==qhi,'ambiguous critical length',p0);q0=qlo
 if gcd(p0,q0)!=1:
  small.append({'p':p0,'q':q0,'status':'noncoprime shape excluded','swaps':0});continue
 c=[((i+1)*p0)//q0-(i*p0)//q0 for i in range(q0)]
 D=2**q0-3**p0;need(D>1,'clock exception',q0,p0)
 rowcount=0;maxg=1
 for i in range(q0):
  j=(i+1)%q0
  if c[i]==c[j]:continue
  w=c.copy();w[i],w[j]=w[j],w[i]
  # Direct affine composition B; denominator doubles after each bit.
  B=0;pow2=1
  for bit in w:
   B=3**bit*B+bit*pow2;pow2*=2
  positions=[a for a,bit in enumerate(w) if bit]
  B2=sum(2**a*3**(p0-1-r) for r,a in enumerate(positions))
  need(B==B2,'carry implementations',p0,q0,i)
  g=gcd(D,B);need(g<D,'integer candidate',p0,q0,i)
  maxg=max(maxg,g);swaps+=1;rowcount+=1
 small.append({'p':p0,'q':q0,'status':'all actual swaps excluded','swaps':rowcount,'maximum_content':maxg})
report={'status':'EXACT certificate for stated one-swap theorem; not general Collatz','script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'checks':CHECKS,'log_series_terms':220,'global_p_cutoff_exclusive':10**12,'legendre_threshold':threshold,'matveev_prefactor_upper':str(Cupper),'continued_fraction_prefix':cf,'terminal_convergent':{'q':q,'p':p},'excluded_large_convergents':rows,'small_shapes':small,'small_swap_occurrences':swaps,'small_positive_coprime_shapes':sum(x['swaps']>0 for x in small),'log_intervals':{'log2_lower':str(L2),'log2_upper':str(U2),'log3_lower':str(L3),'log3_upper':str(U3)}}
print(json.dumps(report,indent=2))
