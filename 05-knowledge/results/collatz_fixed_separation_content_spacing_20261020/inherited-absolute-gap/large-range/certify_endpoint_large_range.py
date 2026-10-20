#!/usr/bin/env python3
"""Exact rational/continued-fraction and direct small-word certificate.
External input: real Matveev bound as stated by Bugeaud--Mignotte--Siksek,
arXiv:math/0403046, Theorem 9.4. Internal input: proved endpoint content bound g < 52**ceil(sqrt(p)).
This file certifies the large-range closure only; small all-clock endpoints are separate.
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
need(F(4*threshold)/exp_lower<L2,'absolute-gap Legendre threshold')
need(L3>F(13,12),'log3 lower for short relation')
need(sum((F(4)**j/factorial(j) for j in range(20)),F(0))>52,'log52 upper')
need(F(13,12)-F(24,25)>F(1,9),'subexponential margin')
need(F(1,9)>F(10**9,10**12),'cutoff monotonicity')
# CF extraction is certified because both exact interval endpoints have same floor.
lo,hi=alpha_lo,alpha_hi;hprev,h=0,1;kprev,k=1,0
cf=[];rows=[]
for i in range(100):
 ai=lo.numerator//lo.denominator
 need(ai==hi.numerator//hi.denominator,'CF ambiguity',i)
 cf.append(ai);hprev,h=h,ai*h+hprev;kprev,k=k,ai*k+kprev
 q,p=h,k;low=q*L2-p*U3;high=q*U2-p*L3
 if threshold<=p<10**12:
  if high<0:
   absolute_lower=-high;sign='negative'
  else:
   need(low>0,'ambiguous sign',q,p)
   absolute_lower=low;sign='positive'
  cap=min(p//9,256);upper=F(2,2**cap)
  need(absolute_lower>upper,'unexcluded absolute gap',q,p)
  status=sign+' gap exceeds twice exponential upper bound'
  rows.append({'q':q,'p':p,'status':status,'gap_lower':str(low),'gap_upper':str(high)})
 if p>=10**12:break
 lo,hi=1/(hi-ai),1/(lo-ai)
else:raise RuntimeError('no terminal convergent')
need(len(rows)==18,'convergent count')
report={'status':'EXACT large-range certificate for the two endpoint chain families, not general two-edit Collatz','script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'checks':CHECKS,'log_series_terms':220,'global_p_cutoff_exclusive':10**12,'legendre_threshold':threshold,'matveev_prefactor_upper':str(Cupper),'continued_fraction_prefix':cf,'terminal_convergent':{'q':q,'p':p},'excluded_large_convergents':rows,'small_range':'Separate all-clock endpoint check required for p<100','log_intervals':{'log2_lower':str(L2),'log2_upper':str(U2),'log3_lower':str(L3),'log3_upper':str(U3)}}
print(json.dumps(report,indent=2))
