#!/usr/bin/env python3
from fractions import Fraction as F
from math import factorial
import json
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def beta(a,b):return F(factorial(a-1)*factorial(b-1),factorial(a+b-1))
rows=[]
for a in range(1,5):
 for b in range(2,7):
  c=F(a,b-1)
  for J in range(0,25):
   ray=sum((beta(j+a,b)/beta(a,b) for j in range(1,J+1)),F(0))
   tail=beta(J+1+a,b-1)/beta(a,b)
   need(ray+tail==c,'full root incoming mass and tail')
  if b>2:
   inv1=F(a+b-1,b-1);inv2=F((a+b-1)*(a+b-2),(b-1)*(b-2));mean1=F(b,a+b)
   bound=(6*inv2+3*inv1-mean1-5)/c
   partial=F(0)
   for n in range(3,1002,6):
    x=n;t=K=clock=0
    while x!=1:
     z=3*x+1;v=0
     while z%2==0:z//=2;v+=1
     x=z;t+=1;K+=(v-1)//2;clock+=v+1
     if t>10000:raise RuntimeError('uncertified finite source')
    weight=beta(K+a,t-1+b)/beta(a,b)/c
    partial+=weight*clock
   need(partial<=bound,'different-law finite clock partial bounded')
   rows.append({'alpha':a,'beta':b,'leaf_mass_before_normalization':str(c),'literal_clock_upper':str(bound),'finite_source_partial':str(partial)})
need(F(1,2)*(6*F(3*2,2*1)+3*F(3,2)-F(3,4)-5)!=F(67,2),'normalization inversion hostile')
need(2*(6*F(3*2,2*1)+3*F(3,2)-F(3,4)-5)==F(67,2),'Beta1,3 bound')
need(F(2,2-1)==2,'Beta2,2 inherited deficit mass')
print(json.dumps({'status':'PASS','checks':checks,'scope':'different normalized beta measures; original native lambda remains Beta1,2','rows':rows},indent=2,sort_keys=True))
