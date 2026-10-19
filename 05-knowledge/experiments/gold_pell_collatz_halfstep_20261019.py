#!/usr/bin/env python3
"""Exact Pell/Gold/Collatz half-step identities and finite carrier."""
from fractions import Fraction as F
import json
checks=0
def req(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def q(m,n):return n*n-m*n-(m-1)*(m-2)
def S(m,n):return 3-m-n,3-m-2*n
def R(m,n):return 2*m+3*n-3,3*m+5*n-6
def u(m,n):return (m+8*n-6)%11
def budget(z):return None if z==1 else (9*z-5)*pow(z-1,-1,11)%11
def g(b):return 2 if b is None else (None if b==1 else (1-2*b)*pow(1-b,-1,11)%11)
# Verify polynomial identities by coefficient-complete small grid.
for m in range(-3,4):
 for n in range(-3,4):
  req(q(*S(m,n))==q(m,n),'quadratic invariance')
  req(S(*S(m,n))==R(m,n),'affine square')
  req((5*m-6)**2-5*(2*n-m)**2==-20*q(m,n)-4,'norm identity')
  req(u(*S(m,n))==2*u(m,n)%11,'linear residue halfstep')
points=[F(x,13) for x in (1,2,4,8,16)]+[F(x) for x in (-5,-14,-7,-20,-10)]
labels=[3,9,5,4,1,8,10,7,6,2]
label=dict(zip(points,labels));inverse=dict(zip(labels,points))
def C(x):return 3*x+1 if x.numerator%2 else x/2
def D(x):return inverse[2*label[x]%11]
conic=[(m,n) for m in range(11) for n in range(11) if q(m,n)%11==0]
req(len(conic)==10 and set(map(lambda p:u(*p),conic))==set(range(1,11)),'conic bijection')
rows=[]
for m,n in conic:
 z=u(m,n);v=7*pow(z,-1,11)%11
 X=(z+v)*pow(2,-1,11)%11;Y=(z-v)*pow(8,-1,11)%11
 mm=(X+6)*pow(5,-1,11)%11;nn=(Y+mm)*pow(2,-1,11)%11
 req((mm,nn)==(m,n),'inverse norm coordinate')
 req(budget(u(*S(m,n)))==g(budget(z)),'budget halfstep')
 req(inverse[u(*R(m,n))]==C(inverse[z]),'standard Collatz recurrence')
 req(D(D(inverse[z]))==C(inverse[z]),'finite square root')
 rows.append({'m':m,'n':n,'label':z,'budget':budget(z),'Collatz':str(inverse[z])})
p=(2,2);orbit=[]
for j in range(12):
 req(q(*p)==0,'integer norm orbit')
 req((p[0]>=2 and p[1]>=p[0]) if j%2==0 else (p[0]<0 and p[1]<0),'alternating chamber')
 orbit.append({'step':j,'m':p[0],'n':p[1],'label':u(*p),'Collatz':str(inverse[u(*p)])});p=S(*p)
req({r['label'] for r in orbit if r['step']%2==0}=={1,3,4,5,9},'positive boundary residues')
print(json.dumps({'checks':checks,'modular_conic':rows,'integer_halfstep_orbit':orbit},indent=2))
