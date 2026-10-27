#!/usr/bin/env python3
from fractions import Fraction as F
from math import comb
import json
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def S(n,j):return ((3*n+1)*4**j-1)//3
def bd(n):
 b=F(3*n,3*n+1);return b,1-b*b
def step(n):
 z=3*n+1;a=0
 while z%2==0:z//=2;a+=1
 return z,a
def word(n):
 out=[]
 while n!=1:
  n,a=step(n);out.append(a)
  if len(out)>10000:raise RuntimeError('unverified finite source')
 return out
def W(n):
 a=word(n);t=len(a);k=sum((x-1)//2 for x in a)
 return F(1) if n==1 else F(2,(t+k+1)*comb(t+k,k))
rows=[]
for h in range(17):
 q=2**h;M=2*q;r=(-pow(3,-1,M))%M
 j=1
 while 2*(j+1)<h+1:j+=3
 m=S(1,j);base=(2*m-1)//3
 need(m%M==r and m%3==2 and word(m)==[2*j+2],'unit target with actual root phase')
 ja=2
 while 2*(ja+1)<h+1:ja+=3
 jz=(-base)%3
 while 2*jz<h+1:jz+=3
 certs=[]
 for k in range(8):
  a=S(1,ja+3*k);z=S(base,jz+3*k)
  need(a%M==r and z%M==r and a%3==0 and z%3==0,'actual same residue lambda sources')
  need(step(a)[0]==1 and step(a)[1]%2==0,'positive native phase')
  need(step(z)[0]==m and step(z)[1]%2==1,'negative native phase shared target')
  ba,da=bd(a);bz,dz=bd(z);A=W(a)/da;B=W(z)/dz;D=B*(1+bz)**2
  lower=A*D/(A*ba*ba+D);opt=A*ba/(A*ba*ba+D)
  need(A*(1-ba*opt)**2+D*opt*opt==lower,'exact relaxed minimizer')
  need(lower>=min(A,D)/2,'growing harmonic lower bound')
  for x in [F(-2),F(-1),F(0),F(1,2),F(1),F(2)]:
   need(A*(1-ba*x)**2+D*x*x>=lower,'signal grid')
  certs.append({'k':k,'root_predecessor':a,'odd_phase_predecessor':z,'lower_energy_per_root_squared':str(lower)})
 need(F(certs[-1]['lower_energy_per_root_squared'])>10**6,'finite high-cap hostile')
 rows.append({'index_modulus':q,'integer_modulus':M,'common_odd_residue':r,'target':m,'target_root_ray_depth':j,'certificates':certs})

def profile(n):return -1 if n%3==2 else 1
partial=F(0)
for n in range(3,1002,6):
 y,a=step(n);b,d=bd(n);force=profile(y)-(-1)**a*b*profile(n)
 need(profile(y)==(-1)**a,'ternary target phase')
 need(force*force/d==F(1,6*n+1),'ternary finite first-stage cost')
 partial+=W(n)*force*force/d
need(partial<F(1,19),'ternary partial below global bound')
hostile=[]
for k in range(13):
 n=S(1,3*k+1);y,a=step(n);b,d=bd(n)
 need(y==1 and a%2==0 and profile(n)==-1,'ternary all-depth root hostile')
 cost=(profile(y)-b*profile(n))**2/d
 need(cost==6*n+1,'all-depth expansive invoice')
 hostile.append({'k':k,'source':n,'occupation_weight':str(W(n)),'occupation_invoice':str(W(n)*cost)})
print(json.dumps({'status':'PASS','checks':checks,'dyadic_moduli':rows,'ternary_first_stage':{'sources':167,'partial_paid':str(partial),'global_bound':'1/19'},'ternary_all_depth_hostile':hostile,'scope':'dyadic periodic first-stage no-go, not all moduli; odd modulus first-stage repair fails occupation energy'},indent=2,sort_keys=True))
