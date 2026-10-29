#!/usr/bin/env python3
"""Exact independent Burnside/core controls for deletion-harmonic series."""
from fractions import Fraction as F
from math import factorial,gcd,comb
from collections import Counter
import json
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def parts(n,minimum=1):
 if not n:yield ();return
 for p in range(minimum,n+1,2):
  for q in parts(n-p,p):yield (p,)+q
def energy(mu):return sum((p-1)//2 for p in mu)+sum(gcd(x,y) for i,x in enumerate(mu) for y in mu[i+1:])
def centralizer(mu):
 z=1
 for p,k in Counter(mu).items():z*=p**k*factorial(k)
 return z
def count(n):
 s=sum((F(2**energy(mu),centralizer(mu)) for mu in parts(n)),F(0))
 need(s.denominator==1,'Burnside integrality');return s.numerator
def cores(N):
 for m in range(N+1):
  for mu in parts(m,3):yield mu
def main():
 N=26;a=[count(n) for n in range(N+1)];h=[a[n]-(a[n-1] if n else 0) for n in range(N+1)]
 need(h[:11]==[1,0,0,1,2,8,44,400,6424,184656,9541520],'initial harmonics')
 beta=[F(factorial(n)*a[n],2**comb(n,2)) for n in range(N+1)]
 eta=[F(factorial(n)*h[n],2**comb(n,2)) for n in range(N+1)]
 for n in range(1,N+1):need(eta[n]==beta[n]-n*F(2)**(1-n)*beta[n-1],'normalized difference')
 core_b=[F(0)]*(N+1);core_h=[F(0)]*(N+1);bank={}
 for mu in cores(N):
  m=sum(mu);c=len(mu);e=energy(mu);A=F(factorial(m),centralizer(mu))*F(2)**(e-comb(m,2));lam=F(2)**(-(m-c));k=(m-c)//2
  bank.setdefault(k,[]).append((m,mu))
  for n in range(m,N+1):
   core_b[n]+=A*comb(n,m)*lam**(n-m)
   core_h[n]+=A*comb(n,m)*lam**(n-m)
   if n>=m+1:core_h[n]-=(m+1)*A*F(2)**(-m)*comb(n,m+1)*(lam/2)**(n-m-1)
  if m<=12:
   for z in [F(1,3),F(2,3),F(-1,3),F(-2)]:
    w=z/2;den=1-lam*w
    val=A*w**m/den**(m+1)
    der=(A*m*w**(m-1)/den**(m+1) if m else 0)+A*(m+1)*lam*w**m/den**(m+2)
    need(z*val+z*z*der/2==(m+1)*A*F(2)**(-m)*z**(m+1)/den**(m+2),'paired rational term')
 need(core_b==beta,'inherited core coefficients independently recover Burnside')
 need(core_h==eta,'paired core coefficients equal harmonic transform')
 for k in range(9):
  b=bank[k];mx=max(m for m,mu in b);top=[mu for m,mu in b if m==mx]
  need(mx==3*k and top==[(3,)*k],'unique maximum core and pole order')
  need(3*k+1==(3*(2*k))//2+1 and 3*k+2==(3*(2*k+1))//2+1,'interlaced orders')
 for n in range(1,N+1):
  tr=sum(F((n-j)*h[j],n) for j in range(n+1))
  need(n*tr==sum(a[:n]),'trace telescoping')
  nexttrace=sum(F((n+1-j)*h[j],n+1) for j in range(n+1))
  need((n+1)*nexttrace-n*tr==a[n],'count from consecutive traces')
 need(eta[:6]==[F(1),F(0),F(0),F(3,4),F(3,4),F(15,16)],'normalized initial terms')
 print(json.dumps({'status':'PASS','checks':checks,'orders_checked':[0,N],'harmonic_dimensions_first11':h[:11],'pole_orders_at_2_power_r':[3*r//2+1 for r in range(12)],'scope':'Exact coefficients and core identities; global continuation inherited, not inferred from samples'},indent=2,sort_keys=True))
if __name__=='__main__':main()
