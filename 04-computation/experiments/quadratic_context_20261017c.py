#!/usr/bin/env python3
"""Exact order/conductor, golden residue, and norm-19 identities."""
import json
from math import gcd,isqrt
import sympy as s
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def forms(D):
 out=[]
 for a in range(1,isqrt(-D//3)+1):
  for b in range(-a,a+1):
   if (b*b-D)%(4*a):continue
   c=(b*b-D)//(4*a)
   if c<a or gcd(gcd(a,b),c)>1 or ((abs(b)==a or a==c) and b<0):continue
   out.append((a,b,c))
 return out
def golden_order(m):
 a,b=1,0
 for k in range(1,6*m*m+1):
  a,b=b%m,(a+b)%m
  if (a,b)==(1,0):return k
 raise RuntimeError('finite unit order not found')
def main():
 rows=[]
 for n in (7,8,11,15,19,21,37):
  d=n
  for p in range(2,isqrt(n)+1):
   while d%(p*p)==0:d//=p*p
  D0=-d if -d%4==1 else -4*d;D=-4*n;f=isqrt(D//D0)
  need(f*f*D0==D,'order discriminant retains the exact conductor')
  fs=forms(D);fund=forms(D0);amb=lambda x:x[1]==0 or abs(x[1])==x[0] or x[0]==x[2]
  rows.append(dict(n=n,order_discriminant=D,fundamental_discriminant=D0,conductor=f,forms=fs,fundamental_forms=fund,idoneal=all(amb(x) for x in fs),golden_unit_order=golden_order(n)))
 need([r['idoneal'] for r in rows]==[True,True,False,True,False,True,True],'11 and19 are the conductor-sensitive hostile controls')
 need(len(forms(-44))==len(forms(-76))==3,'odd class groups are invisible to the genus quotient')
 x,t=s.symbols('x t');f=x**3+2*x*x-4*x-4
 h=s.resultant(f.subs(x,x-t),t*t-t-1,t)
 need(s.Poly(h,x).degree()==6 and h.subs(x,0)==-19,'primitive cubic-golden sum has norm minus19')
 need(s.rem((4*t-3)*(4*(1-t)-3),t*t-t-1,t)==-19,'golden displacement norm')
 need(s.rem((t+2)*(3-t),t*t-t-1,t)==5,'pentagon quadratic radicand norm')
 need(s.discriminant(f,x)==592,'cubic polynomial discriminant 4 times field discriminant148')
 need(s.discriminant(x**3-x*x+x+1,x)==-44,'the conductor11 cubic is not the totally real cubic carrier')
 branches=[]
 for a in (5,15):
  need((a*a-a-1)%19==0,'both golden roots modulo19')
  branches.append(dict(phi=a,order=s.n_order(a,19),delta=(4*a-3)%19,sum=(4+a)%19,pentagon_radicand=(a+2)%19))
 need(pow(3,9,19)==18,'adjoining real sqrt3 exchanges the two primes above every inherited fixed19 prime')
 degree_two=[]
 for p in s.primerange(3,5000):
  if p in (5,7,11,19,37) or p%4!=3 or p%5!=4:continue
  if any(s.legendre_symbol(d,p)!=1 for d in (3,7,11)):continue
  ws=sorted(int(z)%p for z in s.polys.polytools.ground_roots(f,modulus=p))
  phis=sorted(int(z)%p for z in s.polys.polytools.ground_roots(x*x-x-1,modulus=p))
  good=[a for a in phis if s.legendre_symbol((4*a-3)%p,p)==1]
  if ws and good:degree_two.append(dict(p=int(p),w=ws,phi=good))
 need(degree_two[0]['p']==1319,'first qualifying unramified odd degree-two residue prime in the declared search')
 # Literal mixed equilateral edges, with source0 unchanged.
 eta=(1+s.sqrt(3)*s.I)/2
 need(s.expand_complex(eta*s.conjugate(eta))==1 and s.expand_complex((eta-1)*(s.conjugate(eta)-1))==1,'0,v,eta*v is a unit triangle for every unit v')
 return dict(status='FINITE-EXACT forms and identities; no primality or convergence criterion is inferred',rows=rows,prime19=branches,cubic_polynomial=str(h),escaped_field_prime_search=dict(bound_exclusive=5000,excluded_bad_primes=[2,5,7,11,19,37],degree_two_fixed_primes=degree_two),checks=checks)
if __name__=='__main__':print(json.dumps(main(),default=int,sort_keys=True,separators=(',',':')))
