#!/usr/bin/env python3
"""Exact chart, boundary, and integer-cycle controls for the new chain bound."""
from fractions import Fraction as F
from math import gcd,isqrt,lcm
from collections import Counter
from hashlib import sha256
import json,sys
import sympy as s
sys.set_int_max_str_digits(0)
checks=Counter()
def check(b,n):
 checks[n]+=1
 if not b:raise RuntimeError(n)
def fp(t):
 return t*(t**4-8*t**3+2*t*t+1)*(5*t**6-12*t**5+23*t**4-17*t*t+12*t-3)/((2*t**5-t**4-4*t**3-2*t*t+2*t-1)*(t**6+20*t**5-29*t**4+16*t**3+3*t*t-4*t+1))
def fm(t):return 1/fp(1/t)
def f(t):return 4*t*(t*t-1)/(t*t+1)**2
def ppt(t):
 n,m=t.numerator,t.denominator;d=2 if n%2 and m%2 else 1
 a,b,c=(n*n-m*m)//d,2*n*m//d,(n*n+m*m)//d
 check(a*a+b*b==c*c and gcd(a,b)==1,'primitive_ppt')
 if a%2==0:a,b=b,a
 return a,b,c
def square_root(r):
 n,d=isqrt(r.numerator),isqrt(r.denominator)
 check(n*n==r.numerator and d*d==r.denominator,'exact_rational_square')
 return F(n,d)
t,h=s.symbols('t h');symbolic={}
for eps,Fmap in [(1,fp),(-1,fm)]:
 n,d=s.fraction(s.cancel(Fmap(t)))
 check(max(s.degree(n,t),s.degree(d,t))==11,'degree_exactly_11')
 check(s.degree(s.gcd(n,d),t)==0,'coprime_map_forms')
 n,d=s.Poly(n.subs(t,1+h),h),s.Poly(d.subs(t,1+h),h)
 if d.nth(0)<0:n,d=-n,-d
 p1,r1=s.div(n-(1+h)*d,s.Poly(h,h));p2,r2=s.div((1+4*h)*d-n,s.Poly(h,h))
 check(r1.is_zero and r2.is_zero,'growth_divisibility')
 out={}
 for name,p in [('D',d),('growth_lower',p1),('growth_upper',p2),('derivative',n.diff()*d-n*d.diff())]:
  negative=sum((F(-int(p.nth(i)),64**i) for i in range(1,p.degree()+1) if p.nth(i)<0),F(0))
  gap=F(int(p.nth(0)))-negative
  check(gap>0,'uniform_sign_polynomial')
  out[name]={'constant':int(p.nth(0)),'negative_tail':str(negative),'positive_gap':str(gap)}
 symbolic[str(eps)]=out
# Exact odd-prime denominator transport.
NN,DD=s.fraction(s.cancel(fp(t)))
PP,RR=s.div(NN*NN+DD*DD,(t*t+1)**2,t)
check(RR==0,'hypotenuse_square_factor')
check(s.resultant(NN,DD,t)==-2**102,'odd_prime_good_reduction')
check(s.rem(PP,t*t+1,t)==-49152*(t+1),'old_prime_exact_doubling')
# Explicit lower bound f(1+h)>=h throughout 0<h<=1/3.
p=s.Poly(4*(1+h)*(2+h)-((1+h)**2+1)**2,h)
negative=sum((F(-int(p.nth(i)),3**i) for i in range(1,p.degree()+1) if p.nth(i)<0),F(0))
check(F(int(p.nth(0)))-negative>0,'spread_dominates_cusp_parameter')
boundaries=[]
for T in [(9,40,41),(21,220,221),(99,4900,4901)]:
 a,b,c=T;f0=F(2*a*b,c*c)
 check(0<f0<F(1,2),'admissible_boundary_seed')
 for eps in [1,-1]:
  k=F(abs(a+eps*b),c);K=k*k;u=F(eps);z=2*k
  x=2*K+eps;y=2*eps*k*(K+eps)
  check(K==1+eps*f0,'boundary_k_square')
  check(z*z==K*(u*u+1)**2-4*u*(u*u-1),'boundary_on_quartic')
  check(y*y==(x-K)*(x*x-1),'boundary_on_cubic')
  check(x==(k*z+K+2*u+K*u*u)/(2*u*u),'boundary_forward_x')
  check(y==(x-K)*(x*u-1)/k,'boundary_forward_y')
  check(u==(k*y+x-K)/(x*(x-K)),'boundary_inverse_u')
  boundaries.append({'seed':T,'epsilon':eps,'point':[str(x),str(y)]})
cycles=[]
for eps,Fmap in [(1,fp),(-1,fm)]:
 for length in [1,3]:
  t0=1+F(1,64*4**length);ts=[t0]
  for j in range(length):
   old=ts[-1];new=Fmap(old)
   check(old<new<1+4*(old-1),'cusp_growth_probe')
   check(0<f(old)<f(new)<F(1,2),'strict_spread_growth')
   square_root(1+eps*f(old)-eps*f(new))
   ts.append(new)
  triples=[ppt(t) for t in ts]
  for old,new in zip(triples,triples[1:]):
   C,Cp=old[2],new[2]
   check(Cp%(C*C)==0 and gcd(C,Cp//(C*C))==1,'exact_hypotenuse_prime_transport')
  H=lcm(*(T[2] for T in triples));M=2*H*H
  ports=[]
  for a,b,c in triples:
   g=H//c;ports.append((2*g*g*(c*c-4*a*b),2*g*g*(c*c+4*a*b)))
  labels=[M]
  for low,high in ports:labels.extend((low,high) if eps==1 else (high,low))
  check(len(set(labels))==len(labels),'cycle_distinct_vertices')
  check(min(labels)>0 and max(labels)<4*H*H,'cycle_positive_bound')
  roots=[]
  for i,x in enumerate(labels):
   v=x+labels[(i+1)%len(labels)];r=isqrt(v)
   check(r*r==v,'integer_cycle_edge_square');roots.append(r)
  encoded=json.dumps({'labels':labels,'roots':roots},separators=(',',':'))
  rec={'epsilon':eps,'contact_count':length,'petal_count':len(triples),
       'integer_vertices':len(labels),'hypotenuse_bits':[T[2].bit_length() for T in triples],
       'common_scale_bits':H.bit_length(),'largest_label_bits':max(labels).bit_length(),
       'certificate_sha256':sha256(encoded.encode()).hexdigest()}
  if length==1:rec.update(labels=labels,integer_square_roots=roots)
  cycles.append(rec)
print(json.dumps({'status':'FINITE-EXACT controls; universal theorem proved separately',
 'checks':dict(checks),'uniform_chart':symbolic,'boundary_points':boundaries,'cycles':cycles},indent=2))
