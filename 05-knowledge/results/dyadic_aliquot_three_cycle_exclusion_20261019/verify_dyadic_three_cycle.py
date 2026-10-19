#!/usr/bin/env python3
import json
from fractions import Fraction
import sympy as s
checks=0
def need(b,msg):
 global checks
 checks+=1
 if not b: raise RuntimeError(msg)
c,d,e=s.symbols('c d e');p,q,r=s.symbols('p q r')
D=c*d+c*e+d*e-c-d-e+1
Nf=2*c*c*d-2*c*c+2*c*d*d+2*c*d*e-3*c*d-2*c*e+2*c-2*d*e+d+2*e-1
Nr=Nf+(c-d)*(2*(c+d)-1)
sol=s.solve([d*q-(c-1)*p-2*c+1,e*r-(d-1)*q-2*d+1,c*p-(e-1)*r-2*e+1],[p,q,r])
need(s.cancel(sol[r]-Nf/D)==0,'forward elimination')
need(s.cancel(sol[r].xreplace({c:d,d:c})-Nr/D)==0,'reverse elimination')
L=s.expand(Nr-(c-1)*D);U=s.expand((3*c-1)*D-Nf)
ids=[L.coeff(e)-(c-1)*(d-c-1),L.subs(e,2*d)-(d*(4*(c-1)*d-c*(c+1)+3)+c*(c-1)),U.coeff(e)-(3*c*c+c*d-2*c+d-1),U.subs(e,2*d)-(c*((7*c-5)*d-c+2)+2*d*(d-1)),Nf-(2*c-1)*D-(c*(2*d*d-1)-e*(2*c*c-c+d-1)),Nr-(2*c-1)*D-(2*(c-1)*(c+d*d)+d-e*(2*c*c-c+d-1))]
for expr in ids:need(s.expand(expr)==0,'polynomial identity')
w,u,v=s.symbols('w u v',nonnegative=True)
for expr in [L,U]:
 poly=s.Poly(expr.subs(e,2*d+v).subs(d,2*c+u).subs(c,2+w).expand(),w,u,v)
 need(all(co>0 for _,co in poly.terms()),'positive orthant certificate')
need(s.expand((2*(c-1)*(c+d*d)+d).subs(d,2*c)-2*c*c*(4*c-3))==0,'exceptional numerator')
need(s.discriminant(2*c*c-3*c+2,c)==-7,'last impossible equation')

def affine(cores):
 # Compose source-variable maps n -> (1-1/core)n+2core-1.
 slope=Fraction(1);carry=Fraction(0)
 for z in cores:
  slope=Fraction(z-1,z)*slope
  carry=Fraction(z-1,z)*carry+2*z-1
 n0=carry/(1-slope)
 values=[];n=n0
 for z in cores:
  values.append(n/z)
  n=Fraction(z-1,z)*n+2*z-1
 need(n==n0,'source closure')
 return values
words=0
for a in range(1,17):
 for b in range(a+1,25):
  for g in range(b+1,33):
   C,E,F=2**a,2**b,2**g
   for cores in [(C,E,F),(E,C,F)]:
    vals=affine(cores);rr=vals[2]
    need(C-1<rr<3*C-1,'all-size interval control')
    need(rr.denominator!=1,'no integral largest-core variable')
    words+=1
# Non-native repeated-core exception and genuine fixed controls.
need(affine((2,2,16))==[Fraction(23),Fraction(13),Fraction(1)],'ghost retained')
for a in range(1,9):
 C=2**a;need(affine((C,C,C))==[Fraction(2*C-1)]*3,'equal-core fixed template')
print(json.dumps({'status':'PASS','checks':checks,'distinct_oriented_core_words':words,'exponent_control_ranges':'1<=a<=16, a<b<=24, b<g<=32; both orientations','method':'symbolic identities/positive-coefficient proofs plus independent rational source-map composition','scope':'finite controls supplement the all-exponent written proof; no prime search'},indent=2,sort_keys=True))
