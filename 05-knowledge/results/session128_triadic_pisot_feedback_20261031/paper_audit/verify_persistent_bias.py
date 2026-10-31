#!/usr/bin/env python3
"""Exact-rational interval certificate for an Amanda persistent-parity witness.
No floating-point calculation is used in any pass/fail gate. Python 3 stdlib.
The proof of the harmonic zero-count lemma is in PAPER_AUDIT.md.
"""
from fractions import Fraction as F
from math import factorial, isqrt
import json
from pathlib import Path

D = 10**450

def floorf(x): return x.numerator // x.denominator

def ceilf(x): return -floorf(-x)

class I:
    def __init__(self, lo, hi=None):
        lo = F(lo); hi = lo if hi is None else F(hi)
        if lo > hi: raise ValueError('reversed interval')
        self.lo = F(floorf(lo*D), D)
        self.hi = F(ceilf(hi*D), D)
    def __add__(self, other):
        other=as_i(other); return I(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self, other): return self+-as_i(other)
    def __rsub__(self, other): return as_i(other)+-self
    def __mul__(self, other):
        other=as_i(other)
        values=[self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi]
        return I(min(values), max(values))
    __rmul__=__mul__
    def reciprocal(self):
        if self.lo<=0<=self.hi: raise ZeroDivisionError('interval includes zero')
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,other): return self*as_i(other).reciprocal()
    def __rtruediv__(self,other): return as_i(other)*self.reciprocal()
    def __pow__(self,n):
        if n<0:return self.reciprocal()**(-n)
        ans=I(1)
        for _ in range(n):ans=ans*self
        return ans
    def square(self):
        if self.lo<=0<=self.hi:return I(0,max(self.lo*self.lo,self.hi*self.hi))
        return I(min(self.lo*self.lo,self.hi*self.hi),max(self.lo*self.lo,self.hi*self.hi))
    def widen(self,e):return I(self.lo-e,self.hi+e)
    def bracket(self, digits=22):
        d=10**digits
        return [decimal_integer(floorf(self.lo*d),digits),decimal_integer(ceilf(self.hi*d),digits)]

def decimal_integer(n,digits):
    sign='-' if n<0 else '';n=abs(n);return sign+str(n//(10**digits))+'.'+str(n%(10**digits)).zfill(digits)

def as_i(x):return x if isinstance(x,I) else I(x)

def atan_reciprocal(q,N=80):
    a=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(N)),F(0))
    b=a+F((-1)**N,(2*N+1)*q**(2*N+1))
    return I(min(a,b),max(a,b))

def trig(x,kind,N=105):
    # Horner evaluation of the Taylor polynomial, plus a Lagrange remainder.
    x=as_i(x);y=x.square();offset=1 if kind=='sin' else 0
    z=I(F((-1)**N,factorial(2*N+offset)))
    for k in range(N-1,-1,-1):z=z*y+F((-1)**k,factorial(2*k+offset))
    if offset:z=x*z
    degree=2*N+offset
    a=max(abs(x.lo),abs(x.hi))
    remainder=a**(degree+1)/factorial(degree+1)
    return z.widen(remainder)

checks=[]
def require(name,condition):
    if not condition:raise RuntimeError('FAILED: '+name)
    checks.append(name)

# Machin's formula is exact: pi = 16 atan(1/5) - 4 atan(1/239).
pi=16*atan_reciprocal(5)-4*atan_reciprocal(239)
root=isqrt(5*D*D)
require('sqrt5_integer_bracket',root*root<5*D*D<(root+1)*(root+1))
sqrt5=I(F(root,D),F(root+1,D));phi=(1+sqrt5)/2
golden=2*pi/(phi*phi)
h=[trig(j*golden,'sin') for j in (1,2)]
require('h1_above_one_third',h[0].lo>F(1,3))
require('h2_below_minus_one_third',h[1].hi<-F(1,3))
# Canonical normalized golden weights w_k=phi^(-k-2).
A=[-4/(phi**2),9/(phi**3)]
u=sum((4*a*trig(3*j*golden,'cos') for j,a in enumerate(A,1)),I(0))
v=sum((-4*a*trig(3*j*golden,'sin') for j,a in enumerate(A,1)),I(0))
w=sum((-a*trig(6*j*golden,'cos') for j,a in enumerate(A,1)),I(0))
z=sum((a*trig(6*j*golden,'sin') for j,a in enumerate(A,1)),I(0))
R2=u.square()+v.square();S2=w.square()+z.square()
require('two_simple_zeros_harmonic_dominance',(R2-5*S2).lo>0)

def g(x):
    x=I(F(x));return u*trig(x,'sin')+v*trig(x,'cos')+w*trig(2*x,'sin')+z*trig(2*x,'cos')

r1=I(F('3.0450368193027045'),F('3.0450368193027046'))
r2=I(F('5.5300956371054050'),F('5.5300956371054052'))
require('root_brackets_ordered_inside_period',0<r1.lo<r1.hi<r2.lo<r2.hi<2*pi.lo)
values={
  'zero':g('0'),
  'r1_left':g(r1.lo), 'r1_right':g(r1.hi),
  'r2_left':g(r2.lo), 'r2_right':g(r2.hi),
}
require('g0_negative',values['zero'].hi<0)
require('r1_negative_to_positive',values['r1_left'].hi<0<values['r1_right'].lo)
require('r2_positive_to_negative',values['r2_right'].hi<0<values['r2_left'].lo)
# The two-zero lemma plus signs implies positivity exactly between r1 and r2.
parity=(r2-r1)/pi-1
require('persistent_parity_strict_negative',parity.hi<F('-0.20898'))
require('parity_numeric_enclosure',F('-0.208982')<parity.lo<=parity.hi<F('-0.208980'))

# Independent exact Borwein m=7 value verification from simplex-tail probability.
a=[F(1,2*k+1) for k in range(1,8)];delta=sum(a,F(0))-1
prod=F(1)
for x in a:prod*=x
ratio=1-delta**7/(2**6*factorial(7)*prod)
paper_coefficient=F(467807924713440738696537864469,935615849440640907310521750000)
require('borwein_delta_small_simplex',0<delta<2*min(a))
require('borwein_fraction_correct',ratio==2*paper_coefficient)

out={
  'status':'verified','arithmetic':'exact Fraction intervals with outward rounding at 450 decimal places',
  'trigonometry':'degree 210/211 Taylor-Horner polynomials plus rigorous Lagrange remainder',
  'pi':'Machin formula; alternating rational arctan brackets with 80 terms',
  'fixed_schedule':[1,2],'signs':[1,-1], 'backbone_initial_triple':[3,2,1], 'weights':'w_k=phi^(-k-2)',
  'h1':h[0].bracket(),'h2':h[1].bracket(),
  'fundamental_amplitude_squared':R2.bracket(),'five_times_second_amplitude_squared':(5*S2).bracket(),
  'harmonic_gap':(R2-5*S2).bracket(),
  'r1':r1.bracket(),'r2':r2.bracket(),'parity':parity.bracket(),
  'endpoint_values':{k:val.bracket(30) for k,val in values.items()},
  'checks':checks,'native_memory':'b_t=Pi*(1-(1-lambda)^t), for every 0<lambda<=1',
  'conclusion':'Under the explicit native update equations, parity is nonzero forever; signs freeze because |b_t|<=1/3 while h1>1/3 and h2<-1/3. An additional hard acceptance test could reject this trajectory but is not a correction/convergence mechanism.'
}
print(json.dumps(out,indent=2,sort_keys=True))
