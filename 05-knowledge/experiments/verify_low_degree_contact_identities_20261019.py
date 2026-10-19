#!/usr/bin/env python3
"""Elementary low-degree contact identities; Python standard library only."""
from fractions import Fraction as F
from collections import Counter
import json
from rational_contact_certificate_core import Poly,variables
checks=Counter()
def check(b,name):
    checks[name]+=1
    if not b:raise RuntimeError(name)
class Rat:
    def __init__(self,n,d=1):
        self.n=n if isinstance(n,Poly) else Poly.scalar(1,n)
        self.d=d if isinstance(d,Poly) else Poly.scalar(1,d)
        if not self.d:raise ZeroDivisionError()
    def _c(self,x):return x if isinstance(x,Rat) else Rat(x)
    def __add__(self,x):
        x=self._c(x);return Rat(self.n*x.d+x.n*self.d,self.d*x.d)
    __radd__=__add__
    def __neg__(self):return Rat(-self.n,self.d)
    def __sub__(self,x):return self+-self._c(x)
    def __rsub__(self,x):return self._c(x)+-self
    def __mul__(self,x):
        x=self._c(x);return Rat(self.n*x.n,self.d*x.d)
    __rmul__=__mul__
    def __truediv__(self,x):
        x=self._c(x);return Rat(self.n*x.d,self.d*x.n)
    def __rtruediv__(self,x):return self._c(x)/self
    def __pow__(self,k):return Rat(self.n**k,self.d**k)
    def __eq__(self,x):
        x=self._c(x);return self.n*x.d==x.n*self.d
    def derivative(self):
        def der(p):return Poly(1,{(m[0]-1,):c*m[0] for m,c in p.terms.items() if m[0]})
        return Rat(der(self.n)*self.d-self.n*der(self.d),self.d*self.d)
t=Rat(variables(1)[0]);f=lambda t:4*t*(t*t-1)/(t*t+1)**2
for eps in [1,-1]:
    k=(t*t+2*eps*t-1)/(t*t+1)
    check(k*k==1+eps*f(t),'contact_square_chart')
    v=variables(1)[0]
    check(4-4*(1-v)*(-1-v)==4*(2-v*v),'branch_discriminant')
u,c=variables(2)
p=(2*c*u+c*c+1)**2-8*(u+c)**2
r=Poly(2)
for mon,coef in p.terms.items():r+=Poly(2,{(mon[0]%2,mon[1]):coef*2**(mon[0]//2)})
check(r==(c*c-3)*(4*c*u+c*c+5),'branch_permutation_remainder')
x=(t*t-1)/(t*t+1);z=2*t/(t*t+1);g=x*x-z*z
check(x*x+z*z==1,'circle_identity')
check(f(t)==2*x*z,'spread_harmonic')
check(f(t)**2+g*g==1,'spread_critical_values')
check(f(t).derivative()==-4*(t**4-6*t*t+1)/(t*t+1)**3,'spread_derivative')
r,s=variables(2)
check((r*s)**2+(r*r-s*s)**2/4==(r*r+s*s)**2/4,'harmonic_norm_identity')
for h,sign in [(t,1),(-t,-1),(1/t,-1),(-1/t,1),((t-1)/(t+1),-1),((1-t)/(t+1),1),((t+1)/(t-1),1),((-t-1)/(t-1),-1)]:
    check(f(h)==sign*f(t),'trivial_degree_one_symmetry')
for a,b in [(1,0),(0,1),(1,1),(2,1),(3,2),(5,3)]:
    U=F(a*a-b*b,a*a+b*b);V=F(2*a*b,a*a+b*b)
    for sign in [1,-1]:
        h=sign*(a*t+b)/(-b*t+a);A=sign*(U*U-V*V);B=-2*sign*U*V
        check(f(h)==A*f(t)+B*g,'rational_orthogonal_harmonic_identity')
        check(A*A+B*B==1,'unit_harmonic_coefficients')
print(json.dumps({'status':'FINITE-EXACT elementary identity checks; standard library only','checks':dict(checks)},indent=2))
