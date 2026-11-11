#!/usr/bin/env python3
"""Exact finite controls for new intertwiners and rational support theorem."""
import sympy as s
x,t=s.symbols('x t')
a=x*x+6
count=0
def require(ok):
    global count
    if not ok:
        raise AssertionError('Extension control failed')
    count+=1
def L(m,f): return s.expand(-x*f+a*s.diff(f,x)/(2*m))
def A(m,f): return s.expand(a*s.diff(f,x,2)-2*(2*m-1)*x*s.diff(f,x)+2*m*(2*m-1)*f)
def mu(m,j):
    if j%2:return s.Integer(0)
    return s.prod(s.Rational(6*(2*i-1),2*m-2*i+1) for i in range(1,j//2+1))
def ell(m,f): return sum(c*mu(m,j) for (j,),c in s.Poly(f,x).terms())
def weight(p):
    if p==0:return -s.oo
    return max(2*j-r for (r,j),c in s.Poly(p,t,x).terms())
for m in range(1,12):
    for j in range(2*m+1):
        f=x**j
        require(s.expand(L(m+1,a*f)-s.Rational(m,m+1)*a*L(m,f))==0)
        require(s.simplify(ell(m+1,a*f)-s.Rational(12*(m+1),2*m+1)*ell(m,f))==0)
        if m>=2:
            require(s.expand(A(m,L(m,f))-s.Rational(m-1,m)*L(m-1,A(m,f)))==0)
# Lawful rational carriers and their exact finite coefficients.
for P,Q,B in [(1,1-t*t*x,0),(t*x**3,1-t-t*t*x,5),
              (t**5*(1+t*t*x),1-t-t**4*x*x,-5)]:
    require(weight(P)<=B and weight(Q)<=0)
    f=s.series(P/Q,t,0,25).removeO().expand()
    require(weight(f)<=B)
    require(weight(s.diff(f,x))<=B-2)
    require(weight(s.diff(f,t))<=B+1)
# Positive-weight denominator eventually exits every bounded cone.
require(weight(1-t*x)==1)
f=s.series(1/(1-t*x),t,0,13).removeO().expand()
require(weight(f)==12)
# Rationality alone and finite passing rows are insufficient.
N=12
f=s.series(1/(1-t**N*x**N),t,0,N+1).removeO().expand()
require(f==1+t**N*x**N)
require(weight(f)==N)
# Native depth flag does not survive radial raising.
d16=(16-1)//2+2; d17=(17-1)//2+2
require(d16==9 and d17==10)
require(s.degree(a*x**d16,x)>d17)
print(f'PASS: {count} exact row-transport and rational-cone controls.')
