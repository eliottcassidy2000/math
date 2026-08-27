#!/usr/bin/env python3
"""Symbolic denominator probe for (1-omega)f+Tf."""
import sympy as sp

# The common target scaling may be set to one for pole calculations: the
# coefficient-form elliptic addition law is homogeneous of weights (2,3).
A,B,C,D,cc,ss=sp.symbols('A B C D cc ss')
X_scaled=cc*((ss**3*(D-B))/(ss**2*(C-A)))**2-ss**2*A-ss**2*C
X_unit=cc*((D-B)/(C-A))**2-A-C
Y_scaled=(ss**3*(D-B)/(ss**2*(C-A)))*(ss**2*A-X_scaled)-ss**3*B
Y_unit=((D-B)/(C-A))*(A-X_unit)-B
assert sp.cancel(X_scaled-ss**2*X_unit)==0
assert sp.cancel(Y_scaled-ss**3*Y_unit)==0

p,z,t=sp.symbols('p z t')
u=t*t
phi=z**4-z**2+1
omega=sp.rem(sp.Poly(z**4,z),sp.Poly(phi,z)).as_expr()
af=(u-p**2)/(2*t)
bf=(u+p**3)/(2*t)
ag=z**2*(p**2*u-1)/(2*t)
bg=-z**3*(1+p**3*u)/(2*t)
c=(u-1)/(2*t)

def add_x(P,Q):
    a,b=P; aa,bb=Q
    return sp.factor(c*(bb-b)**2-(a+aa)*(aa-a)**2), sp.factor((aa-a)**2)

# First form r=(1-omega)f=f+(-omega f).  Keep X as numerator/denominator.
nr,dr=add_x((af,bf),(omega*af,-bf))
# Y_r = slope*(X_f-X_r)-Y_f with common denominators retained.
slope_num=-2*bf
slope_den=(omega-1)*af
yr=sp.together((slope_num/slope_den)*(af-nr/dr)-bf)
xr=sp.together(nr/dr)

# Add r to g and inspect only the X coordinate.
x=sp.factor(c*(bg-yr)**2-(xr+ag)*(ag-xr)**2)
d=sp.factor((ag-xr)**2)
x_total=sp.cancel(x/d)
raw_num,raw_den=sp.fraction(x_total)
den=sp.Poly(raw_den,t)

sqrt3=2*z-z**3
rels=sp.groebner([p**2-(1+sqrt3)*p+1,phi],p,z,order='lex')
red=sum(rels.reduce(sp.expand(den.coeff_monomial(t**k)))[1]*t**k
        for k in range(den.degree()+1))
print('raw_degree',den.degree())
print('target_scale_homogeneity', 'A:s^2,B:s^3 PASS')
print('reduced_degree',sp.Poly(red,t).degree())
print('factor',sp.factor(red))
lc=sp.factor(sp.Poly(red,t).LC())
print('leading_norm',sp.factor(sp.resultant(sp.resultant(lc,p**2-(1+sqrt3)*p+1,p),phi,z)))
for value in (0,1,-1):
    reduced_value=sp.factor(rels.reduce(sp.expand(raw_num.subs(t,value)))[1])
    value_norm=sp.factor(sp.resultant(
        sp.resultant(reduced_value,p**2-(1+sqrt3)*p+1,p),phi,z))
    print('numerator_at',value,'=',reduced_value,'norm=',value_norm)
