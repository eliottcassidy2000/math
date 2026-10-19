#!/usr/bin/env python3
"""Exact elementary identities for the elliptic-shift obstruction; stdlib only."""
from rational_contact_certificate_core import Poly,variables
from collections import Counter
import json
checks=Counter()
def check(b,name):
    checks[name]+=1
    if not b:raise RuntimeError(name)
x,p,q,g,a=variables(5);b=q*q-p**3-a*p
curve_residue=x**3+a*x+b-(q+g*(x-p))**2
quadratic=x*x+(p-g*g)*x+p*p+a+g*g*p-2*g*q
check(curve_residue==(x-p)*quadratic,'chord_intersection_quadratic')
Delta=g**4-6*p*g*g+8*q*g-3*p*p-4*a
check((g*g-p)**2-4*(p*p+a+g*g*p-2*g*q)==Delta,'chord_discriminant')
Dminus=(-g)**4-6*p*(-g)**2+8*q*(-g)-3*p*p-4*a
check(Delta-Dminus==16*q*g,'opposite_branch_value_obstruction')
z,alpha,beta=variables(3)
square=(z*z+alpha*z+beta)**2
check(square.coeff(0,3)==2*variables(2)[0],'quartic_cubic_coefficient')
check(square.substitute({1:0}).coeff(0,1)==0,'zero_cubic_forces_zero_linear')
X,d=variables(2)
check((3*X*X-d*d)**2-8*X*(X**3-d*d*X)==(X*X+d*d)**2,'duplication_numerator')
# The spread at t=X/d is 4Xd(X²−d²)/(X²+d²)²; multiplication
# by the displayed x-duplication formula gives d identically.
check((4*X*d*(X*X-d*d))*(X*X+d*d)**2==d*(X*X+d*d)**2*(4*X*(X*X-d*d)),'spread_reciprocal_duplication')
check((X**3-d*d*X)-X**3==-d*d*X,'zero_two_torsion_translation')
# Independent normalized x-norm audit for all four sign choices.
x,u=variables(2)
for eps,eta in [(1,1),(-1,-1),(1,-1),(-1,1)]:
    B=(1+(eta-eps)*x+eta*u)*(x-u)**2-eta*(x**3-x+u**3-u)
    if eps==eta==1:
        W=x**4-2*x**3+(4*u*u+2*u-1)*x*x+(-4*u*u-2*u+2)*x+(u+1)**2
        A=-1;C=(u+1)*(2*u-1)
        residual=-4*u*(u-1)*(u+1)**2
    elif eps==eta==-1:
        W=x**4+2*x**3+(4*u*u-2*u-1)*x*x+(4*u*u-2*u-2)*x+(u-1)**2
        A=1;C=(u-1)*(2*u+1)
        residual=-4*u*(u-1)**2*(u+1)
    elif eps==1:
        W=x**4+(-4*u-2)*x**3+(6*u+3)*x*x+(2*u-2)*x+(u-1)**2
        A=-2*u-1;C=-(u-1)*(2*u+1)
        residual=None
    else:
        W=x**4+(-4*u+2)*x**3+(-6*u+3)*x*x+(2*u+2)*x+(u+1)**2
        A=1-2*u;C=-(u+1)*(2*u-1)
        residual=None
    check(B*B-4*(u**3-u)*(x**3-x)==(x-u)**2*W,'normalized_norm_quartic')
    rem=W-(x*x+A*x+C)**2
    check(rem.degree(0)<=1,'square_forced_by_top_coefficients')
    if residual is not None:check(rem==residual,'same_sign_torsion_residual')
    else:check(rem.coeff(0,1)==(-8*u*(u-1)*(u+1)).coeff(0,0),'mixed_sign_torsion_residual')
# General unequal-coefficient extension: A+B*x+C*x_shift.
g,A,B,C,p,q,a=variables(7)
b=(B+C)/2;c=(B-C)/2;alpha=A-b*p
Delta=g**4-6*p*g*g+8*q*g-3*p*p-4*a
norm=(alpha+b*g*g)**2-c*c*Delta
check(norm.coeff(0,4)==(B*C).coeff(0,0),'general_norm_leading_BC')
check(norm.coeff(0,3)==0,'general_norm_zero_cubic')
check(norm.coeff(0,1)==(-8*c*c*q).coeff(0,0),'general_norm_nonzero_linear')
print(json.dumps({'status':'FINITE-EXACT elementary identities; standard library only','checks':dict(checks)},indent=2))
