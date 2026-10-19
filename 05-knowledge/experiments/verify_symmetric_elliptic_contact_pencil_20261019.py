#!/usr/bin/env python3
"""Standard-library exact controls for the symmetric elliptic contact pencil."""
from rational_contact_certificate_core import variables
from fractions import Fraction as F
from collections import Counter
import json
C=Counter()
def check(x,name):
    C[name]+=1
    if not x:raise RuntimeError(name)
t,u,l=variables(3);s=t+u;p=t*u;D=s*s+(p-1)**2;J=p*p+6*p+1-s*s
W=D*D+4*((t**3-t)*(u*u+1)**2+(u**3-u)*(t*t+1)**2)
H=4*J-2*l*D-l*l*s*(p-1);Y=D+l*s*(p-1)
check(W==D*D+4*s*(p-1)*J,'spread_sum_factorization')
check(W-Y*Y==s*(p-1)*H,'pencil_contact_identity')
A=t+2;B=t*t-4*t-1;G=2*t*t-t;eq=A*u*u+B*u+G
check(H.substitute({2:2})==-4*eq,'lambda_two_member')
check(eq==(u+2)*t*t+(u*u-4*u-1)*t+2*u*u-u,'coordinate_swap_symmetry')
quartic=t**4-16*t**3+2*t*t+16*t+1
check(B*B-4*A*G==quartic,'quartic_discriminant')
# Vieta root t' = u(2u-1)/(t(u+2)) in the t-fiber over u.
a=u+2;b=u*u-4*u-1;c=2*u*u-u
check(c*c+b*c*t+a*c*t*t==c*(a*t*t+b*t+c),'vieta_preserves_curve_cleared')
def curve(a,b):return (a+2)*b*b+(a*a-4*a-1)*b+2*a*a-a
def spread(a):return 4*a*(a*a-1)/(a*a+1)**2
orbit=[];a=b=F(1)
for i in range(5):
    check(curve(a,b)==0,'rational_orbit_on_curve')
    y=1+2*(a+b)*(a*b-1)/((1+a*a)*(1+b*b))
    check(1+spread(a)+spread(b)==y*y,'rational_orbit_contact')
    orbit.append((a,b))
    a,b=b,b*(2*b-1)/(a*(b+2))
check(len(set(orbit))==5,'five_distinct_orbit_points')
def trim(a,p):
    a=[x%p for x in a]
    while a and not a[-1]:a.pop()
    return a
def remainder(a,b,p):
    a=trim(a,p);b=trim(b,p)
    while len(a)>=len(b) and a:
        z=a[-1]*pow(b[-1],-1,p)%p;k=len(a)-len(b)
        for j,v in enumerate(b):a[j+k]=(a[j+k]-z*v)%p
        a=trim(a,p)
    return a
counts={}
for p,expected in [(7,12),(13,20)]:
    a=[1,16,2,-16,1];b=[16,4,-48,4]
    while b:a,b=b,remainder(a,b,p)
    check(len(trim(a,p))==1,'good_reduction_squarefree')
    n=2
    for x in range(p):
        z=(x**4-16*x**3+2*x*x+16*x+1)%p
        n+=sum(y*y%p==z for y in range(p))
    check(n==expected,'finite_field_point_count');counts[str(p)]=n
check(spread(F(1,3))==F(-24,25),'explicit_inadmissible_orbit_value')
# Biquadratic model and the exact second-sign branch obstruction.
a,l=variables(2)
Lp=2*(l+2)*(1+a*a)+l*l*a;Lm=2*(2-l)*(1+a*a)-l*l*a
check(Lp+Lm==8*(1+a*a),'biquadratic_circle_relation')
check(l**4-16*(l+2)**2==(l*l-4*l-8)*(l*l+4*l+8),'first_branch_discriminant')
check(l**4-16*(2-l)**2==(l*l+4*l-8)*(l*l-4*l+8),'second_branch_discriminant')
B=(1+a*a)**2-2*l*a*(1+a*a)-l*l*a*a
check(2*(1+a*a)**2-(1+a*a+l*a)**2==B,'second_contact_numerator')
for c,d,L,sign in [(2*(l+2),l*l,Lp,-1),(2*(2-l),-l*l,Lm,1)]:
    check(c*c*B-L*(L-2*(d+c*l)*a)==l*l*(l*l+sign*8*l-16)*a*a,'second_contact_branch_exclusion')
# The pullbacks of (8) reduce exactly to the original H relation.
t,u,l=variables(3);p=t*u;s=t+u;den=(p-1)**2
H=4*(p*p+6*p+1-s*s)-2*l*(s*s+den)-l*l*s*(p-1)
check(8*(p+1)**2-2*(l+2)*(den+s*s)-l*l*s*(p-1)==H,'first_biquadratic_model')
check(8*(t-u)**2-2*(2-l)*(den+s*s)+l*l*s*(p-1)==-H,'second_biquadratic_model')
print(json.dumps({'status':'FINITE-EXACT elliptic contact pencil; standard library only','checks':dict(C),'point_counts':counts,'orbit':[[str(x),str(y)] for x,y in orbit]},indent=2))
