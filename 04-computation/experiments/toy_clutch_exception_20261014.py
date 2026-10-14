#!/usr/bin/env python3
"""Exact certificate for the critical-free toy clutch and its Laurent obstruction."""
import json
import sympy as s

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

def zero(poly, message):
    check(s.expand(poly) == 0, message)

x,z,t,y,b,e,B,J,A,V = s.symbols('x z t y b e B J A V')
R1=2*y*(y+B)*(2*y+B)+V*A
R2=V**3+V**2*y+J*y*y*(y+B)
K=(-A**3*J**3+12*A**2*J**2*V**2-4*A*B**3*J**2*V
   +4*A*B**2*J*V**2+24*A*B*J*V**3-48*A*J*V**4-16*A*V**4
   -8*B**4*J*V**2+8*B**3*J*V**3+32*B**2*V**4-96*B*V**5+64*V**6)
zero(s.resultant(R1,R2,y)-V**3*K,'universal resultant')
d=s.symbols('d')
F=(B**4*(1-2*d)+x*(4*B**3*d*d-2*B**3*d+4*B**2*d-4*B**2)
   +x*x*(-12*B*d*d+12*B*d+4)+8*x**3*(d-1)**3)
zero(K.subs({V:x,A:2*x,J:2*x*d-B})+8*x**3*F,'toy factorization')
Fhalf=s.factor(F.subs({B:b+x/2,d:s.Rational(1,2)}))
zero(Fhalf-x*(-2*b*b+(b+4)*x),'half-slope factor')
Kstar=s.factor((-8*x**3*Fhalf).subs(b,-4))
check(Kstar==256*x**4,'exception resultant')
P=(x*z*z+(b+x/2)*z)**2+2*x*z+x
u=x*t+2*b
Pt=-2*b+u+u*u*(t-1)**2/16
zero(P.subs(z,(t-1)/2)-Pt,'normal form')
Q=1+u*(t-1)**2/8
zero(s.diff(Pt,x)-t*Q,'x derivative')
zero(s.diff(Pt,t)-x*Q-u*u*(t-1)/8,'t derivative')
zero(s.diff(Pt,t).subs(t,0)-(b+4)*x/4+b*b/2,'boundary derivative')
critical={x:2*b*b/(b+4),z:-s.Rational(1,2)}
check(s.factor(s.diff(P,x).subs(critical))==0,'critical x derivative')
check(s.factor(s.diff(P,z).subs(critical))==0,'critical z derivative')
check(s.factor(s.det(s.hessian(P,(x,z))).subs(critical))==-(b+4)**2/4,'Morse determinant')
Faff=s.Poly(s.expand(F.subs({B:b+e*x,d:e})),x)
zero(Faff.coeff_monomial(x**4)-e**4*(2*e-1),'affine leading coefficient')
zero(Faff.coeff_monomial(1)-b**4*(1-2*e),'affine constant coefficient')
check(s.degree(Faff.as_expr().subs(e,0),x)==3,'constant-clutch degree')
Pstar=s.expand(P.subs(b,-4))
checks=0
for i in range(9):
    for j in range(9-i):
        H=x**i*z**j
        jac=s.diff(Pstar,x)*s.diff(H,z)-s.diff(Pstar,z)*s.diff(H,x)
        q=H.subs({x:8/t,z:(t-1)/2})
        zero(jac.subs({x:8/t,z:(t-1)/2})-2*t*s.diff(q,t),'Laurent restriction')
        checks+=1
# Degree estimates from the inherited accessory data deg(V)=16, deg(A)=8.
# n=8 has extra cancellation in J; the listed upper bound still suffices.
def terms(n):
    j=n+15
    return {'-A3J3':24+3*j,'12A2J2V2':16+2*j+32,
            '-4AB3J2V':8+3*n+2*j+16,'4AB2JV2':8+2*n+j+32,
            '24ABJV3':8+n+j+48,'-48AJV4':8+j+64,
            '-16AV4':72,'-8B4JV2':4*n+j+32,'8B3JV3':3*n+j+48,
            '32B2V4':2*n+64,'-96BV5':n+80,'64V6':96}
for n in range(9):
    ds=terms(n)
    check(ds['64V6']>max(v for k,v in ds.items() if k!='64V6'),'low-degree envelope')
for n in range(10,41):
    ds=terms(n)
    check(ds['-4AB3J2V']==5*n+54,'high-degree formula')
    check(ds['-4AB3J2V']>max(v for k,v in ds.items() if k!='-4AB3J2V'),'high-degree envelope')
a,v,beta=s.symbols('a v beta')
leading=-4*a*beta**3*(2*v*beta)**2*v+8*beta**3*(2*v*beta)*v**3
zero(leading-16*beta**4*v**3*(v-a*beta),'degree-nine resonance')
print(json.dumps({
    'status':'PASS','universal_resultant_verified':True,
    'exception_clutch':'B=x/2-4','exception_K':str(Kstar),
    'half_slope_F':str(Fhalf),'critical_free_exception':True,
    'critical_point_if_b_not_minus_4':['2*b^2/(b+4)','-1/2'],
    'hessian_determinant':'-(b+4)^2/4',
    'jacobian_mate_obstruction':'On xt=8, Jac(P,H)=2*t*d(H|C)/dt; constant Laurent coefficient is zero.',
    'laurent_monomials_verified':checks,
    'accessory_degrees':{'n<=8':96,'n=9':'99 unless beta=v/a; resonance degree unsettled','n>=10':'5*n+54'},
    'original_fixed_accessory_HJ':'OPEN; not refuted by the toy exception'
},indent=2,sort_keys=True))
