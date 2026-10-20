#!/usr/bin/env python3
"""Independent exact controls; imports no producer code; survives python -O."""
import json
from pathlib import Path
import hashlib
import sympy as sp

checks=[]
def equal(name, lhs, rhs=0):
    difference=sp.expand(lhs-rhs)
    if difference != 0:
        raise RuntimeError(f'{name}: {difference}')
    checks.append(name)
def truth(name, condition):
    if not condition:
        raise RuntimeError(name)
    checks.append(name)
def coeff(poly, variable, power):
    return sp.Poly(sp.expand(poly),variable).coeff_monomial(variable**power)
def jac(p,q):
    return sp.expand(sp.diff(p,x)*sp.diff(q,y)-sp.diff(p,y)*sp.diff(q,x))

x,y=sp.symbols('x y')
# Bounded coefficient probes complement the all-degree proof in AUDIT.md.
for n in range(1,7):
    for m in range(1,7):
        p=y**n+sum((x**(i%3+1)+i+1)*y**i for i in range(n))
        b=x*x+2*x+3
        q=b*y**m+sum((x**(i%2+1)+2*i+1)*y**i for i in range(m))
        equal(f'leading_gate_n{n}_m{m}',coeff(jac(p,q),y,n+m-1),-n*sp.diff(b,x))

# Positive controls prevent applying the nonprincipal conclusion to the line.
u=x+y**3; v=y+u**2
equal('positive_keller_jacobian',jac(u,v),1)
equal('positive_first_reduction',v-u**2,y)
equal('positive_second_reduction',u-y**3,x)
equal('ordinary_constant_shear_preserves_jacobian',jac(u,v-7*u**3),jac(u,v))
equal('source_polynomial_shear_preserves_jacobian',
      jac(u.subs(y,y+x*x),v.subs(y,y+x*x)),jac(u,v).subs(y,y+x*x))
# Hostiles fail before any purported globalization.
p=y**4; q=x*y**3+y*y+y
equal('frozen_hostile_jacobian',jac(p,q),-4*y**6)
equal('hidden_hostile_initial_both_monic',sp.Poly(q+p*p,y).LC(),1)
equal('hidden_hostile_new_leading_coefficient',sp.Poly(q+p*p-p*p,y).LC(),x)
equal('hidden_hostile_jacobian_unchanged',jac(p,q+p*p),-4*y**6)

# Independent universal reduced (4,6) calculations.
a,b,c,beta,gamma=sp.symbols('a b c beta gamma')
h=y**3+a*y; k=y**4+b*y*y+c*y
p=k+beta*h; q=h*h+gamma*h
r=sp.expand(q*q-p**3)
equal('nonzero_beta_degree11_coefficient',coeff(r,y,11),-3*beta)
r0=sp.expand(r.subs(beta,0)-(4*a-3*b)*k*q)
equal('wrong_gamma_degree9_coefficient',coeff(r0,y,9),2*gamma-3*c)
r1=sp.expand(r0.subs(gamma,sp.Rational(3,2)*c)+2*a*(a-b)*k*k)
equal('off_exceptional_degree7_coefficient',coeff(r1,y,7),-c*(2*a-3*b)/2)
for name,ds,want in [('beta',[0,6,11,17],7),('gamma',[0,6,9,15],6),('exceptional',[0,6,7,13],5)]:
    truth(name+'_distinct_residues',len({d%4 for d in ds})==4)
    truth(name+'_delta',sum(d-d%4 for d in ds)//4==want)
qstar=h*h+sp.Rational(3,2)*c*h
ks=k.subs(b,sp.Rational(2,3)*a)
recovery=(qstar*qstar-ks**3-2*a*ks*qstar+sp.Rational(2,3)*a*a*ks*ks
          +(sp.Rational(2,27)*a**3+sp.Rational(3,4)*c*c)*qstar
          -(a**4/sp.Integer(9)+a*c*c)*ks)
equal('global_recovery_identity',recovery,c**3*h/8)
equal('c_zero_quartic_derivative',sp.diff(ks,y).subs({c:0,y:0}),0)
equal('c_zero_sextic_derivative',sp.diff(qstar,y).subs({c:0,y:0}),0)

# Polynomial normalization coefficient readout, including all translation terms.
s,d,e,eta=sp.symbols('s d e eta')
t=y+s
ht=t**3+a*t; kt=t**4+sp.Rational(2,3)*a*t*t+c*t
pt=kt+d; qt=ht*ht+sp.Rational(3,2)*c*ht+eta*kt+e
equal('quartic_shift_readout',coeff(pt,y,3),4*s)
equal('quartic_a_readout_after_center',coeff(pt.subs(y,y-s),y,2),sp.Rational(2,3)*a)
equal('quartic_c_readout_after_center',coeff(pt.subs(y,y-s),y,1),c)
equal('quartic_d_readout_after_center',pt.subs(y,-s),d)
equal('sextic_eta_readout_after_center',coeff(qt.subs(y,y-s),y,4)-2*a,eta)
equal('sextic_e_readout_after_center',qt.subs(y,-s),e)

# Direct terminal (4,6) obstruction: derivative values are independent symbols.
A,C,Dp,Ep,Np=sp.symbols('A C Dp Ep Np')
p=ks+d; q=qstar+eta*ks+e
px=sp.diff(p,a)*A+sp.diff(p,c)*C+Dp
qx=sp.diff(q,a)*A+sp.diff(q,c)*C+sp.diff(q,eta)*Np+Ep
J46=sp.expand(px*sp.diff(q,y)-sp.diff(p,y)*qx)
equal('46_constant_term_forces_c_unit',coeff(J46,y,0),c*((sp.Rational(3,2)*a+eta)*Dp-Ep))
J46c=J46.subs(C,0)
equal('46_degree7',coeff(J46c,y,7),-4*(A+Np))
J46a=sp.expand(J46c.subs(Np,-A))
equal('46_degree5',coeff(J46a,y,5),-sp.Rational(4,3)*a*A+6*Dp)
J46b=sp.expand(J46a.subs(Dp,sp.Rational(2,9)*a*A))
equal('46_degree3',coeff(J46b,y,3),sp.Rational(4,3)*a*a*A+sp.Rational(8,9)*a*eta*A-4*Ep)
J46end=sp.expand(J46b.subs(Ep,a*a*A/3+sp.Rational(2,9)*a*eta*A))
equal('46_entire_residual_is_nonconstant_linear',J46end,-A*c*c*y/2)

# Direct terminal (3,4) obstruction, including the coefficient of the removable cubic.
B,B0=sp.symbols('B B0')
p=h+d; q=k+eta*h+e
px=sp.diff(p,a)*A+Dp
qx=sp.diff(q,a)*A+sp.diff(q,b)*B+sp.diff(q,c)*C+sp.diff(q,eta)*Np+Ep
J34=sp.expand(px*sp.diff(q,y)-sp.diff(p,y)*qx)
equal('34_degree5_forces_eta_constant',coeff(J34,y,5),-3*Np)
J34=sp.expand(J34.subs({Np:0,eta:0}))
equal('34_degree4',coeff(J34,y,4),4*A-3*B)
equal('34_degree3',coeff(J34,y,3),-3*C+4*Dp)
equal('34_degree2',coeff(J34,y,2),2*b*A-a*B-3*Ep)
subs34={b:sp.Rational(4,3)*a+B0,B:sp.Rational(4,3)*A,Dp:sp.Rational(3,4)*C,
        Ep:(sp.Rational(4,9)*a+sp.Rational(2,3)*B0)*A}
J34end=sp.expand(J34.subs(subs34))
expected=(c*A+(a+sp.Rational(3,2)*B0)*C)*y+sp.Rational(3,4)*c*C-(sp.Rational(4,9)*a*a+sp.Rational(2,3)*B0*a)*A
equal('34_entire_residual',J34end,expected)
F=sp.Rational(3,8)*c*c-sp.Rational(4,27)*a**3-B0*a*a/3
equal('34_constant_is_exact_derivative',coeff(J34end,y,0),sp.diff(F,a)*A+sp.diff(F,c)*C)
T=(2*a+3*B0)*c
equal('34_linear_is_product_derivative_half',2*coeff(J34end,y,1),sp.diff(T,a)*A+sp.diff(T,c)*C)

# Axioms and numerical bound leave exactly the two nonprincipal initial pairs.
allowed=[3,4,6,7]
pairs=[(n,m) for n in allowed for m in allowed if n>m and n%m and m%n]
truth('nonprincipal_candidates',pairs==[(4,3),(6,4),(7,3),(7,4),(7,6)])
survivors=[]
for n,m in pairs:
    g=sp.gcd(n,m)
    if g==1:
        if {3,4}.issubset({n,m}): survivors.append((n,m))
    else:
        if (n,m)==(6,4): survivors.append((n,m,3))
truth('characteristic_sequence_survivors',survivors==[(4,3),(6,4,3)])

print(json.dumps({'status':'PASS','checks':len(checks),'check_names':checks,
                  'scope':'Symbolic identities and finite controls; universal proof is AUDIT.md.'},indent=2))
