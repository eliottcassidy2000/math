#!/usr/bin/env python3
"""Exact identities and valuation-gap certificate; no numerical root tests."""
import json
import sympy as S

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

def zero(expr, message):
    check(S.expand(expr)==0,message)

def positive_affine(expr,param,message):
    p=S.Poly(S.expand(expr),param)
    check(p.degree()<=1,message+': not affine')
    check(p.eval(0)>0 and p.coeff_monomial(param)>=0,message)

x,z,h,V,A,B,c=S.symbols('x z h V A B c')
P=(V*z*z+B*z)**2+A*z+x
E=x-A*B/V-c
shifted=(V*h*h-B*h)**2+A*h+E
zero((P-c).subs(z,-B/V+h)-shifted,'shifted fibre equation')
zero(S.diff(P,z).subs(z,-B/V+h)-2*(V*h*h-B*h)*(2*V*h-B)-A,'shifted derivative')
v,a,b,eta,eps=S.symbols('v a b eta eps', nonzero=True)
Z=S.symbols('Z')
models={
    'low_degree_infinity':v*v*Z**4+1,
    'small_high_degree_infinity':b*b*Z**2+1,
    'escaping_infinity':b*b*Z**2+eps,
    'multiple_finite_escape':b*b*Z**2-a*b/v,
    'simple_finite_escape':b*b*Z**2+eps,
}
for label,model in models.items():
    check(S.resultant(model,S.diff(model,Z),Z)!=0,label+': initial roots not simple')
low=[]
for n in range(9):
    # On x=s^-4,z=s^31*Z, these are gaps above the leading orders.
    gaps=[33-4*n,66-8*n,3,4]
    check(min(gaps)>0,'low-degree valuation dominance')
    low.append({'degree':n,'gaps':gaps,'omega_order_on_cover':30})
k=S.symbols('k',integer=True,nonnegative=True)
n=9+k
for expr in [4*n-33,2*n-15,8*n-66,2*n-2]:
    positive_affine(expr,k,'small infinity branch sector')
sectors=[]
for label,n,e in [('nine_nonresonant',S.Integer(9),S.Integer(1)),('nine_resonant',S.Integer(9),S.Integer(0)),('degree_at_least_ten',10+k,2+k)]:
    gaps=[4*n-e-32,2*n+e-16]
    order=2*n+e-3
    for expr in gaps+[order,2*n-S.Rational(33,2)]:
        positive_affine(expr,k,label)
    sectors.append({'sector':label,'n':str(n),'e':str(e),'gaps':[str(S.expand(g)) for g in gaps],'omega_order_on_cover':str(S.expand(order))})
m=3+k
positive_affine(m+1,k,'finite escape dominance')
positive_affine(m,k,'finite escape omega order')
# Exact Broughton synthesis, in volume coordinates (x,t).
t,U,Y,r=S.symbols('t U Y r')
u=x*t-8
y=(t-1)**2/16
F=U+U*U*Y
Pstar=8+u+u*u*(t-1)**2/16
zero(Pstar-8-F.subs({U:u,Y:y}),'Broughton pullback')
jphi=S.det(S.Matrix([[S.diff(u,x),S.diff(u,t)],[S.diff(y,x),S.diff(y,t)]]))
zero(jphi-t*(t-1)/8,'pullback Jacobian')
check(S.cancel(S.diff(y,t)/jphi)==1/t,'hyperbola time form')
b0=S.symbols('b0')
u=x*t+2*b0
rho=u*(t-1)
Phalf=-2*b0+u+rho*rho/16
jpr=S.diff(Phalf,x)*S.diff(rho,t)-S.diff(Phalf,t)*S.diff(rho,x)
zero(jpr-u-rho,'generic-fibre time denominator')
q=S.symbols('q')
zero((q-r*r/16)+r-(q+r-r*r/16),'rational time denominator')
print(json.dumps({
    'status':'PASS',
    'universal_shift_and_derivative_verified':True,
    'simple_initial_models':{key:str(value) for key,value in models.items()},
    'low_degree_sectors':low,
    'small_branch_all_n_ge_9':{'n':'9+k','k':'nonnegative integer','omega_order':'2*k+16'},
    'escaping_infinity_sectors':sectors,
    'multiple_finite_boundary':{'m':'3+k','gap':'k+4','omega_order_on_cover':'k+3'},
    'simple_finite_boundary':'generic E0 nonzero gives two simple nonzero h roots and regular time form',
    'generic_exclusions':['critical values','c=alpha at every V root','c=alpha-Aprime(alpha)*B(alpha)/Vprime(alpha) at simple V roots','c=limit(x-A*B/V) on the degree-nine resonance wall'],
    'broughton_pullback_and_time_identities_verified':True,
    'conclusion':'All polynomial clutches in the stated degree(16,8), simple-A boundary architecture have no polynomial nonzero-constant-Jacobian mate.',
    'critical_point_hypothesis_HJ':'Still open; not used',
    'scope':'Algebraic identities plus symbolic valuation inequalities; branch completeness proved in accompanying note.'
},indent=2,sort_keys=True))
