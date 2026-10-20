#!/usr/bin/env python3
"""Exact controls for the separately stated one-versus-two appendix."""
import json
import sympy as sp

y,t,s,a,b,c,d=sp.symbols('y t s a b c d')
U,V,lam,mu,nu,e,f=sp.symbols('U V lam mu nu e f')
checks={}
def check(name,value):
    if not bool(value):
        raise RuntimeError(name)
    checks[name]=True
def zero(expr):
    return sp.cancel(sp.expand(expr))==0

h=y**3+a*y
k=y**4+sp.Rational(2,3)*a*y**2+c*y
q=h*h+sp.Rational(3,2)*c*h
L=s**3+sp.Rational(4,3)*a*s-c
quad=y*y-s*y+s*s+a
hrem=sp.rem(h,quad,y)
krem=sp.rem(k,quad,y)
check('collision:h_remainder_constant',zero(hrem+s*(s*s+a)))
check('collision:k_linear_coefficient',zero(sp.Poly(krem,y).nth(1)+L))
check('collision:separation_squared',zero(sp.discriminant(quad,y)+3*s*s+4*a))
check('collision:resultant_minus27c2',sp.factor(sp.resultant(L,-3*s*s-4*a,s))==-27*c*c)
check('collision:common_h_value',zero(sp.rem(hrem+c-a*s/3,L,s)))
T=sp.Matrix([V,U*U+sp.Rational(3,2)*c*U])
det=T.jacobian([U,V]).det()
check('singular:forced_determinant',zero(det+2*U+sp.Rational(3,2)*c))
check('singular:determinant_at_collision',zero(det.subs(U,-c+a*s/3)-c/2+2*a*s/3))
s1,s2=sp.symbols('s1 s2')
check('singular:distinct_determinants',zero((c/2-2*a*s1/3)-(c/2-2*a*s2/3)+2*a*(s1-s2)/3))
r=sp.symbols('r')
triple=sp.Poly((s-r)**3,s)
check('cubic:quadratic_coefficient',triple.nth(2)==-3*r)
check('cubic:triple_root_zero_constant',triple.nth(0).subs(r,0)==0)
aff=sp.Matrix([lam*U+e,mu*V+nu*U+f])
check('intrinsic:affine_lift_determinant',aff.jacobian([U,V]).det()==lam*mu)
H=h.subs({a:0,y:d*t})
K=k.subs({a:0,c:d**3,y:d*t})
Q=q.subs({a:0,c:d**3,y:d*t})
check('scaling:h',zero(H-d**3*t**3))
check('scaling:k',zero(K-d**4*(t**4+t)))
check('scaling:q',zero(Q-d**6*(t**6+sp.Rational(3,2)*t**3)))
check('scaling:second_witness',zero(8*(t**6+sp.Rational(3,2)*t**3)+3-(8*t**6+12*t**3+3)))
print(json.dumps({'status':'PASS: exact controls for the exact-count appendix','sympy_version':sp.__version__,'check_count':len(checks),'checks':checks},sort_keys=True,indent=2))
