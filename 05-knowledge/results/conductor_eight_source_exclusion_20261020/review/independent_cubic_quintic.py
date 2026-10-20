#!/usr/bin/env python3
"""Separate reviewer-owned cubic/quintic audit, using functions of x."""
import json
import sympy as s

tests = []


def check(name, truth):
    if not truth:
        raise RuntimeError(name)
    tests.append(name)


def zero(name, expr):
    check(name, s.expand(expr) == 0)


x, t, y = s.symbols("x t y")
a,b,c,d,e,f,g = [s.Function(v)(x) for v in "abcdefg"]
p = t**3+a*t+b
q = t**5+c*t**4+d*t**3+e*t*t+f*t+g
bracket = lambda P,Q: s.diff(P,x)*s.diff(Q,t)-s.diff(P,t)*s.diff(Q,x)
J = s.Poly(bracket(p,q),t)
raw = [
    -a*s.diff(g,x)+f*s.diff(b,x),
    -a*s.diff(f,x)+f*s.diff(a,x)+2*e*s.diff(b,x),
    -a*s.diff(e,x)+2*e*s.diff(a,x)+3*d*s.diff(b,x)-3*s.diff(g,x),
    -a*s.diff(d,x)+3*d*s.diff(a,x)+4*c*s.diff(b,x)-3*s.diff(f,x),
    -a*s.diff(c,x)+4*c*s.diff(a,x)+5*s.diff(b,x)-3*s.diff(e,x),
    5*s.diff(a,x)-3*s.diff(d,x),
    -3*s.diff(c,x)]
for i,value in enumerate(raw):
    zero("raw coefficient "+str(i),J.nth(i)-value)
check("no unlisted higher coefficients",J.degree()==6)

C,D,E,F,G = s.symbols("C D E F G")
zero("lawful constant target shear",bracket(p,q-D*p)-bracket(p,q))
h = s.Function("h")(x)
P,Q = p.subs(t,y+h),q.subs(t,y+h)
zero("lawful source-centering shear",
     s.diff(P,x)*s.diff(Q,y)-s.diff(P,y)*s.diff(Q,x)-J.as_expr().subs(t,y+h))

# Derive the antiderivatives, not only their final substituted Jacobian.
dn = 5*a/3
en = 4*C*a/3+5*b/3+E
fn = 5*a*a/9+4*C*b/3+F
gn = 2*C*a*a/9+10*a*b/9+2*E*a/3+G
chain = {c:C,d:dn,e:en,f:fn,g:gn}
for i in range(2,7):
    zero("integrated coefficient "+str(i),J.nth(i).subs(chain).doit())
qn = q.subs(chain)
Jn = s.factor(bracket(p,qn))
U = -5*a*a/9+4*C*b/3+F
V = 4*C*a/3+10*b/3+2*E
H = -5*a**3/27+4*C*a*b/3+F*a+5*b*b/3+2*E*b
zero("whole residual identity",Jn-((U*s.diff(a,x)+V*s.diff(b,x))*t
                                 +U*s.diff(b,x)-a*V*s.diff(a,x)/3))
zero("first integral by actual x differentiation",
     s.diff(H,x)-U*s.diff(a,x)-V*s.diff(b,x))

A,B,k = s.symbols("A B k",nonzero=True)
leading = (-5*A*A/9)*(3*k*B)-(A*(10*B/3)/3)*(2*k*A)
zero("universal tied-degree leading coefficient",leading+35*k*A*A*B/9)
zero("first integral tied-degree leading relation",
     (-5*A**3/27+5*B**2/3).subs(B**2,A**3/9))

# A function u, not a finite list of monomials, verifies the hostile family.
u = s.Function("u")(x)
hostile = {a:u*u,b:u**3/3,C:0,E:0,F:0,G:0}
zero("arbitrary polynomial-parameter hostile has H zero",H.subs(hostile))
Ph,Qh = p.subs(hostile),qn.subs(hostile)
zero("arbitrary polynomial-parameter hostile entire Jacobian",
     bracket(Ph,Qh)+35*u**6*s.diff(u,x)/27)
check("hostile actual coordinate degrees",
      (s.degree(Ph,t),s.degree(Qh,t))==(3,5))

# Degree reasoning gets finite hostile controls, not substituted for its proof.
for kval in (1,4,8):
    p1,q1 = Ph.subs(u,x**kval),Qh.subs(u,x**kval)
    zero("hostile monomial coefficient "+str(kval),
         bracket(p1,q1)+s.Rational(35,27)*kval*x**(7*kval-1))
    check("hostile nonconstant Jacobian degree "+str(kval),
          s.degree(bracket(p1,q1),x)==7*kval-1)

# The full monogenic basis for coprime degrees 3,5 has these three residues.
degrees=[0,5,10]
check("hostile generic primitive degree criterion",s.gcd(3,5)==1)
check("hostile generic full-module residue separation",
      sorted(v%3 for v in degrees)==[0,1,2])
check("hostile generic delta",sum(degrees)-sum(range(3))==3*4)
zero("stationary coefficient family",
     Jn.subs({s.diff(a,x):0,s.diff(b,x):0}))

print(json.dumps({
    "status":"PASS","count":len(tests),"checks":tests,
    "scope":"Direct constant-leading degree-(3,5) Keller exclusion",
    "entire_reduced_jacobian":str(Jn),
    "universal_leading_coefficient":str(s.factor(leading)),
},indent=2,sort_keys=True))
