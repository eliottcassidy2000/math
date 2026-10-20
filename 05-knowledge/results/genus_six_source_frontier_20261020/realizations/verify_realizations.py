#!/usr/bin/env python3
"""Exact independent certificates for the displayed polynomial realizations."""
import json
from math import gcd
from functools import reduce
import sympy as s

checks=[]
def check(name, condition):
    if not bool(condition):
        raise RuntimeError('FAILED: '+name)
    checks.append(name)

def zero(name, expression):
    check(name, s.expand(expression)==0)

t,a,b,z,X,Y=s.symbols('t a b z X Y')
R=s.Rational

p4=t**4+t
q10=t**10+R(5,2)*t**7+R(15,8)*t**4+R(5,16)*t
r5=s.expand(q10*q10-p4**5-R(5,64)*p4*p4)
zero('quartic-decisive-remainder',r5-(t**5/R(64)+5*t*t/R(256)))
zero('quartic-derivative-remainder',s.rem(s.diff(q10,t),s.diff(p4,t),t)+R(5,8))
check('quartic-derivative-gcd',s.gcd(s.diff(q10,t),s.diff(p4,t))==1)
check('quartic-field-index-gcd',gcd(s.degree(p4,t),s.degree(r5,t))==1)

pb=t**4+b*t
qb=t**10+R(5,2)*b*t**7+R(15,8)*b*b*t**4+R(5,16)*b**3*t
rb=s.expand(qb*qb-pb**5-R(5,64)*b**4*pb*pb)
zero('quartic-parameter-identity',rb-b**5*t**5/64-5*b**6*t*t/256)
zero('quartic-parameter-derivative-remainder',s.rem(s.diff(qb,t),s.diff(pb,t),t)+5*b**3/8)
check('quartic-nonprimitive-boundary',pb.subs(b,0)==t**4 and qb.subs(b,0)==t**10)

p6=t**6+t**4
q9=t**9+R(3,2)*t**7+R(3,8)*t**5-R(5,32)*t**3-R(3,128)*t
r4=s.expand(q9*q9-p6**3+R(3,16)*p6*p6-R(7,1024)*p6)
zero('sextic-decisive-remainder',r4-t**4/2048-9*t*t/16384)
zero('sextic-derivative-factorization',s.diff(p6,t)-2*t**3*(3*t*t+2))
check('sextic-derivative-at-zero',s.diff(q9,t).subs(t,0)==-R(3,128))
zero('sextic-derivative-quadratic-remainder',s.rem(s.diff(q9,t),3*t*t+2,t)+R(27,128))
check('sextic-derivative-gcd',s.gcd(s.diff(p6,t),s.diff(q9,t))==1)
check('sextic-field-index-gcd',reduce(gcd,[int(s.degree(v,t)) for v in [p6,q9,r4]])==1)

D=a*a-4*b
p=t**6+a*t**4+b*t*t
q=t**9+3*a*t**7/2+(3*b/2+3*a*a/8)*t**5+a*(-5*a*a+36*b)*t**3/32+3*(16*b*b-a**4)*t/128
A=-3*a*D/16
B=D**2*(7*a*a+8*b)/1024
r=s.expand(q*q-p**3-A*p*p-B*p)
zero('sextic-parameter-remainder',r-D**3*t*t*(9*a*a+8*a*t*t-4*b)/16384)
zero('sextic-parameter-specialization-p',p.subs({a:1,b:0})-p6)
zero('sextic-parameter-specialization-q',q.subs({a:1,b:0})-q9)
zero('sextic-parameter-remainder-leading',s.Poly(r,t).coeff_monomial(t**4)-a*D**3/2048)
G=3*z*z+2*a*z+b
Q1=9*z**4+21*a*z**3/2+(15*b/2+15*a*a/8)*z*z+3*a*(-5*a*a+36*b)*z/32+3*(16*b*b-a**4)/128
zero('sextic-parameter-derivative-factorization',s.diff(p,t)-2*t*G.subs(z,t*t))
zero('sextic-parameter-q-derivative',s.diff(q,t)-Q1.subs(z,t*t))
zero('sextic-parameter-q-derivative-at-zero',Q1.subs(z,0)+3*D*(a*a+4*b)/128)
zero('sextic-parameter-q-derivative-remainder',s.rem(Q1,G,z)-3*D*(-a*a+12*a*z+12*b)/128)
zero('sextic-parameter-q-derivative-resultant',s.resultant(Q1,G,z)-R(6561,16384)*D**4)
zero('sextic-genus-four-degeneration',r.subs(a,0)-b**4*t*t/64)
zero('sextic-genus-four-previous-p',p.subs({a:0,b:4})-(t**6+4*t*t))
zero('sextic-genus-four-previous-q',q.subs({a:0,b:4})-(t**9+6*t**5+6*t))
h=t**3+a*t/2
zero('sextic-nonprimitive-boundary-p',p.subs(b,a*a/4)-h*h)
zero('sextic-nonprimitive-boundary-q',q.subs(b,a*a/4)-h**3)
check('sextic-nonimmersed-divisor-control',s.gcd(s.diff(p.subs({a:2,b:-1}),t),s.diff(q.subs({a:2,b:-1}),t))==t)

# Exhaustiveness and exceptional-branch controls in the even/odd ansatz.
u,v=s.symbols('u v')
qans=t**9+3*a*t**7/2+(3*b/2+3*a*a/8)*t**5+u*t**3+v*t
Aans=a**3/8-3*a*b/2+2*u
vans=7*a**4/128-9*a*a*b/16+a*u/2+3*b*b/8
rans=s.Poly(s.expand(qans*qans-p**3-Aans*p*p),t)
zero('parity-ansatz-degree-twelve',rans.coeff_monomial(t**12))
zero('parity-ansatz-degree-ten-elimination',rans.coeff_monomial(t**10).subs(v,vans))
zero('parity-ansatz-degree-eight-obstruction',rans.coeff_monomial(t**8).subs(v,vans)-D*(5*a**3-36*a*b+32*u)/128)
c=u-a**3/8
zero('parity-ansatz-exceptional-composite-q',qans.subs(v,vans).subs(b,a*a/4)-h**3-c*h)
zero('parity-ansatz-exceptional-A',Aans.subs(b,a*a/4)-2*c)
zero('parity-ansatz-exceptional-B',rans.coeff_monomial(t**6).subs(v,vans).subs(b,a*a/4)-c*c)

for n,H,basis in [(4,X**5+R(5,64)*X**2,[1,Y,Y**2-(X**5+R(5,64)*X**2),Y*(Y**2-(X**5+R(5,64)*X**2))]),
                  (6,X**3-R(3,16)*X**2+R(7,1024)*X,[Y**j*(Y**2-(X**3-R(3,16)*X**2+R(7,1024)*X))**i for i in range(3) for j in range(2)])]:
    M=s.Matrix([[s.Poly(v,Y).coeff_monomial(Y**j) for j in range(n)] for v in basis])
    check(f'module-unimodularity-{n}',s.det(M)==1)
    check(f'module-upper-triangular-{n}',all(M[i,j]==0 for i in range(n) for j in range(i+1,n)))

def semigroup(gens,bound):
    S={0}
    for m in range(bound+1):
        if m in S:
            S.update(m+g for g in gens if m+g<=bound)
    return S

for name,n,gens,basis,gaps in [
 ('quartic',4,[4,5],[s.Integer(1),q10,r5,q10*r5],[1,2,3,6,7,11]),
 ('sextic',6,[4,6,9],[s.Integer(1),q9,r4,q9*r4,r4*r4,q9*r4*r4],[1,2,3,5,7,11])]:
    degrees=[int(s.degree(s.expand(v),t)) for v in basis]
    check(name+'-residue-separated',len({d%n for d in degrees})==n)
    check(name+'-apery-codimension',(sum(degrees)-n*(n-1)//2)//n==6)
    bound=120
    S=semigroup(gens,bound)
    union={d+n*j for d in degrees for j in range(bound//n+1) if d+n*j<=bound}
    check(name+'-degree-monoid-bounded-equality',S==union)
    check(name+'-degree-gaps',[i for i in range(bound+1) if i not in union]==gaps)
    check(name+'-conductor-residue-block',all(i in S for i in range(12,12+n)) and 11 not in S)

print(json.dumps({'status':'PASS','check_count':len(checks),'checks':checks,
                  'concrete_degrees':{'quartic':[0,10,5,15],'sextic':[0,9,4,13,8,17]},
                  'bounded_semigroup_check_universe':'integers 0..120; all-degree proof in accompanying note'},indent=2,sort_keys=True))
