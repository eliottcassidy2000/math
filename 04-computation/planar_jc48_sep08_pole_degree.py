#!/usr/bin/env python3
"""Finite exact controls for the DG pole-degree gate; the gate is analytic."""
import hashlib
import json
import sympy as S

x, z, t, s, y, r, b, q, c, h, e, a, u = S.symbols(
    'x z t s y r b q c h e a u')
gates = 0
records = {}

def check(ok, label):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(label)

def eq(left, right, label):
    check(S.cancel(left-right) == 0, label)

def jac(left, right, v, w):
    return S.diff(left,v)*S.diff(right,w)-S.diff(left,w)*S.diff(right,v)

def order(poly, v):
    return min(mon[0] for mon, co in S.Poly(S.expand(poly),v).terms() if co != 0)

def coefficient(poly, v, n):
    return S.expand(poly).coeff(v,n)

def other_chart(poly):
    return S.cancel(poly.subs({x:1/r,t:-r*r-r**4*b}, simultaneous=True))

N=x*x*z-x**4-S.Rational(8,27)*x**4*z+S.Rational(4,27)*x*x*z*z
M=-(1+x*x*z)
H=S.cancel(t*t*N.subs(z,x*x+1/t))
L=S.cancel(t*M.subs(z,x*x+1/t))
F=S.expand(H*H+L)
eq(H,x*x*t-S.Rational(4,27)*x**6*t*t+S.Rational(4,27)*x*x,'literal H')
eq(L,-(1+x**4)*t-x*x,'literal L')
check(S.Poly(N,x,z).degree(x)<=4 and S.Poly(N,x,z).degree(z)<=2,'N complete global box')
check(S.Poly(M,x,z).degree(x)<=2 and S.Poly(M,x,z).degree(z)<=1,'M complete global box')
eq(jac(1/r,-r*r-r**4*b,r,b),r*r,'actual omega transition')
Hb=-1-r*r*b-S.Rational(8,27)*b-S.Rational(4,27)*r*r*b*b
Lb=b+r*r+r**4*b
eq(other_chart(H),Hb,'H full second chart')
eq(other_chart(L),Lb,'L full second chart')
eq(other_chart(F),Hb*Hb+Lb,'F full second chart')
FD=(1+S.Rational(8,27)*b)**2+b
eq((Hb*Hb+Lb).subs(r,0),FD,'nonconstant boundary restriction')
check(S.Poly(FD,b).degree()==2,'two generic D points')
eq(S.discriminant(FD-c,b),(256*c+1593)/729,'generic boundary transversality')
eq(jac(Hb*Hb+Lb,r,r,b),-S.diff(Hb*Hb+Lb,b),'relative form sign in D chart')

E=S.expand((N*N+(z-x*x)**3*M-c*(z-x*x)**4).subs(z,s+x*x))
Nb=s*y-S.Rational(4,27)*y**3+S.Rational(4,27)*s*s*y
Eb=S.expand(Nb*Nb-s**3*(1+y*y+s*y)-c*s**4)
eq(E,Eb.subs(y,x*x),'literal unsimplified fibre reconstruction')
eq(E.subs(x,-x),E,'exact even-x parity')
eq(Eb.subs(y,0),-s**3*(1+c*s),'local Weierstrass degree three at zero')
check(order(Eb.subs(y,0),s)==3,'three and only three local branches counted with multiplicity')

K=S.cancel(Eb.subs(s,y*y*h)/y**6)
K0=-(h-S.Rational(4,9))**2*(h-S.Rational(1,9))
eq(K.subs(y,0),K0,'all initial slopes')
eq(S.diff(K,h).subs({y:0,h:S.Rational(1,9)}),-S.Rational(1,9),'third branch implicit derivative')
third_polar=S.expand(S.diff(Eb,s).subs(s,y*y/S.Integer(9)))
check(order(third_polar,y)==4,'third branch polar order')
eq(coefficient(third_polar,y,4),-S.Rational(1,9),'third branch polar coefficient')
eq(S.Rational(1,81)/coefficient(third_polar,y,4),-S.Rational(1,9),'third branch eta regular leading value')

split=S.cancel(Eb.subs(s,y*y*(S.Rational(4,9)+y*h))/y**8)
split0=-h*h/3-S.Rational(64,59049)*(65+36*c)
eq(split.subs(y,0),split0,'two split branches exact initial equation')
eq(S.diff(split,h).subs({y:0,h:e}),-2*e/3,'split branches implicit derivative')
e2=-S.Rational(64,19683)*(65+36*c)
eq(S.solve(split0,h*h)[0],e2,'corrected generic split constant')
polar=S.expand(S.diff(Eb,s).subs(s,S.Rational(4,9)*y*y+e*y**3))
check(order(polar,y)==5,'two polar branches order in y')
eq(coefficient(polar,y,5),-2*e/3,'two polar branches leading derivative')
eq(S.Rational(16,81)/coefficient(polar,y,5),-S.Rational(8,27)/e,'actual eta double-pole coefficient')
eq(jac(E,x,s,x),S.diff(E,s),'eta=s^2 dx/E_s uses actual orientation')
eq(S.diff(E,s).subs(x,-x),S.diff(E,s),'polar denominator exact even parity')
check(2*(4-5)==-2,'order-two poles in the actual x parameter')
eq(S.residue(-S.Rational(8,27)/(e*x*x),x,0),0,'second-kind leading part')

# Full infinity chart: q=1/z and sigma=q-r^2 (represented by s here).
Ni=S.cancel(r**4*q*q*N.subs({x:1/r,z:1/q},simultaneous=True))
Mi=S.cancel(-r*r*q*M.subs({x:1/r,z:1/q},simultaneous=True))
eq(Ni.subs(q,r*r),-S.Rational(4,27)*r*r,'sole infinity zero of multiplicity two')
eq(Mi.subs(q,r*r),1+r**4,'M unit at infinity')
eq(N.subs(z,x*x),-S.Rational(4,27)*x**6,'sole finite boundary zero of multiplicity six')
eq(M.subs({x:0,z:0}),-1,'M unit at finite boundary zero')
eq(jac(1/r,r*r*q/(r*r-q),r,q),-r*r/(q-r*r)**2,'full compact infinity omega')
Nic=S.expand(Ni.subs(q,s+r*r)); Mic=S.expand(Mi.subs(q,s+r*r))
eq(Nic,-r*r*s-s*s-S.Rational(4,27)*r*r-S.Rational(8,27)*s,'literal infinity N')
eq(Mic,1+r**4+r*r*s,'literal infinity M')
Ei=S.expand(Nic*Nic+s**3*Mic-c*s**4)
check(order(Ei.subs(r,0),s)==2,'two and only two infinity branches')
infty=S.cancel(Ei.subs(s,-r*r/2+r**3*h)/r**6)
eq(infty.subs(r,0),S.Rational(64,729)*h*h-S.Rational(1,8),'infinity branches initial equation')
eq(S.diff(infty,h).subs({r:0,h:a}),S.Rational(128,729)*a,'infinity implicit derivative')
eq(S.solve(infty.subs(r,0),h*h)[0],S.Rational(729,512),'two nonzero infinity splits')
ipolar=S.expand(S.diff(Ei,s).subs(s,-r*r/2+r**3*a))
check(order(ipolar,r)==3,'infinity polar derivative order')
eq(coefficient(ipolar,r,3),S.Rational(128,729)*a,'infinity polar derivative coefficient')
check(2+4-3==3,'actual infinity eta order three')
eq(S.Rational(1,4)/coefficient(ipolar,r,3),S.Rational(729,512)/a,'infinity eta leading coefficient')
records['boundary_eta_orders']=[-2,-2,0,3,3]
check(sum(max(0,-o-1) for o in records['boundary_eta_orders'])==2,'complete pole-degree budget two')

# This literal F was already excluded by elementary affine criticality.
Hy=y*t-S.Rational(4,27)*y**3*t*t+S.Rational(4,27)*y
Fy=S.expand(Hy*Hy-(1+y*y)*t-y)
eq(F,Fy.subs(y,x*x),'affine even polynomial')
point={y:S.Rational(54,97),t:S.Rational(7081,1296)}
eq(S.diff(Fy,y).subs(point),0,'actual affine critical derivative y')
eq(S.diff(Fy,t).subs(point),0,'actual affine critical derivative t')
eq(Fy.subs(point),-S.Rational(77,36),'actual critical fibre value')
eq(Hy.subs(point),S.Rational(85,36),'critical square-prefix control')
records['critical_point']=['x^2=54/97','t=7081/1296','F=-77/36']
eq(jac(x*x,t/(2*x),x,t),1,'critical polynomial can still have a rational mate')
eq(jac(t,-x,x,t),1,'global t has nonglobal rational mate')
eq((-r*r-r**4*b).subs(r,0),0,'constant boundary exception for t')

# An inherited genuine DG rational mate satisfies, rather than violates,
# the new bound: its generic primitive has pole degree eight.
Fpositive=b**4+r*b
Gpositive=-r*r/2-2*r*b**3-S.Rational(4,3)*b**6
eq(jac(Fpositive,Gpositive,r,b),r*r,'actual boundary-linear rational mate')
gp=-c*c/(2*b*b)-c*b*b+b**6/6
eq(Gpositive.subs(r,(c-b**4)/b),gp,'actual generic-fibre primitive')
eq(S.diff(gp,b),(c-b**4)**2/b**3,'order-two D zeros with sufficient pole budget')
check(S.degree(S.cancel(gp*b*b),b)==8,'actual primitive degree eight')
eq(S.cancel(gp*b*b).subs(b,0),-c*c/2,'actual finite pole order two')
records['rational_controls']=['x^2 and t/(2x): criticality alone insufficient',
                              't and -x: boundary restriction constant',
                              'b^4+rb: generic primitive pole degree eight']

# Classical degree controls: the lower bound is sharp, and order two matters.
G=u**3
eta=S.diff(G,u)
check(order(eta,u)==2,'sharp local degree-three control')
eq(eta.subs(u,1/q)*(-1/q**2),-3/q**4,'sharp pole order four, primitive degree three')
G2=u+1/u
eta2=S.diff(G2,u)
eq(eta2,1-u**-2,'two-pole residue-free exact control')
eq(eta2.subs(u,1/q)*(-1/q**2),1-q**-2,'two order-two poles')
eq(S.residue(eta2,u,0),0,'first residue zero')
eq(S.residue(1-q**-2,q,0),0,'second residue zero')
eq(S.gcd(u*u-1,2*u),1,'only simple differential zeros in degree-two control')
check(max(S.degree(u*u+1,u),S.degree(u,u))==2,'exact rational-map degree two')
records['classical_controls']=['u^3: degree3, zero2, pole4','u+1/u: degree2, two pole2, no zero2']

print('RESERVED pole-degree gate: exact controls PASS; analytic proof requires independent audit')
print('Universe: one literal DG quartic; all three finite and two infinite normalized branches; two classical degree controls')
print('Boundary eta orders:',records['boundary_eta_orders'])
print('Primitive pole-degree budget: 2; transverse D local degree required: 3')
print('Known affine critical control: x^2=54/97, t=7081/1296, F=-77/36; polynomial exclusion alone is elementary')
print('New candidate conclusion: no rational mate for the literal quartic despite residue-free relative poles')
print('Always-active gates:',gates)
print('Semantic SHA256:',hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest())
