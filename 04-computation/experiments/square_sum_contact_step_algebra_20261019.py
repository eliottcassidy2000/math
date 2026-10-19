#!/usr/bin/env python3
"""Exact symbolic check of the formulas and finite certificates in the companion note."""
import json
import sympy as s

checks=0
def require(test,label):
 global checks
 checks+=1
 if not test: raise RuntimeError(label)
def fac(x): return s.factor(s.cancel(x))
def eq(a,b,label): require(fac(a-b)==0,label)
def f(v): return 4*v*(v*v-1)/(v*v+1)**2

u,z,k,t=s.symbols('u z k t')
K=k*k;quartic=K*(u*u+1)**2-4*u*(u*u-1)
def red(x):
 n,d=s.fraction(s.cancel(x));n=s.rem(n,z*z-quartic,z);d=s.rem(d,z*z-quartic,z)
 d0=s.Poly(d,z).nth(0);d1=s.Poly(d,z).nth(1)
 n=s.rem(s.expand(n*(d0-d1*z)),z*z-quartic,z)
 return fac(n/(d0*d0-d1*d1*quartic))
X=(k*z+K+2*u+K*u*u)/(2*u*u);Y=(X-K)*(X*u-1)/k
a=K*(u*u-1)-u
require(red(Y*Y-(X-K)*(X*X-1))==0,'quartic to cubic')
require(red(a*a-K*z*z-u*(u-2*K)*(2*K*u+1))==0,'conjugate product')
for sign in [1,-1]:
 l=red((Y-sign*k)/X);xp=red(l*l+K-X);yp=red(-Y+l*(X-xp))
 vp=red((k*yp+xp-K)/(xp*(xp-K)))
 direct=(a+sign*k*z)/(u*(u-2*K));conjugate=(2*K*u+1)/(a-sign*k*z)
 require(red(vp-direct)==0,f'translation {sign}')
 require(red(direct-conjugate)==0,f'conjugate {sign}')
 zp=red((2*vp*vp*xp-K-2*vp-K*vp*vp)/k)
 reverse=sign*(K*(vp*vp-1)-vp-u*vp*(vp-2*K))/k
 require(red(zp-reverse)==0,f'z translation {sign}')

h=lambda v:v*v-2*v-1
F=2+1/u;G=1/(u-2)
eq((2*u+1)/((u*u-1)-u-h(u)),F,'F degeneration')
eq(((u*u-1)-u-h(u))/(u*(u-2)),G,'G degeneration')
eq((F*F-1)-F-u*F*(F-2),h(F),'F z branch')
eq(-((G*G-1)-G-u*G*(G-2)),h(G),'G z branch')
eq(1-f(u),h(u)**2/(u*u+1)**2,'square defect')

kt=(t*t+2*t-1)/(t*t+1);Kt=kt*kt
Q=((t*t+4*t+1)/(t*t+1),8*t*(t+1)/(t*t+1)**2);P=(s.S.Zero,kt)
def add(A,B):
 x,y=A;v,w=B
 ell=fac((3*x*x-2*Kt*x-1)/(2*y)) if A==B else fac((y-w)/(x-v))
 xp=fac(ell*ell+Kt-x-v);return xp,fac(-y+ell*(x-xp))
def inverse(T):
 x,y=T;v=fac((kt*y+x-Kt)/(x*(x-Kt)))
 return v,fac((2*v*v*x-Kt-2*v-Kt*v*v)/kt)
eq(Q[0],X.subs({k:kt,u:t,z:t*t+1}),'anchor X')
eq(Q[1],Y.subs({k:kt,u:t,z:t*t+1}),'anchor Y')
U=(2*t**3+t*t-2*t+1)/(t*(t*t-4*t+1))
V=(t*t+4*t+1)/(t**3+2*t*t+t-2)
ZU=-(t**6+24*t**5-13*t**4-16*t**3+19*t*t-8*t+1)/(t*t*(t*t-4*t+1)**2)
ZV=-(t**6+8*t**5+19*t**4+16*t**3-13*t*t-24*t+1)/(t**3+2*t*t+t-2)**2
for sign,expectedv,expectedz in [(1,U,ZU),(-1,V,ZV)]:
 v,w=inverse(add(Q,(P[0],sign*P[1])))
 eq(v,expectedv,f'q1 u {sign}');eq(w,expectedz,f'q1 z {sign}')
 require(s.limit(v,t,s.oo)==(2 if sign==1 else 0),f'q1 u limit {sign}')
 require(s.limit(w,t,s.oo)==-1,f'q1 z limit {sign}')
A=2*t**3+t*t-2*t+1;B=t*t-4*t+1
C=5*t**7+14*t**6+15*t**5-34*t**4-t**3+18*t*t-11*t+2
D=t**6+28*t**5+19*t**4-40*t**3+19*t*t-4*t+1
U2=B*C/(A*D)
V2=-(t**3+2*t*t+t-2)*(t**6+4*t**5+19*t**4+40*t**3+19*t*t-28*t+1)/((t*t+4*t+1)*(2*t**7+11*t**6+18*t**5+t**4-34*t**3-15*t*t+14*t-5))
P2=add(P,P)
for sign,expected in [(1,U2),(-1,V2)]:
 v,w=inverse(add(Q,(P2[0],sign*P2[1])));eq(v,expected,f'q2 u {sign}')

# q1 inequalities: exact polynomial differences, endpoint values, derivative signs.
eq(U+1,(t-1)*(3*t*t-1)/(t*(t*t-4*t+1)),'q1 lower sign')
eq(3*(t**3+2*t*t+t-2)-(t*t+4*t+1),(t-1)*(3*t*t+8*t+7),'q1 minus upper')
require((-3*t**3-4*t*t+5*t+8).subs(t,s.Rational(4,3))==s.Rational(4,9),'q1 minus lower endpoint')
eq((t**3+2*t*t+t-2)-(t*t+4*t+1),(t+1)*(t*t-3),'q1 minus small')
eq((U-2)*t*(t*t-4*t+1),9*t*t-4*t+1,'q1 plus >2')
eq((5-U)*t*(t*t-4*t+1),3*t**3-21*t*t+7*t-1,'q1 plus <5')
require(f(s.Rational(4,3))==f(s.Integer(7))==s.Rational(336,625),'admissibility gap')
require(f(s.Rational(3,2))==f(s.Integer(5))==s.Rational(120,169),'q1 spread lower bound')

# Patch certificate: positive shifted coefficients prove signs, no sampling.
derivative=fac(s.diff(U2,t));num,den=s.fraction(derivative)
eq(den,A*A*D*D,'derivative denominator')
coeff=s.Poly(num.subs(t,t+2),t).all_coeffs()
require(len(coeff)==17 and all(x>0 for x in coeff),'derivative positive certificate')
require(all(x>0 for x in s.Poly(A.subs(t,t+2),t).all_coeffs()),'A positive')
require(all(x>0 for x in s.Poly(D.subs(t,t+2),t).all_coeffs()),'D positive')
require(U2.subs(t,24)>1 and U2.subs(t,32)<s.Rational(5,4),'patch endpoint bounds')
require(f(s.Rational(5,4))==s.Rational(720,1681)<s.Rational(1,2),'patch spread')
quartics=[7*t**4-16*t**3+6*t*t-1,t**4+8*t**3+2*t*t-8*t+1,3*t**4+8*t**3+6*t*t-8*t+3]
sextics=[t**6-26*t**5+17*t**4+12*t**3-17*t*t+6*t-1,t**6+24*t**5-13*t**4-16*t**3+19*t*t-8*t+1]
for j,g in enumerate(quartics):require(all(x>0 for x in s.Poly(g.subs(t,t+24),t).all_coeffs()),f'positive distinctness quartic {j}')
for j,g in enumerate(sextics):
 require(s.Poly(g,t).LC()==1 and abs(g.subs(t,0))==1,f'rational-root certificate {j}')
 require(g.subs(t,1)!=0 and g.subs(t,-1)!=0,f'no rational root {j}')
E=t**8+4*t**7+104*t**6+44*t**5-110*t**4+12*t**3+48*t*t-28*t+5
H=29*t**8+52*t**7-72*t**6-100*t**5+218*t**4-164*t**3+64*t*t-12*t+1
factor=s.prod(quartics+sextics)
eq(f(U2)-f(t),8*(t*t+1)*(t*t+2*t-1)**2*(5*t*t-4*t+1)*factor/(E*E*H*H),'spread difference factorization')

pell=[s.Integer(0),s.Integer(1)]
for j in range(1,32):pell.append(2*pell[-1]+pell[-2])
limits=[]
for q in range(1,31):
 p0,p1,p2=pell[q-1:q+2];up=p2/p1;um=-p0/p1
 require(p2*p2-2*p2*p1-p1*p1==(-1)**q,f'Pell identity {q}')
 require(2<=up<=s.Rational(5,2),f'plus interval {q}')
 require(um==2-up,f'minus symmetry {q}')
 require(f(up)==1-1/(p2*p2+p1*p1)**2>=s.Rational(24,25),f'plus spread {q}')
 require(f(um)==1-1/(p0*p0+p1*p1)**2,f'minus spread {q}')
 if q>=2:require(f(um)>=s.Rational(24,25),f'minus lower bound {q}')
 if q<=6:limits.append({'q':q,'u_plus':str(up),'u_minus':str(um),'spread_plus':str(f(up))})
# Independent small-t cusp proof: the opposite limiting quartic component.
eq(((u*u-1)-u-h(u))/(u*(u-2)),G,'small-cusp plus map')
eq((2*u+1)/((u*u-1)-u-h(u)),F,'small-cusp minus map')
eq((G*G-1)-G-u*G*(G-2),-h(G),'small-cusp plus z branch')
eq(-((F*F-1)-F-u*F*(F-2)),-h(F),'small-cusp minus z branch')
Xs=(2*u+1)/(u*u);Ys=(-u*u+2*u+1)*(u+1)/u**3
eq(Ys*Ys,(Xs-1)**2*(Xs+1),'nodal component on cubic')
eq((Ys+Xs-1)/(Xs*(Xs-1)),u,'nodal inverse chart')
eq((2*t**3+t*t-2*t+1)/(t*(t*t-4*t+1))+1,(t-1)*(3*t*t-1)/(t*(t*t-4*t+1)),'small-cusp q1 boundary sign')
print(json.dumps({'status':'PASS','exact_checks':checks,'pell_limits':limits,'q2_patch':['24','32'],'q2_patch_spread_upper':'720/1681'},indent=2))
