"""Independent symbolic audit of the fixed-F fourfold unit response.

No producer import. The accompanying proof pays the all-degree lower bound,
all affine components, constant field and canonical connection hypotheses.
"""
from pathlib import Path
from hashlib import sha256
import sys
import sympy as s

sys.stdout.reconfigure(encoding='utf-8',newline='\n')
z,t,a,h,c,u,r,b,w=s.symbols('z t a h c u r b w')
count=0
def need(ok,why):
    global count
    count+=1
    if not ok:raise ArithmeticError(why)
def zero(expr,why):need(s.cancel(expr)==0,why)

# Derive in translated coordinates, with the special value c independent.
K=z**3*t+z-2*h
g=a*z*K
F=c+g
q=1/(3*a*z**3)
def D(f):return s.diff(F,z)*s.diff(f,t)-s.diff(F,t)*s.diff(f,z)
zero(D(q)-1,'rational primitive, direct Hamiltonian derivation')
zero(s.diff(F,t)-a*z**4,'only possible affine critical locus z=0')
zero(s.diff(F,z).subs(z,0)+2*a*h,'nonzero transverse derivative when ah nonzero')
zero((z**3*t+z-K)/(2*h)-1,'comaximal components explicit Bezout identity')

# Rational inverse pays the entire differential-constant field.
t_inverse=(u/a-z*(z-2*h))/z**4
zero(g.subs(t,t_inverse)-u,'k(z,t)=k(F)(z) inverse')
for j in range(5):
    zero(D(z**j)+a*z**4*s.diff(z**j,z),'derivation on fixed-F rational coordinate')

# Prove the full scalar principal part by a denominator-independent remainder.
alpha=-8*a*a*h**3/3
beta=-2*a*h
remaining=s.cancel(q-alpha/g**3-beta/g**2)
num,den=s.fraction(remaining)
need(s.Poly(den,z).eval(0)!=0,'scalar principal-part remainder regular at E1')
zero(s.limit(g**3*q,z,0)-alpha,'independent highest principal coefficient')
zero(s.limit(g**2*(q-alpha/g**3),z,0)-beta,'independent second principal coefficient')
zero(s.limit(g*(q-alpha/g**3-beta/g**2),z,0),'simple scalar coefficient vanishes')
P3=s.cancel(g**3*q)
need(s.denom(P3)==3 or s.denom(P3)==1,'order-three witness polynomial over coefficient field')
zero(D(P3)-g**3,'upper torsion order witness')
zero(s.limit(z*g**2*q,z,0)-4*a*h*h/3,'g squared primitive has a genuine simple E1 pole')

# Full affine E2 is k[z,z^-1]; w=1/z adds its unique DG boundary point.
zero(K.subs({z:1/w,t:2*h*w**3-w*w}),'E2 parametrization')
zero(q.subs(z,1/w)-w**3/(3*a),'primitive regular on E2 and at boundary w=0')
Finf=c-3*a*h*h+4*a*h**3*r-a*h**4*r*r-a*b*(1-h*r)**4
zero(F.subs({z:1/r-h,t:-r*r-r**4*b})-Finf,'global second chart')
zero(s.diff(Finf,b).subs(r,0)+a,'no boundary critical point')
bi=-2*h**5*w**3-7*h**4*w*w-8*h**3*w-3*h*h
zero(Finf.subs({r:w/(1+h*w),b:bi})-c,'complete E2 boundary chart')
zero(Finf.subs({r:0,b:-3*h*h})-c,'actual added E2 point')

# An explicit polynomial Bezout field independently pays the connection input.
Az=-(1+z/h+z*z/h**2+z**3/h**3+2*z**3*t/h)/(2*a*h)
At=s.cancel((1-Az*s.diff(F,z))/(a*z**4))
need(s.denom(At).free_symbols <= {a,h},'second Bezout coefficient polynomial in z,t')
zero(Az*s.diff(F,z)+At*s.diff(F,t)-1,'explicit V(F)=1')
def V(f):return Az*s.diff(f,z)+At*s.diff(f,t)
m=s.diff(Az,z)+s.diff(At,t)
zero(D(V(q))-m,'canonical derivative primitive gives div V')
next_part=-3*alpha/g**4-2*beta/g**3
vrem=s.cancel(V(q)-next_part)
need(s.Poly(s.denom(vrem),z).eval(0)!=0,'connection differentiates entire principal part')
for j in range(8):
    actual=s.diff(alpha/u**3+beta/u**2,u,j)
    expected=(-1)**j*(alpha*s.factorial(j+2)/(2*u**(j+3))+beta*s.factorial(j+1)/u**(j+2))
    zero(actual-expected,'canonical response coefficient formula')
    need(s.limit(u**(j+3)*actual,u,0)!=0,'order j+3 leading coefficient nonzero')

# Hostile boundary: dropping h!=0 changes reducedness and the order bound.
F0=F.subs(h,0);g0=g.subs(h,0);q0=q
P2=s.cancel(g0*g0*q0)
need(s.denom(P2)==3 or s.denom(P2)==1,'singular h=0 boundary already has order at most two')
zero(s.diff(F0,z)*s.diff(P2,t)-s.diff(F0,t)*s.diff(P2,z)-g0*g0,'boundary lower annihilator witness')
zero(s.diff(F0,z).subs(z,0),'boundary gradient actually vanishes')
zero(s.diff(F0,t).subs(z,0),'boundary repeated component actually critical')

print('PASS: fixed-F hypotheses, all affine poles, two component labels, exact primary order three.')
print('Canonical derivative j has exact primary order j+3; full scalar parts retained.')
print('Hostile h=0: repeated critical component admits an order-at-most-two witness; excluded parameter is necessary.')
print(f'Always-active exact gates: {count}')
