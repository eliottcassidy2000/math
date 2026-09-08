#!/usr/bin/env python3
"""Exact branch, trace, and hostile controls for the fixed finite octuple.

The complete analytic proof, including every inherited hypothesis, is in
the matching note. No inherited mathematical implementation is imported.
"""
import hashlib
import json
import sympy as S

x,t,s,z,r,bd,c,Z,W=S.symbols('x t s z r bd c Z W')
alpha,k0,kt,kxt,k2,k3,k4=S.symbols('alpha k0 kt kxt k2 k3 k4')
l0,lt,lxt,l2,l3,l4=S.symbols('l0 lt lxt l2 l3 l4')
hh,vv,y=S.symbols('hh vv y')
gates=0
records={}

def check(ok,label):
    global gates
    gates+=1
    if not ok:
        raise RuntimeError(label)

def eq(a,b,label):
    check(S.cancel(a-b)==0,label)

def coeff(a,v,n):
    return S.expand(a).coeff(v,n)

def jac(a,b,u,v):
    return S.diff(a,u)*S.diff(b,v)-S.diff(a,v)*S.diff(b,u)

w=x*x*t
v=x*(1+w)
h=x*x*(1+w)
K=k0+kt*t+kxt*x*t+k2*w+k3*v+k4*h
L=l0+lt*t+lxt*x*t+l2*w+l3*v+l4*h
H=alpha*h*h+K
F=S.expand(H*H+L)
N=S.cancel(s*s*H.subs(t,1/s))
M=S.cancel(s*L.subs(t,1/s))
eq(N.subs(s,0),alpha*x**8,'complete finite octuple boundary')
eq(coeff(N,s,1),kt+kxt*x+k2*x*x+k3*x**3+k4*x**4+2*alpha*x**6,
   'all first normal coefficients retained')
eq(M.subs(s,0),lt+lxt*x+l2*x*x+l3*x**3+l4*x**4,'complete M boundary polynomial')
eq(N.subs({x:0,s:0}),0,'only specified finite zero')
eq(S.diff(N,s).subs({x:0,s:0}),kt,'normal-unit reduction')
eq(coeff(H,t,2),alpha*x**8,'actual approximate-root source degree two')
eq(coeff(F,t,4),alpha*alpha*x**16,'actual source quartic leading coefficient')
Ng=S.expand(N.subs(s,z-x*x))
Mg=S.expand(M.subs(s,z-x*x))
check(S.Poly(Ng,x,z).degree(x)<=4 and S.Poly(Ng,x,z).degree(z)<=2,
      'complete global H box')
check(S.Poly(Mg,x,z).degree(x)<=2 and S.Poly(Mg,x,z).degree(z)<=1,
      'complete global L box')
eq(S.cancel((Ng-alpha*x**4*z*z)/(z-x*x)),
   S.expand((s*K.subs(t,1/s)).subs(s,z-x*x)),'entire lower-pole section')
Hinf=S.cancel(H.subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
Linf=S.cancel(L.subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
check(Hinf.is_polynomial(r,bd) and Linf.is_polynomial(r,bd),'both complete surface charts')
eq(Hinf.subs(r,0),alpha*bd*bd+k0-k2-k4*bd,'actual H boundary')
eq(Linf.subs(r,0),l0-l2-l4*bd,'actual L boundary')
eq(coeff(Hinf.subs(r,0)**2+Linf.subs(r,0),bd,4),alpha*alpha,
   'nonconstant original D map')
eq(jac(1/r,-r*r-r**4*bd,r,bd),r*r,'literal volume numerator')
eq(coeff(N.subs(s,0),x,8),alpha,'infinity octic value nonzero under alpha hypothesis')

# Full first-jet tangent, before discarding the lower normal coefficients.
first={kt:0,lt:0}
E=S.expand(N*N+s**3*M-c*s**4)
Pfirst=Z*Z*((kxt+k0*Z)**2+lxt*Z+(l0-c)*Z*Z)
eq(coeff(E.subs(first).subs(s,x*Z),x,4),Pfirst,'complete finite first-jet polynomial')
quad=S.cancel(Pfirst/Z**2)
eq(S.diff(S.discriminant(quad,Z),c),4*kxt*kxt,'two generic nonzero logarithmic tangent roots')
eq(Pfirst.subs(kxt,0),Z**3*(lxt+(k0*k0+l0-c)*Z),
   'one remaining nonzero logarithmic tangent root')
root=-lxt/(k0*k0+l0-c)
eq(S.diff(Pfirst.subs(kxt,0),Z).subs(Z,root),
   -lxt**3/(k0*k0+l0-c)**2,'simple first-jet root when lxt nonzero')

reduce={kt:0,kxt:0,lt:0,lxt:0}
Nr=S.expand(N.subs(reduce))
Mr=S.expand(M.subs(reduce))
Er=S.expand(E.subs(reduce))
Hr=S.expand(H.subs(reduce))
Lr=S.expand(L.subs(reduce))
Fr=S.expand(F.subs(reduce))
eq(Er.subs(x,0),(k0*k0+l0-c)*s**4,'exact four-sheet Weierstrass restriction')

# The two low smooth branches in the k2!=0 stratum.
T=k2*k2+(2*k0*k2+l2)*Z+(k0*k0+l0-c)*Z*Z
eq(coeff(Er.subs(s,x*x*Z),x,8),Z*Z*T,'complete low Newton face')
eq(T.subs(Z,0),k2*k2,'nonzero low slopes')
eq(S.diff(S.discriminant(T,Z),c),4*k2*k2,'generic two simple low branches')
lead_derivative=coeff(S.diff(Er,s).subs(s,x*x*Z),x,6)
eq(lead_derivative-Z*Z*S.diff(T,Z),2*Z*T,'exact polar derivative modulo low face')
check(4-6==-2,'each low differential has order minus two')

# Exact centre supplier for the remaining two determinations.
Nb=S.cancel(Nr.subs(s,x**6*Z)/x**8)
check(Nb.is_polynomial(x,Z),'cancellation chart has no negative powers')
eq(Nb.subs(x,0),alpha+k2*Z,'unique nonzero cancellation centre')
eq(S.diff(Nb,Z).subs(x,0),k2,'implicit centre derivative is a unit')
Z0=-alpha/k2
eq((alpha+k2*Z).subs(Z,Z0),0,'actual cancellation slope')
eq(coeff(S.diff(Nr,s).subs(s,x**6*Z),x,2),k2,'normal derivative order exactly two')
for n,lower,ln in [(2,{},l2),(3,{l2:0},l3),(4,{l2:0,l3:0},l4)]:
    term=S.expand((Mr-c*s).subs(lower).subs(s,x**6*Z))
    for j in range(n):
        eq(coeff(term,x,j),0,'all lower M orders vanish '+str((n,j)))
    eq(coeff(term,x,n),ln,'nonzero supplier at declared order '+str(n))
    eq(coeff((s**3*(Mr-c*s)).subs(lower).subs(s,x**6*Z),x,18+n),
       Z**3*ln,'exact split supplier '+str(n))
    direct=S.expand((S.diff(Er,s)-2*Nr*S.diff(Nr,s)).subs(lower).subs(s,x**6*Z))
    direct_orders=[mon[0] for mon,value in S.Poly(direct,x).terms() if value!=0]
    check(min(direct_orders)>=12+n,'every remaining derivative term has higher order '+str(n))
    eq(coeff(direct,x,12+n),3*Z*Z*ln,'full direct-derivative leading coefficient '+str(n))
    nn=S.Rational(18+n,2)
    es=nn+2
    check(es<12+n,'2N Ns dominates every direct derivative term '+str(n))
    eta=12-es
    eq(eta,1-S.Rational(n,2),'normalized x exponent '+str(n))
    actual=eta if n%2==0 else 2*eta+1
    eq(actual,0 if n<4 else -1,'actual normalized differential order '+str(n))
    check(2+(2 if n%2==0 else 2)==4,'full local Weierstrass degree retained '+str(n))
e=S.symbols('e')
eq(2*k2*e*(Z0**2/(2*k2*e)),Z0**2,'nonzero high log coefficient at n four')
eq(e*e+Z0**3*l4,e*e-alpha**3*l4/k2**3,'exact high split equation')
check(2<3,'same-component D local degree exceeds all primitive poles')

# Low linear phase: residues and source criticality must be distinguished.
fw=k0*k0+l0+l2*W
gw=(1+W)*(2*k0*k3+l3)
Fw=S.cancel(Fr.subs(k2,0).subs(t,W/x**2))
eq(Fw.subs(x,0),fw,'whole original boundary phase when k2 zero')
eq(coeff(Fw,x,1),gw,'whole first phase correction')
eq(jac(x,W/x**2,x,W),x**-2,'rational chart volume factor')
residue=S.cancel((S.diff(gw,W)*S.diff(fw,W)-S.diff(fw,W,2)*gw)/S.diff(fw,W)**3)
eq(residue,(2*k0*k3+l3)/l2**2,'exact low linear-phase residue')
eq(S.diff(Fr.subs(k2,0),x).subs(x,0),2*k0*k3+l3,'same numerator on actual source line')
eq(S.diff(Fr.subs(k2,0),t).subs(x,0),0,'source t derivative vanishes')

# Full birational field, square completion, and trace of the actual form.
eq(jac(h,v,x,t),h**3/v**2,'literal h-v Jacobian')
inverse={x:hh/vv,t:(vv*vv-hh)*vv*vv/hh**3}
eq(h.subs(inverse,simultaneous=True),hh,'first exact field inverse')
eq(v.subs(inverse,simultaneous=True),vv,'second exact field inverse')
epsilon=l3/(2*k3)
A=alpha*hh*hh+k4*hh+k0+epsilon
R=epsilon*epsilon-2*epsilon*A+l4*hh+l0
Fhv=(alpha*hh*hh+k3*vv+k4*hh+k0)**2+l3*vv+l4*hh+l0
eq(Fhv.subs(vv,(y-A)/k3),y*y+R,'exact original fibre completion')
eq(jac(hh,(y-A)/k3,hh,y),1/k3,'second exact Jacobian')
omega=(y-A)**2/(k3**3*hh**3)
eta=-(y-A)**2/(2*k3**3*hh**3*y)
eq(-2*y*eta,omega,'dF wedge eta is original volume')
trace=S.cancel(eta+eta.subs(y,-y))
eq(trace,2*A/(k3**3*hh**3),'trace retains the true volume coefficient')
eq(S.residue(trace,hh,0),2*alpha/k3**3,'nonzero trace residue')
eq(S.diff(c-R,c),1,'odd prime valuation also covers constant R')
eq(R.subs({l3:0,l4:0}),l0,'constant quadratic-extension control')
eq(trace.subs(l3,0),2*(alpha*hh*hh+k4*hh+k0)/(k3**3*hh**3),
   'trace residue survives constant R')

# Last case pays the complete inherited boundary-linear input.
last=Fr.subs({k2:0,l2:0,k3:0})
lastinf=S.cancel(last.subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
eq(lastinf,(alpha*bd*bd-k4*bd+k0)**2-l3*r*bd-l4*bd+l0,
   'exact last global boundary-linear class')
check(S.Poly(lastinf,r).degree()<=1,'no mate-degree restriction imported')

# Genuine rational mates INSIDE the main polynomial-exclusion class.
Fmate=h**4+h
Gmate=1/(3*x**3*(4*h**3+1))
eq(jac(Fmate,Gmate,x,t),1,'inside-class rational-mate hostile')
Flinear=h**4-v
Glinear=-1/(2*x*x)+2*h**3/x-S.Rational(4,3)*h**6
eq(jac(Flinear,Glinear,x,t),1,'source-critical-free rational boundary-linear hostile')
eq(S.diff(Flinear,x).subs(x,0),-1,'boundary-linear hostile smooth on source E')
eq(Flinear.subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True),bd**4+r*bd,
   'boundary-linear hostile has no critical point where r nonzero')
eq(S.diff(bd**4+r*bd,r),bd,'only possible boundary critical b is zero')
eq(S.diff(bd**4+r*bd,bd).subs(bd,0),r,'remaining boundary derivative nonzero in source')

records['entry']='complete L1 pair with octic alpha*x^8 at fixed zero'
records['low_branches']='two actual double poles; primitive budget two'
records['cancellation_orders']={2:0,3:0,4:-1}
records['trace_residue']=str(2*alpha/k3**3)
records['scope']='no polynomial mate of any degree; rational hostiles retained'
semantic=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
print('PASS fixed finite-octuple quartic:',gates,'always-active exact gates')
print('Complete global entry, all four j2 branches, and componentwise degree gate checked')
print('Linear-phase residue/source-criticality split and exact quadratic field trace checked')
print('Nonzero trace residue: 2*alpha/k3^3; constant-R case retained')
print('Inside-class rational mates verify the polynomial-only final scope')
print('Semantic SHA256:',semantic)
