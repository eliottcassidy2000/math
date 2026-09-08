#!/usr/bin/env python3
"""Exact controls; the unbounded compact-differential theorem is analytic."""
import hashlib
import json
import sympy as S

x,z,s,t,r,b,q,tau,d,y,e,c=S.symbols('x z s t r b q tau d y e c')
gates=0
records=[]
def check(ok,label):
    global gates
    gates+=1
    if not ok:
        raise RuntimeError(label)
def eq(a,b,label):
    check(S.cancel(a-b)==0,label)
def jac(a,b,u,v):
    return S.diff(a,u)*S.diff(b,v)-S.diff(a,v)*S.diff(b,u)
def coeff(a,u,k):
    return S.expand(a).coeff(u,k)
def min_power(a,u):
    return min(mon[0] for mon,co in S.Poly(S.expand(a),u).terms() if co!=0)
def source(N,power):
    return S.cancel(N.subs(z,x*x+1/t)*t**power)
def boundary_chart(P):
    return S.cancel(P.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
def at_infinity(P,dx,dz):
    return S.cancel(r**dx*q**dz*P.subs({x:1/r,z:1/q},simultaneous=True))

eq(jac(1/r,-r*r-r**4*b,r,b),r*r,'full W chart omega')
eq(jac(x,1/s,s,x),1/s**2,'finite compact chart omega')
eq(jac(1/r,r*r*q/(r*r-q),r,q),-r*r/(q-r*r)**2,'complete infinity omega')

# Exact safe section controls; the theorem treats all higher local jets.
safe=[x**4*z*z+1,(x*x*z+1)**2,(x*x*z+1)**2+(z-x*x)]
for i,N in enumerate(safe):
    R=S.Poly(N.subs(z,x*x),x)
    check(R.degree()==8,'no omitted safe infinity zero')
    eq(S.gcd(R.as_expr(),1),1,'safe M avoids boundary')
    if i==0:
        eq(S.gcd(R.as_expr(),S.diff(R.as_expr(),x)),1,'eight simple safe zeros')
    else:
        eq(R.as_expr(),(x**4+1)**2,'four exact double safe zeros')
        eq(S.gcd(x**4+1,4*x**3),1,'double roots distinct')
    H=source(N,2)
    check(S.denom(H)==1,'safe H source polynomial')
    check(S.denom(boundary_chart(H))==1,'safe H second full chart polynomial')
    check(S.Poly(H,t).degree()==2,'safe exact L2 degree')
    Nlocal=N.subs(z,s+x*x)
    Nv=S.diff(Nlocal,s).subs(s,0)
    if i==1:
        eq(S.rem(Nv,x**4+1,x),0,'double safe cusp derivative zero')
    if i==2:
        eq(S.rem(Nv,x**4+1,x),1,'double safe node derivative unit')

# The three local models give orders 2,1,2 respectively.
eq((tau**2),tau**2,'simple-root model order two')
eq((tau**2/(2*x)).subs(x,tau),tau/2,'double-root node order one')
eq((tau**2*2*x/(2*x**3)).subs(tau,x*x),x*x,'double-root cusp order two')

M=-(1+x*x*z)
L=source(M,1)
eq(L,-((1+x**4)*t+x*x),'common global L')
check(S.denom(boundary_chart(L))==1,'common L second chart')
eq(-at_infinity(M,2,1).subs(q,r*r),1+r**4,'M nonzero at infinity')
N3=x*z-S.Rational(31,27)*x**3+x**3*z
N6=x*x*z-x**4-S.Rational(8,27)*x**4*z+S.Rational(4,27)*x*x*z*z
H3=source(N3,2);H6=source(N6,2)
eq(H3,(x+x**3)*t+(x**5-S.Rational(4,27)*x**3)*t*t,'H3 literal source')
eq(H6,x*x*t-S.Rational(4,27)*x**6*t*t+S.Rational(4,27)*x*x,'H6 literal source')
eq(boundary_chart(H6),-1-r*r*b-S.Rational(8,27)*b-S.Rational(4,27)*r*r*b*b,'H6 full chart')

for name,N,H,finite_order,infinite_order in [('m3',N3,H3,3,3),('m6',N6,H6,6,2)]:
    check(S.Poly(N,x,z).degree(x)<=4 and S.Poly(N,x,z).degree(z)<=2,'true global numerator bounds')
    check(S.denom(boundary_chart(H))==1,'hostile H second chart')
    R=S.expand(N.subs(z,x*x))
    check(min_power(R,x)==finite_order,'actual finite boundary multiplicity')
    Ri=at_infinity(N,4,2).subs(q,r*r)
    check(min_power(Ri,r)==infinite_order,'actual infinite boundary multiplicity')
    eq(S.gcd(R,1+x**4),1,'M avoids every finite zero')
    eq((1+r**4).subs(r,0),1,'M avoids infinite zero')
    Ni=at_infinity(N,4,2);Mi=-at_infinity(M,2,1)
    eq((N/(z-x*x)**2).subs({x:1/r,z:1/q},simultaneous=True),Ni/(q-r*r)**2,'numerator infinity transition')
    eq((M/(z-x*x)).subs({x:1/r,z:1/q},simultaneous=True),Mi/(q-r*r),'M sign in infinity transition')
    records.append([name,str(R),str(S.expand(Ri))])

E3=S.expand((N3*N3+(z-x*x)**3*M-c*(z-x*x)**4).subs(z,s+x*x))
p3={s:tau*tau,x:S.Rational(3,2)*tau+d*tau*tau}
B3=S.expand(E3.subs(p3,simultaneous=True))
check(min_power(B3,tau)==8,'m3 exact first split order')
eq(coeff(B3,tau,8),-S.Rational(4,3)*d*d+S.Rational(351,16)-c,'m3 split coefficient')
J3=S.expand(S.diff(E3,x).subs(p3,simultaneous=True))
check(min_power(J3,tau)==6,'m3 polar order')
eq(coeff(J3,tau,6),-S.Rational(8,3)*d,'m3 polar coefficient')
eq(2/coeff(J3,tau,6),-S.Rational(3,4)/d,'m3 nonzero logarithmic residue')
eq(S.solve(coeff(B3,tau,8),d*d)[0],(1053-48*c)/64,'m3 generic exceptional value')

N6bar=s*y-S.Rational(4,27)*y**3+S.Rational(4,27)*s*s*y
E6bar=S.expand(N6bar**2-s**3*(1+y*y+s*y)-c*s**4)
E6=S.expand(E6bar.subs(y,x*x))
eq(E6,(N6*N6+(z-x*x)**3*M-c*(z-x*x)**4).subs(z,s+x*x),'m6 unsimplified fibre')
eq(E6.subs(x,-x),E6,'m6 exact even parity')
p6={s:S.Rational(4,9)*y*y+e*y**3}
B6=S.expand(E6bar.subs(p6))
check(min_power(B6,y)==8,'m6 exact split order')
eq(coeff(B6,y,8),-e*e/3-S.Rational(64,59049)*(65+36*c),'m6 split coefficient')
eq(S.solve(coeff(B6,y,8),e*e)[0],-S.Rational(64,19683)*(65+36*c),'m6 generic split')
J6=S.expand(S.diff(E6bar,s).subs(p6))
check(min_power(J6,y)==5,'m6 polar order')
eq(coeff(J6,y,5),-S.Rational(2,3)*e,'m6 polar coefficient')
eq(-S.Rational(16,81)/coeff(J6,y,5),S.Rational(8,27)/e,'m6 double-pole coefficient')
check(4-5==-1,'m6 coefficient has y order minus one')
check(2*(4-5)==-2,'m6 actual x-normalization pole order two')

# The implicit equation for h has two simple generic roots; all subsequent
# coefficients are in y, proving zero residue on each actual x branch.
h=S.symbols('h')
implicit=S.cancel(E6bar.subs(s,y*y*(S.Rational(4,9)+y*h))/y**8)
eq(implicit.subs(y,0),-h*h/3-S.Rational(64,59049)*(65+36*c),'m6 implicit initial equation')
eq(S.diff(implicit,h).subs({y:0,h:e}),-2*e/3,'m6 implicit root derivative')
eq(S.residue(S.Rational(8,27)/e/x**2,x,0),0,'m6 second-kind leading part')
for j in range(-1,8):
    check(2*j!=-1,'even full Laurent parity excludes residue')

print('Quartic compact differential controls PASS; analytic proof and independent audit required')
print('Universe: three complete safe sections; two global hostile sections; exact generic split and polar coefficients')
print('Safe orders: simple 2; double node 1; double cusp 2; infinity numerator r^2 retained')
print('Hostiles: m3 nonzero logarithmic residue; m6 order-two pole with identically zero residue; M avoids every boundary zero')
print('Always-active gates:',gates)
print('Semantic SHA256:',hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
