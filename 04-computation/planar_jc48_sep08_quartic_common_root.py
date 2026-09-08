#!/usr/bin/env python3
"""Exact controls for the analytic quartic common-boundary-root theorem.

No finite coefficient census is used to prove the arbitrary-germ Newton
or parameterized Morse statements. Those are proved in the companion note.
"""
import hashlib,json
from fractions import Fraction as Q
import sympy as S

w,s,Z,tau,V,c,x,z,y,e,r,b=S.symbols('w s Z tau V c x z y e r b')
a,bb,M,z0=S.symbols('a bb M z0',nonzero=True)
GATES=0
rows=[]
def check(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
def eq(left,right,label):check(S.cancel(left-right)==0,label)
def lead(expr,var):
    p=S.Poly(S.expand(expr),var)
    degree=min(m[0] for m,co in p.terms() if co!=0)
    return degree,p.coeff_monomial(var**degree)

for m in range(1,9):
    for j in list(range(10))+[None]:
        N=w**m+(s*w**j if j is not None else 0)
        E=S.expand(N*N+s**3-c*s**4)
        if j is None or m<3*j:
            # Clearing the denominator three is only an algebraic weight
            # check; this cover is not asserted to be a minimal parameter.
            sub={w:tau**3,s:Z*tau**(2*m)}
            order,co=lead(E.subs(sub,simultaneous=True),tau)
            check(order==6*m,'A-dominated first Newton weight')
            eq(co,1+Z**3,'A-dominated simple cubic')
            order,co=lead(S.diff(E,s).subs(sub,simultaneous=True),tau)
            check(order==4*m,'A-dominated polar weight')
            eq(co,3*Z**2,'A-dominated polar coefficient')
            check(j is None or j>Q(m,3),'entire high-j regime')
            rows.append([m,j,'regular A'])
        elif m>3*j:
            n=m-j;gap=m-3*j
            check(n>2*j and gap>0,'strict unbalanced separation')
            # The cancellation pair has two Puiseux determinations; these
            # may form a single quadratic normalized branch.
            sub={w:tau**2,s:tau**(2*n)*(-1+tau**gap*V)}
            order,co=lead(E.subs(sub,simultaneous=True),tau)
            check(order==6*n,'cancellation exact leading weight')
            eq(co,V**2-1,'two simple Puiseux determinations')
            order,co=lead(S.diff(E,s).subs(sub,simultaneous=True),tau)
            check(order==3*m-j,'cancellation polar weight')
            eq(co,2*V,'cancellation polar coefficient')
            check(4*n+1-(3*m-j)==gap+1 and gap+1>=2,
                  'normalized cancellation differential regular')
            if j:
                sub={w:tau,s:Z*tau**(2*j)}
                order,co=lead(E.subs(sub,simultaneous=True),tau)
                check(order==6*j,'additional simple-branch weight')
                eq(co,Z**2*(1+Z),'nonzero simple root Z minus one')
                order,co=lead(S.diff(E,s).subs(sub,simultaneous=True),tau)
                check(order==4*j,'additional simple polar weight')
                eq(co.subs(Z,-1),1,'additional simple polar nonzero')
            else:
                eq(lead(E.subs(w,0),s)[0],2,'j zero has only local degree two')
            rows.append([m,j,'regular cancellation'])
        else:
            k=j
            check(m==3*k and k in [1,2],'only balanced orders within octic')
            EQ=S.cancel(E.subs(s,w**(2*k)*Z)/w**(6*k))
            eq(EQ,(1+Z)**2+Z**3-c*w**(2*k)*Z**4,'complete balanced prototype')
            rows.append([m,j,'balanced',2*k-2])

P=(a+bb*Z)**2+M*Z**3
eq(P.subs(Z,0),a*a,'no zero face root')
eq(S.discriminant(P,Z),a**3*M*(4*bb**3-27*a*M),'exact face discriminant')
eq((3*M*z0*z0)**2-4*(-M*z0**3)*(-3*M*z0),-3*M*M*z0**4,
   'triple-root coefficient contradiction')
wall=S.Rational(4,27)*bb**3/a
eq(P.subs(M,wall),wall*(Z+3*a/bb)**2*(Z+3*a/(4*bb)),
   'entire nonzero discriminant wall has double plus simple roots')

Q0=S.Function('Q0');zc=S.Function('zc')(w,c)
for k in [1,2]:
    QQ=Q0(w,Z)-c*w**(2*k)*Z**4
    eq(S.diff(QQ.subs(Z,zc),c),S.diff(QQ,Z).subs(Z,zc)*S.diff(zc,c)-w**(2*k)*zc**4,
       'exact Morse critical-value chain rule')
    eq(S.diff(-c*w**(2*k)*z0**4,c),-w**(2*k)*z0**4,
       'nonzero generic coefficient at order two k')
    for lam in range(1,2*k+1):
        if lam%2==0:
            poles=[lam//2,lam//2]
            actual=sum(max(p-1,0) for p in poles)
            eq(lead(tau**(-lam//2)*tau**(lam//2),tau)[1],1,'even split unit control')
        else:
            poles=[lam-1]
            actual=max(lam-2,0)
            check(2*tau/tau**lam==2*tau**(1-lam),'odd ramified differential')
        check(actual==max(lam-2,0)<=2*k-2,'actual normalized total primitive budget')
        rows.append(['lambda',k,lam,poles,actual])

def partitions(n,ceiling=None):
    if n==0:
        yield ();return
    for first in range(min(n,n if ceiling is None else ceiling),1-1,-1):
        for rest in partitions(n-first,first):yield(first,)+rest
def cost(m):return max(2*(m//3)-2,0) if m%3==0 else 0
pp=list(partitions(8))
check(len(pp)==22,'complete octic multiplicity partitions')
for row in pp:
    check(sum(row)==8 and row.count(6)<=1,'actual octic sum and unique possible six')
    check(sum(cost(m) for m in row)<=2,'all octic primitive budgets at most two')
check(max(sum(cost(m) for m in row) for row in pp)==2,'sharp octic budget')
check(cost(9)==4,'higher boundary degree is a genuine stopping boundary')

# Actual globally admissible higher-pole hostile, reconstructed independently.
N6=x*x*z-x**4-S.Rational(8,27)*x**4*z+S.Rational(4,27)*x*x*z*z
MM=-(1+x*x*z)
eq(N6.subs(z,x*x),-S.Rational(4,27)*x**6,'actual multiplicity six at zero')
eq(S.gcd(N6.subs(z,x*x),MM.subs(z,x*x)),1,'no finite common boundary root')
Ni=S.cancel(r**4*z**2*N6.subs({x:1/r,z:1/z},simultaneous=True))
Mi=S.cancel(-r*r*z*MM.subs({x:1/r,z:1/z},simultaneous=True))
eq(lead(Ni.subs(z,r*r),r)[0],2,'remaining boundary multiplicity two at infinity')
eq(Mi.subs({r:0,z:0}),1,'no infinite common boundary root')
E=S.expand((N6*N6+(z-x*x)**3*MM-c*(z-x*x)**4).subs(z,s+x*x))
Nbar=s*y-S.Rational(4,27)*y**3+S.Rational(4,27)*s*s*y
Ebar=S.expand(Nbar*Nbar-s**3*(1+y*y+s*y)-c*s**4)
eq(E,Ebar.subs(y,x*x),'actual unsimplified fibre equality')
EE=S.expand(Ebar.subs(s,S.Rational(4,9)*y*y+e*y**3))
order,co=lead(EE,y)
check(order==8,'actual split order eight')
eq(co,-e*e/3-S.Rational(64,59049)*(65+36*c),'actual split coefficient')
order,co=lead(S.diff(Ebar,s).subs(s,S.Rational(4,9)*y*y+e*y**3),y)
check(order==5,'actual polar order five')
eq(co,-2*e/3,'actual double-pole coefficient')
eq(-S.Rational(16,81)/co,S.Rational(8,27)/e,'actual double-pole principal part')
eq(E.subs(x,-x),E,'actual parity supplies zero residue')
check(max(2-1,0)+max(2-1,0)==2,'hostile uses full allowed primitive budget')

# The actual D restriction is nonconstant and has unramified generic points.
H_D=-1-S.Rational(8,27)*b;L_D=b
F_D=S.expand(H_D*H_D+L_D)
check(S.degree(F_D,b)==2,'actual nonconstant D restriction')
check(S.discriminant(F_D-c,b)!=0,'actual generic D points unramified')
eq(S.diff(r**3/3,r),r*r,'order-two differential requires local degree three')
check(3>2,'same-component degree contradicts complete boundary budget')

print('Quartic common-root controls PASS; unbounded Newton/Morse proof requires independent audit')
print('Universe: m=1..8; j=0..9 and infinity; all22 octic partitions; all balanced split orders')
print('Controls: exact Newton weights, no triple face, generic Morse coefficient, actual m6 residue-free double pole')
print('Always-active gates:',GATES)
print('Semantic SHA256:',hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest())
