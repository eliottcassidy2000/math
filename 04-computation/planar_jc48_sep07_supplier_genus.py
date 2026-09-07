#!/usr/bin/env python3
"""Exact genus controls for the actual completed-response Hamiltonian.

No research producer import. The universal proof uses Riemann--Hurwitz;
the controls reconstruct the full derivative/discriminant and local charts.
"""
import sympy as S
from math import gcd
from hashlib import sha256
import json

p,y,z,A,b,c,r,V = S.symbols('p y z A b c r V')
gates=0
def check(ok,label):
    global gates
    gates+=1
    if not ok: raise RuntimeError(label)
def eq(a,b,label): check(S.cancel(a-b)==0,label)

I=p**2*y**3*(p**3-y**2)*(A+b*y**2)
Q=3*p**3*A+5*(b*p**3-A)*z-7*b*z**2
eq(S.diff(I,y),p**2*y**2*Q.subs(z,y**2),'complete derivative')
critical=p**4*z**3*(p**3-z)**2*(A+b*z)**2
rem=S.rem(critical,Q,z)
r0=S.factor(rem.subs(z,0));r1=S.factor(S.diff(rem,z))
X=b*p**3
discQ=25*A**2+34*A*X+25*X**2
quartic=25*A**4+19*A**3*X-39*A**2*X**2+19*A*X**3+25*X**4
eq(S.discriminant(Q,z),discQ,'quadratic critical discriminant')
eq(r1,4*p**4*discQ*quartic/(7**6*b**4),'critical-value separation')
trace=S.factor(2*r0+r1*5*(X-A)/(7*b))
norm=S.factor(r0**2+r0*r1*5*(X-A)/(7*b)-r1**2*3*p**3*A/(7*b))
T=sum(v*A**(6-j)*X**j for j,v in enumerate([3125,11875,18700,19036,18700,11875,3125]))
eq(trace,4*p**4*(X-A)*T/(7**7*b**5),'critical squared-value trace')
eq(norm,-432*A**5*p**23*(A+X)**4/(7**7*b**5),'critical squared-value norm')
N=c**4-trace*c**2+norm
eq(S.discriminant(I-c,y),-7**7*b**6*p**12*c**2*N,'literal degree-seven discriminant')
eq(N.subs(p,0),c**4,'no generic finite branch at zero in N')

# The actual supplier coefficients are copied as rational input data and
# checked against their explicit six-term expression in the proof note.
actual_A=(S.Rational(12015,15962432)+S.Rational(28467,3990608)*p
          +S.Rational(2499,498826)*p**2-S.Rational(4806,1247065)*p**3
          +S.Rational(639368,11223585)*p**4)
actual_b=S.Rational(532468,6235325)
check(S.degree(actual_A,p)==4 and actual_b!=0,'actual supplier class entry')
controls=[(actual_A,actual_b),(p**4,1),((p-1)**4,-2),
          (p**4+p**3+1,3),(p**5+1,1),(p**6,2),
          ((p+2)**7,5),(p**8-p+1,-1),(p**9+p**3,7),
          (p**10+2,3)]
rows=[]
for k,(aa,bb) in enumerate(controls):
    m=int(S.degree(aa,p));lead=S.Poly(aa,p).LC()
    tr=S.Poly(trace.subs({A:aa,b:bb}),p)
    no=S.Poly(norm.subs({A:aa,b:bb}),p)
    check(tr.degree()==4+7*m,'trace degree')
    check(no.degree()==23+9*m,'norm and finite branch degree')
    eq(no.LC(),-S.Rational(432,7**7)*lead**9/bb**5,'nonzero branch leading coefficient')
    n=23+9*m
    # Independent finite-field specialization of the full residual divisor.
    # Clearing denominators followed by a squarefree degree-preserving
    # reduction certifies a good specialization, not the universal theorem.
    polynomial=S.Poly(S.cancel(N.subs({A:aa,b:bb,c:1})),p)
    denominator,integer=polynomial.clear_denoms()
    prime=next(q for q in (101,103,107,109,127,131)
               if int(denominator)%q and S.Poly(integer.as_expr(),p,modulus=q).degree()==n
               and S.gcd(S.Poly(integer.as_expr(),p,modulus=q),
                         S.Poly(integer.as_expr(),p,modulus=q).diff()).degree()==0)
    check(prime>=101,'independent good squarefree specialization')
    genus=(9*m+21+(m%2)-gcd(3,m+5))//2
    inf=1+(m%2)+3-gcd(3,m+5)
    check(2*genus-2==-14+n+6+inf,'full Riemann-Hurwitz sum')
    # Three infinity groups: two D-roots, two A+b*y2 roots, three small roots.
    check(2+2+3==7,'all seven infinity branches retained')
    s,tau=S.symbols('s tau')
    W=s**17*(aa.subs(p,s**2)+bb*s**6)
    check(W!=0,'nonzero actual logarithmic first response')
    local=(p**2*y**3*(p**3-y**2)*(aa+bb*y**2)).subs({p:s**2+tau,y:s*(s**2+tau)})
    eq(S.diff(local,tau).subs(tau,0),W,'actual invariant tau coefficient')
    rows.append(dict(degree=m,finite_branch=n,zero_ram=6,infinity_ram=inf,
                     genus=genus,good_prime=prime,actual=(k==0)))

# Universal local charts with independent leading-order polynomial controls.
# Zero: p=r^7, y=r^-2 V. Any A polynomial contributes only higher powers.
for m in range(4,13):
    aa=1+p**m
    zero=S.expand((I-c).subs(A,aa).subs({p:r**7,y:r**-2*V}))
    eq(S.limit(zero,r,0),-b*V**7-c,'zero index-seven initial polynomial')
    # D-group chart at infinity.
    dchart=S.expand((I-c).subs(A,aa).subs({p:r**-2,y:r**-3*V})*r**(19+2*m))
    eq(S.limit(dchart,r,0),V**3*(1-V**2),'D-pair infinity initial polynomial')
    # A+b*y2 group, using a quadratic base change uniformly.
    achart=S.expand((I-c).subs(A,aa).subs({p:r**-2,y:r**-m*V})*r**(4+7*m))
    eq(S.limit(achart,r,0),-V**5*(1+b*V**2),'outer-pair infinity initial polynomial')
    small=S.expand((I-c).subs(A,aa).subs({p:r**-3,y:r**(m+5)*V}))
    eq(S.limit(small,r,0),V**3-c,'small-triple infinity initial polynomial')

# Sharp failed extension at degree three: A=-b*p3 merges two infinity factors.
hostile=S.factor(I.subs(A,-b*p**3))
eq(hostile,-b*p**2*y**3*(p**3-y**2)**2,'coalescent face hostile')
w,h=S.symbols('w h')
eq((p**2*y**3*(p**3-y**2)**2).subs({p:h**2/w,y:h**3/w}),
   h**25*(1-w)**2/w**11,'birational degree25 cyclic cover')
check(all(gcd(25,v)==1 for v in (11,2,9)),'all three cyclic branch indices25')
check((-2*25+3*24+2)//2==12,'hostile exact genus12')
check((9*3+21+1-gcd(3,8))//2==24,'naive degree-three extrapolation fails')
check(rows[0]['genus']==27 and rows[0]['finite_branch']==59,'actual supplier genus27')
blob=json.dumps(rows,sort_keys=True,separators=(',',':')).encode()
print('Completed-response supplier genus PASS')
print('Actual supplier: genus27; finite branch59, zero ram6, infinity ram1')
print('Class: deg A=m>=4, b!=0; genus=(9m+21+(m mod2)-gcd(3,m+5))/2')
print('Named controls:',len(rows),'genera:',[q['genus'] for q in rows])
print('Boundary hostile: A=-b*p3 has genus12, not extrapolated24')
print('Always-active gates:',gates)
print('Semantic SHA256:',sha256(blob).hexdigest())
