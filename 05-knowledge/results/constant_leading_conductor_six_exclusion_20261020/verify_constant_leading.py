#!/usr/bin/env python3
"""Exact controls; proofs and inherited scope are in the accompanying note."""
import hashlib
import json
from pathlib import Path
import sympy as S

x,y,t=S.symbols('x y t')
checks=[]
def check(name, condition, detail=None):
    ok=bool(condition)
    checks.append({'name':name,'passed':ok,**({'detail':detail} if detail is not None else {})})
    if not ok: raise RuntimeError('FAILED: '+name)
def zero(name, expr): check(name,S.expand(expr)==0)
def jac(p,q,var=y): return S.diff(p,x)*S.diff(q,var)-S.diff(p,var)*S.diff(q,x)

# Independent formal coefficient control for all degree pairs in the stated box.
# The all-degree proof is the degree estimate, not this finite universe.
for n in range(1,8):
    for m in range(1,9):
        pp=S.symbols(f'p0:{n}')
        px=S.symbols(f'px0:{n}')
        qq=S.symbols(f'q0:{m+1}')
        qx=S.symbols(f'qx0:{m+1}')
        a=S.Symbol('a')
        p=a*y**n+sum(pp[i]*y**i for i in range(n))
        q=sum(qq[i]*y**i for i in range(m+1))
        p_x=sum(px[i]*y**i for i in range(n))
        q_x=sum(qx[i]*y**i for i in range(m+1))
        top=S.Poly(p_x*S.diff(q,y)-S.diff(p,y)*q_x,y).coeff_monomial(y**(n+m-1))
        zero(f'top coefficient n={n},m={m}',top+n*a*qx[m])

# Target-reduction positive control: a genuine polynomial automorphism.
p=x+y*y
q=y+p*p
zero('positive Keller control',jac(p,q)-1)
q1=S.expand(q-p*p)
zero('first constant target degree reduction',q1-y)
p1=S.expand(p-q1*q1)
zero('second reduction reaches base coordinate',p1-x)
zero('target reduction preserves Jacobian',jac(p1,q1)-1)
check('terminal constant coordinate gives full source',p1==x and q1==y)

# Exact semigroup count in the frozen reduced basis argument.
for degs,delta in [([0,6,11,17],7),([0,6,9,15],6),([0,6,7,13],5),([0,6,3,9],3)]:
    check('residue basis delta '+str(degs),len({v%4 for v in degs})==4 and sum(v-v%4 for v in degs)//4==delta)

# Universal exceptional class identities over Q[a,c,t].
a,c,beta,gamma,b=S.symbols('a c beta gamma b')
h=t**3+a*t
k=t**4+b*t*t+c*t
p=k+beta*h
q=h*h+gamma*h
r=S.expand(q*q-p**3)
zero('nonzero beta gives degree eleven',S.Poly(r,t).coeff_monomial(t**11)+3*beta)
p=k
r=S.expand(q*q-p**3-(4*a-3*b)*p*q)
zero('wrong gamma gives degree nine',S.Poly(r,t).coeff_monomial(t**9)-(2*gamma-3*c))
q=S.expand(h*h+S.Rational(3,2)*c*h)
r1=S.expand(q*q-k**3-(4*a-3*b)*k*q+2*a*(a-b)*k*k)
zero('wrong exceptional constraint gives degree seven',S.Poly(r1,t).coeff_monomial(t**7)+c*(2*a-3*b)/2)
k=S.expand(k.subs(b,2*a/3))
recovery=q*q-k**3-2*a*k*q+S.Rational(2,3)*a*a*k*k+(2*a**3/S.Integer(27)+3*c*c/S.Integer(4))*q-(a**4/S.Integer(9)+a*c*c)*k
zero('universal polynomial cubic recovery',recovery-c**3*h/8)
zero('quartic derivative at origin',S.diff(k,t).subs(t,0)-c)
zero('sextic derivative at origin',S.diff(q,t).subs(t,0)-3*a*c/2)
zero('quartic ramification at c=0',S.diff(k,t).subs({t:0,c:0}))
zero('sextic ramification at c=0',S.diff(q,t).subs({t:0,c:0}))

# Coefficient alignment uses no division by a nonconstant polynomial.
s,d,eta,e=S.symbols('s d eta e')
p_y=S.expand((k+d).subs(t,y+s))
q_y=S.expand((q+eta*k+e).subs(t,y+s))
zero('quartic cubic coefficient pins source translation',S.Poly(p_y,y).coeff_monomial(y**3)-4*s)
p_t=S.expand(p_y.subs(y,t-s))
q_t=S.expand(q_y.subs(y,t-s))
zero('quartic quadratic coefficient determines a',S.Poly(p_t,t).coeff_monomial(t*t)-2*a/3)
zero('quartic linear coefficient determines c',S.Poly(p_t,t).coeff_monomial(t)-c)
zero('quartic constant coefficient determines d',p_t.subs(t,0)-d)
zero('sextic quartic coefficient determines eta',S.Poly(q_t,t).coeff_monomial(t**4)-2*a-eta)
zero('sextic constant coefficient determines e',q_t.subs(t,0)-e)

# The predecessor hostile fails the new whole-identity leading gate.
P=y**4; Q=x*y**3+y*y+y
zero('hostile Jacobian is minus four y-six',jac(P,Q)+4*y**6)
check('hostile cubic leading coefficient is not constant',S.diff(S.Poly(Q,y).LC(),x)==1)

# A first jet still admits Jacobian one at the special node.
A=4*y*y-2*y+1; B=S.Integer(2)
p0=y**4; q0=y*y+y
zero('Bezout for immersed nodal fiber',A*S.diff(q0,y)-S.diff(p0,y)*B-1)
P=p0+x*A; Q=q0+x*B
zero('first jet exact Jacobian residual',jac(P,Q)-(1+x*(4-16*y)))
zero('first jet is Keller on the selected fiber',jac(P,Q).subs(x,0)-1)
check('first jet is not an entire Keller pair',S.Poly(jac(P,Q),x,y).total_degree()>0)

# Recursively construct eight orders of the unrestricted formal lift.
N=8
ps=[p0]; qs=[q0]
for n in range(1,N+1):
    E=S.Integer(0)
    # coefficient x^(n-1) from all already-defined pairs (i,j), i+j=n.
    for i in range(1,n):
        j=n-i
        E += i*ps[i]*S.diff(qs[j],y)-j*S.diff(ps[i],y)*qs[j]
    residual=S.expand((S.Integer(1 if n==1 else 0)-E)/n)
    ps.append(S.expand(A*residual))
    qs.append(S.expand(B*residual))
    zero('formal lift order '+str(n),n*(ps[n]*S.diff(q0,y)-S.diff(p0,y)*qs[n])+E-(1 if n==1 else 0))
P=sum(ps[i]*x**i for i in range(N+1)); Q=sum(qs[i]*x**i for i in range(N+1))
J=S.Poly(S.expand(jac(P,Q)-1),x)
for n in range(N): zero('independent full formal lift coefficient '+str(n),J.coeff_monomial(x**n))
formal_degrees=[{'order':i,'P_y_degree':int(S.degree(ps[i],y)) if ps[i]!=0 else None,'Q_y_degree':int(S.degree(qs[i],y)) if qs[i]!=0 else None} for i in range(N+1)]
check('formal control actually increases y degree',max(z['P_y_degree'] or -1 for z in formal_degrees)>4)

# Independent closed-form construction; no recursion used for its coefficients.
C=S.expand(A*S.diff(B,y)-S.diff(A,y)*B)
zero('nodal formal deformation coefficient',C-(4-16*y))
Fs=[]
for n in range(1,N+1):
    fn=S.expand(S.binomial(S.Rational(1,2),n)*2**n*C**(n-1))
    Fs.append(fn)
    zero('closed-form agrees with recurrence P '+str(n),ps[n]-A*fn)
    zero('closed-form agrees with recurrence Q '+str(n),qs[n]-B*fn)
F=sum(Fs[n-1]*x**n for n in range(1,N+1))
implicit=S.Poly(S.expand(F+C*F*F/2-x),x)
for n in range(N+1): zero('closed-form implicit equation order '+str(n),implicit.coeff_monomial(x**n))
aa,bb,ay,by,py,qy,fx,fy,ff=S.symbols('aa bb ay by py qy fx fy ff')
Jform=aa*fx*(qy+by*ff+bb*fy)-(py+ay*ff+aa*fy)*bb*fx
zero('universal formal Jacobian factorization',Jform-fx*(aa*qy-bb*py+(aa*by-ay*bb)*ff))

# Actual source shear invariant on a moving exceptional family.
ht=y**3+x*y; kt=y**4+S.Rational(2,3)*x*y*y+y; wt=ht*ht+S.Rational(3,2)*ht
shift=x*x+x
zero('polynomial source shear preserves Jacobian',jac(kt.subs(y,y+shift),wt.subs(y,y+shift))-jac(kt,wt).subs(y,y+shift))

result={
 'status':'ALL EXACT CONTROLS PASS',
 'checks':checks,
 'check_count':len(checks),
 'universes':{'top_coefficient_degree_box':{'n':[1,7],'m':[1,8]},'formal_lift_orders':N},
 'formal_lift_degrees':formal_degrees,
 'scope':'Finite symbolic controls support, but do not replace, the universal proofs and the pinned characteristic-sequence input.',
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
print(json.dumps(result,indent=2,sort_keys=True))
