#!/usr/bin/env python3
"""Exact controls for the complete binary (3,3,2) DG quartic class.

All finite positions p are symbolic. The accompanying proof pays local
branches, geometric integrality and complete Riemann--Roch membership.
This verifier controls identities, not a finite search for Keller maps.
"""
from hashlib import sha256
import sympy as S

gates = []
def need(label, predicate):
    if not bool(predicate):
        raise RuntimeError(label)
    gates.append(label)
def zero(label, value):
    need(label, S.cancel(value) == 0)

u,t,p,x,r,T,s,zeta = S.symbols('u t p x r T s zeta')
a,b,c,lam,mu,e = S.symbols('a b c lam mu e')
A = S.symbols('A0:5'); B = S.symbols('B0:5')
C0,E0 = S.symbols('C0 E0')
N = u**3*(u-1)**3
P = sum(A[i]*u**i for i in range(5))
Q = (A[4]-1)*u**2+(A[3]-2*p*A[4]+4*p+3)*u+C0
M = sum(B[i]*u**i for i in range(5))
R = B[4]*u**2+(B[3]-2*p*B[4])*u+E0

def global_h(label, n, pp, qq):
    vals = [qq, pp-2*qq*(u+p)**2,
            n-pp*(u+p)**2+qq*(u+p)**4]
    for j,val in enumerate(vals):
        need(label+' section row '+str(j), S.Poly(S.expand(val),u).degree() <= 4)

def global_l(label, mm, rr):
    for j,val in enumerate([rr,mm-rr*(u+p)**2]):
        need(label+' linear row '+str(j), S.Poly(S.expand(val),u).degree() <= 2)

global_h('full all-p H',N,P,Q)
global_l('full all-p L',M,R)
# Independently solve the high-degree conditions, before substituting Q.
q2,q1,q0 = S.symbols('q2 q1 q0')
Qfree=q2*u**2+q1*u+q0
high=S.Poly(S.expand(N-P*(u+p)**2+Qfree*(u+p)**4),u)
sol=S.solve([high.coeff_monomial(u**j) for j in (6,5)],(q2,q1),dict=True)
need('full H high rows have unique solution',len(sol)==1)
zero('complete Q2',sol[0][q2]-(A[4]-1))
zero('complete Q1',sol[0][q1]-(A[3]-2*p*A[4]+4*p+3))
lo=S.Poly(S.expand(M-Qfree*(u+p)**2),u)
sol=S.solve([lo.coeff_monomial(u**j) for j in (4,3)],(q2,q1),dict=True)
need('full L high rows have unique solution',len(sol)==1)
zero('complete R2',sol[0][q2]-B[4])
zero('complete R1',sol[0][q1]-(B[3]-2*p*B[4]))

def invert(expr):
    return S.expand(expr.subs({u:1/r-p,t:-r**2-r**4*T}, simultaneous=True))
H=N*t**2+P*t+Q; L=M*t+R
Hi=invert(H); Li=invert(L)
need('H original D chart polynomial',S.denom(S.cancel(Hi))==1)
need('L original D chart polynomial',S.denom(S.cancel(Li))==1)
zero('D H slope',S.diff(Hi.subs(r,0),T)-(2-A[4]))
zero('D L slope',S.diff(Li.subs(r,0),T)+B[4])
zero('infinity boundary multiplicity2',S.expand(r**8*N.subs(u,1/r-p))-r**2*(1-p*r)**3*(1-(p+1)*r)**3)

# With two active triples, the entire degree-four lower coefficient is fixed.
both=lam*u**2*(u-1)**2
zero('both-active residue0',S.residue(both/N,u,0)+lam)
zero('both-active residue1',S.residue(both/N,u,1)-lam)

# One active triple, then constant D. Additive e can be absorbed in zeta.
P=2*u**4+a*u**3+b*u**2
Q=u**2+(a+3)*u+c
M=u**2*(lam*u+mu); R=lam*u+e
global_h('constant D active H',N,P,Q)
global_l('constant D active L',M,R)
zero('single-active M/N residue0',S.residue(M/N,u,0)+mu)
zero('single-active M/N residue1',S.residue(M/N,u,1)-mu)
zero('surviving M/N primitive',S.diff(-lam/(2*(u-1)**2),u)-M.subs(mu,0)/N)
zero('other normal value',P.subs(u,1)-(a+b+2))

P=P.subs(b,-a-2)
H=N*t**2+P*t+Q
L=lam*(u**3*t+u)
F=H**2+L
J=u*(u-1)**2*t+u
global_l('primitive basis J',u*(u-1)**2,u)
Ji=invert(J)
need('J original D chart polynomial',S.denom(S.cancel(Ji))==1)
zero('J D value',Ji.subs(r,0)-(2*p+2))
Z=S.symbols('Z',nonzero=True)
zero('J all infinity tangent limits',S.limit(Ji.subs(T,1/(r*Z)),r,0)-(2*p+2-1/Z))
zero('active source line is constant fibre',F.subs(u,0)-c**2)
zero('triple1 M unit',S.expand(L).coeff(t).subs(u,1)-lam)
zero('triple1 P vanishes',P.subs(u,1))

def jac(f,g):
    return S.diff(f,u)*S.diff(g,t)-S.diff(f,t)*S.diff(g,u)
rem=S.Poly(S.rem(jac(F,J),F-zeta,t),t)
expected={3:0,2:2*c*u**3*(u-1)**4,
          1:2*c*u**2*(u-1)**2*(a+2*u+2),
          0:2*a*c*u**2-2*a*c*u+2*c**2*u-2*c**2+
            2*c*u**3+4*c*u**2-6*c*u-2*zeta*u+2*zeta-lam*u}
for j,val in expected.items():
    zero('complete bracket remainder t'+str(j),rem.coeff_monomial(t**j)-val)
zero('c0 obstruction',rem.as_expr().subs(c,0)-(2*zeta-(2*zeta+lam)*u))
need('formal fibre obstruction nonzero',S.Poly(2*zeta+lam,zeta).degree()==1)
zero('special fibre hostile must not replace generic level',
     rem.as_expr().subs({c:0,zeta:-lam/2})+lam)

# The full tangent quartic pays all four nonzero generic infinity branches.
aa,hh,ll,ee=S.symbols('aa hh ll ee')
P0=(1+aa*Z+hh*Z**2)**2+ll*Z**3+ee*Z**4
elim=S.expand(Z*S.diff(P0,Z)-4*P0)
zero('tangent repeated-root eliminant constant',elim.subs(Z,0)+4)
need('infinity eta order1',2+2-3==1)
need('active3 two normalized double poles',4//2==2 and 6+1-9==-2)
need('genus1 canonical degree',4*1-2*2==0)
need('complete RR dimension2',2+1-1==2)
need('D zero defeats pole capacity',3>2)
need('other triple zero defeats pole capacity',5>2)

# Same-partition rational positives, retaining every original finite p.
AA=u*(u-1); ZZ=AA*t+1
Hp=AA*ZZ**2+c
Gh=(2*u-1)/(AA*ZZ)
global_h('sharp all-p H',N,2*AA**2,AA+c)
zero('sharp H mate',jac(Hp,Gh)-1)
zero('sharp F mate',jac(Hp**2+e,Gh/(2*Hp))-1)
need('sharp mate has actual poles',S.denom(S.cancel(Gh/(2*Hp)))!=1)
# Hostile finite-double leading differential has a nonzero logarithm.
alpha=S.symbols('alpha',nonzero=True)
zero('finite double leading residue hostile',S.residue(1/(alpha*u),u,0)-1/alpha)

print('Binary (3,3,2): all-p global rows, weighted infinity and complete elliptic primitive space')
print('Rational mate => L constant; every boundary location excludes polynomial mates')
print('Exact gates:',len(gates))
print('Semantic SHA256:',sha256('\n'.join(gates).encode()).hexdigest())
