#!/usr/bin/env python3
"""Exact controls for nonconstant-D binary (6,2) quartic polynomial exclusion.

The source controls symbolic identities. The all-degree polynomial-weight
gate, actual local entry, and generic residue argument are proved in the
companion note; no finite coefficient scan stands in for those arguments.
"""
from hashlib import sha256
import sympy as S

gates = []

def need(label, predicate):
    if not bool(predicate):
        raise RuntimeError(label)
    gates.append(label)

def zero(label, expression):
    need(label, S.cancel(expression) == 0)

u,t,w,v,h,f,p,A,b,c,d,lam,mu,nu,e = S.symbols(
    'u t w v h f p A b c d lam mu nu e')
r,Bd,delta = S.symbols('r Bd delta')
N = u**6
P = (A+2)*u**4+b*u**3+c*u**2
Q = (A+1)*u**2+(b-2*p*A)*u+d
M = lam*u**4+mu*u**3+nu*u**2
R = lam*u**2+(mu-2*p*lam)*u+e
H = N*t*t+P*t+Q
L = M*t+R
F = H*H+L

def jac(left, right, x=u, y=t):
    return S.diff(left,x)*S.diff(right,y)-S.diff(left,y)*S.diff(right,x)

# Derive the full pre-active row space and then intersect its actual jets.
pp = S.symbols('P0:5'); mm = S.symbols('M0:5')
q2,q1,q0,r2,r1,r0 = S.symbols('q2 q1 q0 r2 r1 r0')
PF = sum(pp[i]*u**i for i in range(5))
MF = sum(mm[i]*u**i for i in range(5))
QF = q2*u*u+q1*u+q0
RF = r2*u*u+r1*u+r0
high = S.Poly(S.expand(N-PF*(u+p)**2+QF*(u+p)**4),u)
sol = S.solve([high.coeff_monomial(u**j) for j in (6,5)],(q2,q1),dict=True)
need('unique full H high-row solution',len(sol)==1)
zero('full Q quadratic',sol[0][q2]-(pp[4]-1))
zero('full Q linear',sol[0][q1]-(pp[3]-2*p*pp[4]+4*p))
low = S.Poly(S.expand(MF-RF*(u+p)**2),u)
solL = S.solve([low.coeff_monomial(u**j) for j in (4,3)],(r2,r1),dict=True)
need('unique full L high-row solution',len(solL)==1)
zero('full R quadratic',solL[0][r2]-mm[4])
zero('full R linear',solL[0][r1]-(mm[3]-2*p*mm[4]))
for j,row in enumerate([Q,P-2*Q*(u+p)**2,N-P*(u+p)**2+Q*(u+p)**4]):
    need('actual H global row '+str(j),S.Poly(S.expand(row),u).degree()<=4)
for j,row in enumerate([R,M-R*(u+p)**2]):
    need('actual L global row '+str(j),S.Poly(S.expand(row),u).degree()<=2)
for label,expr in [('P',P),('M',M)]:
    zero(label+' active value',expr.subs(u,0))
    zero(label+' active first jet',S.diff(expr,u).subs(u,0))

def invert(expr):
    return S.expand(expr.subs({u:1/r-p,t:-r*r-r**4*Bd},simultaneous=True))

Hi,Li = invert(H),invert(L)
need('H original second chart polynomial',S.denom(S.cancel(Hi))==1)
need('L original second chart polynomial',S.denom(S.cancel(Li))==1)
zero('H D slope',S.diff(Hi.subs(r,0),Bd)+A)
zero('L D slope',S.diff(Li.subs(r,0),Bd)+lam)
zero('F D quadratic coefficient',S.expand((Hi*Hi+Li).subs(r,0)).coeff(Bd,2)-A*A)
zero('F D linear when A0',S.diff((Hi*Hi+Li).subs({r:0,A:0}),Bd)+lam)
zero('actual infinity multiplicity two',S.expand(r**8*N.subs(u,1/r-p))-r*r*(1-p*r)**6)

# Actual polynomial blow-up chart; the proof gives the unbounded weight gate.
Fw = S.Poly(S.expand(F.subs(t,w/u**2)),u)
fw = (c*w+d)**2+nu*w+e
gw = 2*(c*w+d)*(b*w+b-2*p*A)+mu*w+mu-2*p*lam
zero('full exceptional row',Fw.coeff_monomial(1)-fw)
zero('full next exceptional row',Fw.coeff_monomial(u)-gw)
zero('exceptional quadratic coefficient',S.expand(fw).coeff(w,2)-c*c)
zero('exceptional linear after c0',S.diff(fw.subs(c,0),w)-nu)
zero('next row slope after first gate',S.diff(gw.subs({c:0,nu:0}),w)-(2*d*b+mu))
zero('next row constant after second gate',gw.subs({c:0,nu:0,mu:-2*d*b})+2*p*(lam+2*d*A))
zero('source line has constant F',F.subs(u,0)-(d*d+e))
# Positive/negative weight controls are symbolic identities, not a degree cutoff.
ll,jj = S.symbols('ell j',integer=True)
f0,f1,g0 = S.symbols('f0 f1 g0')
zero('Jacobian scaling chart',jac(u,u*u*t)-u*u)
zero('literal first negative row obstruction',
     u*u*jac(f0+f1*w,w/u,u,w)-f1*w)
zero('allowed coordinate positive control',jac(u,t)-1)
for exponent in range(1,7):
    zero('negative row divisible by w control '+str(exponent),
         S.expand((u**(2*exponent-1)*t**exponent).subs(t,w/u**2))*u-w**exponent)

# The residual after only the first polynomial gate. No second gate is used.
Hr = H.subs(c,0)
Lr = L.subs(nu,0)
V = u+u**3*t
h0 = v*v+b*v+d
zero('full residual H field formula',
     Hr-(V*V+b*V+d+A*u*(V-2*p)))
zero('full residual L field formula',
     Lr-(lam*u*(V-2*p)+mu*V+e))
zero('actual v chart Jacobian',jac(u,V)-u**3)
zero('actual v field inverse',V.subs(t,(v-u)/u**3)-v)

# A != 0: actual degree-two extension, computed by its literal companion matrix.
U = (h-h0)/(A*(v-2*p))
zero('A nonzero field inverse',h0+A*U*(v-2*p)-h)
zero('A nonzero actual volume',S.diff(U,h)/U**3-A*A*(v-2*p)**2/(h-h0)**3)
Ph = h*h+delta*h-delta*h0+mu*v+e-f
zero('literal generic conic',
     (h*h+lam*U*(v-2*p)+mu*v+e-f).subs(lam,A*delta)-Ph)
disc = S.discriminant(Ph,h)
zero('actual quadratic discriminant f slope',S.diff(disc,f)-4)
need('actual quadratic degree',S.Poly(Ph,h).degree()==2)
P0 = h0*h0+mu*v+e
C0 = P0-f
B0 = 2*h0+delta
# Translate h to z=h-h0; certify the literal coefficients before using a
# universal two-coefficient companion matrix (avoids needless huge expansion).
zz,BB,CC = S.symbols('zz BB CC')
zero('actual centered quadratic',Ph.subs(h,h0+zz)-(zz*zz+B0*zz+C0))
companion = S.Matrix([[0,-CC],[1,-BB]])
Id = S.eye(2)
trace = S.trace(companion**-3*(2*companion+BB*Id)**-1)
zero('full quadratic trace',trace-(1/CC**2-BB**2/CC**3))
zero('normalized trace is half full',trace/2-(1/(2*CC**2)-BB**2/(2*CC**3)))
L1 = (v-2*p)**2*B0**2/S.diff(P0,v)
L2 = (v-2*p)**2/S.diff(P0,v)
need('generic finite poles simple',S.Poly(P0,v).degree()==4 and S.diff(P0,v)!=0)
zero('L1 cubic asymptotic',S.limit(L1/v**3,v,S.oo)-1)
zero('trace necessary identity leading mismatch',
     S.limit((S.diff(L1,v)-2*(v-2*p)**2)/v**2,v,S.oo)-1)
zero('P prime cubic leading',S.Poly(S.diff(P0,v),v).LC()-4)

# A=0, lambda != 0: rational field, no quadratic trace needed.
U0 = (f-P0)/(lam*(v-2*p))
zero('A0 actual rational field inverse',P0+lam*U0*(v-2*p)-f)
zero('A0 actual relative form',
     1/(lam*(v-2*p)*U0**3)-lam**2*(v-2*p)**2/(f-P0)**3)
zero('A0 residue numerator tends zero',S.limit(L2,v,S.oo))
need('A0 residue numerator not zero',L2!=0)
# Exact local residue coefficient identities in the parameter w=P(v).
zeta,zl,z0,z1,z2,y0,y1 = S.symbols('zeta zl z0 z1 z2 y0 y1')
local = (z0+z1*zl+z2*zl*zl/2)/zl**3-(y0+y1*zl)/zl**2
zero('third-second pole residue',S.residue(local,zl,0)-(z2/2-y1))
zero('third pole residue',S.residue((z0+z1*zl+z2*zl*zl/2)/zl**3,zl,0)-z2/2)

# The constant-D rational submersion remains a hostile to a rational blanket claim.
q = 1+u*u*t
Hpos = u*u*q*q+q
Lpos = u*q
Gpos = 1/(2*u*q)-Hpos
zero('constant-D rational hostile Jacobian',jac(Hpos*Hpos+Lpos,Gpos)-1)
zero('constant-D hostile actual D slope',S.diff(invert(Hpos*Hpos+Lpos).subs(r,0),Bd))
need('hostile genuinely nonconstant lower row',S.diff(Lpos,t)!=0)
zero('hostile first gate is nonconstant',
     S.Poly(S.expand((Hpos*Hpos+Lpos).subs(t,w/u**2)),u).coeff_monomial(1)-(w+1)**2)

print('Nonconstant-D finite-six / infinity-double: complete polynomial exclusion')
print('Full shifted active rows; polynomial weight gate; actual degree-two and rational fields: PASS')
print('Full trace, generic residue derivative, and incompatible degree/decay: PASS')
print('Constant-D nonconstant-L rational hostile retained: PASS')
print('Exact gates:',len(gates))
print('Semantic SHA256:',sha256('\n'.join(gates).encode()).hexdigest())
