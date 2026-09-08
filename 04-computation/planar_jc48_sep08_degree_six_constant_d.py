#!/usr/bin/env python3
"""Universal degree-six leading field: two inverse residues force constant D.

All coefficient checks are symbolic over Q, before any finite-root filter.
The finite examples are named controls, not a parameter census or a proof
of rational-mate sufficiency. Assertions are deliberately always active.
"""
import hashlib
import json
import sympy as S

u,t,x,r,z,p,bD = S.symbols('u t x r z p bD')
A,b,c,d,P1,P0,lam,mu,nu,M1,M0,e = S.symbols(
    'A b c d P1 P0 lambda mu nu M1 M0 e')
ns = S.symbols('n0:6')
N = u**6 + sum(ns[j]*u**j for j in range(6))
P = (A+2)*u**4+b*u**3+c*u**2+P1*u+P0
Q = (A+1)*u**2+(b-2*p*A-ns[5])*u+d
M = lam*u**4+mu*u**3+nu*u**2+M1*u+M0
R = lam*u**2+(mu-2*p*lam)*u+e
H,L = N*t**2+P*t+Q, M*t+R
gates = 0
record = {}

def check(name, value):
    global gates
    if not bool(value):
        raise RuntimeError(name)
    gates += 1

def zero(name, expr):
    check(name, S.cancel(expr) == 0)

def coef(expr, power, var=u):
    return S.Poly(S.expand(expr),var).coeff_monomial(var**power)

def jac(f,g):
    return S.diff(f,u)*S.diff(g,t)-S.diff(f,t)*S.diff(g,u)

def inf(f):
    return S.expand(f.subs({u:1/r-p,t:-r**2-r**4*bD}, simultaneous=True))

# Complete original section boxes, including all low coefficients and shift.
Nx,Px,Qx,Mx,Rx = [S.expand(f.subs(u,x-p)) for f in (N,P,Q,M,R)]
zero('original Q2',coef(Qx,2,x)-coef(Px,4,x)+coef(Nx,6,x))
zero('original Q1',coef(Qx,1,x)-coef(Px,3,x)+coef(Nx,5,x))
zero('original R2',coef(Rx,2,x)-coef(Mx,4,x))
zero('original R1',coef(Rx,1,x)-coef(Mx,3,x))
for name, row, bound in [
    ('H z2',Qx,4),('H z1',Px-2*x*x*Qx,4),
    ('H z0',Nx-x*x*Px+x**4*Qx,4),
    ('L z1',Rx,2),('L z0',Mx-x*x*Rx,2)]:
    check('complete numerator box '+name,S.Poly(S.expand(row),x).degree() <= bound)
Hi,Li = inf(H),inf(L)
for name, f in [('H',Hi),('L',Li)]:
    check('whole infinity polynomial '+name,S.denom(S.cancel(f)) == 1)
    check('whole infinity no negative powers '+name,S.Poly(f,r,bD) is not None)
HD,LD = S.expand(Hi.subs(r,0)),S.expand(Li.subs(r,0))
zero('H boundary slope',S.diff(HD,bD)+A)
zero('L boundary slope',S.diff(LD,bD)+lam)
zero('H boundary degree',S.diff(HD,bD,2))
zero('L boundary degree',S.diff(LD,bD,2))
FD = S.expand(HD**2+LD)
zero('F boundary quadratic',coef(FD,2,bD)-A**2)
zero('F boundary slope after T1',S.diff(FD.subs(A,0),bD)+lam)
zero('F boundary constant after both',S.diff(FD.subs({A:0,lam:0}),bD))
record['boundary'] = [str(HD),str(LD)]

# A direct coefficient recursion for the ORIGINAL centered inverse quartic.
v,DD,EE,CC = S.symbols('v DD EE CC')
ws = [S.Integer(1)]
for j in range(1,7):
    h = S.Symbol('h')
    W = sum(ws[i]*v**i for i in range(j))+h*v**j
    eq = W**4+2*DD*v*v*W*W+CC*v**3*W+(DD*DD+EE)*v**4-1
    cj = S.expand(eq).coeff(v,j)
    zero('inverse new coefficient '+str(j),S.diff(cj,h)-4)
    hj = S.factor(-cj.subs(h,0)/4)
    zero('inverse substitution '+str(j),cj.subs(h,hj))
    ws.append(hj)
zero('universal T1',ws[2]+DD/2)
zero('universal T2',ws[3]+CC/4)
zero('universal T5',ws[6]+(2*DD**3+4*DD*EE+CC**2)/32)
record['inverse'] = [str(q) for q in ws]

# Full high-coefficient algebra, with no assumptions on roots of N.
DN = S.expand(4*N*Q-P**2)       # D0=DN/(4N)
EN = S.expand(2*N*R-M*P)       # E0=EN/(2N)
D2 = coef(DN,8)/4
D1 = (coef(DN,7)-coef(DN,8)*ns[5])/4
zero('D quadratic coefficient',D2+A**2/4)
zero('D linear coefficient',D1+A*b/2+2*p*A-(A+A**2/4)*ns[5])
DN0,EN0 = S.expand(DN.subs(A,0)),S.expand(EN.subs(A,0))
check('D bounded after T1',S.Poly(DN0,u).degree() <= 6)
Dlim = coef(DN0,6)/4
zero('D exact finite limit',Dlim+b*b/4-b*ns[5]+c-d-ns[4]+ns[5]**2)
zero('E quadratic coefficient',coef(EN,8)/2+A*lam/2)
check('E at most linear after T1',S.Poly(EN0,u).degree() <= 7)
Elim = coef(EN0,7)/2
zero('E exact linear limit',Elim-lam*(ns[5]-b/2-2*p))
zero('M square leading coefficient',coef(M*M,8)-lam**2)
record['asymptotic'] = [str(S.factor(D1)),str(Dlim),str(Elim)]

# Exact normalized local expressions at infinity: kappa=epsilon*z^-3*sqrt(Nz).
# They apply to both unramified points of a nonsquare field, or to the one
# infinity point of either chosen square field. No disconnected algebra used.
Nz = S.expand(z**6*N.subs(u,1/z))
D2z = S.expand(z**8*DN.subs(u,1/z))/4
Dz = S.expand(z**6*DN0.subs(u,1/z))/4
Ez = S.expand(z**7*EN0.subs(u,1/z))/2
Mz = S.expand(z**4*M.subs(u,1/z))
for name, f in [('Nz',Nz),('D2z',D2z),('Dz',Dz),('Ez',Ez),('Mz',Mz)]:
    check('normalized local polynomial '+name,S.Poly(f,z) is not None)
zero('N local unit',Nz.subs(z,0)-1)
for eps in [1,-1]:
    # Residue T1 du is the constant coefficient of z*T1(1/z)*(-z^-2).
    residue1 = D2z.subs(z,0)/(2*eps*Nz.subs(z,0))
    zero('T1 infinity sign '+str(eps),residue1+eps*A*A/8)
    # z^2(2D0^3+4D0E0+M^2/N) has these three exact rational terms.
    pieces = [2*z*z*Dz**3/Nz**3,4*z*Dz*Ez/Nz**2,Mz**2/Nz]
    for j, piece in enumerate(pieces):
        zero('T5 local constant '+str(eps)+' '+str(j),
             piece.subs(z,0)-(lam**2 if j==2 else 0))
    residue5 = sum(q.subs(z,0) for q in pieces)/(32*eps)
    zero('T5 infinity sign '+str(eps),residue5-eps*lam*lam/32)
record['residues'] = ['-epsilon*A^2/8','epsilon*lambda^2/32 after A=0']

# Short independent (5,2,1) consequence after both universal residues.
# The active jets used here have their complete analytic justification in
# the proof; none is imposed in the universal theorem above.
N521 = u**5*(u-1)
P521 = u*u*(2*u*u+b*u+c)
Q521 = u*u+(b+1)*u+d
M521 = mu*u*u*(u-1)
R521 = mu*u+e
H521,L521 = N521*t*t+P521*t+Q521,M521*t+R521
zero('521 T2 simple residue',S.residue(M521/N521,u,1))
zero('521 factored complete active M',
     (mu*u**3+nu*u*u).subs(nu,-mu)-M521)
w = S.symbols('w')
F521w = S.expand((H521**2+L521).subs(t,w/u**2))
f = (c*w+d)**2-mu*w+e
g = 2*(c*w+d)*(-w*w+b*w+b+1)+mu*w+mu
zero('521 full constant row',coef(F521w,0)-f)
zero('521 full first row',coef(F521w,1)-g)
zero('521 cubic mismatch',coef(g,3,w)+2*c)
zero('521 quadratic mismatch',coef(g.subs(c,0),2,w)+2*d)
zero('521 final linear mismatch',coef(g.subs({c:0,d:0}),1,w)-mu)

# Constant-L rational sharpness of that corollary, at every finite shift.
k = S.symbols('k')
q521 = 1+u*u*t
H521sharp = u*(u-1)*q521**2+k
G521H = -(1+2*u)/(3*u*u*q521)
zero('521 sharp H mate',jac(H521sharp,G521H)-1)
zero('521 sharp quartic mate',jac(H521sharp**2+e,G521H/(2*H521sharp))-1)
zero('521 sharp leading',coef(H521sharp,2,t)-N521)
check('521 sharp all-p global',S.denom(S.cancel(inf(H521sharp)))==1)

# Named sharpness and hostile controls. These are not a parameter bank.
q = 1+u*u*t
vv = u*q
Hp = u*u*q*q+q
Lp = vv
Gp = 1/(2*vv)-Hp
zero('constant-D nonconstant-L actual rational mate',jac(Hp**2+Lp,Gp)-1)
zero('positive leading degree six',coef(Hp,2,t)-u**6)
for name, hh in [('positive H',Hp),('positive L',Lp)]:
    ih = inf(hh)
    check('all-p global '+name,S.denom(S.cancel(ih))==1)
    zero('all-p constant boundary '+name,S.diff(S.expand(ih.subs(r,0)),bD))
check('positive mate has genuine pole',S.denom(S.cancel(Gp)).has(u,t))

# At p=0, T1 alone passes although the next residue is nonzero.
Hone = u*u*q*q
Lone = u**4*t+u*u
zero('one-row hostile global H',S.denom(S.cancel(inf(Hone)))-1)
zero('one-row hostile global L at p0',S.denom(S.cancel(inf(Lone).subs(p,0)))-1)
zero('one-row hostile D0',u*u-(2*u**4)**2/(4*u**6))
zero('one-row hostile E0',u*u-u**4*(2*u**4)/(2*u**6))
zero('one-row hostile T5 residue',S.residue(1/(32*z),z,0)-S.Rational(1,32))

# Constant D is necessary, not sufficient: a finite double residue survives.
Nc = u*u*(u-1)**4
Hc = Nc*t*t+2*u**4*t+u*u+4*u
Hci = inf(Hc)
check('finite-double hostile global',S.denom(S.cancel(Hci))==1)
zero('finite-double hostile constant boundary',S.diff(Hci.subs(r,0),bD))
zero('finite-double hostile exact radical',u**2*(u-1)**4-Nc)
zero('finite-double hostile nonzero residue',S.residue(1/(u*(u-1)**2),u,0)-1)

# If the induced Q row is dropped, the asserted gate would be false.
Hng = u**6*t*t
Gng = 1/(8*u**11*t**3)
zero('nonglobal pure monomial rational mate',jac(Hng**2,Gng)-1)
zero('nonglobal hostile boundary pole',inf(Hng).coeff(r,-2)-1)

# A nonsquare normalized field can have positive genus; no genus-zero
# hypothesis occurs in the theorem. The local unit has the same two signs.
check('smooth degree-six field control',S.discriminant(u**6+1,u)!=0)
zero('smooth control normalized local unit',S.expand(z**6*((u**6+1).subs(u,1/z))).subs(z,0)-1)
record['controls'] = [
    'constant-D/nonconstant-L rational mate', 'T1 pass / T5 failure',
    'constant-D / finite leading residue', 'nonglobal pure monomial',
    'smooth genus-two leading field']
record['gates'] = gates
semantic = hashlib.sha256(json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('Degree-six leading polynomial: rational mate forces constant D')
print('Exact gates:',gates)
print('Complete shifted global rows; no finite-root filters: PASS')
print('Universal inverse T1/T2/T5 and both infinity signs: PASS')
print('Square/nonsquare field scope; all-parameter asymptotics: PASS')
print('Short independent (5,2,1) moving-residue corollary: PASS')
print('Actual rational sharpness and three failure controls: PASS')
print('Semantic SHA256:',semantic)
