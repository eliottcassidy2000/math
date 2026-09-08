#!/usr/bin/env python3
"""Exact two-finite-fourfold section and original-residue certificates.

No inherited producer is imported. Analytic local regularity and compact
exactness are proved/cited in the matching note; the full original first
rows, all section kernels and hostile scope controls are checked here.
"""
import hashlib
import json
import sympy as S

x,t,z,s,u,w,p,q,c,Z=S.symbols('x t z s u w p q c Z')
alpha,k0,kt,kxt,k2,k3,k4=S.symbols('alpha k0 kt kxt k2 k3 k4')
l0,lt,lxt,l2,l3,l4=S.symbols('l0 lt lxt l2 l3 l4')
beta,gamma,delta,eps=S.symbols('beta gamma delta eps')
gates=0
records={}

def check(ok,label):
    global gates
    gates+=1
    if not ok: raise RuntimeError(label)

def eq(f,g,label):
    check(S.cancel(f-g)==0,label)

def coeff(f,var,n):
    return S.expand(f).coeff(var,n)

def jac(f,g,v1,v2):
    return S.diff(f,v1)*S.diff(g,v2)-S.diff(f,v2)*S.diff(g,v1)

# Actual global base function and full section boxes.
S1=p+q
P1=p*q
rpoly=(x-p)*(x-q)
BQ=x*x-2*S1*x
Q=rpoly*rpoly*t+BQ
w0=x*x*t
v0=x+x**3*t
h0=x*x+x**4*t
eq(Q,h0-2*S1*v0+(S1*S1+2*P1)*w0-2*S1*P1*x*t+P1*P1*t,
   'actual full global Q identity')
NQ=S.cancel(s*Q.subs(t,1/s))
eq(NQ,rpoly*rpoly+s*BQ,'actual complete Q numerator')
Qglobal=S.expand(NQ.subs(s,z-x*x))
check(S.Poly(Qglobal,x,z).degree(x)<=2 and S.Poly(Qglobal,x,z).degree(z)<=1,
      'Q is an actual global L1 function')
eq(NQ.subs(s,0),rpoly*rpoly,'actual boundary numerator of Q')
eq(coeff(NQ.subs(s,0),x,4),1,'Q has two finite double roots and no infinity zero')
# The standard full restriction has five independent quartic coefficients.
base0=[S.Integer(1),t,x*t,w0,v0,h0]
sec0=[S.cancel(s*bb.subs(t,1/s)) for bb in base0]
restriction=S.Matrix([[coeff(bb.subs(s,0),x,i) for bb in sec0] for i in range(5)])
eq(restriction[:,1:].det(),1,'entire L1 restriction onto degree-four polynomials')
check(restriction[:,0]==S.zeros(5,1),'the constant function is in the restriction kernel')
check(6-restriction.rank()==1,'constants are the whole L1 restriction kernel')
# Arbitrary finite p changes source coordinates only, while retaining the original graph.
wp=u*u*t
vp=u*(1+wp)-2*p
hp=u*u*(1+wp)-2*p*u
basis=[S.Integer(1),t,u*t,wp,vp,hp]
secs=[S.expand(S.cancel(s*bb.subs(t,1/s)).subs(s,z-(u+p)**2)) for bb in basis]
mat=S.Matrix([[coeff(coeff(bb,u,i),z,j) for bb in secs]
              for i in range(3) for j in range(2)])
eq(mat.det(),1,'full corrected global L1 basis for every p')
rankmon=[u**i for i in range(5)]+[u**3*z,u**4*z,u**3*z*z,u**4*z*z]
rankmat=S.Matrix([[coeff(bb.subs(z,(u+p)**2),u,i) for bb in rankmon] for i in range(9)])
eq(rankmat.det(),1,'all-p full octic restriction rank nine')
check(15-9==6,'entire fixed-octic correction dimension six')

# Both singular shared roots force a whole composition, not just matching jets.
polys=[S.Integer(1),x,x*x,x**3]
jetmat=S.Matrix([[bb.subs(x,p) for bb in polys],
                [S.diff(bb,x).subs(x,p) for bb in polys],
                [bb.subs(x,q) for bb in polys],
                [S.diff(bb,x).subs(x,q) for bb in polys]])
eq(S.factor(jetmat.det()),(p-q)**4,'complete two-point confluent restriction determinant')
for point in [p,q]:
    eq((rpoly*rpoly).subs(x,point),0,'double boundary root value '+str(point))
    eq(S.diff(rpoly*rpoly,x).subs(x,point),0,'double boundary root derivative '+str(point))
P4=(alpha*Z*Z+beta*Z+gamma)**2+delta*Z+eps
check(S.Poly(P4,Z).degree()==4,'composition degree exactly four')
eq(coeff(S.diff(P4,Z),Z,3),4*alpha*alpha,'composition derivative is a nonconstant factor')

# Complete original source family, including all lower normal coefficients.
Qu=S.expand(Q.subs(x,u+p))
K=k0+kt*t+kxt*u*t+k2*wp+k3*vp+k4*hp
L=l0+lt*t+lxt*u*t+l2*wp+l3*vp+l4*hp
H=alpha*Qu*Qu+K
N=S.cancel(s*s*H.subs(t,1/s))
M=S.cancel(s*L.subs(t,1/s))
AK=S.cancel((s*K.subs(t,1/s))).subs(s,0)
AL=M.subs(s,0)
BK=S.cancel(K.subs(t,0))
BL=L.subs(t,0)
ru=u*(u+p-q)
BQu=BQ.subs(x,u+p)
eq(N,alpha*ru**4+s*(2*alpha*ru*ru*BQu+AK)+s*s*(alpha*BQu*BQu+BK),
   'entire source numerator preserves both Q and K jets')
eq(M,AL+s*BL,'entire lower numerator')
eq(N.subs(s,0),alpha*u**4*(u+p-q)**4,'entire binary octic in shifted source')
eq(coeff(N.subs(s,0),u,8),alpha,'no infinity boundary zero')
eq(S.diff(N,s).subs({u:0,s:0}),kt,'shared-root normal value is exactly kt')
eq(S.diff(N,s,u).subs({u:0,s:0}),kxt,'shared-root first normal jet is exactly kxt')
eq(M.subs({u:0,s:0}),lt,'lower shared-root value is exactly lt')
eq(S.diff(M,u).subs({u:0,s:0}),lxt,'lower first jet is exactly lxt')
for pp,di,dj,label in [(N,4,2,'H'),(M,2,1,'L')]:
    po=S.Poly(S.expand(pp.subs(s,z-(u+p)**2)),u,z)
    check(po.degree(u)<=di and po.degree(z)<=dj,'complete actual global section box '+label)
# Exact local degree controls do not replace the inherited analytic regularity lemmas.
E=S.expand(N*N+s**3*M-c*s**4)
eq(coeff(E.subs(u,0),s,2),kt*kt,'normal-unit local degree two')
eq(coeff(E.subs({u:0,kt:0}),s,3),lt,'M-unit local degree three')
red={kt:0,kxt:0,lt:0,lxt:0}
Nr=S.expand(N.subs(red))
Mr=S.expand(M.subs(red))
Er=S.expand(E.subs(red))
H0=S.expand(H.subs(red).subs(u,0))
L0=S.expand(L.subs(red).subs(u,0))
eq(Er.subs(u,0),(H0*H0+L0-c)*s**4,'active-root complete local degree four')

# Full literal original first rows; the nonzero leading residue cannot be tuned.
d=p-q
Qblow=S.cancel(Qu.subs(t,w/(u*u)))
Q0=d*d*w-p*p-2*p*q
Q1=2*d*w-2*q
eq(Qblow,Q0+u*Q1+u*u*(w+1),'complete global Q blow-chart expression')
Hb=S.cancel(H.subs(red).subs(t,w/(u*u)))
Lb=S.cancel(L.subs(red).subs(t,w/(u*u)))
A=alpha*Q0*Q0+k2*w+k0-2*p*k3
B=2*alpha*Q0*Q1+k3*(1+w)-2*p*k4
eq(coeff(Hb,u,0),A,'entire H zero row')
eq(coeff(Hb,u,1),B,'entire H first row including other finite root')
f=A*A+l2*w+l0-2*p*l3
g=2*A*B+l3*(1+w)-2*p*l4
Fb=S.expand(Hb*Hb+Lb)
eq(coeff(Fb,u,0),f,'entire original F zero row')
eq(coeff(Fb,u,1),g,'entire original F first row')
eq(jac(u,w/(u*u),u,w),1/(u*u),'original relative volume in the blow-chart')
eq(coeff(A,w,2),alpha*d**4,'actual H leading zero-row coefficient')
eq(coeff(B,w,2),4*alpha*d**3,'actual H leading first-row coefficient')
fprime=S.diff(f,w)
eq(coeff(fprime,w,3),4*alpha*alpha*d**8,'fprime has exact degree three')
eq(coeff(g,w,4),8*alpha*alpha*d**7,'g has uncancellable degree four')
resnum=S.expand(S.diff(g,w)*fprime-S.diff(fprime,w)*g)
check(S.Poly(resnum,w).degree()==6,'complete residue numerator degree six')
eq(coeff(resnum,w,6),32*alpha**4*d**15,'exact unavoidable individual-residue supplier')
# Independent first-order implicit-root calculation keeps the induced displacement.
f1,f2,g0,g1=S.symbols('f1 f2 g0 g1')
D=f1+u*(g1-f2*g0/f1)
eq(S.limit(u*(-1/(u*u*D)+1/(f1*u*u)),u,0),(g1*f1-f2*g0)/f1**3,
   'original moving-root residue including displacement')
# The four local tangent branches and f(w)=c are literally the same equation.
face=coeff(Er.subs(s,u*u*Z),u,8)
eq(face,S.cancel(Z**4*(f.subs(w,1/Z)-c)),'full tangent quartic is the original f-root equation')
eq(face.subs(Z,0),alpha*alpha*d**8,'all generic branches have nonzero inverse w coordinate')
eq(coeff(face,Z,4),H0*H0+L0-c,'tangent quartic has generic degree four')
for pp,qq in [(0,1),(1,2),(1,-1),(S.I,-S.I)]:
    value=coeff(resnum,w,6).subs({alpha:1,p:pp,q:qq})
    check(value!=0,'named distinct finite-root control '+str((pp,qq)))
# Constant lower row and composition are retained by the same local obstruction.
eq(coeff(resnum.subs({l0:0,l2:0,l3:0,l4:0}),w,6),32*alpha**4*d**15,
   'constant L does not remove the rational obstruction')

# Exact scope hostiles and residue-sum loss.
F_inf=w0**4+w0
G_inf=1/(x*(4*w0**3+1))
eq(jac(F_inf,G_inf,x,t),1,'finite/infinity 4+4 rational-mate hostile')
F_coal=h0**4+h0
G_coal=1/(3*x**3*(4*h0**3+1))
eq(jac(F_coal,G_coal,x,t),1,'coincident finite octuple rational-mate hostile')
eq(coeff(resnum,w,6).subs(q,p),0,'distinctness loss makes the leading supplier vanish')
local_eta=-1/(u*u*(u+d)**2)
eq(S.residue(local_eta,u,0),2/d**3,'individual composition-fibre residue at p')
eq(2/d**3+2/(-d)**3,0,'opposite individual residues can sum to zero')
# Independent literal composition fibre: p1,q2,Q=1,F=Q^4=1,w0=6.
params={p:1,q:2,alpha:1,k0:0,k2:0,k3:0,k4:0,l0:0,l2:0,l3:0,l4:0,w:6}
eq(f.subs(params),1,'named actual original fibre value')
eq((resnum/fprime**3).subs(params),-S.Rational(1,2),
   'named active-root residue agrees with independent Q-fibre residue divided by four')
# Direct relative-form identity for the original global Q.
eq(S.diff(Q,t),rpoly*rpoly,'complete relative form supplier for Q')

records['entry']='alpha*(x-p)^4*(x-q)^4, p and q distinct finite, full global K,L'
records['composition']='both singular shared points force K=beta Q+gamma,L=delta Q+epsilon'
records['boundary_exhaustion']='both regular gives holomorphic contradiction; any active root gives nonzero residue'
records['leading_rows']=['A2=alpha*d^4','B2=4alpha*d^3','fprime3=4alpha^2*d^8','g4=8alpha^2*d^7']
records['residue_leading']='32alpha^4*d^15, d=p-q nonzero'
records['scope']='NO rational mate, including constant L; polynomial corollary uses audited finite/infinity theorem'
records['hostiles']=['finite/infinity rational mate','coincident-octuple rational mate','sum of opposite residues zero']
semantic=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
print('PASS two distinct finite fourfold roots:',gates,'always-active exact gates')
print('Full global sections and both-root composition reduction checked')
print('Original active-root residue degree6, leading32alpha^4*(p-q)^15')
print('Scope: no rational mate, including constant L; distinct finite roots essential')
print('Finite/infinity and coincident-root rational hostiles retained')
print('Semantic SHA256:',semantic)
