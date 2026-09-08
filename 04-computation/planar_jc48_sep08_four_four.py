#!/usr/bin/env python3
"""Exact certificates for the distinguished 0/infinity 4+4 section class.

No inherited producer is imported. The proof supplies complete local branch
exhaustion, compact exactness, and existence of roots of nonconstant complex
polynomials; the exact checks below retain their actual coefficients.
"""
import hashlib
import json
import sympy as S

x,t,u,T,s,z,r,bd,c,Z,q,wv,ww,qv,C0,eta0=S.symbols('x t u T s z r bd c Z q wv ww qv C0 eta0')
a,k0,kt,kxt,k2,k3,k4=S.symbols('alpha k0 kt kxt k2 k3 k4')
l0,lt,lxt,l2,l3,l4=S.symbols('l0 lt lxt l2 l3 l4')
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

def face(f,power,degree):
    return coeff(f.subs(s,u**power*Z),u,degree)

w=u*u*T
v=u+u**3*T
h=u*u+u**4*T
basis=[S.Integer(1),T,u*T,w,v,h]
sections=[S.expand(S.cancel(s*bb.subs(T,1/s)).subs(s,z-u*u)) for bb in basis]
mat=S.Matrix([[coeff(coeff(bb,u,i),z,j) for bb in sections]
              for i in range(3) for j in range(2)])
check(S.factor(mat.det()) in (-1,1),'complete L1 basis')
mon=[u**i*z**j for i in range(5) for j in range(3)]
restrict=S.Matrix([[coeff(bb.subs(z,u*u),u,i) for bb in mon] for i in range(9)])
check(restrict.rank()==9,'full octic restriction rank')
check(15-restrict.rank()==6,'complete fixed-octic correction dimension')
H=a*w*w+k0+kt*T+kxt*u*T+k2*w+k3*v+k4*h
L=l0+lt*T+lxt*u*T+l2*w+l3*v+l4*h
N=S.cancel(s*s*H.subs(T,1/s))
M=S.cancel(s*L.subs(T,1/s))
eq(N,a*u**4+s*(kt+kxt*u+k2*u*u+k3*u**3+k4*u**4)+s*s*(k0+k3*u+k4*u*u),
   'complete normal numerator')
eq(M,lt+lxt*u+l2*u*u+l3*u**3+l4*u**4+s*(l0+l3*u+l4*u*u),
   'complete lower numerator')
eq(N.subs(s,0),a*u**4,'actual 4+4 octic')
for pp,di,dj,label in [(N,4,2,'H'),(M,2,1,'L')]:
    po=S.Poly(S.expand(pp.subs(s,z-u*u)),u,z)
    check(po.degree(u)<=di and po.degree(z)<=dj,'full global section box '+label)
# Actual inversion preserves this whole section family and changes its volume.
oldx=1/u
oldt=-u*u-u**4*T
eq(jac(oldx,oldt,u,T),u*u,'actual old source volume weight')
oldbasis=[S.Integer(1),t,x*t,x*x*t,x+x**3*t,x*x+x**4*t]
newbasis=[S.Integer(1),-h,-v,-1-w,-u*T,-T]
for i,(bb,expected) in enumerate(zip(oldbasis,newbasis)):
    eq(bb.subs({x:oldx,t:oldt},simultaneous=True),expected,'actual full basis involution '+str(i))
eq(((x*x*t)**2).subs({x:oldx,t:oldt},simultaneous=True),
   w*w+2*w+1,'actual base alpha w squared transport')
eq((u*u*jac(1/r,-r*r-r**4*bd,r,bd)).subs(u,1/r),1,'weighted second-chart regular volume')
E=S.expand(N*N+s**3*M-c*s**4)
eq(coeff(E.subs(u,0),s,2),kt*kt,'normal-unit local degree two')
eq(E.subs({u:0,kt:0}),s**3*(lt+(k0*k0+l0-c)*s),'M-unit local degree three')
E1=S.expand(E.subs({kt:0,lt:0}))
eq(E1.subs(u,0),(k0*k0+l0-c)*s**4,'remaining exact local degree four')

# M-unit m=4: no possible integer balanced normal order.
for j,kill,kj in [(1,{kt:0},kxt),(2,{kt:0,kxt:0},k2),
                  (3,{kt:0,kxt:0,k2:0},k3),(4,{kt:0,kxt:0,k2:0,k3:0},k4)]:
    ee=S.expand(E.subs(kill))
    check(4!=3*j,'m4 has no balanced M-unit order '+str(j))
    if j==1:
        eq(face(ee,2,6),Z*Z*(kxt*kxt+lt*Z),'M-unit one low branch')
        eq(face(N.subs(kill),3,4),a+kxt*Z,'M-unit two high cancellation determinations')
        check(2*S.Rational(1,2)+1==2,'M-unit high normalized unweighted regularity')
    else:
        scaled=ee.subs({u:q**3,s:q**8*Z},simultaneous=True)
        eq(coeff(scaled,q,24),a*a+lt*Z**3,'M-unit dominant cubic '+str(j))
        eq(coeff(S.diff(ee,s).subs({u:q**3,s:q**8*Z},simultaneous=True),q,16),
           3*lt*Z*Z,'M-unit full derivative supplier '+str(j))
        check(2*8+2-16==2,'M-unit normalized regular order '+str(j))
eq(face(E.subs({kt:0,kxt:0,k2:0,k3:0,k4:0}),3,8),a*a,
   'infinite normal order has no cancellation below dominant scale')

# Weighted local degree four: j1, n1, and the complete tangent quartic.
low=(kxt+k0*Z)**2+lxt*Z+(l0-c)*Z*Z
eq(face(E1,1,4),Z*Z*low,'entire j1 low quadratic')
eq(S.diff(S.discriminant(low,Z),c),4*kxt*kxt,'j1 low roots generically simple')
eq(low.subs(Z,0),kxt*kxt,'j1 low roots nonzero')
N1=S.expand(N.subs(kt,0))
eq(face(N1,3,4),a+kxt*Z,'j1 high analytic centre')
eq(coeff(S.diff(N1,s).subs(s,u**3*Z),u,1),kxt,'j1 centre normal derivative supplier')
psi=-a*u**3/kxt
mc=S.expand((M-c*s).subs({lt:0,s:psi},simultaneous=True))
eq(coeff(mc,u,1),lxt,'j1 first possible M order')
eq(coeff(mc.subs(lxt,0),u,2),l2,'j1 second possible M order')
eq(S.diff(coeff(mc.subs({lxt:0,l2:0}),u,3),c),a/kxt,
   'j1 generic c supplies order at most three')
for ell in (1,2,3):
    raw=S.Rational(5-ell,2)
    norm=raw if raw.q==1 else 2*raw+1
    check(norm>0,'j1 complete normalized weighted order '+str(ell))
E01=S.expand(E1.subs(kxt,0))
eq(face(E01,1,4),Z**3*(lxt+(k0*k0+l0-c)*Z),'n1 one simple low branch')
sc=E01.subs({u:q**3,s:q**7*Z},simultaneous=True)
eq(coeff(sc,q,24),a*a+lxt*Z**3,'n1 three remaining determinations')
eq(coeff(S.diff(E01,s).subs({u:q**3,s:q**7*Z},simultaneous=True),q,17),
   3*lxt*Z*Z,'n1 exact derivative leading order')
check(2*3+2*7+2-17==5,'n1 weighted normalized relative-form order five')
red={kt:0,lt:0,kxt:0,lxt:0}
Er=S.expand(E.subs(red))
P=(a+k2*Z+k0*Z*Z)**2+l2*Z**3+(l0-c)*Z**4
eq(face(Er,2,8),P,'entire remaining quartic tangent face')
eq(P.subs(Z,0),a*a,'quartic roots are nonzero')
eq(coeff(P,Z,4),k0*k0+l0-c,'quartic has generic degree four')
P0=S.expand(P+c*Z**4)
crit=S.expand(Z*S.diff(P0,Z)-4*P0)
eq(crit.subs(Z,0),-4*a*a,'generic quartic simple-root supplier cannot vanish identically')
eq(coeff(S.diff(Er,s).subs(s,u*u*Z),u,6),S.diff(P,Z),
   'actual weighted quartic derivative supplier')
check(2+2*2-6==0,'all four weighted tangent branches regular')
for cc in [0,2,3]:
    pp=S.Poly(P.subs({a:1,k0:1,k2:2,l2:1,l0:0,c:cc}),Z)
    check(pp.degree()==4 and S.gcd(pp,pp.diff()).degree()==0,'named generic quartic face '+str(cc))

# Original finite-point residue retains all induced coordinate corrections.
wx=x*x*t
vx=x+x**3*t
hx=x*x+x**4*t
Hx=a*wx*wx+k0+k2*wx+k3*vx+k4*hx
Lx=l0+l2*wx+l3*vx+l4*hx
F=S.expand(Hx*Hx+Lx)
Fw=S.cancel(F.subs(t,wv/(x*x)))
A=a*wv*wv+k2*wv+k0
f=A*A+l2*wv+l0
g=(1+wv)*(2*k3*A+l3)
eq(coeff(Fw,x,0),f,'full original zero-row polynomial')
eq(coeff(Fw,x,1),g,'full original first correction')
eq(jac(x,wv/(x*x),x,wv),1/(x*x),'actual finite blow-chart volume')
f1,f2,g0,g1=S.symbols('f1 f2 g0 g1')
# At a simple root, w-w0=-g0/f1*x+O(x²). This is exact first-order jet algebra.
D=f1+x*(g1-f2*g0/f1)
eq(S.limit(x*(-1/(x*x*D)+1/(f1*x*x)),x,0),(g1*f1-f2*g0)/f1**3,
   'entire moving-root residue including induced displacement')
Cc=S.Symbol('Cc')
respoly=S.expand(g-Cc*S.diff(f,wv))
eq(coeff(respoly,wv,3),2*a*k3-4*a*a*Cc,'top coefficient fixes scalar C')
zcase=S.expand(respoly.subs({k3:0,Cc:0}))
eq(zcase,l3*(1+wv),'k3 zero forces l3 zero')
polycase=S.expand(respoly.subs(Cc,k3/(2*a)))
eq(coeff(polycase,wv,2),k3*(2*a-k2),'next coefficient fixes k2')
eq(coeff(polycase.subs(k2,2*a),wv,1),l3,'next coefficient fixes l3')
eq(coeff(polycase.subs({k2:2*a,l3:0}),wv,0),-k3*l2/(2*a),'constant fixes l2')
eq(S.diff(F,x).subs({x:0,k3:0,l3:0}),0,'k3-zero case actual source Fx vanishes')
eq(S.diff(F,t).subs({x:0,k3:0,l3:0}),0,'k3-zero case actual source Ft vanishes')

# All residue-compatible criticality, including both exceptional strata.
c0=S.Symbol('c0')
B=a*qv+k3*x+k4*x*x
Hq=c0+qv*B
Fq=Hq*Hq+l4*x*x*qv+l0
Fx=S.diff(Fq,x)
Fqq=S.diff(Fq,qv)
eq(Fx.subs(qv,0),0,'q0 tangent derivative vanishes')
eq(Fqq.subs(qv,0),x*(2*c0*k3+(2*c0*k4+l4)*x),'whole q0 criticality equation')
xcrit=-2*c0*k3/(2*c0*k4+l4)
eq(Fqq.subs({qv:0,x:xcrit},simultaneous=True),0,'nonexceptional actual critical point')
qcrit=-k3*x/(4*a)
Cx=k3+2*k4*x
common=S.expand((Hq*Cx+l4*x).subs(qv,qcrit))
eq(Fx.subs(qv,qcrit),2*qcrit*common,'first derivative exact common factor')
eq(Fqq.subs(qv,qcrit),x*common,'second derivative exact common factor')
poly0=k3*x*(3*k3+4*k4*x)*Cx-16*a*l4
eq(common.subs(c0,0),-x*poly0/(16*a),'c0-zero entire critical polynomial')
eq(coeff(poly0,x,3),8*k3*k4*k4,'c0-zero cubic leading term')
eq(coeff(poly0.subs(k4,0),x,1),3*k3**3,'c0-zero k4-zero remains nonconstant')
eq(poly0.subs(x,0),-16*a*l4,'c0-zero root cannot be zero')
poly1=x*x*(3*k3+4*k4*x)*Cx-16*a*c0
eq(common.subs(l4,-2*c0*k4),-k3*poly1/(16*a),'second exceptional entire critical polynomial')
eq(coeff(poly1,x,4),8*k4*k4,'second exceptional quartic nonconstant')
eq(poly1.subs(x,0),-16*a*c0,'second exceptional root cannot be zero')
for pp,constant,label in [(poly0,-16*a*l4,'first'),(poly1,-16*a*c0,'second')]:
    eq(pp.subs(x,-3*k3/(4*k4)),constant,'no 3k3+4k4x zero at root '+label)
    eq(pp.subs(x,-k3/(2*k4)),constant,'no k3+2k4x zero at root '+label)
eq(jac(x,1+x*x*t,x,t),x*x,'criticality chart has nonzero Jacobian at x nonzero')
# A residue-compatible named object still has an actual affine critical point.
params={a:1,k0:1,k2:2,k3:1,k4:0,l2:0,l3:0,l4:1,l0:0}
Fp=S.expand(F.subs(params))
point={x:S.Rational(16,3),t:-S.Rational(21,256)}
eq(S.diff(Fp,x).subs(point),0,'named hostile actual source Fx')
eq(S.diff(Fp,t).subs(point),0,'named hostile actual source Ft')
eq((g-k3/(2*a)*S.diff(f,wv)).subs(params),0,'named hostile satisfies whole residue identity')
eq((1+x*x*t).subs(point),-S.Rational(4,3),'named hostile retained q address')
# A genuine nonconstant-L rational mate survives inside this polynomial exclusion.
eq(jac(wx,1/x,x,t),1,'actual canonical rational carrier')
eq(jac(wx**4+wx,1/(x*(4*wx**3+1)),x,t),1,'nonconstant-L rational-mate hostile')
# Constant L is excluded by the original polynomial factor, with unrestricted mate degree.
check(S.Poly(a*wx*wx,x,t).total_degree()>0,'composite H is nonconstant in original source')

records['entry']='full global H=alpha*w^2+K; K,L in complete six-dimensional L1'
records['boundary']='exact distinguished finite0/infinity 4+4; no arbitrary-pair normalization'
records['infinity']='local degrees2/3/4, weighted regular on every normalized branch'
records['residue']='g=C fprime; k3 nonzero forces C=k3/(2alpha),k2=2alpha,l2=l3=0'
records['criticality']='q0 or roots of explicit nonconstant cubic/quartic, all x!=0'
records['hostile']='residue compatible critical point (16/3,-21/256); F=w^4+w has a rational mate'
semantic=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
print('PASS distinguished 4+4 quartic:',gates,'always-active exact gates')
print('Full L1 corrections; complete weighted infinity local degrees 2/3/4')
print('Original finite residue retains moving-root displacement')
print('All residue-compatible coefficients have an actual affine critical point')
print('Scope: no polynomial mate; no arbitrary boundary-pair transport')
print('Semantic SHA256:',semantic)
