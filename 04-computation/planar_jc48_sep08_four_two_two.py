#!/usr/bin/env python3
"""SUPERSEDED RESEARCH DRAFT: active-fourfold/symmetric 4+2+2 controls.

No companion proof or independent audit was completed for this route.
Its proposed analytic consumers below are provenance, not proved dependencies.
The separately proved planar_jc48_sep08_leading_exactness theorem now excludes
the entire 4+2+2 partition even for rational mates at all boundary locations.
This weaker critical-polynomial experiment is retained for its factorization
and forbidden-root mechanism. It is outside the audited checkpoint manifest.
No inherited mathematical producer is imported.
"""
import hashlib
import json
import sympy as S

x,t,s,z,w,q,p,r1,r2,d,alpha,c=S.symbols('x t s z w q p r1 r2 d alpha c')
k0,kt,kxt,k2,k3,k4,l0,lt,lxt,l2,l3,l4=S.symbols('k0 kt kxt k2 k3 k4 l0 lt lxt l2 l3 l4')
beta,c0,C=S.symbols('beta c0 C')
gates=0

def check(value,label):
    global gates
    gates+=1
    if not value:
        raise RuntimeError(label)

def eq(a,b,label):
    check(S.cancel(a-b)==0,label)

def cf(f,var,n):
    return S.expand(f).coeff(var,n)

def jac(a,b,v1,v2):
    return S.diff(a,v1)*S.diff(b,v2)-S.diff(a,v2)*S.diff(b,v1)

# Universal original-coordinate first-row calculation at an active order-four root.
a0,a1,b0,b1,c00,c10,d0,d1,e0,e1=S.symbols('a0 a1 b0 b1 c00 c10 d0 d1 e0 e1')
H0=a0*w*w+b0*w+c00
H1=a1*w*w+b1*w+c10
f=H0*H0+d0*w+e0
g=2*H0*H1+d1*w+e1
rn=S.expand(S.diff(g,w)*S.diff(f,w)-S.diff(f,w,2)*g)
eq(cf(g,w,4),2*a0*a1,'complete general first-row leading term')
eq(cf(S.diff(f,w),w,3),4*a0*a0,'general fprime exact leading term')
eq(cf(rn,w,6),8*a0**3*a1,'general individual-residue leading supplier')
f1,f2,g0,g1=S.symbols('f1 f2 g0 g1')
fw=f1+x*(g1-f2*g0/f1)
eq(S.limit(x*(-1/(x*x*fw)+1/(f1*x*x)),x,0),(g1*f1-f2*g0)/f1**3,
   'implicit displacement retained in the original residue')
Aother=alpha*(x+p-r1)**2*(x+p-r2)**2
eq(S.diff(Aother,x).subs(x,0)/Aother.subs(x,0),2/(p-r1)+2/(p-r2),
   'complete finite 4+2+2 logarithmic derivative')
eq(2/(p-r1)+2/(p-r2),2*(2*p-r1-r2)/((p-r1)*(p-r2)),
   'necessary midpoint equality without a coordinate normalization')
eq(S.diff(Aother,x).subs({x:0,p:(r1+r2)/2}),0,'midpoint leading-jet control')
check(S.diff(Aother,x).subs({x:0,p:0,r1:1,r2:2,alpha:1})!=0,
      'named non-midpoint local residue obstruction')

# Actual global symmetric family, not a translation or nonlinear target shear.
ww=x*x*t
vv=x*(1+ww)
hh=x*x*(1+ww)
Q=hh-d*d*ww
eq(Q,x*x+x*x*(x*x-d*d)*t,'literal global Q')
eq(S.cancel(s*Q.subs(t,1/s)).subs(s,z-x*x),x*x*(z-d*d),'actual Q section box')
basis=[S.Integer(1),t,x*t,ww,vv,hh]
secs=[S.expand(S.cancel(s*v.subs(t,1/s)).subs(s,z-x*x)) for v in basis]
mat=S.Matrix([[cf(cf(v,x,i),z,j) for v in secs] for i in range(3) for j in range(2)])
eq(mat.det(),1,'whole global L1 basis')
rankmon=[x**i for i in range(5)]+[x**3*z,x**4*z,x**3*z*z,x**4*z*z]
rankmat=S.Matrix([[cf(v.subs(z,x*x),x,i) for v in rankmon] for i in range(9)])
eq(rankmat.det(),1,'whole octic restriction rank nine')
check(15-9==6,'fixed-octic correction kernel dimension')
H=alpha*Q*Q+k0+kt*t+kxt*x*t+k2*ww+k3*vv+k4*hh
L=l0+lt*t+lxt*x*t+l2*ww+l3*vv+l4*hh
N=S.cancel(s*s*H.subs(t,1/s))
M=S.cancel(s*L.subs(t,1/s))
eq(N.subs(s,0),alpha*x**4*(x*x-d*d)**2,'full symmetric binary octic')
eq(cf(N.subs(s,0),x,8),alpha,'infinity is not a boundary zero')
eq(S.diff(N,s).subs({x:0,s:0}),kt,'active normal value is actual kt')
eq(S.diff(N,s,x).subs({x:0,s:0}),kxt,'active first normal jet is actual kxt')
eq(M.subs({x:0,s:0}),lt,'active lower value is actual lt')
eq(S.diff(M,x).subs({x:0,s:0}),lxt,'active lower first jet is actual lxt')
for point in (d,-d):
    local=S.expand((alpha*x**4*(x*x-d*d)**2).subs(x,point+s))
    eq(cf(local,s,2),4*alpha*d**6,'each other finite root has exact multiplicity two')
for pp,di,dj in ((N,4,2),(M,2,1)):
    po=S.Poly(S.expand(pp.subs(s,z-x*x)),x,z)
    check(po.degree(x)<=di and po.degree(z)<=dj,'full actual global section bounds')

# All first rows after the necessary finite shared-point jets.
red={kt:0,kxt:0,lt:0,lxt:0}
Fr=S.expand((H*H+L).subs(red))
Fb=S.cancel(Fr.subs(t,w/(x*x)))
Ar=alpha*d**4*w*w+k2*w+k0
fr=Ar*Ar+l2*w+l0
gr=(1+w)*(2*k3*Ar+l3)
eq(cf(Fb,x,0),fr,'complete symmetric zero row')
eq(cf(Fb,x,1),gr,'complete symmetric first row')
eq(jac(x,w/(x*x),x,w),1/(x*x),'actual blow-chart volume')
rp=S.expand(gr-C*S.diff(fr,w))
eq(cf(rp,w,3),2*alpha*d**4*k3-4*alpha**2*d**8*C,'residue scalar leading equation')
rr=S.expand(rp.subs(C,k3/(2*alpha*d**4)))
eq(cf(rr,w,2),k3*(2*alpha*d**4-k2),'residue forces k2')
eq(cf(rr.subs(k2,2*alpha*d**4),w,1),l3,'residue forces l3')
eq(cf(rr.subs({k2:2*alpha*d**4,l3:0}),w,0),-k3*l2/(2*alpha*d**4),
   'residue forces l2')
eq(rr.subs({k2:2*alpha*d**4,l2:0,l3:0}),0,'full residue identity after constraints')
eq(rp.subs({k3:0,C:0}),l3*(1+w),'k3-zero complete residue case')
eq(S.diff(Fr,x).subs({x:0,k3:0,l3:0}),0,'k3-zero actual source Fx')
eq(S.diff(Fr,t).subs({x:0,k3:0,l3:0}),0,'k3-zero actual source Ft')

# Entire polynomial family after the necessary residues.
Hq=c0+q*(alpha*(d*d-x*x)**2*q+k3*x+beta*x*x)
Fq=Hq*Hq+l4*x*x*q+l0
constraints={k2:2*alpha*d**4,l2:0,l3:0,k4:beta-2*alpha*d*d,k0:c0+alpha*d**4}
eq(Fb.subs(w,q-1).subs(constraints),Fq,'complete residual map with every row retained')
eq(S.diff(Fq,x).subs(q,0),0,'q-zero first partial')
eq(S.diff(Fq,q).subs(q,0),x*(2*c0*k3+(2*c0*beta+l4)*x),'q-zero second partial')
xzero=-2*c0*k3/(2*c0*beta+l4)
eq(S.diff(Fq,q).subs({q:0,x:xzero},simultaneous=True),0,'generic q-zero critical supplier')
eq(x*S.diff(Fq,x)-2*q*S.diff(Fq,q),
   -2*q*Hq*(4*alpha*d*d*(d*d-x*x)*q+k3*x),'critical-curve elimination')
qc=-k3*x/(4*alpha*d*d*(d*d-x*x))
L1=S.Rational(3,4)*k3+beta*x+k3*x*x/(4*d*d)
L2=k3/2+beta*x+k3*x*x/(2*d*d)
Hc=c0-k3*x*x*L1/(4*alpha*d*d*(d*d-x*x))
eq(Hq.subs(q,qc),Hc,'complete H on the second critical curve')
eq(S.diff(Fq,x).subs(q,qc),2*qc*(2*Hc*L2+l4*x),'full first partial on curve')
eq(S.diff(Fq,q).subs(q,qc),x*(2*Hc*L2+l4*x),'full second partial on curve')
P5=S.expand(k3*x*L1*L2-2*alpha*d*d*l4*(d*d-x*x))
P6=S.expand(x*x*L1*L2-2*alpha*c0*(d**4-x**4))
eq((2*Hc*L2+l4*x).subs(c0,0),-x*P5/(2*alpha*d*d*(d*d-x*x)),
   'c0-zero exact degree-five critical equation')
eq((2*Hc*L2+l4*x).subs(l4,-2*c0*beta),-k3*P6/(2*alpha*d*d*(d*d-x*x)),
   'exceptional nonzero-c0 exact degree-six critical equation')
eq(cf(P5,x,5),k3**3/(8*d**4),'degree five has nonzero top coefficient')
eq(cf(P6,x,6),k3**2/(8*d**4),'degree six has nonzero top coefficient')
eq(P5.subs(x,0),-2*alpha*d**4*l4,'zero is forbidden but not a P5 root')
eq(P6.subs(x,0),-2*alpha*c0*d**4,'zero is forbidden but not a P6 root')
for sign in (-1,1):
    point=sign*d
    eq(L1.subs(x,point),k3+beta*point,'first factor at forbidden nonzero point')
    eq(L2.subs(x,point),k3+beta*point,'second factor at forbidden nonzero point')
    eq(P5.subs(x,point),k3*point*(k3+beta*point)**2,'P5 forbidden-root condition')
    eq(P6.subs(x,point),d*d*(k3+beta*point)**2,'P6 forbidden-root condition')
    eq(S.diff(P5,x).subs({x:point,beta:-k3/point},simultaneous=True),4*alpha*d*d*l4*point,
       'each actual forbidden P5 root is simple')
    eq(S.diff(P6,x).subs({x:point,beta:-k3/point},simultaneous=True),8*alpha*c0*point**3,
       'each actual forbidden P6 root is simple')
eq(jac(x,1+x*x*t,x,t),x*x,'source chart requires only nonzero x')

# Cheap hostile: some algebraic roots really are forbidden, so choosing any root is invalid.
host5=P5.subs({alpha:1,d:1,k3:1,l4:1,beta:-1})
eq(host5.subs(x,1),0,'P5 has an actual forbidden-root hostile')
eq(S.diff(host5,x).subs(x,1),4,'hostile P5 forbidden root is simple')
host6=P6.subs({alpha:1,d:1,k3:1,c0:1,beta:-1})
eq(host6.subs(x,1),0,'P6 has an actual forbidden-root hostile')
eq(S.diff(host6,x).subs(x,1),8,'hostile P6 forbidden root is simple')

# Original-source critical controls, with all denominators and excluded addresses paid exactly.
sourceF=S.expand(Fq.subs(q,1+x*x*t))
eq(S.diff(sourceF,x).subs({alpha:1,d:1,k3:1,c0:1,beta:0,l4:1,x:-2,t:-S.Rational(1,4)}),0,
   'literal q-zero source Fx')
eq(S.diff(sourceF,t).subs({alpha:1,d:1,k3:1,c0:1,beta:0,l4:1,x:-2,t:-S.Rational(1,4)}),0,
   'literal q-zero source Ft')
controls=[]
for name,params,poly in (
    ('degree5',{alpha:1,d:1,k3:1,c0:0,beta:0,l4:1},P5.subs({alpha:1,d:1,k3:1,beta:0,l4:1})),
    ('degree6',{alpha:1,d:1,k3:1,c0:1,beta:-1,l4:2},S.cancel(host6/(x-1)))):
    poly=S.Poly(poly,x).clear_denoms()[1].as_expr()
    check(S.degree(S.gcd(poly,x*(1-x*x)),x)==0,'control all selected roots have allowed addresses '+name)
    qspec=qc.subs(params);tspec=S.cancel((qspec-1)/(x*x))
    for derivative in (S.diff(sourceF,x),S.diff(sourceF,t)):
        num,den=S.fraction(S.cancel(derivative.subs(params).subs(t,tspec)))
        eq(S.rem(num,poly,x),0,'actual-source derivative vanishes on all control roots '+name)
        check(S.degree(S.gcd(den,poly),x)==0,'control derivative denominator is a unit '+name)
    controls.append([name,str(poly)])

records={'general':'active finite fourfold requires boundary a1=0 in original x coordinate',
         'midpoint':'all-finite 4+2+2 rational mate requires p=(r1+r2)/2',
         'symmetric':'fixed fourfold point0 and doubles+-d; no polynomial mate, d and alpha nonzero',
         'critical_degrees':[5,6],'forbidden':'zero excluded; +-d only simple roots',
         'controls':controls,'hostile':'both critical polynomials may have a simple forbidden root',
         'scope':'no arbitrary-midpoint transport and no symmetric rational exclusion asserted'}
digest=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
print('SUPERSEDED DRAFT controls pass:',gates,'always-active exact gates; analytic consumers unaudited')
print('General residue leading8*a0^3*a1; all-finite 4+2+2 requires midpoint')
print('Fixed-zero symmetric class: exact degree5/6 critical suppliers and simple forbidden roots')
print('Original-source controls:',json.dumps(controls))
print('Proposed scope only: off-midpoint rational exclusion; fixed-zero midpoint polynomial exclusion')
print('Semantic SHA256:',digest)
