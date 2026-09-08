#!/usr/bin/env python3
"""Exact full-section transport and forbidden-root certificates for 4+4.

No inherited mathematical implementation is imported. The finite universe
is all 21 quintic and four cubic forbidden-root multiplicity patterns.
"""
import hashlib
import json
import sympy as S

u,t,s,z,p,x,r,T,c,wv,q,Q,s0=S.symbols('u t s z p x r T c wv q Q s0')
alpha,k0,kt,kxt,k2,k3,k4=S.symbols('alpha k0 kt kxt k2 k3 k4')
l0,lt,lxt,l2,l3,l4=S.symbols('l0 lt lxt l2 l3 l4')
c0,Bc,Cc,Lc,Cscalar=S.symbols('c0 B C L Cscalar')
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

# Full global source family, not an automorphism claim for translation.
w=u*u*t
v=u*(1+w)-2*p
h=u*u*(1+w)-2*p*u
eq(h,u*v,'corrected full-field carrier relation')
basis=[S.Integer(1),t,u*t,w,v,h]
sections=[S.expand(S.cancel(s*bb.subs(t,1/s)).subs(s,z-(u+p)**2)) for bb in basis]
mat=S.Matrix([[coeff(coeff(bb,u,i),z,j) for bb in sections]
              for i in range(3) for j in range(2)])
eq(mat.det(),1,'all-p full L1 section determinant')
for bb in sections:
    check(S.Poly(bb,u,z).degree(u)<=2 and S.Poly(bb,u,z).degree(z)<=1,
          'each corrected section satisfies original global box')
# This triangular minor has determinant one for every p, not merely generically.
rankmon=[u**i for i in range(5)]+[u**3*z,u**4*z,u**3*z*z,u**4*z*z]
rankmat=S.Matrix([[coeff(bb.subs(z,(u+p)**2),u,i) for bb in rankmon] for i in range(9)])
eq(rankmat.det(),1,'all-p octic restriction rank-nine witness')
check(15-9==6,'whole correction kernel dimension')
H=alpha*w*w+k0+kt*t+kxt*u*t+k2*w+k3*v+k4*h
L=l0+lt*t+lxt*u*t+l2*w+l3*v+l4*h
N=S.cancel(s*s*H.subs(t,1/s))
M=S.cancel(s*L.subs(t,1/s))
eq(N.subs(s,0),alpha*u**4,'exact finite-p fourth-order root')
eq(N,alpha*u**4+s*(kt+kxt*u+k2*u*u+k3*u**3+k4*u**4)
   +s*s*(k0-2*p*k3+(k3-2*p*k4)*u+k4*u*u),'all corrected H normal jets')
eq(M,lt+lxt*u+l2*u*u+l3*u**3+l4*u**4
   +s*(l0-2*p*l3+(l3-2*p*l4)*u+l4*u*u),'all corrected L normal jets')
for pp,di,dj,label in [(N,4,2,'H'),(M,2,1,'L')]:
    po=S.Poly(S.expand(pp.subs(s,z-(u+p)**2)),u,z)
    check(po.degree(u)<=di and po.degree(z)<=dj,'full actual section box '+label)
# Literal actual W inversion, with the complete unit factor at infinity.
Hinv=S.cancel(H.subs({u:1/r-p,t:-r*r-r**4*T},simultaneous=True))
Linv=S.cancel(L.subs({u:1/r-p,t:-r*r-r**4*T},simultaneous=True))
check(Hinv.is_polynomial(r,T) and Linv.is_polynomial(r,T),'actual inversion retains both global functions')
Ninv=S.cancel(s*s*Hinv.subs(T,1/s))
Minv=S.cancel(s*Linv.subs(T,1/s))
eq(Ninv.subs(s,0),alpha*r**4*(1-p*r)**4,'entire transformed boundary unit factor')
eq(coeff(Ninv.subs(s,0),r,4),alpha,'actual infinity leading coefficient never vanishes')
eq(jac(1/r,-r*r-r**4*T,r,T),r*r,'actual weighted original volume')

# The m4 infinity lemma with all higher unit terms, not just a pure monomial.
A5,A6,A7,A8,b0,b1,b2,b3,b4,b5,b6,d0,d1,d2,d3,d4,e0,e1,e2,cc0,cc1,cc2,cc3,cc4=S.symbols(
    'A5 A6 A7 A8 b0 b1 b2 b3 b4 b5 b6 d0 d1 d2 d3 d4 e0 e1 e2 cc0 cc1 cc2 cc3 cc4')
Na=alpha*r**4+A5*r**5+A6*r**6+A7*r**7+A8*r**8
Bn=b0+b1*r+b2*r*r+b3*r**3+b4*r**4+b5*r**5+b6*r**6
Cn=cc0+cc1*r+cc2*r*r+cc3*r**3+cc4*r**4
Da=d0+d1*r+d2*r*r+d3*r**3+d4*r**4
Ma=Da+s*(e0+e1*r+e2*r*r)
NN=Na+s*Bn+s*s*Cn
EE=S.expand(NN*NN+s**3*Ma-c*s**4)
eq(coeff(EE.subs(r,0),s,2),b0*b0,'higher-unit normal degree two unchanged')
eq(EE.subs({r:0,b0:0}),s**3*(d0+(cc0*cc0+e0-c)*s),'higher-unit M degree three unchanged')
E4=S.expand(EE.subs({b0:0,d0:0}))
eq(E4.subs(r,0),(cc0*cc0+e0-c)*s**4,'higher-unit remaining degree four')
low=(b1+cc0*Q)**2+d1*Q+(e0-c)*Q*Q
eq(coeff(E4.subs(s,r*Q),r,4),Q*Q*low,'higher-unit j1 low face unchanged')
eq(coeff(NN.subs(b0,0).subs(s,r**3*Q),r,4),alpha+b1*Q,'higher-unit j1 centre unchanged')
mc=S.expand((Ma-c*s).subs({d0:0,d1:0,d2:0,s:-alpha*r**3/b1},simultaneous=True))
eq(S.diff(coeff(mc,r,3),c),alpha/b1,'generic centre order cutoff three survives all higher unit terms')
sc=E4.subs(b1,0).subs({r:s0**3,s:s0**7*Q},simultaneous=True)
eq(coeff(sc,s0,24),alpha*alpha+d1*Q**3,'higher-unit n1 cubic unchanged')
eq(coeff(S.diff(E4.subs(b1,0),s).subs({r:s0**3,s:s0**7*Q},simultaneous=True),s0,17),
   3*d1*Q*Q,'higher-unit n1 derivative leading order unchanged')
face=(alpha+b2*Q+cc0*Q*Q)**2+d2*Q**3+(e0-c)*Q**4
eq(coeff(E4.subs({b1:0,d1:0}).subs(s,r*r*Q),r,8),face,
   'higher-unit complete tangent quartic unchanged')
face0=S.expand(face+c*Q**4)
eq((Q*S.diff(face0,Q)-4*face0).subs(Q,0),-4*alpha*alpha,'generic quartic simple-root supplier')

# Full finite moving-root residue and every shifted coefficient.
red={kt:0,kxt:0,lt:0,lxt:0}
Hr=S.expand(H.subs(red))
Lr=S.expand(L.subs(red))
F=S.expand(Hr*Hr+Lr)
Fuw=S.cancel(F.subs(t,wv/(u*u)))
A=alpha*wv*wv+k2*wv+k0-2*p*k3
f=A*A+l2*wv+l0-2*p*l3
g=2*A*(k3*(1+wv)-2*p*k4)+l3*(1+wv)-2*p*l4
eq(coeff(Fuw,u,0),f,'complete shifted zero-row')
eq(coeff(Fuw,u,1),g,'complete shifted first correction')
eq(jac(u,wv/(u*u),u,wv),1/(u*u),'actual finite-point blow-chart volume')
respoly=S.expand(g-Cscalar*S.diff(f,wv))
eq(coeff(respoly,wv,3),2*alpha*k3-4*alpha*alpha*Cscalar,'top coefficient fixes residue scalar')
k2forced=2*alpha-4*alpha*p*k4/k3
l2forced=-4*alpha*p*l4/k3
pp=S.expand(respoly.subs(Cscalar,k3/(2*alpha)))
eq(coeff(pp,wv,2),k3*(2*alpha-k2)-4*alpha*p*k4,'shifted k2 constraint')
eq(coeff(pp.subs(k2,k2forced),wv,1),l3,'shifted l3 constraint')
eq(coeff(pp.subs({k2:k2forced,l3:0}),wv,0),-2*p*l4-k3*l2/(2*alpha),'shifted l2 constraint')
eq(pp.subs({k2:k2forced,l3:0,l2:l2forced}),0,'entire residue identity after all constraints')
eq(S.diff(F,t).subs(u,0),0,'source Ft at finite exceptional line')
eq(S.diff(F,u).subs({u:0,k3:0}),g.subs({wv:0,k3:0}),
   'k3-zero g0 is exactly source Fu at exceptional line')
eq(coeff(respoly.subs(k3,0),wv,3),-4*alpha*alpha*Cscalar,
   'k3-zero forces scalar zero independently of other coefficients')

# The full residual polynomial map and its actual critical fold.
aa=4*alpha/k3
zz=(u*u-p*aa)*q-2*p*u
PP=alpha*q*q+k3*u*q
c0actual=k0-alpha-2*p*k3+p*aa*k4
Hq=PP+k4*zz+c0
Fq=Hq*Hq+l4*zz+l0-l2forced
Fres=S.cancel(Fuw.subs({wv:q-1,k2:k2forced,l3:0,l2:l2forced}))
eq(Fres,Fq.subs(c0,c0actual),'no higher induced row dropped from residual map')
eq(jac(PP,zz,u,q),u*(2*p*k3-q*(k3*u+4*alpha*q)),'entire critical-map Jacobian')
uc=2*p/q-aa*q
Pc=2*p*k3-3*alpha*q*q
zc=aa*aa*q**3-3*p*aa*q
eq(PP.subs(u,uc),Pc,'critical-curve P value')
eq(zz.subs(u,uc),zc,'critical-curve z value')
Hc=Pc+k4*zc+c0
R=S.expand(Hc*(k3*q+2*k4*(p-aa*q*q))+l4*(p-aa*q*q))
eq(S.diff(Fq,u).subs(u,uc),2*R,'actual Fu equals twice critical polynomial')
eq(S.diff(Fq,q).subs(u,uc),(4*p/q**2-aa)*R,'actual Fq exact critical polynomial multiplier')
eq(jac(u,1+u*u*t,u,t),u*u,'source inverse chart regular exactly at u nonzero')

# Normalization retains every nonzero p, alpha, k3, l4 value by choosing a square root.
ps=aa*s0*s0/2
normalize={p:ps,q:s0*Q,k4:Bc*k3/(2*aa*s0),l4:Lc*k3*k3*s0/2,
           c0:ps*k3*(Cc-4)/2}
Rbar=S.expand((Cc-3*Q*Q+Bc*(2*Q**3-3*Q))*(Q+Bc*(1-2*Q*Q)/2)+Lc*(1-2*Q*Q))
eq(R.subs(normalize,simultaneous=True),ps*k3*k3*s0*Rbar/2,'whole normalized polynomial including scale')
eq(uc.subs({p:ps,q:s0*Q},simultaneous=True),aa*s0*(1/Q-Q),'forbidden source addresses exactly Q0 and Qplusminus1')
eq(coeff(Rbar,Q,5),-2*Bc*Bc,'quintic leading coefficient')
eq(coeff(Rbar,Q,4),5*Bc,'quintic next coefficient')
eq(coeff(Rbar,Q,3),4*Bc*Bc-3,'quintic third coefficient')

# Exhaustive forbidden-root multiplicity universe; no numerical root solver.
patterns=[]
for e in range(6):
    for rp in range(6-e):
        sm=5-e-rp
        dd=rp-sm
        vv=rp+sm
        target=-2*Bc*Bc*Q**e*(Q-1)**rp*(Q+1)**sm
        eq(coeff(target,Q,4),2*Bc*Bc*dd,'quintic Q4 pattern '+str((e,rp,sm)))
        eq(coeff(target,Q,3),-Bc*Bc*(dd*dd-vv),'quintic Q3 pattern '+str((e,rp,sm)))
        if dd==0:
            eq(coeff(Rbar-target,Q,4),5*Bc,'d-zero forces forbidden B-zero '+str((e,rp,sm)))
            reason='nonzero 5B'
        else:
            residual=S.cancel(4*dd*dd*coeff(Rbar-target,Q,3).subs(Bc,S.Rational(5,2*dd)))
            eq(residual,100+13*dd*dd-25*vv,'entire integer obstruction '+str((e,rp,sm)))
            check(residual!=0,'nonzero obstruction excludes pattern '+str((e,rp,sm)))
            reason=str(residual)
        patterns.append([e,rp,sm,reason])
check(len(patterns)==21,'all and only 21 quintic multiplicity patterns')
R0=S.expand(Rbar.subs(Bc,0))
eq(R0,-3*Q**3-2*Lc*Q*Q+Cc*Q+Lc,'complete B-zero cubic')
cubic=[]
for rp in range(4):
    sm=3-rp
    target=-3*(Q-1)**rp*(Q+1)**sm
    forced=-3*(-1)**rp
    eq(target.subs(Q,0),forced,'cubic constant fixes nonzero L '+str(rp))
    mismatch=S.expand(R0-target).subs(Lc,forced)
    eq(coeff(mismatch,Q,2),6*(-1)**rp-3*(rp-sm),'cubic parity obstruction '+str(rp))
    check(coeff(mismatch,Q,2)!=0,'all forbidden cubic patterns excluded '+str(rp))
    cubic.append([rp,sm,int(coeff(mismatch,Q,2))])
check(len(cubic)==4,'all four cubic patterns')
# Controls pay the necessary saturation rather than pretending L=0 was covered.
control=S.Poly(Rbar.subs({Bc:0,Cc:0,Lc:1}),Q)
check(control.degree()==3,'named nondegenerate cubic')
check(all(control.eval(vv)!=0 for vv in [-1,0,1]),'named cubic all roots have allowed chart addresses')
eq(Rbar.subs({Bc:0,Cc:3,Lc:0}),-3*Q*(Q-1)*(Q+1),'L-zero hostile only has forbidden roots')

records['entry']='complete corrected all-p L1; alpha*(x-p)^4; nonzero p proof'
records['infinity']='unit alpha*r^4*(1-p*r)^4 retained with actual volume r^2'
records['residue']=['k2=2alpha-4alpha*p*k4/k3','l3=0','l2=-4alpha*p*l4/k3']
records['critical_map']='F=(P+k4*z+c0)^2+l4*z+const; every allowed root gives Fu=Fq=0'
records['quintic_patterns']=patterns
records['cubic_patterns']=cubic
records['hostile']='B=0,C=3,L=0 gives only forbidden roots; excluded before normalization'
semantic=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
print('PASS shifted finite/infinity 4+4:',gates,'always-active exact gates')
print('All-p full global basis and actual critical-map transport checked')
print('Quintic forbidden-root patterns:21; cubic patterns:4; all excluded')
print('Integer obstruction:13(r-s)^2=25(r+s)-100 has no allowed solution')
print('Scope: nonzero finite p polynomial-mate exclusion; L-zero hostile retained')
print('Semantic SHA256:',semantic)
