#!/usr/bin/env python3
"""Finite exact controls for the all-m moving-index constant-D theorem.
The theorem for every m>=2 is analytic; this file has named finite ranges,
full symbolic numerator boxes, positive rational and negative controls.
"""
import hashlib
import json
import sympy as S

x,t,r,b,z,v=S.symbols('x t r b z v')
D,E,C,lam,e=S.symbols('D E C lambda e')
gates=0
records={}
def check(name,value):
    global gates
    if not bool(value):raise RuntimeError(name)
    gates+=1
def zero(name,expr):check(name,S.cancel(expr)==0)
def jac(f,g):return S.diff(f,x)*S.diff(g,t)-S.diff(f,t)*S.diff(g,x)

# Universal identities, with arbitrary scalar numerator values.
A,B,Q,R,SS,XX=S.symbols('A B Q R SS XX')
N=Q*XX**2+B*XX+A;P=2*Q*XX+B;M=R*XX+SS
zero('D numerator',4*N*Q-P*P-(4*A*Q-B*B))
zero('E numerator',2*N*R-M*P-(B*R*XX+2*A*R-2*Q*SS*XX-B*SS))

# Independent inverse recursion versus the exact Lagrange coefficient formula.
# W=kappa*v*y; the source recurrence has constant new coefficient four.
W=S.Integer(1)
rows={}
for n in range(1,10):
    res=S.expand(W**4+2*D*v*v*W*W+C*v**3*W+(D*D+E)*v**4-1).coeff(v,n)
    wn=-res/4;W=S.expand(W+wn*v**n);rows[n]=S.expand(wn)
    # Truncate only after the n-th exact coefficient has been paid.
    W=sum(S.expand(W).coeff(v,j)*v**j for j in range(n+1))
    zero('inverse recurrence '+str(n),S.expand(W**4+2*D*v*v*W*W+C*v**3*W+(D*D+E)*v**4-1).coeff(v,n))

def lagrange(j):
    total=0;h=S.Rational(j,4)
    for aa in range((j+1)//2+1):
        for bb in range((j+1)//3+1):
            remain=j+1-2*aa-3*bb
            if remain<0 or remain%4:continue
            cc=remain//4;n=aa+bb+cc
            coeff=S.prod(h-i for i in range(n))/(S.factorial(aa)*S.factorial(bb)*S.factorial(cc))
            total+=coeff*(2*D)**aa*C**bb*(D*D+E)**cc
    return S.expand(-total/j)
for j in range(1,9):zero('Lagrange independent '+str(j),lagrange(j)-rows[j+1])
zero('T1',rows[2]+D/2)
zero('T5',rows[6]+(2*D**3+4*D*E+C*C)/32)

# All m=2..9: exhaustive monomial weight checks at both claimed indices.
# A monomial D^k C^q E^ell has z-degree2k+3q+4ell.
weight_records=[]
for m in range(2,10):
    j1=2*m-3;j2=6*m-7;top=2*m-2
    first=[];second=[]
    for k in range((j1+1)//2+1):
        for q in range((j1+1)//3+1):
            rest=j1+1-2*k-3*q
            if rest<0 or rest%4:continue
            ll=rest//4;weight=2*k+q+2*ll
            check('first weight bound '+str((m,k,q,ll)),weight<=top)
            if weight==top:first.append((k,q,ll))
    check('first unique pure D '+str(m),first==[(m-1,0,0)])
    for k in range((j2+1)//2+1):
        for q in range((j2+1)//3+1):
            rest=j2+1-2*k-3*q
            if rest<0 or rest%4:continue
            ll=rest//4;weight=q+ll
            check('second weight bound '+str((m,k,q,ll)),weight<=top)
            if weight==top:second.append((k,q,ll))
    check('second unique pure C '+str(m),second==[(0,top,0)])
    dcoef=-S.binomial(S.Rational(j1,2),m-1)/j1
    exact=S.binomial(S.Rational(1,2),m-1)*(-1)**(m-1)
    zero('pure D closed coefficient '+str(m),dcoef-exact)
    ccoef=-S.binomial(S.Rational(j2,4),top)/j2
    check('first nonzero coefficient '+str(m),dcoef!=0)
    check('second nonzero coefficient '+str(m),ccoef!=0)
    weight_records.append([m,j1,j2,str(dcoef),str(ccoef)])
records['weights']=weight_records

# All original global coefficient parameters, m=2..5, with exact leading N.
# Monicity is imposed only on the one leading coefficient that it determines.
for m in range(2,6):
    avec=S.symbols('a0:'+str(2*m+1));bvec=S.symbols('b0:'+str(2*m+1))
    qvec=S.symbols('q0:'+str(2*m-1));rvec=S.symbols('r0:'+str(m+1));svec=S.symbols('s0:'+str(m+1))
    ap=sum(avec[i]*x**i for i in range(2*m+1));bp=sum(bvec[i]*x**i for i in range(2*m+1))
    qp=sum(qvec[i]*x**i for i in range(2*m-1))
    # For m=2 the Bx^m term ties the leading Qx^(2m) term.
    qp=qp.subs(qvec[-1],1-bvec[-1] if m==2 else 1)
    rp=sum(rvec[i]*x**i for i in range(m+1));sp=sum(svec[i]*x**i for i in range(m+1))
    np=S.expand(qp*x**(2*m)+bp*x**m+ap);pp=S.expand(2*qp*x**m+bp);mp=S.expand(rp*x**m+sp)
    hp=np*t*t+pp*t+qp;lp=mp*t+rp
    check('full N degree '+str(m),S.Poly(np,x).degree()==4*m-2)
    zero('full N leading '+str(m),np.coeff(x,4*m-2)-1)
    ht=S.expand(hp.subs({x:1/r,t:-r**m-r**(2*m)*b},simultaneous=True))
    lt=S.expand(lp.subs({x:1/r,t:-r**m-r**(2*m)*b},simultaneous=True))
    check('H whole chart '+str(m),S.denom(S.cancel(ht))==1)
    check('L whole chart '+str(m),S.denom(S.cancel(lt))==1)
    zero('H exact boundary '+str(m),ht.subs(r,0)-bvec[-1]*b-avec[-1])
    zero('L exact boundary '+str(m),lt.subs(r,0)+rvec[-1]*b+svec[-1])
    zero('actual volume '+str(m),S.det(S.Matrix([[S.diff(1/r,r),0],[S.diff(-r**m-r**(2*m)*b,r),-r**(2*m)]]))-r**(2*m-2))
    dn=S.expand(4*ap*qp-bp*bp);en=S.expand(bp*rp*x**m+2*ap*rp-2*qp*sp*x**m-bp*sp)
    zero('original D numerator '+str(m),4*np*qp-pp*pp-dn)
    zero('original E numerator '+str(m),2*np*rp-mp*pp-en)
    check('first D degree '+str(m),S.Poly(dn,x).degree()<=4*m)
    check('first E degree '+str(m),S.Poly(en,x).degree()<=4*m)
    zero('D leading square '+str(m),dn.coeff(x,4*m)+bvec[-1]**2)
    dnext=S.expand(dn.subs(bvec[-1],0));enext=S.expand(en.subs(bvec[-1],0))
    check('bounded D '+str(m),S.Poly(dnext,x).degree()<=4*m-2)
    check('linear E '+str(m),S.Poly(enext,x).degree()<=4*m-1)
    check('linear C '+str(m),S.Poly(mp,x).degree()<=2*m)
    efin=S.expand(enext.subs(rvec[-1],0));mfin=S.expand(mp.subs(rvec[-1],0))
    check('final bounded E '+str(m),S.Poly(efin,x).degree()<=4*m-2)
    check('final bounded C '+str(m),S.Poly(mfin,x).degree()<=2*m-1)

# Sharp actual global controls, including nonconstant L, at every named m.
for m in range(2,10):
    y=x**(m-1)*(1+x**m*t)
    gy=1/((2*m-2)*x**(2*m-2))
    zero('primitive y '+str(m),jac(y,gy)-1)
    f=y**4+lam*y+e;g=gy/(4*y**3+lam)
    # Product/chain path keeps exact expressions small and verifies source J.
    zero('rational mate '+str(m),jac(f,g)-1)
    yt=S.expand(y.subs({x:1/r,t:-r**m-r**(2*m)*b},simultaneous=True))
    zero('sharp y global '+str(m),yt+r*b)
    nh=x**(4*m-2)-1
    badh=y*y-t*t;badl=t
    zero('negative true leading '+str(m),S.expand(badh).coeff(t,2)-nh)
    zero('negative simple root '+str(m),nh.subs(x,1))
    zero('negative residue '+str(m),-1/(4*S.diff(nh,x).subs(x,1))+S.Rational(1,4*(4*m-2)))
    check('negative residue nonzero '+str(m),S.Rational(-1,4*(4*m-2))!=0)
records['range']={'weights':[2,9],'full_boxes':[2,5],'rational_and_negative':[2,9]}
records['scope']='Analytic all-m>=2 degree-(4m-2) gate; finite controls do not assert finite sufficiency.'
blob=json.dumps(records,sort_keys=True,separators=(',',':')).encode()
print('all_m_constant_d exact controls: PASS')
print('gates:',gates)
print('moving indices: 2m-3 and 6m-7; finite weight controls m=2..9')
print('full global boxes m=2..5; original-volume and square-field controls: PASS')
print('all-m rational positive and constant-D negative families: PASS')
print('semantic sha256:',hashlib.sha256(blob).hexdigest())
