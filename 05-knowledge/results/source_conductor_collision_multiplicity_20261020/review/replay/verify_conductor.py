#!/usr/bin/env python3
"""Exact controls for the actual-source conductor and norm-square gate.

Not a search for noninvertible Keller maps. Every positive Keller control is
an automorphism; every collision-bearing example is explicitly non-Keller.
"""
from pathlib import Path
import json
import sympy as S
x,y,u,v,z,t=S.symbols('x y u v z t')
labels=[]
records=[]
def check(value, label):
    if not bool(value):
        raise RuntimeError(label)
    labels.append(label)
def zero(expr):
    return S.cancel(expr)==0
def sub(expr,P,Q):
    return S.expand(expr.subs({u:P,v:Q,z:x}, simultaneous=True))
def matrix_for(P,Q):
    n=S.degree(P,y)
    columns=[]
    for j in range(n):
        rem=S.rem(S.expand(Q**j),P-u,y)
        columns.append([S.expand(rem).coeff(y,i) for i in range(n)])
    return S.Matrix.hstack(*(S.Matrix(c) for c in columns))
def norm(P,c):
    return S.factor(S.resultant(P-u,c,y))
def d_resultant(P,Q):
    U=S.cancel((P.subs(y,t)-P)/(t-y)).expand()
    V=S.cancel((Q.subs(y,t)-Q)/(t-y)).expand()
    if U==1:
        return S.Integer(1)
    # Monic U permits polynomial reduction without changing the resultant.
    return S.factor(S.resultant(U,S.rem(V,U,t),t))
# Positive controls: actual Keller automorphisms, including J != 1.
for r in [1,2,3,4]:
    for k in [1,2,-1]:
        for m in [2,3]:
            P=x+y**r
            Q=k*y+P**m
            R=(v-u**m)**r+k**r*(z-u)
            J=S.expand(S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x))
            c=S.cancel(sub(S.diff(R,z),P,Q)/J)
            M=matrix_for(P,Q)
            index=S.factor(M.det())
            tag=f'Keller r={r}, k={k}, m={m}'
            check(J==k,tag+' Jacobian')
            check(zero(sub(R,P,Q)),tag+' graph relation')
            check(zero(sub(S.diff(R,v),P,Q)-c*S.diff(P,y)),tag+' Rv=conductor*Py')
            check(zero(sub(S.diff(R,u),P,Q)+c*S.diff(Q,y)),tag+' Ru=-conductor*Qy')
            check(zero(c-k**(r-1)),tag+' unit conductor')
            check(zero(norm(P,c)-index**2),tag+' global norm square')
            check(zero(d_resultant(P,Q)-c),tag+' D-resultant and collision determinant')
            records.append({'family':'Keller automorphism','r':r,'k':k,'m':m,'conductor':str(c),'index':str(index)})
# Immersed rational node/tacnode families; same two source parameters collide.
for k in range(1,7):
    P=y**2
    Q=x+y*(y*y-1)**k
    R=(v-z)**2-u*(u-1)**(2*k)
    J=S.expand(S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x))
    c=S.cancel(sub(S.diff(R,v),P,Q)/S.diff(P,y))
    pulled_Rz=sub(S.diff(R,z),P,Q)
    M=matrix_for(P,Q)
    index=S.factor(M.det())
    tag=f'immersed non-Keller contact order {k}'
    check(zero(sub(R,P,Q)),tag+' finite graph relation')
    check(J==-2*y,tag+' explicitly not Keller')
    check(S.gcd(S.diff(P,y),S.diff(Q,y))==1,tag+' every slice immersed')
    check(zero(c-(y*y-1)**k),tag+' exact conductor')
    check(zero(pulled_Rz-J*c),tag+' derivative has extra Jacobian factor')
    check(zero(norm(P,c)-index**2),tag+' norm equals square index')
    check(zero(d_resultant(P,Q)-c),tag+' D-resultant and collision determinant')
    check(S.degree(c,y)==2*S.degree(index,u)==2*k,tag+' conductor degree twice delta')
    check(S.degree(pulled_Rz,y)==2*k+1,tag+' odd Rz degree rejects Keller')
    check(P.subs(y,1)==P.subs(y,-1) and Q.subs(y,1)==Q.subs(y,-1),tag+' two colliding parameters')
    check(S.Poly(c,y).eval(1)==0,tag+' positive conductor support')
    records.append({'family':'immersed non-Keller','contact_order':k,'conductor':str(S.factor(c)),'index':str(index),'degree_Rz':int(S.degree(pulled_Rz,y))})
# A cusp keeps degree parity and norm-square, but violates immersion.
P=y**2; Q=x+y**3; R=(v-z)**2-u**3
c=S.cancel(sub(S.diff(R,v),P,Q)/S.diff(P,y))
check(c==y*y,'cusp conductor has repeated root')
check(S.diff(P,y).subs(y,0)==S.diff(Q,y).subs(y,0)==0,'cusp violates immersion')
check(norm(P,c)==matrix_for(P,Q).det()**2,'cusp passes norm-square: not sufficient')
# Node-to-cusp degeneration: every conductor degree is 2, but discriminant vanishes.
P=y*y; Q=y*(y*y-x); R=v*v-u*(u-z)**2
c=S.cancel(sub(S.diff(R,v),P,Q)/S.diff(P,y))
check(zero(sub(R,P,Q)),'node-to-cusp graph relation')
check(c==y*y-x,'node-to-cusp conductor')
check(S.discriminant(c,y)==4*x,'node-to-cusp discriminant not a unit')
check(S.diff(P,y).subs({x:0,y:0})==S.diff(Q,y).subs({x:0,y:0})==0,'node-to-cusp degeneration violates immersion')
check(zero(norm(P,c)-matrix_for(P,Q).det()**2),'node-to-cusp still satisfies index square')
# Numerical coefficient filters are necessary, not realization certificates.
odd_drop=x*y*y+y+1
check(S.degree(odd_drop,y)==2 and S.degree(odd_drop.subs(x,0),y)==1,'odd specialization drop forbidden')
even_drop=x*x*y**4+y*y-1
check(S.degree(even_drop,y)==4 and S.degree(even_drop.subs(x,0),y)==2,'even specialization drop passes parity only')
check(norm(y*y,y)==-u,'non-square norm excludes a putative conductor y')
# Exact nodal equalizer division under arbitrary polynomial shears.
for q in [0,x,x*x+1]:
    h=(y+q)**2-1
    for N in range(1,9):
        f=(x+2)*y**N+x*y+3
        rem=S.rem(f,h,y)
        r0=S.expand(rem).coeff(y,0)
        r1=S.expand(rem).coeff(y,1)
        diff=S.expand(f.subs(y,-q+1)-f.subs(y,-q-1))
        check(zero(diff-2*r1),f'equalizer shear {q}, degree {N}')
        check(zero((f-rem)-h*S.div(f,h,y)[0]),f'equalizer division shear {q}, degree {N}')
output={'status':'PASS','checks':len(labels),'labels':labels,'records':records,
        'scope':'Exact controls and symbolic identities only. No noninvertible Keller map, global degree exclusion, or Jacobian-conjecture proof.'}
Path(__file__).with_name('results.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','checks':len(labels),'families':len(records)},sort_keys=True))
