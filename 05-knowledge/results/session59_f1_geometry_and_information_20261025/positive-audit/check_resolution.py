#!/usr/bin/env python3
"""Independent resolution controls; no producer source imports."""
from pathlib import Path
from itertools import product
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
checks=0
def check(v,m):
    global checks
    checks+=1
    if not v:raise RuntimeError(m)

s,t=S.symbols('s t')
c=s*s+t*t;a=s*s-t*t+2*s*t;b=s*s-t*t-2*s*t
check(S.expand(a*a+b*b-2*c*c)==0,'conic parameterization')
V=S.Matrix([[1,0,1],[1,2,-1],[1,-2,-1]])
check(abs(V.det())==8,'conic quadratic basis')

# Resolved fan: primitive cyclic rays, self-intersection, conjugation, Picard quotient.
rays=[(1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1)]
dets=[];selfints=[]
for j,v in enumerate(rays):
    vp=rays[(j-1)%8];vn=rays[(j+1)%8]
    det=v[0]*vn[1]-v[1]*vn[0]
    check(det==1,'smooth cone')
    k=S.Matrix([vp[0]+vn[0],vp[1]+vn[1]])
    factor=next(k[i]/v[i] for i in range(2) if v[i])
    check(k==factor*S.Matrix(v),'boundary self intersection relation')
    selfints.append(-int(factor));dets.append(det)
check(selfints==[-2,-1,-2,-1,-2,-1,-2,-1],'exceptional and boundary classes')
M=S.Matrix(rays)
P=S.zeros(8)
for j in range(8):P[(j+4)%8,j]=1
check(P*M==-M,'character lattice central inversion')
check(P*P==S.eye(8),'Galois involution')
SNF=smith_normal_form(M,domain=ZZ)
check([SNF[i,i] for i in range(2)]==[1,1],'Picard quotient torsion free')
trivial=(8+S.trace(P))//2
sign=(8-S.trace(P))//2
check((trivial,sign-2)==(4,2),'Picard representation')

# Blowup presentation Picard basis H1,H2,E00,E01,E10,E11.
Q=S.diag(1,1,0,0,0,0)
for j in range(4):Q[2+(3-j),2+j]=1
check(Q*Q==S.eye(6),'blowup Picard action')
check(S.trace(Q)==2,'nonsplit Picard trace')
check(S.trace(Q*Q)==6,'split extension Picard trace')
intersection=S.diag(0,0,-1,-1,-1,-1)
intersection[0,1]=intersection[1,0]=1
minusK=S.Matrix([2,2,-1,-1,-1,-1])
check((minusK.T*intersection*minusK)[0]==4,'anticanonical degree')
fibers=[S.Matrix(v) for v in [(1,0,-1,-1,0,0),(1,0,0,0,-1,-1),(0,1,-1,0,-1,0),(0,1,0,-1,0,-1)]]
for j,f in enumerate(fibers):
    check((f.T*intersection*f)[0]==-2,'strict fiber self intersection')
    check((minusK.T*intersection*f)[0]==0,'anticanonical contraction')
    for g in fibers[j+1:]:check((f.T*intersection*g)[0]==0,'disjoint exceptional fibers')

# Direct conic and infinity enumeration in prime fields, then actual blowup fiber count.
rows=[]
for q in (3,5,7,11,13,17,19,23,29,31):
    conic=[]
    for first in range(3):
        for tail in product(range(q),repeat=2-first):
            c,a,b=(0,)*first+(1,)+tail
            if (a*a+b*b-2*c*c)%q==0:conic.append((c,a,b))
    infinity=[v for v in conic if v[0]==0]
    centers=list(product(infinity,repeat=2))
    count=len(conic)**2-len(centers)+len(centers)*(q+1)
    epsilon=1 if (q-1) in {x*x%q for x in range(q)} else -1
    check(len(conic)==q+1,'conic point count')
    check(len(infinity)==1+epsilon,'infinity splitting')
    check(count==q*q+(4+2*epsilon)*q+1,'blowup count')
    for n in range(1,7):
        Qn=q**n;ep=epsilon**n
        node_count=2*(1+ep)
        from_resolution=(Qn+1)**2+Qn*node_count
        from_picard=Qn*Qn+(4+2*ep)*Qn+1
        check(from_resolution==from_picard,'extension and Picard consistency')
    rows.append(dict(q=q,epsilon=epsilon,conic_points=len(conic),infinity_points=len(infinity),rational_blowup_centers=len(centers),resolution_points=count))

T,q=S.symbols('T q')
for epsilon in (-1,1):
    denominator=(1-T)*(1-q*q*T)*(1-q*T)**4*(1-epsilon*q*T)**2
    log_derivative=-T*S.diff(denominator,T)/denominator
    series=S.series(log_derivative,T,0,6).removeO().expand()
    for n in range(1,6):check(S.expand(series.coeff(T,n)-(q**(2*n)+(4+2*epsilon**n)*q**n+1))==0,'zeta power sums')

out=dict(status='PASS',always_active_checks=checks,conic_basis_determinant=int(V.det()),fan_rays=rays,fan_self_intersections=selfints,picard_trivial_rank=int(trivial),picard_sign_rank=int(sign-2),direct_prime_field_blowup_counts=rows)
Path(__file__).with_name('resolution_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,sort_keys=True))
