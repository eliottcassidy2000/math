#!/usr/bin/env python3
"""Exact regression controls; finite matrix checks do not replace PROOF.md."""
import json
from pathlib import Path
import sympy as S
x,y,T=S.symbols('x y T')

def require(v,msg):
    if not v: raise RuntimeError(msg)
def bracket(P,Q):
    return S.expand(S.diff(P,x)*S.diff(Q,y)-S.diff(P,y)*S.diff(Q,x))
def zero(f): return S.cancel(f)==0

identities=0
for m in range(2,10):
    for s in range(1,9):
        # Translated/scaled cases ensure this is not just monomial-at-origin code.
        a,b,beta=3,2,5
        z=y-b; P=a*z**m*x+beta+z**s
        e=(m-2+s)//s
        for k in (e,e+1):
            R=(P-beta)**k
            Q=S.cancel(-R/(a*(m-1)*z**(m-1)))
            require(S.denom(Q)==1 or not S.denom(Q).has(x,y),'polynomial descent')
            require(zero(bracket(P,Q)-R),'Poisson identity')
            identities+=1
        if e>1:
            bad=S.cancel(-(P-beta)**(e-1)/(a*(m-1)*z**(m-1)))
            require(S.denom(bad).has(y),'sharp threshold failed')
for m in range(2,10):
    P=2*(y-3)**m*x+7
    Q=S.cancel(-(P-7)/(2*(m-1)*(y-3)**(m-1)))
    require(zero(bracket(P,Q)-(P-7)),'constant B identity')
    identities+=1
for B in [0,y,y**3+2*y+1]:
    P=3*x+B
    for R in [S.Integer(1),P,P**2+1]:
        Q=y*R/3+P**3
        require(zero(bracket(P,Q)-R),'constant A classification')
        identities+=1

# Independent direct coefficient-matrix test of image intersection.
# No integration or candidate companion formula enters construction of M.
# rank(M)+k-rank([M|N]) is dimension of span(N) intersect image(M).
cases=[
    (S.Integer(2),y**3+1,3,'constant'),
    (y,y,0,'simple root'),
    (y**2,y,2,'one root m2'),
    (y**3,y,1,'one root m3 sharp'),
    (y**5,y,0,'one root m5 missing powers'),
    (y**5,y**2,1,'one root m5 B order2'),
    (y**5,y**4,2,'one root m5 B order4'),
    (y**5,S.Integer(7),2,'constant B'),
    (y**2*(y-1)**2,y,0,'two double roots'),
    (y**2*(y-1)**3,y**2,0,'mixed multiple roots'),
    ((y-2)**3,5+(y-2)**2,2,'translated value'),
]
matrices=[]
for A,B,expected,label in cases:
    P=S.expand(A*x+B)
    dx,dy=2,12
    columns=[bracket(P,x**i*y**j) for i in range(dx+1) for j in range(dy+1)]
    target=[S.expand(P**k) for k in range(3)]
    polys=[S.Poly(c,x,y) for c in columns+target]
    mons=sorted(set().union(*(set(p.monoms()) for p in polys)))
    M=S.polys.matrices.DomainMatrix.from_Matrix(S.Matrix([[p.coeff_monomial(mon) for p in polys[:len(columns)]] for mon in mons]))
    MN=S.polys.matrices.DomainMatrix.from_Matrix(S.Matrix([[p.coeff_monomial(mon) for p in polys] for mon in mons]))
    rank=M.rank(); joint=MN.rank(); actual=rank+3-joint
    require(actual==expected,f'image dimension {label}: {actual} != {expected}')
    matrices.append({'case':label,'A':str(A),'B':str(B),'Q_box':[dx,dy], 'image_intersection_dimension':actual,'matrix_rank':rank,'joint_rank':joint})

# Exact residue obstruction on several all-multiple-root examples.
# At least one nonzero residue is mandatory by the all-parameter lemma.
residues=[]
for A,roots in [(y**2*(y-1)**2,[0,1]),(y**3*(y-1)**3,[0,1]),(y**2*(y-1)**2*(y+1)**2,[-1,0,1])]:
    rr=[S.residue(1/A,y,a) for a in roots]
    require(any(r!=0 for r in rr),'multiple-root rational primitive hostile')
    residues.append({'A':str(A),'roots':roots,'residues':[str(r) for r in rr]})
# Degeneration A=0: independent critical-value radical controls.
boundary=[]
for B,C in [(y**2,T),(y**5,T),((y**2-1)**2,T*(T-1)),(y**3-3*y,(T-2)*(T+2)),(2*y+3,S.Integer(1))]:
    R=C.subs(T,B)
    Q=S.cancel(-x*R/S.diff(B,y))
    require(not S.denom(Q).has(x,y),'critical-value descent')
    require(zero(bracket(B,Q)-R),'critical-value bracket')
    resultant=S.Poly(S.resultant(B-T,S.diff(B,y),y),T)
    radical=resultant.sqf_part().monic().as_expr()
    require(zero(radical-C),'resultant radical')
    for factor,multiplicity in S.factor_list(C,T)[1]:
        candidate=S.cancel(-x*(C/factor).subs(T,B)/S.diff(B,y))
        require(S.denom(candidate).has(y),'critical-value missing-factor hostile')
    boundary.append({'B':str(B),'image_generator':str(C),'companion':str(Q)})
receipt={'zero_coefficient_boundary':boundary,'status':'FINITE-EXACT; all-parameter theorem proved separately','symbolic_identity_checks':identities,'direct_matrix_tests':matrices,'multiple_root_residue_controls':residues,'assert_gates':False}
print(json.dumps(receipt,sort_keys=True,indent=2))
