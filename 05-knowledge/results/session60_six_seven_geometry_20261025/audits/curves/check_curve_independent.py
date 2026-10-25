#!/usr/bin/env python3
"""Independent homogeneous Richelot and finite-field controls, no source import."""
from pathlib import Path
from itertools import combinations
from math import comb
import json
import sympy as S

def require(ok,label):
    if not ok: raise ValueError(label)

s,t,x,z=S.symbols('s t x z')
forms=[-s-6*t,s,s+t,s-t,3*s+5*t,-s-5*t,3*s+6*t]
source_json=json.loads((Path(__file__).resolve().parents[2] / 'recovered-bridges' / 'normal.json').read_text())
source_rows={row['omitted_label']:row for row in source_json['richelot']}

def det3(rows):
    a,b,c=rows
    return S.expand(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))

def polynomial(F):
    p=S.Poly(F.subs({s:1,t:x}),x)
    return [int(p.nth(i)) for i in range(7)]

def smooth(coeff,p):
    f=S.Poly(sum(c*x**i for i,c in enumerate(coeff)),x, modulus=p)
    return f.degree()>=5 and S.gcd(f,f.diff()).degree()==0

def legendre(a,p):
    a%=p
    return 0 if not a else (1 if pow(a,(p-1)//2,p)==1 else -1)

def count_fp(coeff,p):
    total=1+legendre(coeff[6],p)
    for u in range(p):
        v=0
        for c in reversed(coeff): v=(v*u+c)%p
        total+=1+legendre(v,p)
    return total

def count_fp2(coeff,p):
    d=next(i for i in range(2,p) if legendre(i,p)==-1)
    # Character is the F_p character of the norm in F_p[v]/(v²-d).
    total=1+int(coeff[6]%p!=0)
    for a in range(p):
        for b in range(p):
            h,k=0,0
            for c in reversed(coeff): h,k=(h*a+k*b*d+c)%p,(h*b+k*a)%p
            total+=1+legendre(h*h-d*k*k,p)
    return total

geometric=[]; finite=[]; finite2=[]; raw_twist_counterexamples=[]
for omitted in range(7):
    pairs=[((omitted+q)%7,(omitted+3*q)%7) for q in (1,2,4)]
    F=[S.expand(forms[i]*forms[j]) for i,j in pairs]
    # Homogeneous Poisson bracket divided by two, independent of affine differentiation.
    B=[S.expand((S.diff(F[(i+1)%3],t)*S.diff(F[(i+2)%3],s)-S.diff(F[(i+1)%3],s)*S.diff(F[(i+2)%3],t))/2) for i in range(3)]
    rows=[[f.coeff(s,2),f.coeff(s,1).coeff(t,1),f.coeff(t,2)] for f in F]
    delta=det3(rows)
    require(delta!=0,('determinant',omitted))
    require(int(delta)==source_rows[omitted]['determinant'],('source determinant pin',omitted))
    target=S.prod(B).expand(); src=S.prod(F).expand()
    # Binary smoothness: each quadratic is separable and pairwise disjoint.
    discriminants=[int(S.discriminant(b.subs(s,1),t)) for b in B]
    resultants=[int(S.resultant(B[i].subs(s,1),B[j].subs(s,1),t)) for i,j in combinations(range(3),2)]
    require(all(discriminants) and all(resultants),('target squarefreeness',omitted))
    require(target.subs({s:1,t:x})==S.sympify(source_rows[omitted]['geometric_target_polynomial'],locals={'u':x}),('source target pin',omitted))
    # Universal Richelot correspondence identity with ascending-coefficient determinant.
    identity=S.expand(sum(F[i].subs({s:1,t:x})*B[i].subs({s:1,t:z}) for i in range(3)))
    require(S.expand(identity+delta*(x-z)**2)==0,('signed correspondence identity',omitted))
    cf=polynomial(src); raw=polynomial(target); normalized=[int(delta)*c for c in raw]
    geometric.append(dict(omitted=omitted,delta=int(delta),bracket_discriminants=discriminants,pairwise_resultants=resultants,identity='sum F_i(x) B_i(z) = -Delta (x-z)^2'))
    for p in S.primerange(7,102):
        # Conservative full-seven branch-good condition, plus nonsingular isogeny.
        points=[None if int(f.coeff(t))%p==0 else (-int(f.coeff(s))*pow(int(f.coeff(t)),-1,p))%p for f in forms]
        if len(set(points))!=7 or int(delta)%p==0 or not smooth(cf,p) or not smooth(normalized,p): continue
        lhs,rhs=count_fp(cf,p),count_fp(normalized,p)
        require(lhs==rhs,('Q-normalized Fp trace',omitted,p,lhs,rhs))
        finite.append(dict(omitted=omitted,p=p,points=lhs))
        unnormalized=count_fp(raw,p)
        if unnormalized!=lhs: raw_twist_counterexamples.append(dict(omitted=omitted,p=p,source=lhs,raw_target=unnormalized))
        if p<=31:
            lhs2,rhs2=count_fp2(cf,p),count_fp2(normalized,p)
            require(lhs2==rhs2,('Q-normalized Fp2 trace',omitted,p,lhs2,rhs2))
            finite2.append(dict(omitted=omitted,p=p,points_over_p2=lhs2))
require(raw_twist_counterexamples,'raw target twist hostile must fail somewhere')
for n in (6,7):
    g=1+(n-4)*2**(n-3)
    require(g==comb(n,4)+2*comb(n,6),('character genus',n))
    # A transverse BASE line has class 2H, not H.
    require(2*g-2==2*(n-4)*2**(n-3),('adjunction base-line class',n))
out=dict(status='PASS',geometric_targets=geometric,normalized_Fp_cases=finite,normalized_Fp2_cases=finite2,raw_twist_counterexamples=raw_twist_counterexamples,scope='Independent audit; Q-normalized finite checks are supplementary to the C-only producer packet, not evidence that arbitrary rescalings preserve the Q quotient.')
Path(__file__).with_name('curve_independent_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(f'PASS: seven homogeneous squarefree targets; seven signed identities; {len(finite)} Fp and {len(finite2)} Fp2 normalized trace controls; {len(raw_twist_counterexamples)} raw-twist failures.')
