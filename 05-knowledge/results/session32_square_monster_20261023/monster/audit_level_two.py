#!/usr/bin/env python3
"""Exact F_p polynomial audit of the Monster/Baby prime split.

Cited inputs: supersingular Legendre Hasse polynomial; X_0(2) coordinate
u=(eta(tau)/eta(2*tau))^24 has j=(u+256)^3/u^2 and w_2:u->4096/u.
We obtain supersingular j-values by a resultant, pull back to u, then quotient
by w_2 using s=u+4096/u. Every polynomial is squarefreed: these are coarse
moduli POINT counts, not ramification or stack-weighted counts.

Run: python3 audit_level_two.py > level_two_exact.json
Dependencies: Python 3, SymPy. No network or input files are used.
"""
from math import comb, prod
import json
import sympy as sp
l,j,u,s=sp.symbols('l j u s')
MONSTER={2:46,3:20,5:9,7:6,11:2,13:3,17:1,19:1,23:1,29:1,31:1,41:1,47:1,59:1,71:1}
BABY={2:41,3:13,5:6,7:2,11:1,13:1,17:1,19:1,23:1,31:1,47:1}

def squarefree(poly,x,p):
    return sp.Poly(poly,x,modulus=p).sqf_part().monic()

def degree_profile(poly):
    return sorted((f.degree(),e) for f,e in sp.factor_list(poly)[1])

def character(a,p):
    z=pow(int(a)%p,(p-1)//2,p)
    return -1 if z==p-1 else z

def ss_j_from_counts(p):
    out=[]
    for a_j in range(p):
        if a_j==0: a,b=0,1
        elif a_j==1728%p: a,b=1,0
        else:
            k=a_j*pow((1728-a_j)%p,-1,p)%p
            a,b=3*k%p,2*k%p
        trace=-sum(character(x**3+a*x+b,p) for x in range(p))
        if trace==0: out.append(a_j)
    return out

def row(p):
    m=(p-1)//2
    H=sum((comb(m,k)**2%p)*l**k for k in range(m+1))
    jrel=256*(l*l-l+1)**3-j*l*l*(l-1)**2
    # Compute an integer resultant using a lift of H, then reduce modulo p.
    hj=sp.resultant(sp.Poly(H,l,modulus=p).as_expr(),jrel,l)
    J=squarefree(hj,j,p)
    expected=p//12+{1:0,5:1,7:1,11:2}[p%12]
    assert J.degree()==expected
    assert sorted(a for a in range(p) if J.eval(a)%p==0)==ss_j_from_counts(p)
    # Form numerator of J((u+256)^3/u^2) exactly over F_p.
    A=sp.Poly((u+256)**3,u,modulus=p)
    B=sp.Poly(u*u,u,modulus=p)
    d=J.degree()
    P=sp.Poly(0,u,modulus=p)
    for (e,),coef in J.terms(): P += int(coef)*(A**e)*(B**(d-e))
    U=squarefree(P.as_expr(),u,p)
    # Its root set must be w_2 invariant.
    rev=sum(int(c)*pow(4096,e,p)*u**(U.degree()-e) for (e,),c in U.terms())
    assert squarefree(rev,u,p)==U
    us=sp.resultant(U.as_expr(),u*u-s*u+4096,u)
    S=squarefree(us,s,p)
    jp=degree_profile(J);up=degree_profile(U);sp_=degree_profile(S)
    assert all(e==1 for d,e in jp+up+sp_)
    assert max(d for d,e in jp+up+sp_)<=2
    nonrational=sum(d for d,e in sp_ if d>1)
    assert nonrational%2==0
    return {'p':p,'ss_j_count':J.degree(),'j_degrees':jp,'ss_X0_2_count':U.degree(),
       'X0_2_degrees':up,'ss_X0plus_2_count':S.degree(),'X0plus_2_degrees':sp_,
       'nonrational_X0plus_2_points':nonrational,
       'genus_X0star_2p_by_Aricheta':nonrational//2,
       'X0plus_2_polynomial_coefficients_high_to_low':[int(c)%p for c in S.all_coeffs()],
       'X0plus_2_factorization':[(str(f.as_expr()),e) for f,e in sp.factor_list(S)[1]]}

if __name__=='__main__':
    M=prod(p**e for p,e in MONSTER.items());B=prod(p**e for p,e in BABY.items())
    out={'method':'Legendre Hasse resultant; X0(2) j-map; Fricke invariant s=u+4096/u; squarefree mod p',
       'source_urls':['https://arxiv.org/html/1809.07421','https://brauer.maths.qmul.ac.uk/Atlas/spor/M/','https://brauer.maths.qmul.ac.uk/Atlas/spor/B/'],
       'orders':{'Monster':M,'Baby':B,'M_over_B':M//B,'M_over_2B':M//(2*B)},
       'valuation_strata':{'equal_positive':[p for p in MONSTER if MONSTER[p]==BABY.get(p,0)],
         'Monster_only':[p for p in MONSTER if p not in BABY],
         'increased_common':[p for p in BABY if MONSTER[p]>BABY[p]]},
       'rows':[row(p) for p in sorted(set(MONSTER)-{2,3}|{37,43})]}
    assert out['valuation_strata']=={'equal_positive':[17,19,23,31,47],'Monster_only':[29,41,59,71],'increased_common':[2,3,5,7,11,13]}
    for r in out['rows']:
        if r['p'] in BABY: assert r['nonrational_X0plus_2_points']==0
        if r['p'] in {29,41,59,71}: assert r['nonrational_X0plus_2_points']==2
    print(json.dumps(out,indent=2))
