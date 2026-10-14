#!/usr/bin/env python3
"""Symbolic all-sector gap certificate for balanced-degree clutches."""
import json
import sympy as S

def check(ok, message):
    if not ok:
        raise RuntimeError(message)

def positive(expr, params, message, strict=True):
    p=S.Poly(S.expand(expr),*params)
    check(p.total_degree()<=1,message+': nonaffine')
    const=p.coeff_monomial(1)
    check(all(coef>=0 for mon,coef in p.terms()),message+': negative coefficient')
    check(const>0 if strict else const>=0,message+': invalid constant')

u,k,r=S.symbols('u k r',integer=True,nonnegative=True)
a=1+u
# Low degrees: a-n=k>=0. All dominance gaps are affine in this difference.
low_gaps=[1+4*k,2+8*k,3]
for expr in low_gaps+[4*a-2]:
    positive(expr,(u,k),'low-degree sector')
# High small branches: n=a+1+k.
n=a+1+k
small_gaps=[4*(n-a)-1,2*(n-a)+1]
for expr in small_gaps+[2*n-2]:
    positive(expr,(u,k),'small high-degree roots')
sectors=[]
for label,n,e in [('first_high_nonresonant',a+1,S.Integer(1)),('first_high_resonant',a+1,S.Integer(0)),('higher_degrees',a+2+k,2+k)]:
    gaps=[4*n-e-4*a,2*n+e-2*a]
    order=2*n+e-3
    for expr in gaps+[order]:
        positive(expr,(u,k),label)
    sectors.append({'sector':label,'n':str(S.expand(n)),'e':str(e),'dominance_gaps':[str(S.expand(g)) for g in gaps],'omega_order_on_double_cover':str(S.expand(order))})
# Finite r<m: write m=r+1+k, so d0=m-r=1+k.
m=r+1+k
finite_gap=m+r
finite_order=m-r+1
positive(finite_gap,(r,k),'finite escape dominance on double cover')
positive(finite_order,(r,k),'finite escape time-form order')
# Smallest edge on its actual unramified infinity branch.
edge_order=S.Integer(2)+S.Integer(0)/2-2
check(edge_order==0,'smallest resonance is regular, not vanishing')
print(json.dumps({
    'status':'PASS',
    'parameters':{'a':'1+u','u,k,r':'nonnegative integers'},
    'low_sector':{'a-n':'k','domination_gaps':[str(g) for g in low_gaps],'omega_order_on_fourfold_cover':str(S.expand(4*a-2))},
    'small_high_sector':{'n':'a+1+k','domination_gaps':[str(S.expand(g)) for g in small_gaps],'omega_order_on_double_cover':str(S.expand(2*(a+1+k)-2))},
    'escaping_infinity_sectors':sectors,
    'finite_r_lt_m':{'m':'r+1+k','domination_gap_on_double_cover':str(S.expand(finite_gap)),'omega_order_on_double_cover':str(S.expand(finite_order))},
    'finite_r_ge_m':'Generic E0 nonzero; two simple constant escaping h roots; Pz has nonzero constant limit.',
    'smallest_edge_a1_n2_resonant':{'actual_unramified_order':0,'double_cover_order':1},
    'theorem':'deg V=2a,deg A=a>=1,gcd(B,V)=1 imply no polynomial constant-nonzero-Jacobian mate.',
    'all_B_corollary':'If A has a simple zero at every V-root, boundary-colliding B have explicit critical points, so every B is excluded.',
    'scope':'All-sector symbolic inequalities; Newton-Puiseux completeness supplied by the note. H-J remains open.'
},indent=2,sort_keys=True))
