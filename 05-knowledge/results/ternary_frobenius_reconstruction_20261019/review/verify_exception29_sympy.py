#!/usr/bin/env python3
"""Independent exact SymPy reconstruction from the 23-row exception note."""
import argparse,json,re
from pathlib import Path
from functools import reduce
from math import gcd
from sympy import symbols,Matrix,Poly,expand
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--base',type=Path);args=parser.parse_args()
HERE=Path(__file__).resolve().parent
candidates=[args.base]if args.base is not None else[HERE.parent,HERE.parent/'session6-periods']
BASE=next((q for q in candidates if(q/'EXCEPTIONAL_PRIME_29.md').is_file()),None)
if BASE is None:raise FileNotFoundError('Pass --base with the exception note directory')
rows=[tuple(map(int,m.groups()))for m in re.finditer(r'^\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(-?\d+)\s*\|', (BASE/'EXCEPTIONAL_PRIME_29.md').read_text(),re.M)]
x,y,z=symbols('x y z');v=(x,y,z)
if len(rows)!=23:raise RuntimeError('table length')
E=expand(sum(k*Matrix([[t**a,t**b,t**c]for t in v]).det()for a,b,c,k in rows))
P=Poly(E,*v);A=E.xreplace({x:-y,y:x-y});A2=E.xreplace({x:y-x,y:-x})
N=Poly(expand(E+A+A2),*v);g=reduce(gcd,(abs(int(c))for c in N.coeffs()),0)
if g!=29:raise RuntimeError(('integer norm gcd',g))
if len(P.terms())!=138 or P.total_degree()!=42:raise RuntimeError('support size/degree')
if P.coeff_monomial(x**2*y**12*z**28)!=1:raise RuntimeError('nonzero coefficient')
if expand(E.xreplace({x:y,y:z,z:x})-E)!=0:raise RuntimeError('cyclicity')
if expand(E.xreplace({x:y,y:x})+E)!=0:raise RuntimeError('alternation')
if any(min(t)<2 or any(a%2 for a in t)for t,c in P.terms()):raise RuntimeError('interior/even')
V=x*x*y*y*z*z*(x*x-y*y)*(y*y-z*z)*(z*z-x*x)
Q=P.exquo(Poly(V,*v));Qexpr=Q.as_expr()
if any(not c.is_Integer for c in Q.coeffs()):raise RuntimeError('integral quotient')
if Q.total_degree()!=30 or any(any(a%2 for a in t)for t,c in Q.terms()):raise RuntimeError('quotient degree')
if expand(Qexpr.xreplace({x:y,y:x})-Qexpr)!=0 or expand(Qexpr.xreplace({x:y,y:z,z:x})-Qexpr)!=0:raise RuntimeError('quotient symmetry')
f={tuple(map(int,k.split(','))):v%29 for row in json.loads((BASE/'next_thresholds.json').read_text())if row['p']==29 for k,v in row['witness'].items()}
if {t:int(c)%29 for t,c in P.terms()}!=f:raise RuntimeError('saved witness')
# Independent cyclic interior lift obstruction, retaining all monomial rows.
import numpy as np
from verify_period_ranks import integer_columns,monomials,rank_forward
cols=integer_columns(42,'interior');mons=monomials(42);ix={t:i for i,t in enumerate(mons)}
M=np.zeros((len(mons),len(cols)+1),dtype=np.int64)
for j,col in enumerate(cols):
    for t,c in col.items():M[ix[t],j]=c%29
for t,c in N.terms():M[ix[t],-1]=(-int(c)//29)%29
r=rank_forward(M[:,:-1],29);aug=rank_forward(M,29)
if (len(mons),len(cols),r,aug)!=(946,274,273,274):raise RuntimeError(('interior lift obstruction',M.shape,r,aug))
print(json.dumps({'status':'PASS','alternants':len(rows),'E_terms':len(P.terms()),'norm_terms':len(N.terms()),'norm_coefficient_gcd':g,'quotient_terms':len(Q.terms()),'quotient_degree':Q.total_degree(),'interior_norm_shape':[len(mons),len(cols)],'rank_mod29':r,'augmented_rank_mod29':aug},indent=2))
