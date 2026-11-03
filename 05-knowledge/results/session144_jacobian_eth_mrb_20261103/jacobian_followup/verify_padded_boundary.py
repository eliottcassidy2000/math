#!/usr/bin/env python3
"""Exact higher-dimensional positive controls for the general boundary corollary.
The universal theorem is proved in prose. These are finite tests, not inference.
"""
import json
import sympy as S
v,t,r,s=S.symbols('v t r s'); l=2-3*v-t
K=t*(1+v)**2+v**2*(4+3*v); H1=l*(v+3*K); H2=l*l*(1+v)*K
rows=[]
for extra in [(1,),(2,),(3,),(5,),(1,1),(1,3),(2,5)]:
    variables=[v,t]+[r,s][:len(extra)]
    H=S.Matrix([H1,H2]+[l**w*q for w,q in zip(extra,variables[2:])])
    actual=S.factor(H.jacobian(variables).det())
    expected=2*l**(2+sum(extra))
    if S.expand(actual-expected)!=0:raise RuntimeError('padded quotient determinant failed: '+str(extra))
    rows.append({'weights':[1,-1,-2]+[-w for w in extra],'jacobian_exponent':2+sum(extra),'status':'PASS'})
print(json.dumps({'status':'PASS','finite_controls':rows,'source_map':'untwisted Alpoge map times identity','theorem_scope':'One positive weight equal to 1; remaining weights strictly negative; no claim for multiple positive or zero weights.'},sort_keys=True,indent=2))
