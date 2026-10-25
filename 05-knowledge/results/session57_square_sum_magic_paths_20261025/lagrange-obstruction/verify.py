#!/usr/bin/env python3
"""Exact, exception-safe verifier: no checks depend on Python assert."""
import json
from pathlib import Path
import sympy as s
OUT=Path(__file__).resolve().parent
u,v=s.symbols('u v')
A=[9*(18*u*u+50*u*v+37*v*v)**2,36*(u*u+13*u*v+16*v*v)**2,36*(9*u*u+53*u*v+64*v*v)**2,4*(3*u*u-25*u*v-48*v*v)**2]
B=[s.Integer(0),20*(u-4*v)*(11*u+19*v)*(u+4*v)*(u+v),36*(3*u+11*v)*(7*u+12*v)*(3*u+4*v)*(u+v),756*(u-v)*(3*u+2*v)*(u+3*v)*(u+2*v)]
F=[16608*u**3+102220*u*u*v+205548*u*v*v+135135*v**3,5760*u**4+50208*u**3*v+168500*u*u*v*v+252948*u*v**3+141345*v**4,19680*u**3+113996*u*u*v+216300*u*v*v+135135*v**3]
def require(condition,message):
 if not condition: raise ValueError(message)
def balances(X): return [s.expand(X[0]+X[i]-sum(X[j] for j in range(1,4) if j!=i)) for i in range(1,4)]
for actual,expected in zip(balances(A),[-v*F[0],F[1],-v*F[2]]): require(s.expand(actual-expected)==0,'factor identity failed')
rows=[[t]+[int(f.subs({u:t,v:1})%47) for f in F] for t in range(47)]
infinity=[int(f.subs({u:1,v:0})%47) for f in F]
require(all(x for row in rows for x in row[1:]) and all(infinity),'projective zero found')
expected_A=[2916,36,2916,36];expected_B=[0,220,2268,2268]
for X,expected in [(A,expected_A),(B,expected_B)]:
 require(all(s.expand(x.subs(v,0)-n*u**4)==0 for x,n in zip(X,expected)),'boundary tuple mismatch')
boundary_balances=[int(x.subs({u:1,v:0})) for x in balances(B)]
require(boundary_balances==[-4316,-220,-220],'boundary B balances mismatch')
# Controls: a bad prime may have projective roots, hence one must check both charts.
controls={}
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43]:
 controls[str(p)]={'affine_roots':[[t for t in range(p) if f.subs({u:t,v:1})%p==0] for f in F],'infinity_roots':[i+1 for i,f in enumerate(F) if f.subs({u:1,v:0})%p==0]}
cert={'source':'https://arxiv.org/html/2508.07806v1#S1.E1','source_equation':'(1.1)','attribution':'Jean Lagrange family, reproduced by Ajai Choudhry (2025), with attribution through Guy, p. 267','prime':47,'columns':['t','F12(t,1)','F13(t,1)','F14(t,1)'],'residue_rows':rows,'infinity_residues':infinity,'A_at_v0_coefficients':expected_A,'B_at_v0_coefficients':expected_B,'B_boundary_balances':boundary_balances,'other_prime_controls':controls,'checks_passed':True}
(OUT/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
print(json.dumps({'checks_passed':True,'projective_evaluations':144,'infinity':infinity,'boundary_B_balances':boundary_balances}))
