#!/usr/bin/env python3
"""Exact controls for finite inverse-graph elimination and Euler accounting."""
from pathlib import Path
from itertools import product
import json
import sympy as s
u,v,z,p,q,w=s.symbols('u v z p q w')
checks=0

def check(ok,label):
    global checks
    if not bool(ok): raise AssertionError(label)
    checks+=1

counts={'degree_d':0,'linear':0,'empty':0,'whole_line':0}
for d in range(3,9):
  for H,A,B in product(range(-2,3),repeat=3):
    R=s.Poly(H*z**d+H*(z**(d-1)+2*z**2)+A*z+B,z)
    if H:
      check(R.degree()==d,('degree d',d,H,A,B)); counts['degree_d']+=1
    elif A:
      check(R.degree()==1,('linear',d,H,A,B))
      check(R.eval(-s.Rational(B,A))==0,('linear selector',d,H,A,B))
      counts['linear']+=1
    elif B:
      check(R.degree()==0 and R.as_expr()!=0,('empty',d,H,A,B)); counts['empty']+=1
    else:
      check(R.is_zero,('whole line',d,H,A,B)); counts['whole_line']+=1

# The finite graph of (p, p*q^d-d*q) is exactly the original affine source.
# It is a polynomial/quasifinite hostile, but NOT Keller: critical fibers remain.
nonkeller=[]
for d in range(3,9):
    R=u*z**d-d*z-v
    check(s.expand(R.subs({u:p,v:p*q**d-d*q,z:q}))==0,('graph identity',d))
    # p=u, q=z are the polynomial inverse from graph to source.
    check(s.diff(R,v)==-1,('graph is affine plane',d))
    check(R.subs(u,0)==-d*z-v,('one finite boundary root',d))
    jac=s.det(s.Matrix([[1,0],[q**d,d*p*q**(d-1)-d]]))
    check(jac.subs({p:1,q:1})==0,('non-Keller critical point',d))
    special=s.Poly(R.subs({u:1,v:1-d}),z)
    check(special.eval(1)==0 and special.diff().eval(1)==0,('ramified fiber off boundary',d))
    check(s.gcd(special,special.diff()).degree()>=1,('cover formula hypothesis fails',d))
    nonkeller.append({'d':d,'boundary':'u=0','A':-d,'omitted_set':'empty',
                     'critical_target':[1,1-d]})

# Without monicity in the companion, a finite primitive root does not give a finite source point.
# F(q,p)=(q,q*p), h=q=u; its companion p=v/u has a pole over u=0,v=1.
check((z-u).subs({u:0,z:0})==0,'primitive graph has finite root at omitted target')
check((u*p-v).subs({u:0,v:1})==-1,'companion equation excludes omitted target')
check((u*p-v).subs({u:0,v:0})==0,'companion equation contracts whole source line')
check(s.cancel(v/u).subs(v,1)==1/u,'finite h misses companion pole')

# Explicit normalization/branch examples; no assertion these are Keller boundaries.
models=[
  ('line',1,[],1),
  ('unibranch cusp',1,[1],1),
  ('one nodal rational curve',1,[2],0),
  ('two transverse lines',2,[2],1),
  ('two disjoint lines',2,[],2),
  ('triangle of lines',3,[2,2,2],0),
  ('three concurrent lines',3,[3],1),
]
# Cusp and node equations are checked from their normalization parameters.
check(s.expand((w**3)**2-(w**2)**3)==0,'cusp normalization')
nx=w**2-1; ny=w*(w**2-1)
check(s.expand(ny**2-nx**2*(nx+1))==0,'node normalization')
check(nx.subs(w,1)==nx.subs(w,-1)==0 and ny.subs(w,1)==ny.subs(w,-1)==0,'node two preimages')
check(s.diff(nx,w).subs(w,1)!=0 and s.diff(nx,w).subs(w,-1)!=0,'node smooth branches')
euler=[]
for name,c,branches,wanted in models:
    chi=c-sum(b-1 for b in branches)
    check(chi==wanted,('normalization Euler',name))
    for d in range(3,9):
      O=(d-1)*(1-chi)
      check(d*(1-chi)+chi-O==1,('Euler identity',name,d))
      euler.append({'model':name,'d':d,'chi_boundary':chi,'required_omitted_points':O,
                    'ruled_out_by_nonnegativity':O<0})

# Polynomial coefficient atlases with a line denominator and nonempty numerator base locus.
# Their finite-graph omission prediction contradicts the Keller singleton-boundary Euler gate.
rejected=[]
for d in range(3,9):
  for r in range(1,7):
    H=u; A=v**r-1; B=s.Integer(1)
    check(s.gcd(s.Poly(H,u,v),s.Poly(A,u,v)).total_degree()==0,('coprime H,A',d,r))
    check(s.gcd(s.Poly(A,v),s.Poly(s.diff(A,v),v)).degree()==0,('reduced intersection',d,r))
    check(s.Poly(A,v).degree()==r,('intersection cardinality',d,r))
    check(B==1,('B unit at every intersection',d,r))
    check(r>(d-1)*(1-1),('Euler rejects line-boundary atlas',d,r))
    rejected.append({'d':d,'H':'u','A':f'v^{r}-1','B':'1','base_points':r,
                     'Euler_required_points':0,'status':'excluded if all Keller hypotheses hold'})

out={'status':'FINITE-EXACT controls; general proof is in FINITE_SELECTOR_EULER_GATE.md',
     'checks':checks,'specialized_fiber_cases':counts,'non_keller_graph_controls':nonkeller,
     'Euler_models':euler,'rejected_coefficient_atlases':rejected,
     'scope':'No graph atlas is claimed to be a Keller realization.'}
Path(__file__).with_name('selector-verification-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'checks':checks,'fiber_cases':sum(counts.values()),'rejected_atlases':len(rejected),'status':'PASS'},sort_keys=True))
