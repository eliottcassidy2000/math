#!/usr/bin/env python3
"""Finite arithmetic controls, not a substitute for the equality proof."""
from itertools import product
import json
from pathlib import Path
pairs=[(i,j) for i in range(-2,3) for j in range(-2,3) if -2<=i+j<=2]
solutions=[]
# Explicit finite universe of integral diagonal exponents.
for a,b,c,d in product(range(-6,7),repeat=4):
 w={-2:a,-1:b,0:0,1:c,2:d}
 if a+d or b+c or d!=2*c: continue
 assert all(w[i]+w[j]-w[i+j]==0 for i,j in pairs)
 solutions.append(w)
assert len(solutions)==7
assert all(w[j]==j*w[1] for w in solutions for j in range(-2,3))
classical=[]
for n in range(1,249):
 for family,dim in [('A',n*(n+2)),('BC',n*(2*n+1)),('D',n*(2*n-1))]:
  if dim==248: classical.append((family,n))
assert not classical
result={'diagonal_exponent_range':[-6,6],'fixed_metric_and_positive_bracket_solutions':solutions,'all_solutions_are_grading_automorphisms':True,'classical_simple_rank_solutions_dim248':classical,'exceptional_dimensions':{'G2':14,'F4':52,'E6':78,'E7':133,'E8':248},'scope':'Finite diagonal probe and dimension arithmetic only; equality rigidity proved separately.'}
Path(__file__).with_name('equality_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
