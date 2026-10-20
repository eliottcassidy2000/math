#!/usr/bin/env python3
"""Exact arithmetic and order-five controls for source-backed rigidity/cutoff proof."""
import itertools as it
import math
import json
from pathlib import Path
OUT=Path(__file__).resolve().parent
# All circulant step sets in order five, one representative of each opposite pair.
sets=[frozenset(s) for s in it.product((1,4),(2,3))]
base=frozenset((1,2))
multiplier_orbit={frozenset((a*x)%5 for x in base) for a in range(1,5)}
assert set(sets)==multiplier_orbit
P=(0,2,3,4,1)
assert all((y-x)%5 in base for x,y in zip(P,P[1:]))
assert (P[-1]-P[0])%5 in base
assert (P[0]-P[-1])%5 not in base
# Independent direct HP, Aut and bypass control at the regular order-five class.
A=[[int(i!=j and (j-i)%5 in base) for j in range(5)] for i in range(5)]
H=Aut=nonclosing=HC=0
for p in it.permutations(range(5)):
 Aut+=all(A[i][j]==A[p[i]][p[j]] for i in range(5) for j in range(5))
 if all(A[x][y] for x,y in zip(p,p[1:])):
  H+=1;nonclosing+=A[p[0]][p[-1]]
  HC+=int(p[0]==0 and A[p[-1]][p[0]])
assert (H,Aut,nonclosing,HC)==(15,5,5,2)
# Explicit cutoff: bases cover the three odd residue classes modulo six.
bases=[]
for n in [133,135,137]:
 k=(n-7)//6
 L=2*4**k*math.factorial(k);U=3**((n-1)//2)
 assert L>U
 assert 4*(k+1)>27
 bases.append({'n':n,'k':k,'double_hamilton_cycle_lower_bound':L,'automorphism_upper_bound':U,'difference':L-U})
# n=131 is not excluded by this pair of numerical bounds.
n=131;k=(n-7)//6
assert 2*4**k*math.factorial(k)<=3**((n-1)//2)
rows={'n5_circulant_step_sets':list(map(sorted,sets)),'n5_multiplier_orbit_of_1_2':list(map(sorted,sorted(multiplier_orbit,key=lambda x:sorted(x)))),'n5_nonclosing_HP':P,'n5_HP':H,'n5_Aut':Aut,'n5_nonclosing_HP_count':nonclosing,'n5_unrooted_directed_HC':HC,'cutoff_base_checks':bases,'propagation':'For n -> n+6, k -> k+1 and the lower/upper ratio is multiplied by 4(k+1)/27 > 1.','largest_not_numerically_excluded_order':131}
(OUT/'independent_rigidity_cutoff.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows),flush=True)
