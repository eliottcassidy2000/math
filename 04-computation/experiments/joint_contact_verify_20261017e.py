#!/usr/bin/env python3
"""Independent exhaustive local table and characteristic-zero pair witnesses.

The finite ring is encoded as integers 0..840, not producer coordinate pairs.
It checks the selected obstruction, not the search over all 888 parameters.
"""
import json
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from contact_catalog_verify_20261017e import add as qadd, mul as qmul, norm as qnorm


def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def pack(a,b):return a%29+29*(b%29)
def unpack(x):return x%29,x//29
def add(x,y):
 a,b=unpack(x);c,d=unpack(y);return pack(a+c,b+d)
def mul(x,y):
 a,b=unpack(x);c,d=unpack(y);return pack(a*c+2*b*d,a*d+b*c)
def power(x,n):
 r=1
 while n:
  if n&1:r=mul(r,x)
  x=mul(x,x);n//=2
 return r
def conj(x):
 a,b=unpack(x);return pack(a,-b)
def norm(x):return mul(x,conj(x))
def reduce(v):
 roots=[1,pack(0,3),pack(0,10),2];out=0
 for c,r in zip(v,roots):
  c=Q(c);need(c.denominator%29!=0,'integral coefficient at29')
  out=add(out,mul(c.numerator*pow(c.denominator,-1,29)%29,r))
 return out


def main():
 base=Path('04-computation/experiments/marked_contact_current_20261017e')
 atlas=json.loads((base/'catalog.json').read_text());receipt=json.loads((base/'joint_place.json').read_text())
 one=(Q(1),Q(0),Q(0),Q(0));units=[tuple(Q(c) for c in u) for u in atlas['units']]
 need(power(2,14)==28,'2 is nonsquare modulo29')
 h,j,z=(pack(*receipt['base_place'][k]) for k in ('sqrt_minus3','sqrt_minus11','zeta15'))
 need(mul(h,h)==26 and mul(j,j)==18 and mul(h,j)==2,'base defining relations')
 need(power(z,15)==1 and power(z,3)!=1 and power(z,5)!=1,'primitive fifteenth root')
 need(add(1,mul(2,power(z,5)))==h,'sqrt(-3) matches the cyclotomic embedding')
 need(all(power(x,29)==conj(x) for x in (h,j,z)),'conjugation fixes the chosen place')
 U=[x for x in range(841) if norm(x)==1]
 need(len(U)==30 and sorted(unpack(x) for x in U)==[tuple(x) for x in receipt['complete_unit_group']],'entire norm-one group')
 masks=[];equation_checks=0
 for row in receipt['selected']:
  original=atlas['rows'][row['row']]
  for key in ('d','e'):
   need(row[key]==original[key] and row[key+'_coefficients']==atlas['differences'][row[key]],'marked source coefficients')
  d,e=reduce(row['d_coefficients']),reduce(row['e_coefficients']);allowed=set()
  for r in U:
   tests=[norm(add(add(d,mul(r,e)),t))==1 for t in U];equation_checks+=len(tests)
   if any(tests):allowed.add(unpack(r))
  need(allowed==set(map(tuple,row['allowed_rotations'])),'complete positive and negative local table')
  masks.append(allowed)
 need(not set.intersection(*masks),'no simultaneous local rotation')
 need(all(a&b for a,b in combinations(masks,2)),'positive local pair controls')
 pairs=set()
 for item in receipt['pairwise_characteristic_zero_witnesses']:
  pairs.add(tuple(item['rows']));ri=item['rotation'];r=units[ri]
  need(qnorm(r)==one,'shared exact unit rotation')
  for index,witness in zip(item['rows'],item['witnesses']):
   row=atlas['rows'][index];d,e=[tuple(Q(c) for c in atlas['differences'][row[k]]) for k in ('d','e')]
   wr,ti,wi=witness;t=units[ti];out=qadd(d,qadd(qmul(r,e),t))
   need(wr==ri and qnorm(t)==one and qnorm(out)==one,'pairwise characteristic-zero equations')
   if wi is not None:need(out==units[wi],'exact retained output')
 expected=set(combinations([r['row'] for r in receipt['selected']],2))
 need(pairs==expected,'all three pair witnesses')
 return dict(status='FINITE-EXACT independent complete local obstruction and all exact pair witnesses',
  selected=[r['row'] for r in receipt['selected']],norm_one_elements=len(U),local_equation_checks=equation_checks,
  allowed_sizes=[len(s) for s in masks],pairwise_intersection_sizes=[len(a&b) for a,b in combinations(masks,2)],
  common_size=0,characteristic_zero_pair_witnesses=len(pairs),
  scope='displayed contact equations; K-linear full graph retractions are ruled out, six-colour forcing remains open')

if __name__=='__main__':print(json.dumps(main(),sort_keys=True))
