#!/usr/bin/env python3
"""Independent occupancy invariance plus exact rational separating rows."""
from fractions import Fraction
from pathlib import Path
import itertools,json,copy
from verify_farkas import states,orders,D,row_step,rows_for
ROOT=Path(__file__).parent

def occupancy():
 occ=[[[0]*D for mask in range(16)] for o in range(6)]
 for j,(o,p,q,r,s) in enumerate(states):
  a,b,c=orders[o]
  mask=int(p!=0)|(int(q!=0)<<(a+1))|(int(r!=0)<<(b+1))|(int(s!=0)<<(c+1))
  occ[o][mask][j]=1
 return occ
OCC=occupancy()
RELATIONS=[[x-y for x,y in zip(OCC[o][mask],OCC[0][mask])] for o in range(1,6) for mask in range(16)]

def verify_quotient():
 checks=0
 for letter in range(8):
  N=[[0]*16 for _ in range(16)]
  for mask in range(16):
   N[mask|1][mask]+=1
   for a in range(3):
    if mask&(1<<(a+1)) or letter&(1<<a):N[mask|(1<<(a+1))][mask]+=1
  for o in range(6):
   for target in range(16):
    actual=row_step(OCC[o][target],letter)
    expected=[sum(N[target][source]*OCC[o][source][j] for source in range(16)) for j in range(D)]
    if actual!=expected:raise RuntimeError(('occupancy quotient',letter,o,target))
    checks+=1
 # At the empty spine every block has its sole mass in occupancy mask zero.
 initial=[int(p==q==r==s==0) for o,p,q,r,s in states]
 if any(sum(a*b for a,b in zip(r,initial))!=0 for r in RELATIONS):raise RuntimeError('initial occupancy relation')
 return checks

def verify(cert):
 rows,rhs=rows_for(cert);rows+=RELATIONS;rhs+=[0]*len(RELATIONS)
 y=[Fraction(x) for x in cert['multipliers']]
 if len(y)!=len(rows):return False
 if sum(a*b for a,b in zip(y,rhs))>=0:return False
 return all(sum(a*r[j] for a,r in zip(y,rows))>=0 for j in range(D))

def main():
 quotient_checks=verify_quotient()
 certs=json.loads((ROOT/'farkas_occupancy_certificates.json').read_text())
 expected={(core,a,b) for core in range(2) for a,b in itertools.product(range(8),repeat=2) if not(a&b) and (a|b).bit_count()==2}
 if {(c['core'],*c['word']) for c in certs}!=expected or len(certs)!=24:raise RuntimeError('wrong universe')
 for c in certs:
  if not verify(c):raise RuntimeError(('bad occupancy dual',c['core'],c['word']))
 mutations=[]
 c=copy.deepcopy(certs[0]);c['multipliers']=['0']*len(c['multipliers']);mutations.append(c)
 c=copy.deepcopy(certs[0]);c['multipliers']=[str(-Fraction(y)) for y in c['multipliers']];mutations.append(c)
 c=copy.deepcopy(certs[0]);c['word']=[0,0];mutations.append(c)
 c=copy.deepcopy(certs[0]);c['multipliers']=c['multipliers'][:-1];mutations.append(c)
 for c in mutations:
  if verify(c):raise RuntimeError('mutation accepted')
 print(json.dumps(dict(occupancy_quotient_row_identities=quotient_checks,occupancy_coordinate_equalities=quotient_checks*D,cross_order_relations=len(RELATIONS),exact_additional_exclusions=len(certs),certificate_coordinate_inequalities=len(certs)*D,remaining_growth_eligible_template='00',rejected_mutations=len(mutations),standard_library_only=True),sort_keys=True))
if __name__=='__main__':main()
