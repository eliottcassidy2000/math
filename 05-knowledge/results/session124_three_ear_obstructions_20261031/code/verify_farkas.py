#!/usr/bin/env python3
"""Standard-library exact verifier, independently rebuilding ordered-run transfer."""
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path
import itertools,json,copy
ROOT=Path(__file__).parent
orders=list(itertools.permutations(range(3)))
states=list(itertools.product(range(6),range(3),range(3),range(3),range(2)))
index={s:i for i,s in enumerate(states)}
D=len(states)

def row_step(row,letter):
 result=[0]*D
 for j,(o,p,q,r,s) in enumerate(states):
  a,b,c=orders[o];aa=(letter>>a)&1;bb=(letter>>b)&1;cc=(letter>>c)&1
  result[j]+=row[index[o,2 if aa else 1,q,r,s]]
  if q!=0 or aa:result[j]+=row[index[o,p,2 if bb else 1,r,s]]
  if r!=0 or bb:result[j]+=row[index[o,p,q,2 if cc else 1,s]]
  if s!=0 or cc:result[j]+=row[index[o,p,q,r,1]]
 return result

def output(core):
 arcs={(0,1),(1,2),(2,0)} if core else {(0,1),(0,2),(1,2)}
 result=[]
 for o,p,q,r,s in states:
  a,b,c=orders[o]
  result.append(int(p!=2 and (q==1 or (q==0 and (a,b) in arcs)) and (r==1 or (r==0 and (b,c) in arcs))))
 return result

def primitive(row):
 g=reduce(gcd,row,0)
 return [a//g for a in row] if g else row

def rows_for(cert):
 out=output(cert['core']);word=cert['word'];rank=cert['rank']
 def mul(row):
  for letter in reversed(word):row=row_step(row,letter)
  return row
 once=mul(out);twice=mul(once)
 R=primitive([x-17*y+16*z for x,y,z in zip(twice,once,out)])
 rows=[]
 for _ in range(rank):rows.append(R);R=primitive(mul(R))
 rows.extend([[int(o==oo and q==0 and r==0 and s==0) for oo,p,q,r,s in states] for o in range(6)])
 return rows,[0]*rank+[1]*6

def verify(cert):
 if cert['core'] not in (0,1) or len(cert['word'])!=2 or any(c not in range(8) for c in cert['word']):return False
 if cert['rank'] not in range(1,17):return False
 ys=[Fraction(x) for x in cert['multipliers']]
 rows,rhs=rows_for(cert)
 if len(ys)!=len(rows):return False
 if sum(y*b for y,b in zip(ys,rhs))>=0:return False
 return all(sum(y*r[j] for y,r in zip(ys,rows))>=0 for j in range(D))

def main():
 certs=json.loads((ROOT/'farkas_certificates.json').read_text())
 expected={(core,a,b) for core in range(2) for a,b in itertools.product(range(8),repeat=2) if not(a&b) and ((a|b).bit_count()%2)==1}
 actual={(c['core'],*c['word']) for c in certs}
 if actual!=expected or len(certs)!=28:raise RuntimeError('wrong certificate universe')
 for cert in certs:
  if not verify(cert):raise RuntimeError(('invalid exact dual',cert['core'],cert['word']))
 mutations=[]
 z=copy.deepcopy(certs[0]);z['multipliers']=['0']*len(z['multipliers']);mutations.append(z)
 z=copy.deepcopy(certs[0]);z['multipliers']=[str(-Fraction(y)) for y in z['multipliers']];mutations.append(z)
 z=copy.deepcopy(certs[0]);z['word']=[0,0];mutations.append(z)
 z=copy.deepcopy(certs[0]);z['multipliers']=z['multipliers'][:-1];mutations.append(z)
 for z in mutations:
  if verify(z):raise RuntimeError('mutation accepted')
 print(json.dumps(dict(exact_certificates=28,excluded_per_core=14,remaining_growth_eligible_per_core=13,required_distinct_templates=16,coordinate_inequalities=28*D,rejected_mutations=len(mutations),standard_library_only=True),sort_keys=True))
if __name__=='__main__':main()
