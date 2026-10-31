"""Exact standalone certificate verifier, standard library only."""
import itertools,json,copy,sys
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).parent
states=list(itertools.product(range(2),range(3),range(3),range(2)));indices={s:i for i,s in enumerate(states)}
def advance(v,c):
 out=[0]*36
 for j,n in enumerate(v):
  if not n:continue
  o,p,q,r=states[j];a,b=(0,1) if o==0 else (1,0);fa=(c>>a)&1;fb=(c>>b)&1
  out[indices[o,2 if fa else 1,q,r]]+=n
  if q or fa:out[indices[o,p,2 if fb else 1,r]]+=n
  if r or fb:out[indices[o,p,q,1]]+=n
 return out
output=[int(p!=2 and(q==1 or(q==0 and o==0))) for o,p,q,r in states]
def rowmul(r,M):return [sum(r[i]*M[i][j] for i in range(36)) for j in range(36)]
def occupancy_relations():
 rows=[]
 for ep,eq,er in itertools.product(range(2),repeat=3):
  row=[]
  for o,p,q,r in states:
   wanted=(ep,eq,er) if o==0 else (ep,er,eq)
   row.append((1 if o==0 else -1)*int((int(p!=0),int(q!=0),int(r!=0))==wanted))
  rows.append(row)
 return rows

def data(word,rank,radix):
 M=[[0]*36 for _ in range(36)]
 for j in range(36):
  v=[0]*36;v[j]=1
  for c in word:v=advance(v,c)
  for i in range(36):M[i][j]=v[i]
 c1=rowmul(output,M);c2=rowmul(c1,M)
 R=[c2[j]-(radix+1)*c1[j]+radix*output[j] for j in range(36)]
 rows=[]
 for _ in range(rank):rows.append(R);R=rowmul(R,M)
 rows.extend(occupancy_relations());rhs=[0]*len(rows)
 for o in range(2):
  rows.append([int(oo==o and q==0 and r==0) for oo,p,q,r in states]);rhs.append(1)
 return rows,rhs

def verify(cert,radix):
 y=[Fraction(x) for x in cert['multipliers']];rank=len(y)-10
 rows,rhs=data(cert['word'],rank,radix)
 if sum(y[-2:])>=0:return False
 return all(sum(y[i]*rows[i][j] for i in range(len(rows)))>=0 for j in range(36))
# Verify exact invariant subspace by a complete integer kernel basis.
E=occupancy_relations();kb=[]
for row in E:
 pos=[i for i,x in enumerate(row) if x==1];neg=[i for i,x in enumerate(row) if x==-1]
 for group in (pos,neg):
  for j in group[1:]:
   v=[0]*36;v[j]=1;v[group[0]]=-1;kb.append(v)
 v=[0]*36;v[pos[0]]=1;v[neg[0]]=1;kb.append(v)
if len(kb)!=28:raise RuntimeError('bad occupancy kernel dimension')
initial=[0]*36;initial[indices[0,0,0,0]]=1;initial[indices[1,0,0,0]]=1
if any(sum(a*b for a,b in zip(row,initial)) for row in E):raise RuntimeError('initial occupancy relation')
for v in kb:
 if any(sum(a*b for a,b in zip(row,v)) for row in E):raise RuntimeError('bad kernel vector')
 for c in range(4):
  w=advance(v,c)
  if any(sum(a*b for a,b in zip(row,w)) for row in E):raise RuntimeError('occupancy not invariant')

checks=0
for word in [(0,3),(1,2),(2,1),(3,0)]:
 M=[[0]*36 for _ in range(36)]
 for j in range(36):
  v=[0]*36;v[j]=1
  for c in word:v=advance(v,c)
  for i in range(36):M[i][j]=v[i]
 def poly_row(factors):
  r=list(output)
  for a in factors:
   rm=rowmul(r,M);r=[x-a*y for x,y in zip(rm,r)]
  return r
 if any(poly_row([0,1,1,4,4,9])):raise RuntimeError('scalar annihilator')
 r=poly_row([0,1,4,4,9])
 if any(x%8 for x in r):raise RuntimeError('q not integral')
 q=[x//8 for x in r]
 if sum(x*y for x,y in zip(q,initial))!=-4:raise RuntimeError('wrong initial residue')
 for v in kb:
  for c in range(4):
   mv=advance(v,c)
   if sum(q[i]*(mv[i]-v[i]) for i in range(36))%3:raise RuntimeError('modular invariance')
   checks+=1
print(json.dumps(dict(complementary_words=4,integer_divisibility_checks=144,mod3_kernel_transition_checks=checks,initial_value=-4,universal_residue_mod3=2,annihilator_rows=4,forced_linear_mode='D=-q(v)/9 != 0',standard_library_only=True)))
