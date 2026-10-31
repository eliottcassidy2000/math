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
length=int(sys.argv[1])
record=json.loads((ROOT/f'periodic_reachable_length{length}_certificates.json').read_text())
if record['length']!=length or record['radix']!=3**length:raise RuntimeError('wrong parameters')
certs=record['certificates'];words={tuple(c['word']) for c in certs}
expected={w for w in itertools.product(range(4),repeat=length) if any(not(c&1) for c in w) and any(not(c&2) for c in w)}
remaining={tuple(x[0]) for x in record['uncertified']}
if words&remaining or words|remaining!=expected or len(words)!=len(certs):raise RuntimeError('wrong word universe')
for c in certs:
 if not verify(c,3**length):raise RuntimeError(('invalid rational certificate',c['word']))
bad=copy.deepcopy(certs[0]);bad['multipliers']=['0']*len(bad['multipliers'])
if verify(bad,3**length):raise RuntimeError('failed zero mutation')
bad=copy.deepcopy(certs[0]);bad['multipliers'][-2:]=['1','1']
if verify(bad,3**length):raise RuntimeError('failed separating-sign mutation')
print(json.dumps(dict(length=length,radix=3**length,exact_certificates=len(certs),coordinate_inequalities=36*len(certs),remaining=sorted(remaining),rejected_mutations=2,standard_library_only=True,occupancy_kernel_transition_checks=112)))
