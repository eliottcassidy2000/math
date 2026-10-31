"""Exact standalone certificate verifier, standard library only."""
import itertools,json,copy
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
def data(word,rank):
 M=[[0]*36 for _ in range(36)]
 for j in range(36):
  v=[0]*36;v[j]=1
  for c in word:v=advance(v,c)
  for i in range(36):M[i][j]=v[i]
 c1=rowmul(output,M);c2=rowmul(c1,M)
 R=[c2[j]-28*c1[j]+27*output[j] for j in range(36)]
 rows=[]
 for _ in range(rank):rows.append(R);R=rowmul(R,M)
 rhs=[0]*rank
 for o in range(2):
  rows.append([int(oo==o and q==0 and r==0) for oo,p,q,r in states]);rhs.append(1)
 defect=[c1[j]-27*output[j] for j in range(36)]
 rows.extend([defect,[-x for x in defect],[-x for x in output]]);rhs.extend([26,26,-3])
 return rows,rhs

def verify(cert):
 y=[Fraction(x) for x in cert['multipliers']];rank=len(y)-5
 rows,rhs=data(cert['word'],rank)
 if any(x!=0 for x in y[-3:]):return False
 constant=sum(x*a for x,a in zip(y,rhs))
 if constant>=0:return False
 coefs=[sum(y[i]*rows[i][j] for i in range(len(rows))) for j in range(36)]
 return all(x>=0 for x in coefs)
certs=json.loads((ROOT/'radix27_farkas_certificates.json').read_text())
words={tuple(c['word']) for c in certs}
expected={w for w in itertools.product(range(4),repeat=3) if any(not(c&1) for c in w) and any(not(c&2) for c in w)}-{(0,0,0)}
if words!=expected or len(certs)!=48:raise RuntimeError('wrong certified template universe')
for c in certs:
 if not verify(c):raise RuntimeError(('invalid rational dual',c['word']))
# All-zero coefficients necessarily destroy the strict separating constant.
bad=copy.deepcopy(certs[0]);bad['multipliers']=['0']*len(bad['multipliers'])
if verify(bad):raise RuntimeError('failed zero mutation')
# Strong theorem uses no digit/H inequalities; any such nonzero multiplier is invalid.
bad=copy.deepcopy(certs[0]);bad['multipliers'][-1]='-1'
if verify(bad):raise RuntimeError('failed sign mutation')
print(json.dumps(dict(exact_template_exclusions=len(certs),coordinate_inequalities=36*len(certs),strict_negative_constants=len(certs),rejected_mutations=2,remaining_eligible_template=[0,0,0],standard_library_only=True,digit_and_H_bound_multipliers_all_zero=True)))
