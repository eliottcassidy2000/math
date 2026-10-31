"""Independent occupancy/modular/Farkas replay using audited tensor-sum matrices."""
from pathlib import Path
import runpy,itertools,json,hashlib
from fractions import Fraction
D=runpy.run_path(str(Path(__file__).with_name('independent_audit.py')))
s=D['s']; matrices=D['matrices']; states=D['states'];c=D['output'];v0=s.Matrix(D['initial']);root=D['ROOT'];wm=D['word_matrix']
E=s.Matrix([[ (1 if o==0 else -1)*int((int(p!=0),int(q!=0),int(r!=0))==((ep,eq,er) if o==0 else (ep,er,eq))) for o,p,q,r in states] for ep,eq,er in itertools.product(range(2),repeat=3)])
if E.rank()!=8 or E*v0!=s.zeros(8,1):raise RuntimeError('occupancy initial/rank')
# Each row has disjoint support and one unit pivot; nullspace denominators1 imply
# free-variable integer basis is a complete Z-kernel basis, not just Q-spanning.
K=s.Matrix.hstack(*E.nullspace())
if K.shape!=(36,28) or any(x.q!=1 for x in K):raise RuntimeError('integral kernel')
for M in matrices:
 if E*M*K!=s.zeros(8,28):raise RuntimeError('occupancy preservation')
reachable=v0
while True:
 A=reachable.row_join(s.Matrix.hstack(*(M*reachable for M in matrices)))
 ids=A.rref()[1];N=A[:,list(ids)]
 if N.cols==reachable.cols:break
 reachable=N
if reachable.cols!=28:raise RuntimeError('reachable rank')
mod_checks=0
for w in ((0,3),(1,2),(2,1),(3,0)):
 M=wm(w);I=s.eye(36)
 ann=c*M*(M-I)**2*(M-4*I)**2*(M-9*I)
 if ann!=s.zeros(1,36):raise RuntimeError('annihilator')
 q=c*M*(M-I)*(M-4*I)**2*(M-9*I)/8
 if any(x.q!=1 for x in q) or (q*v0)[0]!=-4:raise RuntimeError('integer q')
 for A in matrices:
  if any(x%3 for x in q*(A-I)*K):raise RuntimeError('mod3 invariance')
  mod_checks+=28
 # On a formal scalar sequence D*k, p(shift)/8 has value -9*D.
 z=s.Symbol('z');p=z*(z-1)*(z-4)**2*(z-9)/8
 if s.diff(p,z).subs(z,1)!=-9:raise RuntimeError('mode normalization')
# Integer row multiplication instead of producer Krylov rows or floating LP.
def multiply(r,M):return [sum(r[i]*M[i][j] for i in range(36)) for j in range(36)]
units=[[int(oo==o and q==r==0) for oo,p,q,r in states] for o in range(2)]
Erows=[[int(x) for x in E.row(i)] for i in range(8)]
counts=[];digests={}
for name,length,kind in [('radix27_farkas_certificates.json',3,'old'),('periodic_length4_certificates.json',4,'plain'),('periodic_reachable_length2_certificates.json',2,'occupancy'),('periodic_reachable_length4_certificates.json',4,'occupancy')]:
 raw=(root/name).read_bytes();digests[name]=hashlib.sha256(raw).hexdigest();record=json.loads(raw)
 certs=record if kind=='old' else record['certificates'];B=3**length
 universe={w for w in itertools.product(range(4),repeat=length) if any(not(x&1) for x in w) and any(not(x&2) for x in w)}
 certified={tuple(x['word']) for x in certs}
 expected_remaining={(0,)*length}
 if kind=='plain':expected_remaining|={(0,3,0,3),(1,2,1,2),(2,1,2,1),(3,0,3,0)}
 if certified!=universe-expected_remaining or len(certs)!=len(certified):raise RuntimeError(('universe',name))
 for cert in certs:
  M=[[int(x) for x in wm(cert['word']).row(i)] for i in range(36)]
  y=[Fraction(x) for x in cert['multipliers']]
  if kind=='old':
   if any(y[-3:]):raise RuntimeError('legacy bounds used')
   y=y[:-3]
  offset=10 if kind=='occupancy' else 2
  rank=len(y)-offset
  row=list(map(int,c));c1=multiply(row,M);c2=multiply(c1,M)
  r=[c2[j]-(B+1)*c1[j]+B*row[j] for j in range(36)]
  rows=[]
  for j in range(rank):rows.append(r);r=multiply(r,M)
  if kind=='occupancy':rows.extend(Erows)
  rows.extend(units)
  if y[-1]+y[-2]>=0:raise RuntimeError('nonseparating constant')
  if any(sum(a*r[j] for a,r in zip(y,rows))<0 for j in range(36)):raise RuntimeError(('negative coordinate',cert['word']))
 counts.append(dict(file=name,exclusions=len(certs),remaining=sorted(expected_remaining),coordinate_inequalities=len(certs)*36))
print(json.dumps(dict(occupancy_rank=8,reachable_rank=28,integer_kernel_rank=28,mod3_kernel_transition_checks=mod_checks,forced_linear_mode='D=-q/9; q=2 mod3',certificate_families=counts,input_sha256=digests),sort_keys=True))
