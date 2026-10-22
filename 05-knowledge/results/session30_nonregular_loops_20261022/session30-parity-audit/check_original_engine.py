"""Compare independently contracted parity JSON to the original optimized engine."""
from pathlib import Path
import sys,json,numpy as np,time
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE/'session30-decorated-engine'))
from decorated import evaluate
cores=json.loads((BASE/'session30-parity-basis/parity-kernels.json').read_text())
def parity(M,z):
 Q=z*M;Q2=Q@Q;I=np.eye(len(M));im=np.linalg.inv(I-Q2);ip=np.linalg.inv(I+Q2)
 B=[-Q2@ip,Q2@im,Q@im,Q@ip]
 out=[]
 for core in cores:
  spec=','.join('abcd'[i]+'abcd'[j] for i,j in core['edges'])+'->'
  result=sum(term['numerator']*np.einsum(spec,*[B[k] for k in term['indices']],optimize=True) for term in core['terms'])/core['denominator']
  out.append(result)
 return out
fixtures=[]
for n in (4,6):
 # Deletion of one vertex from a regular cyclic tournament gives an almost-regular fixture.
 N=n+1;S=np.array([[0 if i==j else 1 if (j-i)%N in range(1,(N+1)//2) else -1 for j in range(n)] for i in range(n)])
 fixtures.append((n,S/n))
results=[]
for n,M in fixtures:
 for z in (1.,.4+.2j):
  mine=parity(M,z);reference=[v for w in (1,2) for d,v in evaluate(M,w,z,detail=True,reduced=False)]
  err=max(abs(a-b) for a,b in zip(mine,reference))
  assert err<2e-11
  results.append({'n':n,'z':str(z),'cores':len(mine),'max_absolute_difference':float(err)})
print(json.dumps({'status':'PASS','mode':'original five-basis optimized engine, reduced=False','results':results},indent=2))
