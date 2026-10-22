import pathlib,sys,json,time,statistics,functools
import numpy as np
import sympy as sp
ROOT=pathlib.Path(__file__).resolve().parent;ENGINE=ROOT.parent/'session30-decorated-engine'
sys.path.insert(0,str(ENGINE))
from decorated import evaluate as reference,compiled_reduced
x=sp.symbols('x');G=-x*x/(1+x*x);H=x*x/(1-x*x);O=x/(1-x*x);P=x/(1+x*x)
rp=x/(1-x);rm=-x/(1+x);K=(1-x*x)/(1+x*x)
identities=[rp-H-O,rm-H+O,-rp*rp*K+H-P+O,-rm*rm*K+H+P-O,-rp*rm*K+G]
if any(sp.cancel(v)!=0 for v in identities):raise RuntimeError('identity')
cores=json.loads((ROOT/'quotient-kernels.json').read_text())

def basis(M,z=1.):
 Q=z*M;I=np.eye(len(M));Q2=Q@Q;minus=np.linalg.inv(I-Q2);plus=np.linalg.inv(I+Q2)
 return np.array([-Q2@plus,Q2@minus,Q@minus,Q@plus])

@functools.lru_cache(None)
def spec(edges):
 letters='abcd';return ','.join(letters[a] if a==b else letters[a]+letters[b] for a,b in edges)+'->'

@functools.lru_cache(None)
def path(edges,n):
 dummy=[np.empty((n,)) if a==b else np.empty((n,n)) for a,b in edges]
 return np.einsum_path(spec(edges),*dummy,optimize='greedy')[0]

def contract(core,B):
 edges=tuple(tuple(x) for x in core['edges']);subs=spec(edges);pa=path(edges,B.shape[1]);diag=np.diagonal(B,axis1=1,axis2=2)
 val=0
 for t in core['terms']:
  ops=[diag[i] if a==b else B[i] for (a,b),i in zip(edges,t['indices'])]
  val+=t['numerator']*np.einsum(subs,*ops,optimize=pa)
 return val/core['denominator']

def evaluate(M,w,z=1.):
 B=basis(M,z);return sum(contract(c,B) for c in cores if c['weight']==w)

if __name__=='__main__':
 fixtures=json.loads((ENGINE/'fixtures.json').read_text());worst=0;checks=0
 for f in fixtures:
  M=np.array(f['S'])/len(f['S'])
  for z in [1.,.7+.1j]:
   B=basis(M,z)
   for w in [1,2]:
    expected=reference(M,w,z,detail=True)
    for c,(_,v) in zip([c for c in cores if c['weight']==w],expected):
     got=contract(c,B);err=abs(got-v);worst=max(worst,err);checks+=1
     if err>1e-10:raise RuntimeError((f['name'],w,c['edges'],v,got))
 # Benchmarks isolate repeated evaluations after compiling/caching both routes.
 bench=[]
 for n in [10,20,40]:
  q=n+1;S=np.array([[0 if i==j else 1 if 0<(j-i)%q<=q//2 else -1 for j in range(n)] for i in range(n)],float);M=S/n
  reference(M,2);evaluate(M,2)
  rr=[];ss=[]
  for _ in range(3):
   start=time.perf_counter();a=reference(M,2);rr.append(time.perf_counter()-start)
   start=time.perf_counter();b=evaluate(M,2);ss.append(time.perf_counter()-start)
   if abs(a-b)>1e-10:raise RuntimeError('benchmark mismatch')
  bench.append({'n':n,'repeats':3,'dense_four_basis_median_seconds':statistics.median(rr),'sparse_parity_quotient_median_seconds':statistics.median(ss),'ratio_dense_over_sparse':statistics.median(rr)/statistics.median(ss)})
 report={'status':'PASS','exact_scalar_identities':5,'core_contraction_checks':checks,'maximum_absolute_error':float(worst),'benchmark':bench}
 (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
