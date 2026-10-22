"""Runtime-only equal-cap path planning, without editing either engine."""
from pathlib import Path
import sys,json,time,statistics,numpy as np
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE/'session30-parity-basis'))
import parity_engine as pe
from decorated import evaluate,compiled_reduced
from loop_engine_fast import basis
out=[]
for n in (10,20,40):
 N=n+1;M=np.array([[0 if i==j else 1 if 0<(j-i)%N<=N//2 else -1 for j in range(n)] for i in range(n)])/n
 B=basis(M)[:4]
 for d,c in compiled_reduced(2):
  ops=[np.einsum('kii->ki',B) if a==b else B for a,b in c.edges]
  c.pathb=np.einsum_path(c.spec,c.Wb,*ops,optimize=('greedy',2**24))[0]
 pe.PATHS.clear()
 BP=pe.basis(M)
 for d,c,p,W in pe.compiled(2):
  # Match memory cap also for the all-channel fallback branch.
  ops=[np.einsum('kii->ki',BP) if a==b else BP for a,b in c.edges]
  c.pathb=np.einsum_path(c.spec,c.Wb,*ops,optimize=('greedy',2**24))[0]
 evaluate(M,2);pe.evaluate(M,2)
 ts={'old':[],'hybrid':[]}
 for i in range(7):
  order=[('old',evaluate),('hybrid',pe.evaluate)]
  if i%2:order.reverse()
  vals={}
  for name,f in order:
   start=time.perf_counter();vals[name]=f(M,2);ts[name].append(time.perf_counter()-start)
  assert abs(vals['old']-vals['hybrid'])<1e-10
 out.append({'n':n,'repeats':7,'old_median_seconds':statistics.median(ts['old']),'hybrid_median_seconds':statistics.median(ts['hybrid']),'ratio_old_over_hybrid':statistics.median(ts['old'])/statistics.median(ts['hybrid'])})
print(json.dumps({'status':'PASS','planning':'Both generic engines replanned at each n with greedy intermediate cap2**24; runtime cache changes only; K4/doubled-square special routines retained','benchmarks':out},indent=2))
