"""Addendum: independently check compressed tensor axes and hybrid execution.
Only this addendum and its receipt are written; prior audit files stay frozen.
"""
from pathlib import Path
import sys,json,itertools,time,statistics,os
import numpy as np
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE/'session30-parity-basis'))
import parity_engine as candidate
from decorated import evaluate as reference
source=json.loads((BASE/'session30-parity-basis/parity-kernels.json').read_text())
def key(c):return (c['weight'],tuple(map(tuple,c['edges'])))
compiled={key(c):(c,ob,p,W) for w in (1,2) for c,ob,p,W in candidate.compiled(w)}
assert set(compiled)=={key(c) for c in source}
checks=0;dense_size=compressed_size=0;fallback=[]
for c in source:
 cc,ob,p,W=compiled[key(c)];E=len(c['edges'])
 assert list(map(tuple,c['edges']))==list(map(tuple,ob.edges))
 assert c['V']==ob.V and c['denominator']==ob.G
 expected_spec='ijklmnopqrstuvwxyz'[:E]+','+','.join('ijklmnopqrstuvwxyz'[j]+('abcd'[a] if a==b else 'abcd'[a]+'abcd'[b]) for j,(a,b) in enumerate(c['edges']))+'->'
 assert ob.spec==expected_spec
 terms={tuple(t['indices']):t['numerator'] for t in c['terms']}
 for j in range(E):assert p[j]==sorted({idx[j] for idx in terms})
 assert W.shape==tuple(len(q) for q in p)
 for idx in itertools.product(range(4),repeat=E):
  value=W[tuple(p[j].index(k) for j,k in enumerate(idx))] if all(k in p[j] for j,k in enumerate(idx)) else 0
  assert value==terms.get(idx,0)==ob.Wb[idx];checks+=1
 dense_size+=ob.Wb.size;compressed_size+=W.size
 if all(len(q)==4 for q in p):fallback.append(c['edges'])
# Completely arbitrary dense matrices, no transpose or rational-function constraints.
# This explicitly checks contraction axes in both compressed and fallback branches.
rng=np.random.default_rng(4403);B=rng.integers(-2,3,size=(4,3,3),dtype=np.int64)
axis_checks=0;axis_error=0
for c in source:
 numerator=0
 for labels in itertools.product(range(3),repeat=c['V']):
  for t in c['terms']:
   product=t['numerator']
   for k,(a,b) in zip(t['indices'],c['edges']):product*=int(B[k,labels[a],labels[b]])
   numerator+=product
 actual=candidate.core_value(*compiled[key(c)],B)
 error=abs(actual-numerator/c['denominator']);axis_error=max(axis_error,error)
 assert error<1e-10;axis_checks+=1
# General and tournament fixtures; compare every core to old five-basis implementation.
fixtures=[('weighted3',np.array([[0,.1,.2],[-.1,0,-.15],[-.2,.15,0]])),('transitive4',np.sign(np.arange(4)[:,None]-np.arange(4)[None,:])/4)]
for n in (6,10):
 N=n+1;S=np.array([[0 if i==j else 1 if 0<(j-i)%N<=N//2 else -1 for j in range(n)] for i in range(n)])
 fixtures.append(('almostregular'+str(n),S/n))
num_checks=0;worst=0;rows=[]
for name,M in fixtures:
 for z in (1.,.37+.19j):
  for w in (1,2):
   got=candidate.evaluate(M,w,z,detail=True);expected=reference(M,w,z,detail=True,reduced=False)
   assert len(got)==len(expected)
   err=0
   for (a,v),(b,u) in zip(got,expected):
    assert key(a)==key(b)
    err=max(err,abs(v-u));num_checks+=1
   assert err<1e-10;worst=max(worst,err)
   rows.append({'fixture':name,'z':str(z),'weight':w,'core_count':len(got),'max_absolute_difference':float(err)})
# Warm single-thread local timings, alternating evaluation order each repetition.
bench=[]
for n in (10,20,40):
 N=n+1;M=np.array([[0 if i==j else 1 if 0<(j-i)%N<=N//2 else -1 for j in range(n)] for i in range(n)])/n
 reference(M,2);candidate.evaluate(M,2);times={'reference':[],'hybrid':[]}
 for repeat in range(7):
  order=[('reference',reference),('hybrid',candidate.evaluate)]
  if repeat%2:order.reverse()
  vals={}
  for name,fn in order:
   t=time.perf_counter();vals[name]=fn(M,2);times[name].append(time.perf_counter()-t)
  assert abs(vals['reference']-vals['hybrid'])<1e-10
 bench.append({'n':n,'repeats':7,'weight':2,'reference_seconds':statistics.median(times['reference']),'hybrid_seconds':statistics.median(times['hybrid']),'ratio_reference_over_hybrid':statistics.median(times['reference'])/statistics.median(times['hybrid'])})
print(json.dumps({'status':'PASS','all_axis_entries_checked':checks,'dense_entries':dense_size,'compressed_entries':compressed_size,'arbitrary_matrix_core_checks':axis_checks,'arbitrary_matrix_max_error':axis_error,'all_channel_fallback_cores':fallback,'original_engine_core_checks':num_checks,'original_engine_max_error':float(worst),'comparisons':rows,'benchmarks':bench,'thread_environment':{k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']},'timing_scope':'Warm single-thread local observations; includes basis inversion but excludes compilation, not a guaranteed speedup'},indent=2))
