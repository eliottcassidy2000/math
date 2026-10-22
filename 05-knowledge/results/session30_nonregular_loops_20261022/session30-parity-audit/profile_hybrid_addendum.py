from pathlib import Path
import sys,time,json,numpy as np
BASE=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(BASE/'session30-parity-basis'))
import parity_engine as pe
from decorated import compiled_reduced,evaluate
from loop_engine_fast import basis
out=[]
# Match the producer benchmark's initial reference path planning at n=10.
for n in (10,20):
 N=n+1;M=np.array([[0 if i==j else 1 if 0<(j-i)%N<=N//2 else -1 for j in range(n)] for i in range(n)])/n
 B=basis(M)[:4]
 for d,c in compiled_reduced(2):
  t=time.perf_counter();v=c.value_b(B);dt=time.perf_counter()-t
  out.append({'n':n,'edges':d['edges'],'seconds':dt,'path':c.pathb})
 t=time.perf_counter();a=evaluate(M,2);dt=time.perf_counter()-t
 t=time.perf_counter();b=pe.evaluate(M,2);dh=time.perf_counter()-t
 print(json.dumps({'n':n,'reference_seconds':dt,'hybrid_seconds':dh,'error':abs(a-b)}),flush=True)
(Path(__file__).parent/'hybrid-profile-receipt.json').write_text(json.dumps(out,indent=2))
