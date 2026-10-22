import json,time,statistics,pathlib,numpy as np
from parity_engine import evaluate,compiled,CORES,ROOT
from decorated import evaluate as reference
worst=0;checks=0
for f in json.loads((ROOT.parent/'session30-decorated-engine/fixtures.json').read_text()):
 M=np.array(f['S'])/len(f['S'])
 for z in [1.,.7+.1j]:
  for w in [1,2]:
   for (c,a),(_,b) in zip(evaluate(M,w,z,True),reference(M,w,z,True)):
    worst=max(worst,abs(a-b));checks+=1
    assert abs(a-b)<1e-10
bench=[]
for n in [10,20,40]:
 q=n+1;S=np.array([[0 if i==j else 1 if 0<(j-i)%q<=q//2 else -1 for j in range(n)] for i in range(n)],float);M=S/n
 reference(M,2);evaluate(M,2)
 rr=[];ss=[]
 for _ in range(7):
  t=time.perf_counter();a=reference(M,2);rr.append(time.perf_counter()-t)
  t=time.perf_counter();b=evaluate(M,2);ss.append(time.perf_counter()-t)
  assert abs(a-b)<1e-10
 bench.append({'n':n,'repeats':7,'reference_seconds':statistics.median(rr),'parity_compressed_seconds':statistics.median(ss),'ratio_reference_over_parity':statistics.median(rr)/statistics.median(ss)})
shapes=[]
for w in [1,2]:
 for c,ob,p,W in compiled(w):shapes.append({'weight':w,'edges':c['edges'],'palettes':p,'dense_size':ob.Wb.size,'compressed_size':W.size})
res={'status':'PASS','checks':checks,'max_error':worst,'benchmarks':bench,'dense_size':sum(c['dense_size'] for c in shapes),'compressed_size':sum(c['compressed_size'] for c in shapes),'cores':shapes}
(ROOT/'compressed-results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='cores'},indent=2))
