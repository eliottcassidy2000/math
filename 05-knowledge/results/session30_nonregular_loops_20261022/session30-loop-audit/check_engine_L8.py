import json,sys,numpy as np
from fractions import Fraction
from collections import defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'session30-decorated-engine'))
import decorated
refs=[json.loads(x)for x in(HERE/'partition-reference-L8.jsonl').read_text().splitlines()]
controls=[np.array([[0,2,-1],[-2,0,3],[1,-3,0]]),np.array([[0,1,1,1],[-1,0,1,1],[-1,-1,0,1],[-1,-1,-1,0]])]
report=[]
for ref,S in zip(refs,controls):
 M=S/ref['scale'];K=48;r=.45;vals=defaultdict(complex)
 for k in range(K):
  phase=np.exp(2j*np.pi*k/K);z=r*phase
  for w in(1,2):
   for d,value in decorated.evaluate(M,w,z,detail=True):vals[d['excess'],d['leaves']]+=value/phase**8/K/r**8
 rows=[]
 for row in ref['rows']:
  key=row['excess'],row['leaves'];exact=Fraction(row['numerator'],row['denominator']);observed=vals[key];err=abs(observed-float(exact))
  if err>2e-8:raise RuntimeError((ref['n'],key,str(exact),observed,err))
  rows.append(dict(row,observed_real=observed.real,observed_imag=observed.imag,absolute_error=err))
 report.append({'n':ref['n'],'L':8,'connected_partitions':ref['connected'],'rows':rows});print(ref['n'],'PASS',max(x['absolute_error']for x in rows))
(HERE/'engine-L8-results.json').write_text(json.dumps(report,indent=2)+'\n')
