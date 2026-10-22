import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
DEFAULT=HERE.parent/"session30-decorated-engine"/"fixtures.json"
import numpy as np
from decorated_engine import symbolic_kernel,evaluate
raw=json.load(open(sys.argv[1] if len(sys.argv)>1 else DEFAULT));cache={};out=[]
for row in raw:
 S=np.array(row['S']);M=S/len(S)
 for w,terms in row['terms'].items():
  err=0
  for term in terms:
   ed=tuple(tuple(x) for x in term['edges']);a=term['automorphisms'];key=(ed,a)
   if key not in cache:cache[key]=symbolic_kernel(ed,a)
   val=evaluate(ed,cache[key],M);err=max(err,abs(val-term['value']))
  out.append(dict(name=row['name'],weight=w,max_absolute_error=err,cores=len(terms)));print(out[-1],flush=True)
  if err>1e-10:raise RuntimeError('independent engines disagree')
json.dump(out,open(HERE/'engine_comparison.json','w'),indent=2)
