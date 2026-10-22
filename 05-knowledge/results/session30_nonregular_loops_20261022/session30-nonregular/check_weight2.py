import json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
import numpy as np
from decorated_engine import cores,symbolic_kernel,evaluate
cache=[]
for e,l in [(-1,3),(0,2),(1,1),(2,0)]:
 for ed,a in cores(e,l):cache.append((e,l,ed,a,symbolic_kernel(ed,a)))
res=[]
for name,n in [('transitive4',4),('transitive6',6),('R7',7),('R7_delete0',6)]:
 if name.startswith('transitive'):S=np.triu(np.ones((n,n)),1);S=S-S.T
 else:
  N=7;S=np.array([[0 if i==j else (1 if (j-i)%N<=N//2 else -1) for j in range(N)]for i in range(N)],dtype=float)
  if n==6:S=S[1:,1:]
 M=S/n;groups={}
 for e,l,ed,a,sy in cache:
  groups[str((e,l))]=groups.get(str((e,l)),0)+evaluate(ed,sy,M)
 rec={'name':name,'groups':groups,'total':sum(groups.values())};res.append(rec);print(rec,flush=True)
with open(HERE/'weight2_fixtures.json','w')as f:json.dump(res,f,indent=2)
