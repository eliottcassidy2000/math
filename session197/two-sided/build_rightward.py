from itertools import combinations
from collections import deque
from pathlib import Path
import json
from rightward import brute,step,H
T=[122,112,67,102,108,4,32];X=151;cases=[];lines=[]
for a,b in combinations(range(7),2):
 exact=brute(T,[a,b]);initial=tuple(min(X+1,x) for x in exact);seen={initial};queue=deque([initial]);support=set()
 while queue:
  z=queue.popleft();support.add(H(z,2))
  for S in range(4):
   for d in range(3):
    nz=tuple(min(X+1,x) for x in step(z,2,S,d))
    if H(nz,2)<=X and nz not in seen:seen.add(nz);queue.append(nz)
 cases.append({'ports':[a,b],'uncapped_seed':exact,'capped_seed':initial,'states':sorted(seen),'support':sorted(support)})
 lines.append(' '.join(map(str,[a,b]+list(initial))))
 print(a,b,len(seen),sorted(support))
Path('rightward_certificates.json').write_text(json.dumps({'rows':T,'bound':X,'cases':cases},indent=2)+'\n')
Path('rightward_seeds.txt').write_text('\n'.join(lines)+'\n')
