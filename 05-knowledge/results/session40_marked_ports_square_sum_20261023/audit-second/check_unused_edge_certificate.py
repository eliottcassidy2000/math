#!/usr/bin/env python3
import json,itertools,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
source=next(p for p in (ROOT/'session40-square-sum/receipt.json',ROOT/'square-sum/receipt.json') if p.exists())
r=json.load(open(source));cert=r['edge_1_15_exclusion_at25'];E={(a,b) for a in range(1,26) for b in range(a+1,26) if math.isqrt(a+b)**2==a+b};steps=0
expected={tuple(sorted((18,v))) for v in range(1,26) if v!=18}
if {tuple(x['endpoints']) for x in cert}!=expected:raise RuntimeError('endpoint coverage')
for c in cert:
 ends=set(c['endpoints']);target={v:1 if v in ends else 2 for v in range(1,26)};A=set(E);F={(1,15)}
 for operation,v,edges in c['trace']:
  av={e for e in A if v in e};fv=F&av;edges={tuple(e) for e in edges};steps+=1
  if edges!=av-F:raise RuntimeError('trace incidence mismatch')
  if operation=='delete':
   if len(fv)!=target[v]:raise RuntimeError('unsound deletion')
   A-=edges
  elif operation=='force':
   if len(av)!=target[v]:raise RuntimeError('unsound forcing')
   F|=edges
  else:raise RuntimeError('unknown action')
  if not F<=A:raise RuntimeError('forced edge removed')
 reason=c['contradiction'];kind=reason[0]
 if kind=='excess':
  v=reason[1]
  if not sum(v in e for e in F)>target[v]:raise RuntimeError('false excess')
 elif kind=='deficit':
  v=reason[1]
  if not sum(v in e for e in A)<target[v]:raise RuntimeError('false deficit')
 elif kind=='forced_cycle':
  S=set(reason[1])
  if sum(set(e)<=S for e in F)<len(S):raise RuntimeError('false cycle')
 else:raise RuntimeError('unproved branch')
out={'status':'PASS','endpoint_cases':len(cert),'validated_propagation_steps':steps,'forced_edge':[1,15],'n':25,'branching_used':False}
Path(__file__).with_name('unused-edge-certificate-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
