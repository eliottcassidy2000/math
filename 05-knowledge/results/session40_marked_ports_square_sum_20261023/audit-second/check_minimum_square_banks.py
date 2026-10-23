#!/usr/bin/env python3
"""Independent exhaustive combinations lower bound, no hitting-set branch solver."""
import json,itertools,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
source=next((p for p in (ROOT/'session40-square-sum/receipt.json', ROOT/'square-sum/receipt.json') if p.exists()),None)
if source is None: raise RuntimeError('Missing square-sum receipt.json beside published audit packet')
r=json.load(open(source));out=[]
for row in r['rows'][24:]:
 n=row['n'];P=row['paths'];sets=[{tuple(sorted(e)) for e in zip(p,p[1:])} for p in P];ends=[tuple(sorted((p[0],p[-1]))) for p in P]
 pairs=[(i,j) for i,j in itertools.combinations(range(len(P)),2) if ends[i]==ends[j]]
 E=[(a,b) for a in range(1,n+1) for b in range(a+1,n+1) if math.isqrt(a+b)**2==a+b]
 masks=set(sum(1<<k for k,(i,j) in enumerate(pairs) if (e in sets[i])!=(e in sets[j])) for e in E);masks.discard(0);masks=sorted(masks)
 d=row['endpoint_sensor']['minimum'];target=(1<<len(pairs))-1;bank=[tuple(e) for e in row['endpoint_sensor']['bank']]
 sig={(ends[i],tuple(e in sets[i] for e in bank)) for i in range(len(P))}
 if len(sig)!=len(P):raise RuntimeError('bank upper bound')
 tested=0
 if d:
  for chosen in itertools.combinations(masks,d-1):
   tested+=1;covered=0
   for m in chosen:covered|=m
   if covered==target:raise RuntimeError(f'lower bound false n{n}')
 else:
  if pairs:raise RuntimeError('zero bank failure')
 out.append({'n':n,'minimum':d,'distinct_nonzero_edge_tests':len(masks),'smaller_banks_exhaustively_rejected':tested,'endpoint_pair_ambiguities':len(pairs)})
print(out);Path(__file__).with_name('minimum-square-banks-receipt.json').write_text(json.dumps({'status':'PASS','rows':out},indent=2)+'\n')
