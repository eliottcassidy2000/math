#!/usr/bin/env python3
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
source=next(p for p in (ROOT/'session40-square-sum/receipt.json',ROOT/'square-sum/receipt.json') if p.exists())
a=json.load(open(source))['rows'];b=json.load(open(HERE/'square-edge-subsets-receipt.json'))['rows'];checks=0
for x,y in zip(a,b):
 if x['n']!=y['n'] or x['count']!=y['paths_up_to_reversal']:raise RuntimeError('path total')
 checks+=1
 ed={':'.join(map(str,z['edge'])):z['count'] for z in x['edge_counts']}
 if ed!=y['edge_counts']:raise RuntimeError('edge counts')
 checks+=len(ed)
 if x['n']>1:
  ep={':'.join(map(str,z['pair'])):z['count'] for z in x['endpoints']}
  if ep!=y['endpoint_counts']:raise RuntimeError('endpoint counts')
  checks+=len(ep)
out={'status':'PASS','max_order':25,'checks':checks,'methods':'vertex-order DFS versus acyclic degree-constrained edge subsets'}
(HERE/'square-census-agreement-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
