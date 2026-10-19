#!/usr/bin/env python3
"""Replay the universal 2N-1 pair superset, derived by projection."""
import json,argparse
from word_projected_search import search
from independent_search import factor
p=argparse.ArgumentParser();p.add_argument('--output',default='session7-cycle-review/universal_pair_results.json');args=p.parse_args()
bs=[b for b in range(4,1001,4) if len(factor(b))>=3]
results=[]
for b in bs:
 N=2**len(factor(b))+2
 pairs=[(0,k) for k in range(N)]+[(k,k-1) for k in range(1,N)]
 r=search(b,pairs);results.append(r);print(json.dumps(r),flush=True)
with open(args.output,'w') as f:json.dump(dict(method='Independent exact search with universal projected pair superset',results=results),f,indent=2)
