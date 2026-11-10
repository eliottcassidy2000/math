#!/usr/bin/env python3
"""Direct permutation counts, exhaustive chronological forests through order6."""
import csv,itertools,json,sys
from pathlib import Path
folder=Path(sys.argv[1]);summary=[]
for n in range(1,7):
    values={};histories=0
    for choices in itertools.product(*(range(v+1) for v in range(1,n))):
        a=[[False]*n for _ in range(n)];depth=[0]*n;code=0
        for v,c in enumerate(choices,1):
            code=code*(v+1)+c
            depth[v]=depth[c-1]+1 if c else 0
            for u in range(v):a[u][v]=c!=u+1;a[v][u]=c==u+1
        h=sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in itertools.permutations(range(n)))
        values[h]=min(values.get(h,99),max(depth));histories+=1
    with (folder/f'n{n}_witnesses.csv').open() as f:rows=list(csv.DictReader(f))
    observed={int(row['H']):int(row['min_forest_depth']) for row in rows}
    if observed!=values:raise RuntimeError(('forest census mismatch',n))
    summary.append({'n':n,'histories':histories,'distinct_values':len(values),'minimum_depths_match':True})
print(json.dumps({'status':'PASS','method':'all histories independently reconstructed; Hamilton paths by full permutation enumeration','census':summary},indent=2))
