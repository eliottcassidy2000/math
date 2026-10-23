import json, hashlib, math
from fractions import Fraction
from collections import defaultdict
from pathlib import Path

def rank(rows):
    piv={}
    for row in rows:
        a=list(map(Fraction,row))
        for j,p in sorted(piv.items()):
            if a[j]:
                c=a[j];a=[x-c*y for x,y in zip(a,p)]
        k=next((j for j,x in enumerate(a) if x),None)
        if k is not None:
            v=a[k];piv[k]=[x/v for x in a]
    return len(piv)

base=Path(__file__).resolve().parent
source=base/'square-sum'/'receipt.json'
if not source.exists(): source=base.parent/'session40-square-sum'/'receipt.json'
out=[]
for r in json.loads(source.read_text())['rows']:
 if not r['paths'] or r['n']<15:continue
 es=[tuple(e['edge']) for e in r['edge_counts']]; rows=[];gs=defaultdict(list)
 for p in r['paths']:
    if sorted(p)!=list(range(1,r['n']+1)): raise ValueError('not spanning')
    if any(math.isqrt(a+b)**2!=a+b for a,b in zip(p,p[1:])): raise ValueError('not square sum')
    used={tuple(sorted(e)) for e in zip(p,p[1:])};v=[int(e in used) for e in es];rows.append(v);gs[tuple(sorted([p[0],p[-1]]))].append(v)
 diff=[[x-y for x,y in zip(v,rows[0])] for v in rows[1:]]
 wd=[[x-y for x,y in zip(v,g[0])] for g in gs.values() for v in g[1:]]
 # Every fixed-endpoint trade conserves every vertex degree exactly.
 for d in wd:
    for a in range(1,r['n']+1):
        if sum(d[j] for j,e in enumerate(es) if a in e)!=0:raise ValueError('degree conservation')
 B=[[int(a in e) for e in es] for a in range(1,r['n']+1)]
 if r['n']==32:
    extra=[[int(e==f) for e in es] for f in [(5,31),(9,27)]]
    if any(sum(a*b for a,b in zip(d,v)) for d in wd for v in extra): raise ValueError('conditional charge')
    if rank(B)!=32 or rank(B+extra)!=34 or rank(wd)!=12: raise ValueError('kernel basis extension')
 out.append(dict(n=r['n'],paths=len(rows),edges=len(es),endpoint_pairs=len(gs),max_endpoint_multiplicity=max(map(len,gs.values())),total_scatter_rank=rank(diff),within_endpoint_scatter_rank=rank(wd),signless_incidence_rank=rank([[int(a in e) for e in es] for a in range(1,r['n']+1)])))
print(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'ranks':out,'checks':'path validation, all degree trades, exact ranks, Q32 kernel basis extension'},indent=2))
