"""Constant-arc anti-involution obstructions for the context lifts."""
import json
from pathlib import Path
from five_atom_lift import inflate
HERE=Path(__file__).resolve().parent

def symbolic(m,row):
    n=m+5
    # None denotes a free context arc. Known path/source arcs are retained.
    a=[[None]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1,n):
            if u>=m:value=bool(row['atom'][u-m]>>(v-m)&1)
            elif v>=m or u==0 or v==u+1:value=True
            else:continue
            a[u][v]=value;a[v][u]=not value
    records=[]
    for raw in row['paths']:
        p=inflate(raw,m); rev=dict(zip(p,p[::-1])); bad=[]
        for u in range(n):
            for v in range(u+1,n):
                r,s=rev[v],rev[u]
                if a[u][v] is not None and a[r][s] is not None and a[u][v]!=a[r][s]:
                    bad.append((u,v,r,s))
        def name(v):return f'b{v-m}' if v>=m else ('u' if v==0 else f'A{v}')
        records.append([[name(v) for v in q] for q in bad[:2]])
        for u,v in zip(p,p[1:]):
            if a[u][v] is not None:a[u][v]=not a[u][v];a[v][u]=not a[v][u]
    return records

if __name__=='__main__':
    rows=json.loads(HERE.joinpath('five_atom_probe.json').read_text())
    for i,row in enumerate(rows):
        print(json.dumps(dict(atom=i,results={m:symbolic(m,row) for m in (1,2,3,4,5,6,7,8,10,20)})))
