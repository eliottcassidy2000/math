#!/usr/bin/env python3
"""Complete small tournament banks by explicit permutation orbits, n<=6."""
from itertools import combinations,permutations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent


def bank(n):
    pairs=list(combinations(range(n),2));idx={p:i for i,p in enumerate(pairs)}
    relabel=[]
    for p in permutations(range(n)):
        relabel.append([(idx[tuple(sorted((p[i],p[j])))],int(p[i]<p[j])) for i,j in pairs])
    seen=set();out=[]
    for mask in range(1<<len(pairs)):
        if mask in seen:continue
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            a,b=(i,j)if mask>>k&1 else (j,i);rows[a]|=1<<b
        out.append(rows)
        for mapping in relabel:
            image=0
            for k,(dest,forward) in enumerate(mapping):
                if bool(mask>>k&1)==bool(forward):image|=1<<dest
            seen.add(image)
    if len(seen)!=(1<<len(pairs)):raise RuntimeError('orbit bank incomplete')
    return out


if __name__=='__main__':
    result={n:bank(n)for n in range(1,7)}
    expected=[1,1,2,4,12,56]
    if [len(result[n])for n in range(1,7)]!=expected:raise RuntimeError('small tournament census')
    (HERE/'small_banks.json').write_text(json.dumps(result,indent=2)+'\n')
    for n,rows in result.items():
        (HERE/f'bank{n}.txt').write_text(str(len(rows))+'\n'+'\n'.join(' '.join(map(str,t))for t in rows)+'\n')
    print({n:len(r)for n,r in result.items()})
