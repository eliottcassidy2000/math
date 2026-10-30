#!/usr/bin/env python3
"""Independent finite check by literal affine-set enumeration; no external packages."""
import json
from pathlib import Path

def demand(ok, message):
    if not ok:
        raise RuntimeError(message)

def gaussian(n,d):
    if d < 0 or d > n:
        return 0
    a=b=1
    for i in range(d):
        a *= (1 << (n-i))-1
        b *= (1 << (d-i))-1
    return a//b

def subspaces(n):
    levels=[{frozenset([0])}]
    for d in range(n):
        nxt=set()
        for u in levels[-1]:
            for x in range(1<<n):
                if x not in u:
                    nxt.add(u | frozenset(x^y for y in u))
        demand(len(nxt)==gaussian(n,d+1),'subspace enumeration count')
        levels.append(nxt)
    return levels

def check(n):
    N=1<<n
    spaces=subspaces(n)
    flats=[]
    for d,level in enumerate(spaces):
        for u in level:
            remaining=set(range(N))
            while remaining:
                a=min(remaining)
                w=frozenset(a^x for x in u)
                remaining.difference_update(w)
                flats.append((d,u,w))
    checks=0
    # All linear forms; three independently evaluated translations, both signs.
    for r in range(N):
        for v in sorted({0,1,N-1}):
            for eps in range(2):
                counts=[0]*(n+1)
                c=eps^((r&v).bit_count()%2)
                colors=[[0,0] for _ in range(n+1)]
                for d,u,w in flats:
                    vals={eps^((r&x).bit_count()%2)^(x==v) for x in w}
                    actual=len(vals)==1
                    expected=(d==0 or (v not in w and all((r&x).bit_count()%2==0 for x in u)) or (d==1 and v in w and any((r&x).bit_count()%2 for x in u)))
                    demand(actual==expected, f'classification n={n},r={r},v={v},d={d}')
                    if actual:
                        counts[d]+=1
                        colors[d][next(iter(vals))]+=1
                    checks+=1
                demand(counts[0]==N,'singleton count')
                for d in range(1,n+1):
                    expected=((1<<(n-d))-1)*gaussian(n if r==0 else n-1,d)
                    if r and d==1:
                        expected+=1<<(n-1)
                    demand(counts[d]==expected,'total count')
                    if r==0:
                        demand(colors[d][eps]==expected and colors[d][1-eps]==0,'zero-form colors')
                    elif d<n:
                        g=gaussian(n-1,d)
                        low=g*((1<<(n-d-1))-1)
                        high=g*(1<<(n-d-1)) + ((1<<(n-1)) if d==1 else 0)
                        demand(colors[d][c]==low and colors[d][1-c]==high,'nonzero-form colors')
    return dict(n=n,linear_subspaces=sum(map(len,spaces)),affine_flats=len(flats),literal_classifications=checks)

if __name__=='__main__':
    results=[check(n) for n in range(3,6)]
    receipt={'status':'PASS','method':'exhaustive affine sets; all linear forms; defect points 0,1,2^n-1; both output signs','results':results}
    print(json.dumps(receipt,indent=2))
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
