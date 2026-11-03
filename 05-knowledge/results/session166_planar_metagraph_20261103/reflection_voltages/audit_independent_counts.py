from itertools import permutations, product
from collections import defaultdict, Counter
import json

# Independent isomorphism convention: order vertices by their score, then
# take the minimum upper-triangle word over permutations within score cells.
def counts(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    chords=[(i,j) for i in range(n) for j in range(i+2,n)]
    m=len(chords); full=(1<<m)-1; ci={p:i for i,p in enumerate(chords)}
    key_cache={}; order_cache={}
    def canon(rows):
        rows=tuple(rows)
        if rows in key_cache:return key_cache[rows]
        groups=defaultdict(list)
        for i,row in enumerate(rows): groups[row.bit_count()].append(i)
        cells=tuple(tuple(groups[d]) for d in sorted(groups))
        if cells not in order_cache:
            order_cache[cells]=[sum(blocks,()) for blocks in product(*(tuple(permutations(g)) for g in cells))]
        orders=order_cache[cells]
        best=None
        for p in orders:
            value=0
            for i,j in pairs:value=(value<<1)|((rows[p[i]]>>p[j])&1)
            if best is None or value<best:best=value
        key_cache[rows]=best
        return best
    c=[]; opp={}
    for t in range(1<<m):
        rows=[0]*n
        for i,j in pairs:
            backward=j>=i+2 and ((t>>ci[i,j])&1)
            a,b=(j,i) if backward else (i,j)
            rows[a]|=1<<b
        key=canon(rows)
        c.append(key)
        if key not in opp:
            rev=[(((1<<n)-1)^(1<<i)^rows[i]) for i in range(n)]
            opp[key]=canon(rev)
    classes=set(c); fixed={v for v in classes if opp[v]==v}
    free={min(v,opp[v]) for v in classes-fixed}
    sbit=[ci[n-1-j,n-1-i] for i,j in chords]
    def refl(t):
        return sum(((t>>i)&1)<<j for i,j in enumerate(sbit))
    blue=0;black=0;nse=0;loops=Counter();seen=set();adj=defaultdict(set)
    for t in range(1<<(m-1)):
        f=t^full;s=refl(t)
        if s==t:blue+=1;continue
        edge=min(t,f); reflected=min(s,s^full)
        orbit=min(edge,reflected)
        if orbit in seen:continue
        seen.add(orbit);black+=1
        a,b=c[t],c[f]
        u,v=min(a,opp[a]),min(b,opp[b])
        if a not in fixed and b not in fixed:
            nse+=1;adj[u].add(v);adj[v].add(u)
            if u==v:loops['positive' if a==b else 'negative']+=1
    unseen=set(free);components=0
    while unseen:
        root=unseen.pop();stack=[root];components+=1
        while stack:
            for v in adj[stack.pop()]:
                if v in unseen:unseen.remove(v);stack.append(v)
    return dict(n=n,classes=len(classes),SC=len(fixed),NS=len(free),blue=blue,black_orbits=black,NS_edges=nse,NS_components=components,beta=nse-len(free)+components,NS_loops=dict(loops))
print(json.dumps([counts(n) for n in (6,7)],indent=2))
