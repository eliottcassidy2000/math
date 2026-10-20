#!/usr/bin/env python3
"""Independent literal shared-edge hub check on every source class through n7.
Helper code belongs to this independent audit, not to the producer.
"""
from pathlib import Path
import json
from collections import defaultdict
from independent_controls import check, decode, encode, reflect, edge_key, paths, triangle, connected
ROOT=Path(__file__).resolve().parent

def K(rows,path):
    full=(1<<len(rows))-1
    out=[full^r^(1<<i) for i,r in enumerate(rows)]
    for a,b in zip(path,path[1:]):out[a]^=1<<b;out[b]^=1<<a
    return out

def verify(n,data):
    result={'n':n,'nonstrong_sources':0,'qualifying_shared_edge_groups':0,'literal_alternate_routes':0,'first_second_last_routes':0}
    for node in data['nodes']:
        rows=decode(node['representative'],n)
        if connected(rows):continue
        result['nonstrong_sources']+=1
        degrees=[r.bit_count() for r in rows];groups=defaultdict(dict)
        for path in paths(rows):
            t=encode(rows,path)
            if reflect(t,n)==t:continue
            e=edge_key(t,n)
            for u,v in zip(path,path[1:]):groups[path[0],path[-1],u,v].setdefault(e,path)
        for (a,z,u,v),group in groups.items():
            if len(group)<2:continue
            endpoint_sum=degrees[u]+degrees[v]-(u==a)+(v==z)
            if endpoint_sum==n-1:continue
            q=[0]*n;q[a]-=1;q[z]+=1;q[u]-=1;q[v]+=1
            twice_delta=-2*sum(d*c for d,c in zip(degrees,q))-sum(c*c for c in q)
            check(twice_delta%2==0,('integral delta',n,node['node']))
            delta=twice_delta//2
            if not delta:continue
            expected=rows[:]
            for s,t in ((a,z),(u,v)):expected[s]^=1<<t;expected[t]^=1<<s
            check(triangle(expected)-triangle(rows)==delta,('literal hub potential',n,node['node']))
            hubs={};edges={};partners={}
            for e,path in group.items():
                out=K(rows,path);idx=path.index(v);rot=path[idx:]+path[:idx]
                check(all((out[x]>>y)&1 for x,y in zip(rot,rot[1:])),('rotated path',n,node['node']))
                rt=encode(out,rot)
                check(reflect(rt,n)!=rt,('rotation black',n,node['node']))
                final=K(out,rot)
                check(final==expected,('literal common hub',n,node['node']))
                check(connected(out),('strong complement partner',n,node['node']))
                hubs[e]=final;edges[e]=edge_key(rt,n);partners[e]=out
            for e in group:
                f=next(f for f in group if f!=e)
                check(e not in (f,edges[e],edges[f]),('avoids designated edge',n,node['node'],e,f))
                result['literal_alternate_routes']+=1
                if u==a:result['first_second_last_routes']+=1
            result['qualifying_shared_edge_groups']+=1
    return result

out=[]
for n in range(3,8):
    data=json.loads((ROOT/'inputs'/f'census_n{n}_representatives.json').read_text())
    row=verify(n,data);out.append(row);print(row,flush=True)
result={'status':'FINITE-EXACT; all literal checks passed','orders':out}
check(result==json.loads((ROOT/'results'/'independent_common_hubs.json').read_text()),'Bundled result mismatch')
print('Matched bundled independent common-hub result.',flush=True)
