#!/usr/bin/env python3
"""Independent combinatorial audit, using explicit simple paths, not producer DP."""
from itertools import permutations
from collections import Counter
from pathlib import Path
from hashlib import sha256
import json

HERE = Path(__file__).resolve().parent
A = (6,28,8,1,109,15,47)
B = (38,28,8,1,13,94,31)
ROOT_A, ROOT_B = 2, 3

def need(ok, message="independent audit check failed"):
    if not ok: raise RuntimeError(message)

def valid(rows):
    n=len(rows)
    need(all(0 <= r < (1<<n) and not (r>>i)&1 for i,r in enumerate(rows)))
    need(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i)))

def strong(rows):
    n=len(rows)
    for root in range(n):
        reached={root}; todo=[root]
        while todo:
            i=todo.pop()
            for j in range(n):
                if (rows[i]>>j)&1 and j not in reached:
                    reached.add(j); todo.append(j)
        if len(reached)!=n:return False
    return True

_cache={}
def profile(rows,root):
    rows=tuple(rows); key=(rows,root)
    if key in _cache:return _cache[key]
    valid(rows)
    n=len(rows);N=1<<n;h=[0]*N;e=[0]*N;starts=[0]*N; h[0]=1
    # Visit each actual simple directed path exactly once. Retain its exact vertex mask.
    def visit(path,mask):
        h[mask]+=1
        if path[-1]==root:e[mask]+=1
        if path[0]==root:starts[mask]+=1
        for j in range(n):
            if not (mask>>j)&1 and (rows[path[-1]]>>j)&1:
                visit(path+(j,),mask|(1<<j))
    for i in range(n):visit((i,),1<<i)
    full=N-1
    ans=dict(H=h[full],s=starts[full],e=e[full],F=sum(e[M]*h[full^M] for M in range(N)),
             q=h[full^(1<<root)],strong=strong(rows),
             S2=sum(starts[M]*h[full^M] for M in range(N)),P2=sum(h[M]*h[full^M] for M in range(N)),
             F_by_size=[sum(e[M]*h[full^M] for M in range(N) if M.bit_count()==k) for k in range(n+1)],
             H_by_size=[sum(h[M] for M in range(N) if M.bit_count()==k) for k in range(n+1)])
    _cache[key]=ans
    return ans

def relabel(rows,p):
    return tuple(sum(((rows[p[i]]>>p[j])&1)<<j for j in range(len(rows))) for i in range(len(rows)))

def edges_between(a,b):
    return tuple((i,j) for i in range(len(a)) for j in range(i+1,len(a)) if ((a[i]^b[i])>>j)&1)

def flip(rows,edge):
    i,j=edge;out=list(rows);out[i]^=1<<j;out[j]^=1<<i;return tuple(out)

def main():
    dist=Counter(); candidates=[]
    for p in permutations(range(7)):
        if p[ROOT_A]!=ROOT_B:continue
        b=relabel(B,p);edges=edges_between(A,b);d=len(edges);dist[d]+=1
        candidates.append((d,p,b,edges))
    mind=min(dist); minimizers=[c for c in candidates if c[0]==mind]
    paths=[]; cubes=[]
    for d,p,b,edges in minimizers:
        cube={}
        for mask in range(1<<d):
            rows=A
            for k,edge in enumerate(edges):
                if (mask>>k)&1:rows=flip(rows,edge)
            cube[mask]=dict(mask=mask,rows=rows,**profile(rows,ROOT_A))
        cubes.append(dict(map=p,edges=edges,vertices=list(cube.values())))
        for order in permutations(range(d)):
            masks=[0]
            for k in order:masks.append(masks[-1]|(1<<k))
            records=[cube[M] for M in masks]
            paths.append(dict(map=p,edges=edges,order=order,masks=masks,
                              all_strong=all(r['strong'] for r in records),
                              unit_start=all(r['s']==1 for r in records),
                              peak_start=max(r['s'] for r in records),
                              peak_H=max(r['H'] for r in records),
                              trough_H=min(r['H'] for r in records)))
    st=[p for p in paths if p['all_strong']]
    target_map=(2,0,3,1,5,4,6);target_order=(1,0,2,3)
    selected=next(p for p in paths if p['map']==target_map and p['order']==target_order)
    cube=next(c for c in cubes if c['map']==target_map)
    vertices={r['mask']:r for r in cube['vertices']}
    selected['records']=[vertices[M] for M in selected['masks']]
    need(sum(dist.values())==720)
    need(mind==4 and len(minimizers)==5 and len(st)==38)
    need(min(p['peak_start'] for p in st)==2)
    need(not any(p['unit_start'] for p in paths))
    need([(r['H'],r['s'],r['F']) for r in selected['records']]==[(29,1,157),(41,2,134),(29,1,89),(25,1,121),(31,1,157)])
    out=dict(status='PASS: finite local construction only',method='Explicit enumeration of all actual simple directed paths on each induced vertex set; exhaustive root-compatible relabelings and shortest edge orders.',
             conventions={'rows':'bit j of row i means i -> j','map':'p[i] is the original B label assigned to new label i','root_A':ROOT_A,'root_B':ROOT_B},
             endpoints={'A':dict(rows=A,root=ROOT_A,**profile(A,ROOT_A)),'B':dict(rows=B,root=ROOT_B,**profile(B,ROOT_B))},
             root_compatible_relabelings=720,distance_distribution=dict(sorted(dist.items())),minimum_edge_distance=mind,
             minimum_maps=len(minimizers),shortest_edge_orders=len(paths),all_strong_shortest_paths=len(st),
             unit_start_shortest_paths=sum(p['unit_start'] for p in paths),
             minimum_peak_start_among_all_strong=min(p['peak_start'] for p in st),
             minimum_peak_H_among_all_strong=min(p['peak_H'] for p in st),
             selected=selected,cubes=cubes,paths=paths)
    (HERE/'local_audit_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('cubes','paths','selected')},indent=2))
    print('Selected path:',[(r['H'],r['s'],r['F']) for r in selected['records']])
if __name__=='__main__':main()
