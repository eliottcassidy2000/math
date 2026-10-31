#!/usr/bin/env python3
"""Exact local cyclic matchings and the global connected-wiring selector."""
from itertools import product
from math import factorial
from collections import Counter
import json
def need(b,msg):
    if not b:raise RuntimeError(msg)
def matchings(vertices):
    if not vertices:yield ();return
    a=vertices[0]
    for i in range(1,len(vertices)):
        b=vertices[i]
        for rest in matchings(vertices[1:i]+vertices[i+1:]):yield ((a,b),)+rest
def components(P,M,n):
    adj={i:[] for i in range(n)}
    for a,b in P+M:adj[a].append(b);adj[b].append(a)
    need(all(len(x)==2 for x in adj.values()),'two matching incidences at every slot')
    seen=set();count=0
    for i in range(n):
        if i in seen:continue
        count+=1;todo=[i];seen.add(i)
        while todo:
            a=todo.pop()
            for b in adj[a]:
                if b not in seen:seen.add(b);todo.append(b)
    return count
def main():
    counts=[]
    for r in range(1,7):
        original=tuple((2*i,2*i+1) for i in range(r));total=0;connected=0
        for P in matchings(tuple(range(2*r))):
            total+=1;connected+=components(P,original,2*r)==1
        need(connected==2**(r-1)*factorial(r-1),'cyclic local matching formula')
        counts.append({'ports':r,'all_matchings':total,'single_source_cycle_matchings':connected})
    cyclic=[]
    P=tuple((4*i+a,4*i+b) for i in range(3) for a,b in [(0,3),(1,2)])
    for bits in product((0,1),repeat=3):
        M=[]
        for i,flip in enumerate(bits):
            a=[4*i+2,4*i+3];j=(i+1)%3;b=[4*j,4*j+1]
            if flip:b.reverse()
            M+=list(zip(a,b))
        count=components(P,tuple(M),12)
        need((count==1)==(sum(bits)%2==0),'retained zero-length path matches native incidence hostile')
        cyclic.append({'pairing_bits':bits,'cycles':count})
    # Three-leaf tree: every cyclic matching of the six center slots works.
    original=((0,1),(2,3),(4,5));leafpaths=((6,7),(8,9),(10,11))
    wires=((0,6),(1,7),(2,8),(3,9),(4,10),(5,11));valid=0
    for center in matchings(tuple(range(6))):
        if components(center,original,6)!=1:continue
        need(components(center+leafpaths,wires,12)==1,'all local matching states glue on a tree');valid+=1
    need(valid==8,'complete three-port center state space')
    original=((0,1),(2,3),(4,5),(6,7));atlas=Counter();selected=Counter();table=Counter()
    for P in matchings(tuple(range(8))):
        directions=len({a^b for a,b in P});k=components(P,original,8)
        atlas[directions]+=1;selected[directions]+=k==1;table[directions,k]+=1
    need([atlas[d] for d in [1,2,4]]==[7,42,56],'inherited complete atlas partition')
    need([selected[d] for d in [1,2,4]]==[0,24,24],'new connected-cycle atlas restriction')
    print(json.dumps({'status':'passed','local_matching_counts':counts,'cyclic_selector':cyclic,
          'tree_center_states':valid,'atlas_counts':dict(atlas),'connected_atlas_counts':dict(selected),
          'atlas_component_table':[{'directions':d,'components':k,'count':n} for (d,k),n in sorted(table.items())],
          'scope':'Exact boundary incidence matching and marked restriction of the prior105-atlas domain; repeated physical vertices have separate port slots and zero-length paths. No cyclic arithmetic or geometric-selector realization is asserted.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
