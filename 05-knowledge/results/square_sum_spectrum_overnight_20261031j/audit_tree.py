#!/usr/bin/env python3
"""Independent literal Hamiltonian-cycle consumer for port-linear trees."""
from itertools import product
import json

def need(b,msg):
    if not b:raise RuntimeError(msg)
def edge(a,b):return tuple(sorted((a,b)))
def cycles(E):
    V={v for e in E for v in e};adj={v:set() for v in V}
    for a,b in E:adj[a].add(b);adj[b].add(a)
    root=min(V);out=set()
    def visit(path,used):
        v=path[-1]
        if len(used)==len(V):
            if root in adj[v] and path[1]<path[-1]:
                out.add(frozenset(edge(a,b) for a,b in zip(path,path[1:]+[root])))
            return
        for w in adj[v]-used:visit(path+[w],used|{w})
    visit([root],{root});return out

def main():
    # K4 with edge0--3 subdivided three times: two cycles, four forced ports.
    base={edge(a,b) for a in range(4) for b in range(a+1,4)}-{(0,3)}
    ports=[(0,4),(4,5),(5,6),(3,6)];base|=set(ports)
    local=cycles(base);need(len(local)==2 and all(set(ports)<=c for c in local),'complete local source language')
    checks=0;total=0
    for shape in ['path','star']:
        wiring=[(0,0,1,0),(1,1,2,0),(2,1,3,0)] if shape=='path' else [(0,j,j+1,0) for j in range(3)]
        for twists in product((0,1),repeat=3):
            E={edge(7*i+a,7*i+b) for i in range(4) for a,b in base};spent={i:set() for i in range(4)}
            for (i,p,j,q),twist in zip(wiring,twists):
                a=[7*i+v for v in ports[p]];b=[7*j+v for v in ports[q]]
                E.remove(edge(*a));E.remove(edge(*b));spent[i].add(ports[p]);spent[j].add(ports[q])
                if twist:b.reverse()
                E|={edge(x,y) for x,y in zip(a,b)}
            actual=cycles(E);need(len(actual)==16,'literal tree product 2^4')
            recovered=set()
            for c in actual:
                factors=[]
                for i in range(4):
                    factor=frozenset({edge(a-7*i,b-7*i) for a,b in c if a//7==b//7==i}|spent[i])
                    need(factor in local,'local source recovered without a quotient');factors.append(factor)
                recovered.add(tuple(factors))
            need(len(recovered)==16,'literal product bijection')
            checks+=1;total+=len(actual)
    print(json.dumps({'status':'passed','local_cycles':2,'source_vertices':7,'tree_sources':4,
          'all_pairing_contexts':checks,'complete_output_cycles':total,'scope':'Independent exhaustive graph-cycle audit for path and star incidence, including ports sharing vertices. Arithmetic seams are audited separately.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
