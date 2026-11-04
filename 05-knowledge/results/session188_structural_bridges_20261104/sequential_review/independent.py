#!/usr/bin/env python3
from itertools import combinations_with_replacement,product
from collections import Counter
import json


def need(ok,why):
    if not ok:raise RuntimeError(why)


def child_sets(graph):
    c={v:set() for v in graph}
    for v,ps in graph.items():
        for p in set(ps):
            if p in graph:c[p].add(v)
    return c


def phi(graph):return sum(max(0,len(cs)-1) for cs in child_sets(graph).values())


def replace(graph,g,p):
    return {v:tuple(p if x==g else x for x in ps) for v,ps in graph.items() if v!=g}


def main():
    presentations=reroutes=constant_deletions=exact_transfers=0
    for n in range(1,5):
        options=[]
        for i in range(n):
            vs=range(i+2)
            options.append([()]+[(v,) for v in vs]+list(combinations_with_replacement(vs,2)))
        for ps in product(*options):
            graph={i+2:p for i,p in enumerate(ps)}
            degrees=child_sets(graph)
            before=phi(graph)
            presentations+=1
            for g,parents in graph.items():
                constant=replace(graph,g,-1)
                need(phi(constant)<=before,'constant deletion increased potential')
                constant_deletions+=1
                for p in set(parents):
                    after=replace(graph,g,p)
                    need(phi(after)<=before,'alias increased potential')
                    if p in graph and degrees[g]:
                        overlap=len((degrees[p]-{g})&degrees[g])
                        parent_loss=sum(q in graph and len(degrees[q])>=2 for q in set(parents)-{p})
                        need(phi(after)-before==-overlap-parent_loss,'exact transfer')
                        exact_transfers+=1
                    reroutes+=1
    ladder=[]
    for k in range(1,129):
        graph={'G0':('a','b')}
        for i in range(1,k+1):
            graph['A'+str(i)]=('a'+str(i),'b'+str(i))
            graph['G'+str(i)]=('G'+str(i-1),'A'+str(i))
        graph.update(L=('G'+str(k),'y'),R=('G'+str(k),'z'),O=('L','R'))
        need(phi(graph)==1 and len(graph)==2*k+4,'initial ladder')
        paid=0
        for i in range(k,0,-1):
            A,G,P='A'+str(i),'G'+str(i),'G'+str(i-1)
            c=child_sets(graph)
            need(c[A]=={G} and c[G]=={'L','R'},'current ladder incidences')
            paid+=len(c[G])-1
            graph=replace(graph,A,-1)
            graph=replace(graph,G,P)
            need(phi(graph)==1 and len(graph)==2*i+2,'alias ladder transfer')
        need(paid==k,'boundary excess sum')
        graph=replace(graph,'G0',-1)
        graph=replace(graph,'L','y')
        graph=replace(graph,'R','z')
        need(graph=={'O':('y','z')} and phi(graph)==0,'final consumed surplus')
        ladder.append(dict(k=k,boundary_excess=paid,constant_excess=1,total_deletions=2*k+3))
    return dict(status='PASS',at_most_two_pin_DAG_presentations=presentations,
                literal_reroutes=reroutes,constant_deletions=constant_deletions,
                exact_gate_transfer_checks=exact_transfers,ladder_max_k=128,
                last_ladder=ladder[-1])


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
