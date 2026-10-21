#!/usr/bin/env python3
"""Explicit odd-girth witnesses in the fully labelled complement graph."""
from itertools import combinations,permutations
from collections import deque

def path_mask(p,idx):return sum(1<<idx[tuple(sorted(e))] for e in zip(p,p[1:]))
def walecki(k):
    n=2*k
    return [tuple([i]+[z%n for j in range(1,k) for z in (i-j,i+j)]+[(i-k)%n]) for i in range(k)]
def build(ps,n):
    pairs=list(combinations(range(n),2));idx={e:i for i,e in enumerate(pairs)};full=(1<<len(pairs))-1
    start=0;assigned={}
    # Constraints on initial orientations from all requested paths and preceding flips.
    flips=0
    for p in ps:
        for u,v in zip(p,p[1:]):
            e=tuple(sorted((u,v)));j=idx[e];bit=int(u<v)^((flips>>j)&1)
            assert j not in assigned or assigned[j]==bit
            assigned[j]=bit
        flips^=full^path_mask(p,idx)
    assert flips==0
    for j,b in assigned.items():start|=b<<j
    states=[start];t=start
    for p in ps:
        assert all(((t>>idx[tuple(sorted((u,v)))])&1)==(u<v) for u,v in zip(p,p[1:]))
        t^=full^path_mask(p,idx);states.append(t)
    assert t==start and len(set(states[:-1]))==len(ps)
    return states

def tiny_bipartite(n):
    pairs=list(combinations(range(n),2));idx={e:i for i,e in enumerate(pairs)};full=(1<<len(pairs))-1
    ps=list(permutations(range(n)));color={};bip=True;edges=0
    for t in range(1<<len(pairs)):
      for p in ps:
        if all(((t>>idx[tuple(sorted((u,v)))])&1)==(u<v) for u,v in zip(p,p[1:])):
            s=t^(full^path_mask(p,idx));assert s!=t;edges+=1
            if ((full^path_mask(p,idx)).bit_count()%2): assert t.bit_count()%2!=s.bit_count()%2
    print('LABELLED parity check',n,'vertices',1<<len(pairs),'directed incidences',edges)

if __name__=='__main__':
    for k in (3,5,7,9,11):
        ps=walecki(k);edges=[tuple(sorted(e)) for p in ps for e in zip(p,p[1:])]
        assert len(edges)==len(set(edges))==k*(2*k-1)
        ss=build(ps,2*k)
        print('ODD_GIRTH exact',2*k,k,'Hamilton decomposition and directed cycle verified')
    # n=5 triangle: common arc 0->1, all remaining edges occur once.
    ps=[(2+i,0,1,2+(i+1)%3,2+(i+2)%3) for i in range(3)]
    ss=build(ps,5);print('ORDER5 odd-girth3 states',ss,'paths',ps)
    ps9=[(5,0,8,6,1,2,3,7,4),(4,2,3,6,0,1,5,7,8),(5,4,8,3,0,1,7,2,6),(6,5,2,8,1,3,4,0,7),(7,6,4,1,0,2,3,5,8)]
    ss9=build(ps9,9);print('ORDER9 odd-girth5 states',ss9,'paths',ps9)
    for n in (3,4):tiny_bipartite(n)
    print('ALL HOLONOMY CHECKS PASSED')
