#!/usr/bin/env python3
"""Two independent tight-Hall tests and exhaustive unit-port reconfiguration."""
from collections import Counter, deque
from itertools import permutations
import json

CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(label)

def tight_union(E):
    n=len(E); tight=set()
    for mask in range(1,1<<n):
        slots={i for i in range(n) if mask>>i&1}
        ports=set().union(*(E[i] for i in slots))
        if len(ports)==len(slots): tight|=slots
    return frozenset(tight)

def inaccessible(E,M,R):
    reached_ports=set(R)-set(M); reached_slots=set()
    todo=deque(reached_ports)
    while todo:
        p=todo.popleft()
        for i,e in enumerate(E):
            if p in e and i not in reached_slots:
                reached_slots.add(i)
                if M[i] not in reached_ports:
                    reached_ports.add(M[i]);todo.append(M[i])
    return frozenset(set(range(len(E)))-reached_slots)

def universe(s,t):
    R=set(range(t)); tuples=list(permutations(range(t),s)); eligible=0; total_states=0; pairs=0
    components_hist=Counter(); maxdiam=0; strict=0
    for bits in range(1<<(s*t)):
        E=[{j for j in range(t) if bits>>(i*t+j)&1} for i in range(s)]
        states=[M for M in tuples if all(M[i] in E[i] for i in range(s))]
        if not states: continue
        eligible+=1;total_states+=len(states); T=tight_union(E)
        need(len(set().union(*(E[i] for i in T)))==len(T),'union of all tight sets is tight')
        strict+=not T
        graph={}
        for M in states:
            need(inaccessible(E,M,R)==T,'Hall subset enumeration equals independent vacancy BFS')
            neighbours=[]
            for i in range(s):
                for p in E[i]-set(M):
                    N=list(M);N[i]=p;neighbours.append(tuple(N))
            graph[M]=neighbours
        labels={tuple(M[i] for i in sorted(T)) for M in states}
        components_hist[len(labels)]+=1
        for M in states:
            distance={M:0};q=deque([M])
            while q:
                U=q.popleft()
                for V in graph[U]:
                    if V not in distance:
                        distance[V]=distance[U]+1;q.append(V)
            expected={N for N in states if all(N[i]==M[i] for i in T)}
            need(set(distance)==expected,'complete component classification by actual tight restriction')
            need(max(distance.values())<=2*s*s+2*s,'constructive polynomial diameter upper bound')
            maxdiam=max(maxdiam,max(distance.values()));pairs+=len(distance)
    return dict(slots=s,ports=t,graphs=1<<(s*t),feasible_graphs=eligible,
                actual_assignments=total_states,reachable_ordered_pairs=pairs,
                strict_Hall_graphs=strict,component_histogram=dict(components_hist),maxdiam=maxdiam)

def main():
    cases=[universe(s,t) for s,t in ((1,2),(2,2),(2,3),(3,3),(3,4))]
    E=[{0,1},{1,2},{2,0}];M=(0,1,2)
    need(inaccessible(E,M,set(range(4)))==frozenset(range(3)), 'unconnected spare port repairs nothing')
    return dict(status='FINITE-EXACT PASS; classical matching mechanism reproved in port language',cases=cases,checks=CHECKS)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
