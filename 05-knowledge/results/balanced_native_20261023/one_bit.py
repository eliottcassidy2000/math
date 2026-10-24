#!/usr/bin/env python3
"""Constructive two-layer completion for arbitrary connected local port graphs."""
from collections import deque
from itertools import product, permutations, combinations_with_replacement, combinations
import json,random

CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def shortest(G,a,b):
    todo=deque([a]);parent={a:None}
    while todo:
        x=todo.popleft()
        if x==b:break
        for y in sorted(G[x]):
            if y not in parent:parent[y]=x;todo.append(y)
    if b not in parent:raise ValueError('local endpoints disconnected')
    p=[b]
    while p[-1]!=a:p.append(parent[p[-1]])
    return p[::-1]

def complete(graphs,a,b):
    s=len(a);state=[(x,0) for x in a];history=[tuple(state)];events=[]
    need(len(set(a))==s and len(set(b))==s,'actual endpoint injections')
    routes=[shortest(G,x,y) for G,x,y in zip(graphs,a,b)]
    lengths=sum(len(p)-1 for p in routes);done=set()
    def move(i,dest,kind):
        old=state[i]
        need(dest not in state,'actual enlarged port is vacant')
        need((dest[0]==old[0] and dest[1]==1-old[1]) or
             (dest[1]==old[1] and dest[0] in graphs[i][old[0]]),'legal one-slot native edge')
        state[i]=dest;history.append(tuple(state));events.append((i,kind))
        need(len(set(state))==s,'injection at each integer clock')
    while len(done)<s:
        i=min(set(range(s))-done)
        while True:
            need(state[i][0]==a[i], 'unprocessed slot retains original base port')
            for v in routes[i][1:]:
                bit=state[i][1]
                blockers=[j for j,(x,l) in enumerate(state) if x==v]
                need(len(blockers)<=1,'only the travelling slot can create a double base port')
                if blockers and state[blockers[0]][1]==bit:
                    j=blockers[0];move(j,(v,1-bit),'resident layer change')
                move(i,(v,bit),'horizontal')
            done.add(i)
            need(all(state[j][0]==b[j] for j in done),'completed base targets remain fixed')
            residents=[j for j in range(s) if j!=i and state[j][0]==b[i]]
            if not residents:break
            need(len(residents)==1 and residents[0] not in done,'unique unprocessed displaced resident')
            i=residents[0]
    for i in range(s):
        if state[i][1]:move(i,(b[i],0),'final layer reset')
    need(tuple(state)==tuple((x,0) for x in b),'full marked endpoint restored')
    need(sum(k=='horizontal' for i,k in events)==lengths,'each supplied local route used exactly once')
    need(len(events)<=2*lengths+s,'uniform one-bit integer cost')
    return history,events,lengths

def assignments(Es):
    out=[]
    def visit(i,x,used):
        if i==len(Es):out.append(tuple(x));return
        for v in sorted(Es[i]-used):visit(i+1,x+[v],used|{v})
    visit(0,[],set());return out

def faces(r):
    out=[]
    for pattern in product(range(3),repeat=r):
        E={x for x in range(1<<r) if all(v==2 or (x>>j)&1==v for j,v in enumerate(pattern))}
        out.append({x:{x^(1<<j) for j in range(r) if x^(1<<j) in E} for x in E})
    return out

def main():
    cases=0;pairs=0;events=0;maxlen=0
    # All labelled tuples, no slot quotient, for r=2 and up to three slots.
    F=faces(2)
    for s in range(1,4):
        for G in product(F,repeat=s):
            A=assignments([set(g) for g in G])
            if not A:continue
            cases+=1
            for a,b in product(A,repeat=2):
                h,e,L=complete(G,a,b);pairs+=1;events+=len(e);maxlen=max(maxlen,len(e))
    # Complete all local connected graphs on every nonempty subset of four ports.
    graphs=[]
    for mask in range(1,16):
        E=[v for v in range(4) if mask>>v&1];edges=list(combinations(E,2))
        for bits in range(1<<len(edges)):
            G={v:set() for v in E}
            for k,(u,v) in enumerate(edges):
                if bits>>k&1:G[u].add(v);G[v].add(u)
            try:
                for v in E:shortest(G,E[0],v)
            except ValueError:continue
            graphs.append(G)
    general_cases=0;general_pairs=0
    for G in combinations_with_replacement(graphs,2):
        A=assignments([set(g) for g in G])
        if not A:continue
        general_cases+=1
        for a,b in product(A,repeat=2):
            complete(G,a,b);general_pairs+=1
    rng=random.Random(20261023);random_cases=0;full=[]
    for r in range(1,9):
        N=1<<r;cube={x:{x^(1<<j) for j in range(r)} for x in range(N)}
        a=tuple(range(N));b=tuple(x^(N-1) for x in a)
        h,e,L=complete([cube]*N,a,b)
        full.append(dict(r=r,slots=N,horizontal=L,events=len(e),bound=2*L+N))
        for _ in range(16):
            s=min(12,N);a=rng.sample(range(N),s);b=rng.sample(range(N),s);G=[]
            for u,v in zip(a,b):
                fixed=(N-1)^(u^v);E={x for x in range(N) if x&fixed==u&fixed}
                G.append({x:{x^(1<<j) for j in range(r) if x^(1<<j) in E} for x in E})
            complete(G,a,b);random_cases+=1
    return dict(status='FINITE-EXACT PASS; general one-bit theorem proved separately',
                face_menu_tuples=cases,face_endpoint_pairs=pairs,face_events=events,largest_face_return=maxlen,
                connected_local_graphs= len(graphs),general_two_slot_menus=general_cases,
                general_endpoint_pairs=general_pairs,seeded_cases=random_cases,full_occupancy=full,checks=CHECKS)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
