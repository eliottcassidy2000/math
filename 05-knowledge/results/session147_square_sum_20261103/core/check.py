#!/usr/bin/env python3
"""Finite exact certificate, Python 3 standard library; works with and without -O."""
from collections import deque, Counter
from fractions import Fraction
from itertools import combinations
from math import isqrt


def require(condition, message):
    if not condition: raise AssertionError(message)

def graph(vertices):
    v=tuple(vertices)
    if any(type(x) is not int or x <= 0 for x in v): raise ValueError('vertices must be positive integers, not booleans')
    if len(v)!=len(set(v)): raise ValueError('vertices must be distinct')
    return {tuple(sorted((x,y))) for x,y in combinations(v,2) if isqrt(x+y)**2 == x+y}

def components(vertices, edges):
    vertices=set(vertices)
    if any(len(e)!=2 or not set(e)<=vertices or e[0]>=e[1] for e in edges): raise ValueError('invalid edge')
    todo=set(vertices); out=[]
    while todo:
        x=todo.pop(); c={x}; stack=[x]
        while stack:
            x=stack.pop()
            for a,b in edges:
                y=b if a==x else a if b==x else None
                if y in todo: todo.remove(y); c.add(y); stack.append(y)
        out.append(frozenset(c))
    return frozenset(out)

def remove_vertices(V,E,D):
    W=V-set(D)
    return components(W,{e for e in E if set(e)<=W})

def edge_betweenness(V,E):
    # Enumerate shortest-path counts with breadth-first dynamic programming,
    # then sum sigma_s,t(edge)/sigma_s,t over unordered pairs.
    adj={v:set() for v in V}
    for a,b in E: adj[a].add(b); adj[b].add(a)
    result={e:Fraction(0) for e in E}
    for s,t in combinations(sorted(V),2):
        def bfs(root):
            d={root:0}; sig={root:1}; q=deque([root])
            while q:
                x=q.popleft()
                for y in adj[x]:
                    if y not in d:d[y]=d[x]+1;sig[y]=0;q.append(y)
                    if d[y]==d[x]+1:sig[y]+=sig[x]
            return d,sig
        ds,ns=bfs(s);dt,nt=bfs(t)
        if t not in ds: continue
        for a,b in E:
            for x,y in ((a,b),(b,a)):
                if x in ds and y in dt and ds[x]+1+dt[y]==ds[t]:
                    result[(a,b)]+=Fraction(ns[x]*nt[y],ns[t])
    return result

def main():
    A={10,6,3,1,8,17}; B={12,4,5,11}; C={2,7,9,18,16}
    S=A|B|C; V=set(range(1,19)); E=graph(V)
    ES={(1,3),(1,8),(2,7),(3,6),(4,5),(4,12),(5,11),(6,10),(7,9),(7,18),(8,17),(9,16)}
    require(len(S)==15 and graph(S)==ES, 'induced graph mismatch')
    require(components(S,graph(S))==frozenset(map(frozenset,[A,B,C])), 'three components mismatch')
    require(remove_vertices(V,E,{13,14,15})==frozenset(map(frozenset,[A,B,C])), 'vertex deletion mismatch')
    require(len(E)==18 and len(components(V,E))==1,'G18 not connected unicyclic')
    cycle={(1,3),(3,6),(6,10),(10,15),(1,15)}
    bridges={e for e in E if len(components(V,E-{e}))==2}
    require(len(bridges)==13 and E-bridges==cycle,'wrong bridges or cycle')
    require({x for x in V if sum(x in e for e in E)==1}=={16,17,18},'wrong leaves')
    scores=edge_betweenness(V,E)
    require(scores[(4,5)]==80 and scores[(12,13)]==80 and scores[(4,12)]==81,'wrong betweenness')
    require(max(scores.values())==81 and [e for e in scores if scores[e]==81]==[(4,12)],'wrong maximum')
    cuts=[(x,y) for x,y in combinations(sorted(bridges),2) if len(components(V,E-{x,y}))==3]
    require(len(cuts)==78,'wrong edge cut count')
    # Degree sum / induced-edge / surviving-cycle identity, for every triple.
    cyc_vertices={1,3,6,10,15}
    degrees={v:sum(v in e for e in E) for v in V}
    signature_counts=Counter()
    for D in combinations(sorted(V),3):
        D=set(D); beta=int(not D & cyc_vertices)
        s=sum(degrees[x] for x in D); internal=sum(set(e)<=D for e in E)
        number=len(remove_vertices(V,E,D))
        require(number==s-internal-3+beta, 'separator identity')
        if number==3:signature_counts[(s,internal,beta)]+=1
    require(signature_counts==Counter({(5,0,1):75,(6,0,0):162,(6,1,1):43,(7,1,0):63,(7,2,1):2,(8,2,0):4}), 'signature counts')
    equal=[D for D in combinations(sorted(V),3) if sorted(map(len,remove_vertices(V,E,D)))==[5,5,5]]
    require(equal==[(3,6,14),(3,14,17)], 'equal 5 partition mismatch')
    balanced_cuts=[cut for cut in cuts if sorted(map(len,components(V,E-set(cut))))==[5,6,7]]
    require(balanced_cuts==[((2,14),(3,13)),((3,13),(11,14))], 'balanced edge cuts')
    require(not any(sorted(map(len,components(V,E-set(cut))))==[6,6,6] for cut in cuts), 'equal edge cut exists')
    deletion_counts=Counter(tuple(sorted(map(len,remove_vertices(V,E,D)))) for D in combinations(sorted(V),3) if len(remove_vertices(V,E,D))==3)
    require(remove_vertices(V,E,{3,6,14})==frozenset(map(frozenset,[{1,8,10,15,17},{2,7,9,16,18},{4,5,11,12,13}])), 'five-five-five certificate')
    balanced=components(V,E-{(3,13),(11,14)})
    require(sorted(map(len,balanced))==[5,6,7], 'balanced edge-cut certificate')
    # Search smallest maximum vertex in any 3-cycle or 4-cycle.
    N=46; W=set(range(1,N+1)); F=graph(W); adj={x:set() for x in W}
    for a,b in F:adj[a].add(b);adj[b].add(a)
    triangles=[t for t in combinations(range(1,N+1),3) if all(e in F for e in combinations(t,2))]
    quads=set()
    for a,c in combinations(range(1,N+1),2):
        for b,d in combinations(sorted(adj[a]&adj[c]),2): quads.add(tuple(sorted((a,b,c,d))))
    require(min(map(max,triangles))==30 and (6,19,30) in triangles,'triangle threshold')
    require(min(map(max,quads))==46 and (1,3,35,46) in quads,'4-cycle threshold')
    for bad in ([1,1],[0,1],[-1,3],[True,3],[1,3.0]):
        try:graph(bad)
        except ValueError:pass
        else:raise AssertionError('invalid vertices accepted: '+repr(bad))
    for bad in ({(1,1)},{(1,4)},{(3,1)}):
        try:components({1,3},bad)
        except ValueError:pass
        else:raise AssertionError('invalid edges accepted: '+repr(bad))
    print('G18 edges:',sorted(E));print('G18 bridges:',sorted(bridges))
    print('betweenness:',sorted(scores.items(),key=lambda t:(-t[1],t[0])))
    print('3-vertex-deletions giving 3 components:',sum(deletion_counts.values()),'profiles:',sorted(deletion_counts.items()))
    print('2-bridge decompositions:',len(cuts));print('triangle threshold:',min(map(max,triangles)),'4-cycle threshold:',min(map(max,quads)))
    print('All certificates and hostile-input checks passed.')
if __name__=='__main__':main()
