from collections import deque
from itertools import combinations
from math import isqrt

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def graph(n):
    adj = {v:set() for v in range(1,n+1)}
    for x,y in combinations(adj,2):
        q=isqrt(x+y)
        if q*q==x+y:
            adj[x].add(y); adj[y].add(x)
    return adj

def components(adj, removed_vertices=(), removed_edges=()):
    forbidden={frozenset(e) for e in removed_edges}
    left=set(adj)-set(removed_vertices); result=[]
    while left:
        s=min(left); left.remove(s); q=[s]; comp={s}
        for u in q:
            for v in adj[u]&left:
                if frozenset((u,v)) not in forbidden:
                    left.remove(v); comp.add(v); q.append(v)
        result.append(tuple(sorted(comp)))
    return sorted(result)

def edges(adj):
    return [(a,b) for a in adj for b in sorted(adj[a]) if a<b]

def triples(adj):
    return [t for t in combinations(adj,3) if all(y in adj[x] for x,y in combinations(t,2))]

def path_without(adj, start, end, forbidden):
    q=deque([start]); prev={start:None}
    while q:
        u=q.popleft()
        for v in sorted(adj[u]):
            if frozenset((u,v))==frozenset(forbidden): continue
            if v not in prev: prev[v]=u; q.append(v)
    if end not in prev: return None
    out=[]; x=end
    while x is not None: out.append(x); x=prev[x]
    return out[::-1]

def edge_betweenness_unnormalized(adj):
    # count shortest-path fractions for unordered pairs using BFS distances
    E=edges(adj); totals={e:0.0 for e in E}
    for source in adj:
        dist={source:0}; sigma={source:1}; q=deque([source]); stack=[]
        pred={v:[] for v in adj}
        while q:
            v=q.popleft(); stack.append(v)
            for w in adj[v]:
                if w not in dist: dist[w]=dist[v]+1; sigma[w]=0; q.append(w)
                if dist[w]==dist[v]+1:
                    sigma[w]+=sigma[v]; pred[w].append(v)
        dep={v:0.0 for v in adj}
        while stack:
            w=stack.pop()
            for v in pred[w]:
                contribution=(sigma[v]/sigma[w])*(1+dep[w]); totals[tuple(sorted((v,w)))]+=contribution; dep[v]+=contribution
    return {e:v/2 for e,v in totals.items()}

G=graph(18); E=edges(G)
A={10,6,3,1,8,17}; B={12,4,5,11}; C={2,7,9,18,16}
expectedA={(1,3),(1,8),(3,6),(6,10),(8,17)}
expectedB={(4,5),(4,12),(5,11)}
expectedC={(2,7),(7,9),(7,18),(9,16)}
for S, exp in ((A,expectedA),(B,expectedB),(C,expectedC)):
    check({e for e in E if set(e)<=S}==exp,'wrong internal edges')
check(not any(any(x in S and y in T or y in S and x in T for S,T in combinations((A,B,C),2)) for x,y in E),'cross edges')
check(sorted(map(len,components(graph(12))))==[3,4,5],'G12 components')
check(len(E)==18 and len(components(G))==1,'G18 connectivity or edge count')
check([len(components(graph(n),removed_edges=[(4,5)])) for n in range(15,22)]==[2]*7,'4-5 bridge G15..21')
path=path_without(graph(22),4,5,(4,5))
check(path==[4,12,13,3,22,14,11,5],'G22 alternative path')
check(all(isqrt(x+y)**2==x+y for x,y in zip(path,path[1:])),'path edges')
check(len(triples(graph(29)))==0 and triples(graph(30))==[(6,19,30)],'triangle threshold')
squares=[4,9,16,25,36,49]
sum_triples=[]
for a,b,c in combinations(squares,3):
    if (a+b+c)%2==0 and a+b>c:
        sum_triples.append(((a,b,c),((a+b-c)//2,(a+c-b)//2,(b+c-a)//2)))
check(sum_triples==[((25,36,49),(6,19,30))],'triangle sum inverse')
vertex_cuts=[]
for remove in combinations(G,3):
    comps=components(G,removed_vertices=remove)
    if len(comps)==3 and sorted(map(len,comps))==[5,5,5]:vertex_cuts.append(remove)
check(vertex_cuts==[(3,6,14),(3,14,17)],'5/5/5 cuts')
edge_cuts={}
for removed in combinations(E,2):
    comps=components(G,removed_edges=removed)
    if len(comps)==3:edge_cuts[removed]=tuple(sorted(map(len,comps)))
check((6,6,6) not in edge_cuts.values(),'6/6/6 edge cut')
best=[k for k,v in edge_cuts.items() if v==(5,6,7)]
check(best==[((2,14),(3,13)),((3,13),(11,14))],'5/6/7 cuts')
bet=edge_betweenness_unnormalized(G)
check(abs(bet[(4,12)]-81)<1e-9 and abs(bet[(4,5)]-80)<1e-9 and max(bet.values())==81,'betweenness')
check([e for e,v in bet.items() if v==81]==[(4,12)],'betweenness uniqueness')
check(len(path_without(graph(14),2,4,(-1,-1)))==5,'G14 path 2 to 4')
check(path_without(graph(13),2,4,(-1,-1)) is None,'G13 disconnection')
check(len(path_without(graph(76),2,4,(-1,-1)))==5,'G76 distance 4')
check(path_without(graph(77),2,4,(-1,-1))==[2,23,77,4],'G77 distance 3')
sized_cuts=[r for r in combinations(G,3) if sorted(map(len,components(G,removed_vertices=r)))==[4,5,6]]
check(len(sized_cuts)==20,'G18 count 4/5/6 vertex cuts')
print('verified G12 and G18 components, G15–22 bridge threshold, G30 first triangle, square-sum inverse, 5/5/5 vertex cuts, edge cuts and shortest-path edge betweenness')
print('G18 equal vertex cuts:', vertex_cuts)
print('G18 5/6/7 edge cuts:', best)
print('G22 alternative path:', path)
print('G30 triangles:', triples(graph(30)))
print('G18 4-12 and 4-5 betweenness:', bet[(4,12)], bet[(4,5)])
