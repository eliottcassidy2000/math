#!/usr/bin/env python3
"""Center-induced square-sum graph verification on G_n, n=10,...,14,18."""
from collections import deque
from itertools import permutations
from check import graph, components, require

def centers(vertices,edges):
    result=set()
    for component in components(vertices,edges):
        ecc={}
        for source in component:
            dist={source:0};q=deque([source])
            while q:
                x=q.popleft()
                for a,b in edges:
                    y=b if a==x else a if b==x else None
                    if y is not None and y not in dist:dist[y]=dist[x]+1;q.append(y)
            require(len(dist)==len(component), 'disconnected component')
            ecc[source]=max(dist.values())
        result.update(x for x in component if ecc[x]==min(ecc.values()))
    return result

def isomorphisms(n):
    V=set(range(1,n+1)); E=graph(V);C=centers(V,E);F=graph(C)
    k=len(C);base=set(range(1,k+1));B=graph(base);out=[]
    for p in permutations(sorted(C)):
        f=dict(zip(sorted(base),p));image={tuple(sorted((f[a],f[b]))) for a,b in B}
        if image==F:out.append(f)
    return C,F,out

def main():
    expected={10:({3,4,5,7},{(4,5)},4),11:({3,5,7},set(),0),12:({3,4,5,7},{(4,5)},4),13:({7,12,13},{(12,13)},2),14:({4,5},{(4,5)},0),15:({4,5},{(4,5)},0),16:({5},set(),1),17:({4,5},{(4,5)},0),18:({4,5},{(4,5)},0),19:({4,5},{(4,5)},0),20:({4},set(),1),21:({4,5},{(4,5)},0),22:({14,22},{(14,22)},0)}
    for n,(C,F,count) in expected.items():
        c,f,isos=isomorphisms(n)
        require((c,f,len(isos))==(C,F,count),'bad center/isomorphism result n='+str(n))
        print('G'+str(n),'centers',sorted(c),'induced edges',sorted(f),'isomorphisms to G'+str(len(c)),len(isos))
    c,e,isos=isomorphisms(12)
    require({1:4,3:5,2:3,4:7} in isos,'stated iso absent')
    require(not any(list(m.values())==sorted(m.values()) for m in isos),'order preserving iso')
    require(centers(c,e)==c,'center graph not fixed')
    require(all(len(components(set(range(1,n+1)),graph(range(1,n+1))-{(4,5)}))==2 for n in range(14,22)), 'early bridge loss')
    require(len(components(set(range(1,23)),graph(range(1,23))-{(4,5)}))==1, 'n22 bridge not lost')
    # Distances from the center 14 in G22 give radius 5.
    E=graph(range(1,23));d={14:0};q=deque([14])
    while q:
        x=q.popleft()
        for a,b in E:
            y=b if a==x else a if b==x else None
            if y is not None and y not in d:d[y]=d[x]+1;q.append(y)
    require(len(d)==22 and max(d.values())==5, 'G22 radius')
    print('G12 center graph is fixed; all four isomorphisms are non-order-preserving.')
    print('4–5 remains a bridge for n=14..21 and first ceases at n=22; G22 radius 5.')
if __name__=='__main__':main()
