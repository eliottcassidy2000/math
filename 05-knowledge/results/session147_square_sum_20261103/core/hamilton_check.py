#!/usr/bin/env python3
"""Exhaustive ordered Hamilton paths in G_n for n <= 18."""
from check import graph,require

def paths(n):
    v=set(range(1,n+1));E=graph(v)
    neighbors={x:set() for x in v}
    for a,b in E:neighbors[a].add(b);neighbors[b].add(a)
    def dfs(path,remaining):
        if not remaining:
            yield tuple(path);return
        for y in sorted(neighbors[path[-1]] & remaining):
            yield from dfs(path+[y],remaining-{y})
    return [p for s in sorted(v) for p in dfs([s],v-{s})]

def main():
    positive={1,15,16,17}
    for n in range(1,19):
        ps=paths(n);require(len(ps)==(1 if n==1 else 2 if n in positive else 0),'Hamilton count n='+str(n))
        print('G%d ordered Hamilton paths: %d'%(n,len(ps)),ps[:1])
    p=(8,1,15,10,6,3,13,12,4,5,11,14,2,7,9)
    require(p in paths(15) and (p+(16,)) in paths(16) and ((17,)+p+(16,)) in paths(17),'path witness')
    for n in (15,16,17):
        E=graph(range(1,n+1));ps=paths(n)
        require({(4,5),(1,3)}<=E and (2,4) not in E, 'edge presence')
        for path in ps:
            es={tuple(sorted(e)) for e in zip(path,path[1:])}
            require((4,5) in es and (1,3) not in es and (2,4) not in es, 'forced/forbidden/absent')
        print('G%d: 4–5 forced; 1–3 present but forbidden; 2–4 absent.'%n)
    print('All counts passed.')
if __name__=='__main__':main()
