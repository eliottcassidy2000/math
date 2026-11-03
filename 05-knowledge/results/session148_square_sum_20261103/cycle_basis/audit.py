#!/usr/bin/env python3
"""Independent bitmask enumeration; does not import the source checker."""
from math import isqrt
from itertools import combinations
from collections import Counter

def need(x, reason):
    if not x: raise RuntimeError(reason)

def graph(n):
    return [(i,j) for i in range(n) for j in range(i+1,n) if isqrt(i+j+2)**2==i+j+2]

def hpaths(n, edges):
    adjacency=[0]*n
    for a,b in edges:
        adjacency[a]|=1<<b;adjacency[b]|=1<<a
    out=[]; full=(1<<n)-1
    def visit(v,seen,seq):
        if seen==full:
            if seq[0]<v: out.append(tuple(seq))
            return
        possibilities=adjacency[v]&~seen
        while possibilities:
            bit=possibilities&-possibilities;possibilities-=bit
            w=bit.bit_length()-1
            visit(w,seen|bit,seq+[w])
    # A Hamilton path must start at either leaf when a leaf exists.
    leaves=[v for v in range(n) if adjacency[v].bit_count()==1]
    if len(leaves)>2:return out
    starts=leaves if len(leaves)==2 else range(n)
    for v in starts:visit(v,1<<v,[v])
    return out

def orientations(n, edges, paths):
    index={e:i for i,e in enumerate(edges)}
    hist=Counter()
    for p in paths:
        for seq in (p,p[::-1]):
            fixed=bits=0
            for a,b in zip(seq,seq[1:]):
                e=(min(a,b),max(a,b));i=index[e];fixed|=1<<i
                if a<b:bits|=1<<i
            free=[i for i in range(len(edges)) if not (fixed>>i)&1]
            for v in range(1<<len(free)):
                orient=bits
                for j,i in enumerate(free):
                    if (v>>j)&1:orient|=1<<i
                hist[orient]+=1
    counts=Counter(hist.values())
    counts[0]=(1<<len(edges))-len(hist)
    return counts

def basis(n,edges,path):
    idx={e:i for i,e in enumerate(edges)}
    path_edges={tuple(sorted(pair)) for pair in zip(path,path[1:])}
    chords=set(edges)-path_edges
    beta=len(edges)-n+1
    need(len(path_edges)==n-1 and len(chords)==beta, ('beta',n))
    positions={v:i for i,v in enumerate(path)}
    vectors=[]
    for a,b in sorted(chords):
        lo,hi=sorted((positions[a],positions[b]));vec=1<<idx[(a,b)]
        for x,y in zip(path[lo:hi],path[lo+1:hi+1]):vec |= 1<<idx[tuple(sorted((x,y)))]
        chord_projection={e for e in chords if vec>>idx[e]&1}
        need(chord_projection=={(a,b)},('projection',n,a,b))
        parity=[0]*n
        for i,(x,y) in enumerate(edges):
            if vec>>i&1:parity[x]^=1;parity[y]^=1
        need(not any(parity),('not even degree',n,a,b))
        vectors.append(vec)
    span={0}
    for v in vectors:span|={x^v for x in span}
    need(len(span)==1<<beta,('rank',n))
    return span

expect={15:(15,1,1),16:(16,1,1),17:(17,1,1),23:(28,6,3),25:(32,8,10),26:(34,9,12),27:(36,10,35)}
hists={23:{0:268435072,1:384},25:{0:4294962432,1:4608,2:256}}
for n,(m,beta,u) in expect.items():
    edges=graph(n); paths=hpaths(n,edges)
    need((len(edges),len(edges)-n+1,len(paths))==(m,beta,u),('path count',n,len(edges),len(paths)))
    for p in paths:
        span=basis(n,edges,p)
        need(len(span)==2**beta,('span',n))
    counts=orientations(n,edges,paths)
    need(sum(k*v for k,v in counts.items())==u*2**(beta+1),('incidence',n,counts))
    if n in (15,16,17):
        p=paths[0]; pedges={tuple(sorted(x)) for x in zip(p,p[1:])}
        need({tuple(i+1 for i in e) for e in set(edges)-pedges}=={(1,3)},('unique chord',n))
        allcycles=span-{0}; cycle=next(iter(allcycles))
        cedges={tuple(v+1 for v in e) for i,e in enumerate(edges) if cycle>>i&1}
        need(cedges=={(1,3),(3,6),(6,10),(10,15),(1,15)},('G15-17 cycle',n,cedges))
        need(counts=={0:(1<<m)-4,1:4},('unique path orientation',n,counts))
    if n in hists:need(counts==hists[n],('histogram',n,counts))
    print('G%d: edges=%d beta=%d undirected_paths=%d histogram=%s'%(n,m,beta,u,dict(sorted(counts.items()))))

# Independent complete finite cycle-space enumeration for all simple graphs
# containing the specified Hamilton path on 2--5 vertices.
for n in range(2,6):
    all_edges=list(combinations(range(n),2)); fixed={(i,i+1) for i in range(n-1)}
    optional=[e for e in all_edges if e not in fixed]
    for mask in range(1<<len(optional)):
        edges=sorted(fixed|{e for i,e in enumerate(optional) if mask>>i&1})
        span=basis(n,edges,tuple(range(n)))
        eulerian=set()
        for v in range(1<<len(edges)):
            parity=0
            for i,(a,b) in enumerate(edges):
                if v>>i&1:parity^=(1<<a)|(1<<b)
            if not parity:eulerian.add(v)
        need(span==eulerian,('failed spanning',n,mask))
# Hostile overlap: all K3 orientations are successful but 12 incidences.
edges=list(combinations(range(3),2)); counts=orientations(3,edges,hpaths(3,edges))
need(counts=={0:0,1:6,3:2},('overlap example',counts))
need(6+2==8 and 6+3*2==12,'overlap counts')
print('PASS: independent DFS and orientation-extension census, cycle coordinates, all path-containing simple graphs <=5; K3 overlap hostile; explicit exceptions')
