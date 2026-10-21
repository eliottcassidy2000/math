#!/usr/bin/env python3
"""Independent agreement-graph and pinned-path census; no repository imports."""
from itertools import permutations, combinations, product
from collections import Counter
from functools import lru_cache

def setup(n):
    pairs=list(combinations(range(n),2)); ix={p:k for k,p in enumerate(pairs)}
    def arc(t,i,j): return ((t>>ix[min(i,j),max(i,j)])&1) ^ (i>j)
    def rel(t,p): return sum(arc(t,p[i],p[j])<<k for k,(i,j) in enumerate(pairs))
    ps=list(permutations(range(n))); full=(1<<len(pairs))-1
    @lru_cache(None)
    def canon(t): return min(rel(t,p) for p in ps)
    def matrix(t): return [[0 if i==j else arc(t,i,j) for j in range(n)] for i in range(n)]
    def directed_path(edges):
        if len(edges)!=n-1:return None
        out={}; indeg=Counter()
        for a,b in edges:
            if a in out:return None
            out[a]=b;indeg[b]+=1
            if indeg[b]>1:return None
        starts=[v for v in range(n) if indeg[v]==0]
        if len(starts)!=1:return None
        p=[starts[0]]
        while p[-1] in out:
            v=out[p[-1]]
            if v in p:return None
            p.append(v)
        return tuple(p) if len(p)==n else None
    def blue(t,p):return all(arc(t,p[i],p[j])==arc(t,p[n-1-j],p[n-1-i]) for i,j in pairs)
    return pairs,ix,arc,rel,ps,full,canon,matrix,directed_path,blue

def census(n):
    pairs,ix,arc,rel,ps,full,canon,mat,path,blue=setup(n)
    pin=sum(1<<ix[i,i+1] for i in range(n-1)); free=[k for k,(i,j) in enumerate(pairs) if j>i+1]
    reps=sorted({canon(pin+sum(((bits>>j)&1)<<k for j,k in enumerate(free))) for bits in range(1<<len(free))})
    aut={a:sum(rel(a,p)==a for p in ps) for a in reps}
    counts=Counter()
    # independent method A: agreement graphs on arbitrary class representatives
    tested=0
    for a in reps:
      for b in reps:
        for f in ps:
            bf=rel(b,f); es=[(i,j) if arc(a,i,j) else (j,i) for i,j in pairs if arc(a,i,j)==arc(bf,i,j)]
            p=path(es);tested+=1
            if p:counts[a,b,blue(a,p)]+=1
        for c in (False,True):
            assert counts[a,b,c]%(aut[a]*aut[b])==0
            counts[a,b,c]//=aut[a]*aut[b]
    # independent B: all pinned standard Hamilton-path presentations
    pin=sum(1<<ix[i,i+1] for i in range(n-1)); free=[k for k,(i,j) in enumerate(pairs) if j>i+1]
    direct=Counter(); edges=Counter()
    for bits in range(1<<len(free)):
        t=pin+sum(((bits>>j)&1)<<k for j,k in enumerate(free)); kappa=t^(full^pin)
        a,b=canon(t),canon(kappa);c=blue(t,tuple(range(n)))
        direct[a,b,c]+=1
        assert canon(rel(t,tuple(reversed(range(n))))^full)==canon(t^full)
        if t<kappa:
            u=min(a,canon(a^full));v=min(b,canon(b^full))
            edges[min(u,v),max(u,v),c]+=1
    assert +counts==+direct
    fundamental={k:(v if k[2] else v//2) for k,v in edges.items()}
    assert all(c or v%2==0 for (u,w,c),v in edges.items())
    nodes={min(a,canon(a^full)) for a in reps}; deg=Counter()
    for (u,v,c),m in fundamental.items():
        if not c:deg[u]+=m;deg[v]+=m
    bridges=[]
    for (u,v,c),m in fundamental.items():
        if c or u==v or m!=1:continue
        seen={u};stack=[u]
        while stack:
            x=stack.pop()
            for (y,z,d),mult in fundamental.items():
                if d or (y,z)==(u,v):continue
                w=z if x==y else y if x==z else None
                if w is not None and w not in seen:seen.add(w);stack.append(w)
        if v not in seen:bridges.append((u,v,deg[u],deg[v]))
    print('CENSUS',n,'classes',len(reps),'merged',len(nodes),'bijections',tested,'lines',sum(edges.values()),'blue',sum(m for (u,v,c),m in edges.items() if c),'fundamental_black',sum(m for (u,v,c),m in fundamental.items() if not c),'raw_black_loops',sum(m for (u,v,c),m in edges.items() if not c and u==v),'bridges',bridges)
    return fundamental

def port_tests():
    tested=0
    # all 3-vertex bases, triangle or transitive two-vertex modules in combinations
    for k in (2,3):
      bp,bix,ba,*_=setup(k)
      for base in range(1<<len(bp)):
       for sizes in product((1,2,3),repeat=k):
        n=sum(sizes)
        if n>7:continue
        blocks=[];labels=[]
        for i,s in enumerate(sizes):
            block=list(range(len(labels),len(labels)+s));blocks.append(block);labels += [i]*s
        M=[[0]*n for _ in range(n)]
        for u,v in combinations(range(n),2):
            i,j=labels[u],labels[v]
            b=ba(base,i,j) if i!=j else (not (len(blocks[i])==3 and (u,v)==(blocks[i][0],blocks[i][-1])))
            M[u][v]=int(b);M[v][u]=1-int(b)
        def hpaths(p,unused):
            if not unused:yield p;return
            for v in unused:
                if not p or M[p[-1]][v]:yield from hpaths(p+[v],unused-{v})
        for p in hpaths([],set(range(n))):
            # direct nonpath complement
            arcs=set(zip(p,p[1:]));D=[[0]*n for _ in range(n)]
            for u,v in combinations(range(n),2):
                b=M[u][v] if (u,v) in arcs or (v,u) in arcs else 1-M[u][v]
                D[u][v]=b;D[v][u]=1-b
            # maximal runs, reversed blocks with only actual cover/boundary edges restored
            runs=[]
            for v in p:
                if not runs or labels[runs[-1][-1]]!=labels[v]:runs.append([v])
                else:runs[-1].append(v)
            internal={e for r in runs for e in zip(r,r[1:])}
            ports={(r[-1],s[0]) for r,s in zip(runs,runs[1:])}
            assert internal|ports==arcs
            R=[[0 if u==v else 1-M[u][v] for v in range(n)] for u in range(n)]
            for u,v in internal|ports:R[u][v]=1;R[v][u]=0
            assert R==D;tested+=1
    print('PORT_TRANSFER paths',tested,'all exact')

def hostile_and_groups():
    pairs,ix,arc,rel,ps,full,can,mat,path,blue=setup(4)
    a=sum(1<<ix[e] for e in [(0,1),(1,2),(0,3),(1,3),(2,3)])
    b=a ^ sum(1<<ix[i,3] for i in range(3))
    es=[(i,j) if arc(a,i,j) else (j,i) for i,j in pairs if arc(a,i,j)==arc(b,i,j)]
    assert len(es)==3 and max(Counter(u for u,v in es).values())==1 and path(es) is None
    print('HOSTILE degree-only agreement:',es,'cycle plus isolated vertex rejected')
    p3=list(permutations(range(3)));sgn=lambda p:(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
    frames=[(sgn(p),p) for p in p3]
    for g in p3:
        image={(sgn(g)*e,tuple(g[v] for v in p)) for e,p in frames}
        assert image==set(frames)
    assert len({(sgn(g),tuple(g[v] for v in (0,1,2))) for g in p3})==6
    # Two triangle blocks: all block-preserving signed permutations are compared directly.
    M=[[0]*6 for _ in range(6)]
    for u,v in combinations(range(6),2):
        z=1 if u//3!=v//3 else int((v-u)%3==1)
        M[u][v]=z;M[v][u]=1-z
    actual=pred=0
    for swap in (0,1):
      for g,h in product(p3,repeat=2):
        f=tuple(3*swap+g[i] for i in range(3))+tuple(3*(1-swap)+h[i] for i in range(3))
        signed=any(all(M[f[u]][f[v]]==(M[u][v] if e==1 else 1-M[u][v]) for u,v in combinations(range(6),2)) for e in (1,-1))
        expected=sgn(g)==sgn(h)==(-1)**swap
        assert signed==expected;actual+=signed;pred+=expected
    print('S3 signed frames regular; synchronized two-block signed group',actual,'of 72 block permutations')

if __name__=='__main__':
    for n in (3,4,5,6):census(n)
    hostile_and_groups();port_tests();print('ALL CHECKS PASSED')
