#!/usr/bin/env python3
"""Davenport form of the THM-469 seam: (i) max m with FS(A) inside one v_p-level L_v (A an m-set of positive
integers) equals p-1; (ii) the same with FS(A) u FP(A) in one level forces v=0 and m<=p-1;
(iii) clique number of the 1-D Cayley graph {x,y}: v_p(|x-y|)=v equals p.  Brute force on small ranges."""
import itertools
def vp(n,p):
    c=0
    while n%p==0: n//=p; c+=1
    return c
def max_fs(p,v,Nmax,mmax):
    L=[x for x in range(1,Nmax+1) if vp(x,p)==v]
    best=1; ex=None
    # extend greedily-searchable by DFS over increasing elements
    def dfs(A,sums,start):
        nonlocal best,ex
        if len(A)>best: best=len(A); ex=list(A)
        if len(A)>=mmax: return
        for i in range(start,len(L)):
            a=L[i]
            new={s+a for s in sums}|{a}
            if all(vp(s,p)==v for s in new):
                dfs(A+[a],sums|new,i+1)
    dfs([],set(),0)
    return best,ex
def max_fsfp(p,Nmax,mmax):
    best=0; ex=None; lev=None
    for v in range(0,3):
        L=[x for x in range(2,Nmax+1) if vp(x,p)==v]
        def dfs(A,sums,prods,start):
            nonlocal best,ex,lev
            if len(A)>best: best=len(A); ex=list(A); lev=v
            if len(A)>=mmax: return
            for i in range(start,len(L)):
                a=L[i]
                ns={s+a for s in sums}|{a}; npr={q*a for q in prods}|{a}
                if all(vp(s,p)==v for s in ns) and all(vp(q,p)==v for q in npr):
                    dfs(A+[a],sums|ns,prods|npr,i+1)
        dfs([],set(),set(),0)
    return best,ex,lev
for p in [2,3,5,7]:
    for v in [0,1]:
        b,e=max_fs(p,v,Nmax=60*p**v,mmax=p+1)
        print(f"p={p} v={v}: max |A| with FS(A) in L_v = {b} (p-1={p-1}); example {e}")
    b,e,l=max_fsfp(p,Nmax=60,mmax=p+1)
    print(f"p={p}: max |A| (elements>=2) with FS(A)uFP(A) in one v_p level = {b}, level v={l}, example {e}")
    # clique number of Cayley graph on [0,40]
    Nn=40; V=range(Nn)
    best=1
    for v in [0,1]:
        # greedy exact: cliques correspond to sets with pairwise v_p(diff)==v
        import networkx as nx
        G=nx.Graph(); G.add_nodes_from(V)
        for x,y in itertools.combinations(V,2):
            if vp(y-x,p)==v: G.add_edge(x,y)
        w=max(len(c) for c in nx.find_cliques(G))
        print(f"   clique number of Cay([0,{Nn}), v_p(gap)={v}) = {w}")
