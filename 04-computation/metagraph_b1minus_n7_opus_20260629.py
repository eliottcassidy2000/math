"""
Compute b1^-(n) = (b1 + SC - trC1 - 1)/2 for n=6,7 via the metagraph, with score+WL-pruned canon.
trC1 = #(SC-SC edges) - #(NS-pair direct edges). Also report the cylinder ring sizes (shared SC nbrs
of each NS pair) to test the polycylinder decomposition.
"""
from itertools import permutations, groupby
import sys
def refine(n, adj, radj):
    color=[bin(adj[v]).count("1") for v in range(n)]
    for _ in range(n):
        sig=[]
        for v in range(n):
            outc=tuple(sorted(color[w] for w in range(n) if (adj[v]>>w)&1))
            inc=tuple(sorted(color[w] for w in range(n) if (radj[v]>>w)&1))
            sig.append((color[v],outc,inc))
        # recolor
        order=sorted(range(n), key=lambda v: sig[v])
        newc={}; nc=0; prev=None
        col=[0]*n
        for v in order:
            if sig[v]!=prev: nc+=1; prev=sig[v]
            col[v]=nc
        if col==color: break
        color=col
    return color
def canon(n, adj):
    radj=[0]*n
    for v in range(n):
        for w in range(n):
            if (adj[v]>>w)&1: radj[w]|=1<<v
    color=refine(n,adj,radj)
    order=sorted(range(n), key=lambda v: color[v])
    blocks=[list(g) for k,g in groupby(order, key=lambda v: color[v])]
    def gen(bl):
        if not bl: yield []; return
        for perm in permutations(bl[0]):
            for rest in gen(bl[1:]): yield list(perm)+rest
    best=None
    for arr in gen(blocks):
        code=0
        for i in range(n):
            oi=arr[i]; ai=adj[oi]
            for j in range(i+1,n):
                if (ai>>arr[j])&1: code|=1<<(i*n+j)
        if best is None or code<best: best=code
    return best
def metagraph(n):
    base=[(k,k-1) for k in range(n-1,0,-1)]
    tiles=sorted([(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in base and (j,i) not in base])
    m=len(tiles)
    cl={}; tc=[]
    for bits in range(1<<m):
        adj=[0]*n
        for (i,j) in base: adj[i]|=1<<j
        for idx,(i,j) in enumerate(tiles):
            if (bits>>idx)&1: adj[i]|=1<<j
            else: adj[j]|=1<<i
        c=canon(n,adj)
        if c not in cl: cl[c]=len(cl)
        tc.append(cl[c])
    V=len(cl)
    edgeset=set()
    for bits in range(1<<m):
        a=tc[bits]
        for idx in range(m):
            b=tc[bits^(1<<idx)]
            if a!=b: edgeset.add((min(a,b),max(a,b)))
    # R = complement
    def compl(code):
        adj=[0]*n
        for i in range(n):
            for j in range(i+1,n):
                if (code>>(i*n+j))&1: adj[i]|=1<<j
                else: adj[j]|=1<<i
        cadj=[0]*n
        for i in range(n):
            for j in range(n):
                if (adj[i]>>j)&1: cadj[j]|=1<<i
        return canon(n,cadj)
    R={i:cl[compl(c)] for c,i in cl.items()}
    SC=sum(1 for i in range(V) if R[i]==i)
    E=len(edgeset)
    adjG={v:set() for v in range(V)}
    for (a,b) in edgeset: adjG[a].add(b); adjG[b].add(a)
    SCset={i for i in range(V) if R[i]==i}
    scsc=sum(1 for (a,b) in edgeset if a in SCset and b in SCset)
    nspair=sum(1 for (a,b) in edgeset if R[a]==b)
    trC1=scsc-nspair
    # components
    seen=set(); comps=0
    for v in range(V):
        if v in seen: continue
        comps+=1; st=[v]
        while st:
            x=st.pop()
            if x in seen: continue
            seen.add(x); st.extend(adjG[x]-seen)
    b1=E-V+comps
    b1m=(b1 + SC - trC1 - comps)//2
    # cylinder rings
    rings=[]; used=set()
    for i in range(V):
        if R[i]!=i and i not in used:
            j=R[i]; used.add(i); used.add(j)
            ring=adjG[i]&adjG[j]&SCset
            rings.append(len(ring))
    return V,E,comps,b1,SC,scsc,nspair,trC1,b1m,sorted(rings,reverse=True)
for n in [6,7]:
    V,E,comps,b1,SC,scsc,nspair,trC1,b1m,rings=metagraph(n)
    print(f"n={n}: V={V} E={E} comps={comps} b1={b1} SC={SC} SC-SC_edges={scsc} NS-pair_edges={nspair} trC1={trC1}")
    print(f"      => b1^- = {b1m}    apex-7 divides: {b1m%7==0} ({b1m}={'7*'+str(b1m//7) if b1m%7==0 else b1m})")
    print(f"      cylinder ring sizes (top 12): {rings[:12]}  (#cylinders={len(rings)}, sum(ring-1)={sum(r-1 for r in rings)})")
