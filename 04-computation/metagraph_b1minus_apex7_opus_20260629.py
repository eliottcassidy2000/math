"""
Recompute the metagraph R-odd Betti b1^- (the SECONDARY obstruction) directly, n=4,5,6.
Metagraph: vertices=iso-classes, edges=wiggly (flip one tile, base path n->...->1). R=complement.
b1 = E - V + components; b1^- = dim of (-1)-eigenspace of R acting on H_1 = (b1 - tr(R|H1))/2,
tr(R|H1) via Lefschetz = tr(R|C0) - tr(R|C1) (1-complex, with tr on H0 added back).
"""
from itertools import permutations, combinations
def canon(n, adjset):  # adjset: set of directed arcs (i,j); return canonical frozenset
    best=None
    for p in permutations(range(n)):
        t=frozenset((p[i],p[j]) for (i,j) in adjset)
        key=tuple(sorted(t))
        if best is None or key<best: best=key
    return best
def tournament_from_bits(n, tiles, bits):
    arcs=set()
    for k in range(n-1, 0, -1): arcs.add((k, k-1))  # base path k->k-1 (0-indexed vertices)
    for idx,(i,j) in enumerate(tiles):
        if (bits>>idx)&1: arcs.add((i,j))
        else: arcs.add((j,i))
    return arcs
def complement(arcs): return set((j,i) for (i,j) in arcs)
def build(n):
    base=set((k,k-1) for k in range(n-1,0,-1))
    allpairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    tiles=[(i,j) for (i,j) in [(a,b) for b in range(n) for a in range(n) if a>b] if (i,j) not in base and (j,i) not in base]
    tiles=sorted(set(tiles))
    m=len(tiles)
    # map each tiling bits -> canonical class
    cl={}; tilingclass=[]
    for bits in range(1<<m):
        arcs=tournament_from_bits(n,tiles,bits)
        c=canon(n,arcs)
        if c not in cl: cl[c]=len(cl)
        tilingclass.append(cl[c])
    V=len(cl)
    # wiggly edges between classes
    edges=set()
    for bits in range(1<<m):
        a=tilingclass[bits]
        for idx in range(m):
            b=tilingclass[bits^(1<<idx)]
            if a!=b: edges.add((min(a,b),max(a,b)))
    E=len(edges)
    # components
    import sys; sys.setrecursionlimit(10000)
    adj={v:set() for v in range(V)}
    for (a,b) in edges: adj[a].add(b); adj[b].add(a)
    seen=set(); comps=0
    for v in range(V):
        if v in seen: continue
        comps+=1; stack=[v]
        while stack:
            x=stack.pop()
            if x in seen: continue
            seen.add(x)
            stack.extend(adj[x]-seen)
    b1=E-V+comps
    # R = complement: class index -> complement class index
    invcl={i:c for c,i in cl.items()}
    Rmap={}
    for c,i in cl.items():
        comp=complement(set(c)); cc=canon(n,comp); Rmap[i]=cl[cc]
    SC=sum(1 for i in range(V) if Rmap[i]==i)
    # tr(R|C1): fixed edges, signed (+1 preserve orientation, -1 reversed)
    trC1=0
    for (a,b) in edges:
        ra,rb=Rmap[a],Rmap[b]
        if (min(ra,rb),max(ra,rb)) in edges:
            if {ra,rb}=={a,b}:
                trC1 += 1 if (ra,rb)==(a,b) else -1
    # tr(R|H0): # R-fixed components. components as sets:
    # assign comp id
    compid={}; seen=set(); cid=0
    for v in range(V):
        if v in seen: continue
        stack=[v]; members=[]
        while stack:
            x=stack.pop()
            if x in seen: continue
            seen.add(x); members.append(x); stack.extend(adj[x]-seen)
        for mem in members: compid[mem]=cid
        cid+=1
    trH0=sum(1 for c0 in range(cid) if all(compid.get(Rmap[v])==c0 for v in range(V) if compid[v]==c0) and any(compid[v]==c0 for v in range(V)))
    # simpler: a component is R-fixed if R maps its vertex set to itself
    compsets={}
    for v in range(V): compsets.setdefault(compid[v],set()).add(v)
    trH0=sum(1 for s in compsets.values() if set(Rmap[v] for v in s)==s)
    trH1 = trC1 - SC + trH0   # Lefschetz: tr H0 - tr H1 = tr C0 - tr C1 ; tr C0=SC
    # tr H0 - tr H1 = SC - trC1  => trH1 = trH0 - (SC - trC1) = trH0 - SC + trC1
    b1minus=(b1 - trH1)//2
    return V,E,comps,b1,SC,trC1,trH0,trH1,b1minus
for n in [6]:
    V,E,comps,b1,SC,trC1,trH0,trH1,b1m=build(n)
    print(f"n={n}: V={V} E={E} comps={comps} b1={b1}  SC={SC} trC1={trC1} trH0={trH0} trH1={trH1}  => b1^- = {b1m}")
