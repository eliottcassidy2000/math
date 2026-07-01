#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE MERGED-METAGRAPH TILING-COUNT PAIRING PROCESS: blue/black lines as a degree-constrained assignment.

kind-pasteur-2026-07-01-S12. Replicates tournament-tiling-explorer.html EXACTLY and grounds the owner's
description.  A LINE = an unordered pair {t, flip(t)} (flip = complement TILING = flip all m tiles, the d=m
waggly layer).  It is BLUE iff isGridSym(t) (grid symmetry (x,y)->(n-y+1,n-x+1)), else BLACK.  MERGE = the
transpose involution transMask (grid-transpose of a tiling); SC node = transpose-self, NS-merged node = a
transpose pair.  Each tiling is EXACTLY ONE line-endpoint => tiling count of a node = its DEGREE in the line
multigraph (a self-loop, flip(t) in the same node, counts as +2 to one node; else +1 to each of two nodes).
Node types: PURE-BLUE (all tilings grid-sym), PURE-BLACK (none), MIXED (both).  We tabulate per-node
(blue/black)x(self/other) degrees and their PARITIES to test the owner's 3-category claim and the self-loop
conjecture (self-loops only on mixed nodes).
"""
import sys, itertools
from collections import defaultdict
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def build(n):
    VERTS=[n-i for i in range(n)]                       # [n, n-1, ..., 1]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); tileIdx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileIdx[(n-y+1,n-x+1)] for (x,y) in TILES]   # grid-symmetry permutation of tiles
    def isGridSym(bits):
        return all(TRANS[i]==i or bits[i]==bits[TRANS[i]] for i in range(m))
    def bitsToAdj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=VERTS.index(xL); yi=VERTS.index(yL)
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    perms=list(itertools.permutations(range(n)))
    def canon(A):
        best=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        flipMask=mask ^ ((1<<m)-1)
        transBits=[0]*m
        for i in range(m): transBits[TRANS[i]]=bits[i]
        transMask=sum(b<<k for k,b in enumerate(transBits))
        T.append(dict(mask=mask,bits=bits,canon=canon(bitsToAdj(bits)),
                      gridSym=isGridSym(bits),flipMask=flipMask,transMask=transMask))
    # classes
    sigs=sorted(set(t['canon'] for t in T)); cidx={s:i for i,s in enumerate(sigs)}
    for t in T: t['ci']=cidx[t['canon']]
    byMask={t['mask']:t for t in T}
    for t in T: t['transCi']=byMask[t['transMask']]['ci']
    # transpose target per class
    transTgt={}
    for ci in range(len(sigs)):
        reps=[t for t in T if t['ci']==ci]
        transTgt[ci]=reps[0]['transCi']
    # merged nodes (union ci with transTgt[ci])
    parent=list(range(len(sigs)))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for ci in range(len(sigs)):
        a,b=find(ci),find(transTgt[ci])
        if a!=b: parent[max(a,b)]=min(a,b)
    node_of_class=[find(ci) for ci in range(len(sigs))]
    nodes=sorted(set(node_of_class))
    isSC={nd: (transTgt[nd]==nd) for nd in nodes}     # transpose-self
    node_of_tiling={t['mask']: node_of_class[t['ci']] for t in T}
    return dict(n=n,m=m,T=T,nsig=len(sigs),transTgt=transTgt,node_of_class=node_of_class,
                nodes=nodes,isSC=isSC,node_of_tiling=node_of_tiling)

def analyze(D):
    n=D['n']; T=D['T']; nodes=D['nodes']; nom=D['node_of_tiling']; isSC=D['isSC']
    # per-node degrees split blue/black x self/other
    deg=defaultdict(lambda: dict(blueSelf=0,blackSelf=0,blueOther=0,blackOther=0,tilings=0,gsym=0))
    processed=set()
    nlines_blue=nlines_black=0
    for t in T:
        nd=nom[t['mask']]; deg[nd]['tilings']+=1
        if t['gridSym']: deg[nd]['gsym']+=1
        # each line once (mask<flipMask); attribute endpoints
        pair=frozenset((t['mask'],t['flipMask']))
        if pair in processed: continue
        processed.add(pair)
        a=nom[t['mask']]; b=nom[t['flipMask']]; blue=t['gridSym']
        if blue: nlines_blue+=1
        else: nlines_black+=1
        if a==b:
            deg[a]['blueSelf' if blue else 'blackSelf']+=1   # self-loop: +2 tilings already counted
        else:
            deg[a]['blueOther' if blue else 'blackOther']+=1
            deg[b]['blueOther' if blue else 'blackOther']+=1
    # node type
    def ntype(nd):
        g=deg[nd]['gsym']; tot=deg[nd]['tilings']
        return 'pure-blue' if g==tot else ('pure-black' if g==0 else 'mixed')
    print("="*104)
    print(f" n={n}: m={D['m']} tiles, 2^m={1<<D['m']} tilings, {D['nsig']} iso classes, {len(nodes)} merged nodes; "
          f"lines: blue={nlines_blue} black={nlines_black} total={nlines_blue+nlines_black}=2^(m-1)")
    print("="*104)
    print(f"  {'node':>5} {'SC?':>4} {'type':>10} {'#til':>5} | {'bkOth':>5} {'blOth':>5} {'bkSelf':>6} {'blSelf':>6} "
          f"| {'par(bkOth)':>10} {'par(blOth)':>10} {'selfloop?':>9}")
    counts=defaultdict(int)
    for nd in nodes:
        d=deg[nd]; ty=ntype(nd); counts[(('SC' if isSC[nd] else 'NS'),ty)]+=1
        par=lambda x: 'even' if x%2==0 else 'ODD'
        sl = (d['blueSelf']+d['blackSelf'])>0
        # verify tiling count = degree (self-loops count x2)
        degsum=d['blueOther']+d['blackOther']+2*(d['blueSelf']+d['blackSelf'])
        chk='' if degsum==d['tilings'] else f' !!DEG {degsum}!={d["tilings"]}'
        print(f"  {nd:>5} {('SC' if isSC[nd] else 'NS'):>4} {ty:>10} {d['tilings']:>5} | "
              f"{d['blackOther']:>5} {d['blueOther']:>5} {d['blackSelf']:>6} {d['blueSelf']:>6} | "
              f"{par(d['blackOther']):>10} {par(d['blueOther']):>10} {str(sl):>9}{chk}")
    print("  node census (SC/NS x type):")
    for k in sorted(counts): print(f"    {k[0]:>2} {k[1]:>10}: {counts[k]}")
    # self-loop conjecture
    slnodes=[nd for nd in nodes if (deg[nd]['blueSelf']+deg[nd]['blackSelf'])>0]
    sl_types=set(ntype(nd) for nd in slnodes)
    print(f"  SELF-LOOP nodes: {len(slnodes)}; their types = {sl_types or '{}'}  "
          f"=> owner's 'self-loops only on MIXED' : {sl_types<= {'mixed'} if slnodes else 'n/a (none)'}")
    # ---- TYPE x TYPE incidence for blue & black lines (who connects to whom) ----
    inc={'blue':defaultdict(int),'black':defaultdict(int)}
    sl={'blue':defaultdict(int),'black':defaultdict(int)}
    seen2=set()
    for t in T:
        pair=frozenset((t['mask'],t['flipMask']))
        if pair in seen2: continue
        seen2.add(pair)
        a=nom[t['mask']]; b=nom[t['flipMask']]; col='blue' if t['gridSym'] else 'black'
        ta,tb=ntype(a),ntype(b)
        if a==b: sl[col][ta]+=1
        else: inc[col][tuple(sorted((ta,tb)))]+=1
    print("  BLUE line type-pairs (between distinct nodes):  ", dict(inc['blue']))
    print("  BLACK line type-pairs (between distinct nodes): ", dict(inc['black']))
    print("  BLUE self-loops by type:  ", dict(sl['blue']), "   BLACK self-loops by type: ", dict(sl['black']))
    return deg,ntype

for n in [4,5,6]:
    D=build(n); analyze(D)
    print()
print("KEY IDENTITY: tiling count of a merged node = its DEGREE in the blue/black line multigraph")
print("  (each tiling = one line-endpoint; self-loop {t,flip(t)} in one node contributes +2).")
print("DONE.")
