#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
DOES H CLOSE RECONSTRUCTION? + realization-degeneracy metrics for the merged-metagraph blue/black multigraph.

kind-pasteur-2026-07-01-S14. Last session: (category, blue-deg, black-deg) does NOT determine the metagraph
for n>=5 (92 legal degree-preserving 2-swaps at n=5). H helps partially. Two honest tests of "H closes it":
 (i)  H-EDGE-STRUCTURE swaps: 2-swaps preserving degrees AND the multiset of endpoint-H-pairs per colour.
      If 0 while baseline>0, then (degrees + which-H-levels-connect) pins the graph -- H-STRUCTURE closes it.
 (ii) 1-WL COLOUR REFINEMENT seeded by an invariant I: iterate node-colour = (I, sorted multiset of
      (nbr-colour, edge-colour)). If the stable colouring is DISCRETE, I identifies every node ("closes it").
      Ladder I in {category, +degree, +H, +H+c3}.
DEGENERACY METRICS: baseline swap number sigma; H-refined swap number sigma_H; WL-identification level;
twin count (#node-pairs with identical (cat,deg,H)); WL stable-class count.
"""
import sys, itertools
from collections import defaultdict, Counter
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def build(n):
    VERTS=[n-i for i in range(n)]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); tileIdx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileIdx[(n-y+1,n-x+1)] for (x,y) in TILES]
    def gsym(bits): return all(TRANS[i]==i or bits[i]==bits[TRANS[i]] for i in range(m))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=VERTS.index(xL); yi=VERTS.index(yL)
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    P=list(itertools.permutations(range(n)))
    def canon(A):
        best=None
        for p in P:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    def H(A): return sum(1 for p in P if all(A[p[t]][p[t+1]] for t in range(n-1)))
    def c3(A): return sum(1 for i,j,k in itertools.combinations(range(n),3)
                          if (A[i][j]and A[j][k]and A[k][i]) or (A[i][k]and A[k][j]and A[j][i]))
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]; A=adj(bits)
        T.append(dict(mask=mask,canon=canon(A),g=gsym(bits),fl=mask^((1<<m)-1),H=H(A),c3=c3(A)))
    sigs=sorted(set(t['canon'] for t in T)); ci={s:i for i,s in enumerate(sigs)}
    for t in T: t['ci']=ci[t['canon']]
    bym={t['mask']:t for t in T}
    for t in T: t['tci']=bym[t['fl']]['ci'] if False else t['ci']  # placeholder
    # transpose target
    def tmask(t):
        tb=[0]*m; bits=[(t>>k)&1 for k in range(m)]
        for i in range(m): tb[TRANS[i]]=bits[i]
        return sum(b<<k for k,b in enumerate(tb))
    for t in T: t['tci']=bym[tmask(t['mask'])]['ci']
    tgt={c:[t for t in T if t['ci']==c][0]['tci'] for c in range(len(sigs))}
    par=list(range(len(sigs)))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for c in range(len(sigs)):
        a,b=find(c),find(tgt[c])
        if a!=b: par[max(a,b)]=min(a,b)
    noc=[find(c) for c in range(len(sigs))]
    Hcl={t['ci']:t['H'] for t in T}; c3cl={t['ci']:t['c3'] for t in T}
    return dict(n=n,m=m,T=T,noc=noc,tgt=tgt,Hcl=Hcl,c3cl=c3cl,bym=bym)

def metagraph(D):
    T=D['T']; noc=D['noc']; not_=lambda t: noc[t['ci']]
    nodes=sorted(set(noc))
    deg=defaultdict(lambda: dict(bl=0,bk=0,til=0,g=0))
    inc=defaultdict(list)  # node -> list of (nbr, color)  (nbr=node; self-loop nbr=node)
    seen=set()
    for t in T:
        nd=not_(t); deg[nd]['til']+=1; deg[nd]['g']+=(1 if t['g'] else 0)
        pr=frozenset((t['mask'],t['fl']))
        if pr in seen: continue
        seen.add(pr)
        a=not_(t); b=not_(D['bym'][t['fl']]); col=0 if t['g'] else 1  # 0=blue,1=black
        if a==b:
            deg[a]['bl' if col==0 else 'bk']+=2  # self-loop degree +2
            inc[a].append((a,col)); inc[a].append((a,col))
        else:
            deg[a]['bl' if col==0 else 'bk']+=1; deg[b]['bl' if col==0 else 'bk']+=1
            inc[a].append((b,col)); inc[b].append((a,col))
    def ntype(nd):
        g=deg[nd]['g']; tot=deg[nd]['til']
        return 0 if g==tot else (2 if g==0 else 1)  # 0=pure-blue,1=mixed,2=pure-black
    H={nd:D['Hcl'][[c for c in range(len(D['noc'])) if D['noc'][c]==nd][0]] for nd in nodes}
    c3={nd:D['c3cl'][[c for c in range(len(D['noc'])) if D['noc'][c]==nd][0]] for nd in nodes}
    return nodes,deg,inc,ntype,H,c3

def wl_refine(nodes,inc,init):
    color={nd:init[nd] for nd in nodes}
    for _ in range(len(nodes)+2):
        newc={}
        for nd in nodes:
            sig=(color[nd], tuple(sorted((color[x],c) for (x,c) in inc[nd])))
            newc[nd]=sig
        # renumber
        order={s:i for i,s in enumerate(sorted(set(newc.values()), key=lambda z:str(z)))}
        newc={nd:order[newc[nd]] for nd in nodes}
        if len(set(newc.values()))==len(set(color.values())):
            color=newc; break
        color=newc
    return color

def count_swaps(edges, allowed, keyfn=None):
    """count 2-swaps (a,b),(c,d)->(a,d),(c,b) legal under 'allowed'; if keyfn given, require it preserved
       (edge endpoint-key multiset unchanged: {key(a),key(b)},{key(c),key(d)} == {key(a),key(d)},{key(c),key(b)})."""
    cnt=0
    E=edges
    for i in range(len(E)):
        a,b=E[i]
        for j in range(i+1,len(E)):
            c,d=E[j]
            for (x,y),(z,w) in [((a,d),(c,b)),((a,c),(b,d))]:
                if x==y or z==w: continue
                if {frozenset((x,y)),frozenset((z,w))}=={frozenset((a,b)),frozenset((c,d))}: continue
                if not(allowed(x) and allowed(y) and allowed(z) and allowed(w)): continue
                if keyfn is not None:
                    before=Counter([frozenset((keyfn(a),keyfn(b))),frozenset((keyfn(c),keyfn(d)))])
                    after =Counter([frozenset((keyfn(x),keyfn(y))),frozenset((keyfn(z),keyfn(w)))])
                    if before!=after: continue
                cnt+=1; break
    return cnt

for n in [4,5,6]:
    D=build(n); nodes,deg,inc,ntype,H,c3=metagraph(D)
    tname=['PB','Mx','Bk']
    # edge lists per colour (distinct-node edges only; self-loops excluded from swaps)
    bluE=[]; blkE=[]
    seen=set()
    for nd in nodes:
        for (x,col) in inc[nd]:
            if x==nd: continue
            key=frozenset((nd,x))
            # keep each undirected edge once per its multiplicity: approximate by processing pairs
    # rebuild edge lists cleanly from inc
    from collections import Counter as C2
    pe=C2()
    for nd in nodes:
        for (x,col) in inc[nd]:
            if x==nd: continue
            pe[(min(nd,x),max(nd,x),col)]+=1
    for (a,b,col),mu in pe.items():
        mu//=2  # each edge counted from both endpoints
        for _ in range(mu):
            (bluE if col==0 else blkE).append((a,b))
    blkAllow=lambda nd: ntype(nd) in (1,2)
    bluAllow=lambda nd: ntype(nd) in (0,1)
    # baseline swaps
    sb=count_swaps(blkE,blkAllow); su=count_swaps(bluE,bluAllow)
    # H-edge-structure-preserving swaps
    sbH=count_swaps(blkE,blkAllow,keyfn=lambda nd:H[nd]); suH=count_swaps(bluE,bluAllow,keyfn=lambda nd:H[nd])
    # twin count under (cat,deg,H)
    lab=lambda nd:(ntype(nd),deg[nd]['bl'],deg[nd]['bk'],H[nd])
    labs=Counter(lab(nd) for nd in nodes)
    twins=sum(v-1 for v in labs.values() if v>1)
    twingroups=defaultdict(list)
    for nd in nodes:
        if labs[lab(nd)]>1: twingroups[lab(nd)].append(nd)
    if twingroups:
        print("  TWIN GROUPS (type,blue-deg,black-deg,H) -> nodes [+ tiling-count, c3]:")
        for L,nds in twingroups.items():
            det=[(nd,deg[nd]['til'],c3[nd]) for nd in nds]
            print(f"    ({tname[L[0]]},bl={L[1]},bk={L[2]},H={L[3]}) -> nodes {det}")
    # WL identification ladder
    def wl_classes(initfn):
        col=wl_refine(nodes,inc,{nd:initfn(nd) for nd in nodes})
        return len(set(col.values()))
    wl_cat  =wl_classes(lambda nd: ntype(nd))
    wl_deg  =wl_classes(lambda nd: (ntype(nd),deg[nd]['bl'],deg[nd]['bk']))
    wl_H    =wl_classes(lambda nd: (ntype(nd),deg[nd]['bl'],deg[nd]['bk'],H[nd]))
    wl_Hc3  =wl_classes(lambda nd: (ntype(nd),deg[nd]['bl'],deg[nd]['bk'],H[nd],c3[nd]))
    N=len(nodes)
    print(f"n={n}: {N} nodes | black {len(blkE)} edges, blue {len(bluE)} edges")
    print(f"  SWAP numbers: baseline (deg+colour) black={sb} blue={su};  H-edge-structure black={sbH} blue={suH}")
    print(f"    => H-structure closes reconstruction (swaps->0)? {sbH==0 and suH==0}  (baseline unique? {sb==0 and su==0})")
    print(f"  1-WL identification (#stable classes / {N} nodes; DISCRETE=closed):")
    print(f"    category={wl_cat}  +degree={wl_deg}  +H={wl_H}  +H+c3={wl_Hc3}   "
          f"[+H discrete? {wl_H==N}; +H+c3 discrete? {wl_Hc3==N}]")
    print(f"  twin nodes under (cat,deg,H): {twins}  (0 => (cat,deg,H) injective on nodes)")
    print()
print("VERDICT printed per n. Metrics: sigma (baseline swaps), sigma_H (H-structure swaps), WL-level, twins.")
print("DONE.")
