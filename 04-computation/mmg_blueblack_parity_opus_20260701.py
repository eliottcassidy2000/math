"""Merged-metagraph blue/black flip-line structure (replicates tournament-tiling-explorer.html).
Verify the owner's parity structure among pure-blue / pure-black / mixed nodes."""
import numpy as np
from itertools import permutations
def build(n):
    VERTS=[n-i for i in range(n)]                    # [n,n-1,...,1]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]  # x>=y+2
    m=len(TILES); tileidx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileidx[(n-y+1,n-x+1)] for (x,y) in TILES]
    perms=list(permutations(range(n)))
    def bits_adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=n-xL; yi=n-yL
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    def canon(A):
        best=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    def gridsym(bits): return all(bits[i]==bits[TRANS[i]] for i in range(m) if TRANS[i]!=i)
    tilings=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        A=bits_adj(bits); c=canon(A)
        Aop=[[A[j][i] for j in range(n)] for i in range(n)]; cop=canon(Aop)
        tilings.append(dict(mask=mask,bits=bits,canon=c,cop=cop,gs=gridsym(bits),flip=mask^((1<<m)-1)))
    # class ids
    classes=sorted(set(t['canon'] for t in tilings)); cid={c:i for i,c in enumerate(classes)}
    for t in tilings: t['ci']=cid[t['canon']]; t['copi']=cid[t['cop']]
    # merged node id = min(ci, transpose ci)
    merged={}; 
    for c in classes:
        ci=cid[c]; A=None
        # find cop id: pick a tiling with this canon
    # compute transpose target per class
    tpt={}
    for t in tilings:
        tpt[t['ci']]=t['copi']
    def mnode(ci): return min(ci,tpt[ci])
    for t in tilings: t['mn']=mnode(t['ci'])
    return n,m,TILES,tilings,cid,tpt
for n in [4,5,6]:
    n,m,TILES,tilings,cid,tpt=build(n)
    nodes=sorted(set(t['mn'] for t in tilings))
    # per merged node: tiling count, gridsym counts (category), self/edge blue/black
    tau={v:0 for v in nodes}; gs_cnt={v:0 for v in nodes}
    for t in tilings:
        tau[t['mn']]+=1; gs_cnt[t['mn']]+= (1 if t['gs'] else 0)
    # category
    cat={}
    for v in nodes:
        if gs_cnt[v]==0: cat[v]='pureBLACK'
        elif gs_cnt[v]==tau[v]: cat[v]='pureBLUE'
        else: cat[v]='mixed'
    # lines (flip pairs, deduped)
    bymask={t['mask']:t for t in tilings}
    selfB={v:0 for v in nodes}; selfK={v:0 for v in nodes}; edgeB={v:0 for v in nodes}; edgeK={v:0 for v in nodes}
    seen=set(); L=0; blines=0; klines=0
    for t in tilings:
        fm=t['flip']
        key=(min(t['mask'],fm),max(t['mask'],fm))
        if key in seen: continue
        seen.add(key); L+=1
        u=t; w=bymask[fm]; col='blue' if t['gs'] else 'black'
        if col=='blue': blines+=1
        else: klines+=1
        a,b=u['mn'],w['mn']
        if a==b:
            if col=='blue': selfB[a]+=1
            else: selfK[a]+=1
        else:
            if col=='blue': edgeB[a]+=1; edgeB[b]+=1
            else: edgeK[a]+=1; edgeK[b]+=1
    from collections import Counter
    catcnt=Counter(cat.values())
    print(f"\n===== n={n}: m={m}, tilings=2^{m}={2**m}, merged nodes={len(nodes)}, lines=2^{m-1}={L} (blue={blines},black={klines}) =====")
    print(f"  categories: {dict(catcnt)}")
    # verify pairing identity tau = 2*self + edges
    ok=all(tau[v]==2*(selfB[v]+selfK[v])+(edgeB[v]+edgeK[v]) for v in nodes)
    print(f"  tau[v] == 2*selfloops + edges for all v? {ok}")
    # verify self-loops only on mixed?
    sl_nonmixed=[v for v in nodes if (selfB[v]+selfK[v])>0 and cat[v]!='mixed']
    print(f"  self-loops only on MIXED nodes? {len(sl_nonmixed)==0}  (violations: {[(v,cat[v],selfB[v],selfK[v]) for v in sl_nonmixed][:4]})")
    # parity per category
    for c in ['pureBLACK','mixed','pureBLUE']:
        vs=[v for v in nodes if cat[v]==c]
        if not vs: continue
        eb=[edgeB[v] for v in vs]; ek=[edgeK[v] for v in vs]
        print(f"  {c}: n={len(vs)}; edgeBLACK parity all even? {all(x%2==0 for x in ek)}; edgeBLUE parity all odd? {all(x%2==1 for x in eb)}; (blue range {min(eb)}-{max(eb)}, black {min(ek)}-{max(ek)})")
