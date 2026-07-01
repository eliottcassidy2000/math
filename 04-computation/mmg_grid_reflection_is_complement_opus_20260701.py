"""Deep structure of the blue/black flip-line merged metagraph. (1) grid-reflection R=transMask action on class;
(2) category-adjacency (which categories connect via blue/black); (3) tau parity; (4) self-loop distribution."""
from itertools import permutations
from collections import Counter, defaultdict
def build(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES)
    tileidx={t:i for i,t in enumerate(TILES)}; TRANS=[tileidx[(n-y+1,n-x+1)] for (x,y) in TILES]
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
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        A=bits_adj(bits); c=canon(A)
        Aop=[[A[j][i] for j in range(n)] for i in range(n)]; cop=canon(Aop)
        transbits=[0]*m
        for i in range(m): transbits[TRANS[i]]=bits[i]
        Rt=canon(bits_adj(transbits))
        T.append(dict(mask=mask,bits=bits,c=c,cop=cop,Rc=Rt,gs=all(bits[i]==bits[TRANS[i]] for i in range(m) if TRANS[i]!=i),flip=mask^((1<<m)-1)))
    cls=sorted(set(t['c'] for t in T)); cid={c:i for i,c in enumerate(cls)}
    for t in T: t['ci']=cid[t['c']]; t['copi']=cid[t['cop']]; t['Rci']=cid[t['Rc']]
    tpt={t['ci']:t['copi'] for t in T}
    for t in T: t['mn']=min(t['ci'],tpt[t['ci']])
    return n,m,T,cid,tpt
for n in [4,5,6]:
    n,m,T,cid,tpt=build(n)
    nodes=sorted(set(t['mn'] for t in T))
    tau=Counter(t['mn'] for t in T); gs=Counter(t['mn'] for t in T if t['gs'])
    SC={v:(v==tpt[v]) for v in nodes}
    cat={v:('pureBLUE' if gs[v]==tau[v] else 'pureBLACK' if gs[v]==0 else 'mixed') for v in nodes}
    # (1) R action: does transMask preserve class or complement it?
    Rpres=sum(1 for t in T if t['Rci']==t['ci']); Rcomp=sum(1 for t in T if t['Rci']==t['copi'])
    print(f"\n== n={n} == m={m}, {len(nodes)} nodes; R(grid-refl): class-preserving {Rpres}/{len(T)}, ->complement {Rcomp}/{len(T)}")
    # (2) category adjacency + self-loops
    bym={t['mask']:t for t in T}; seen=set()
    adjB=defaultdict(int); adjK=defaultdict(int); selfloop=defaultdict(lambda:[0,0])  # per cat: [blue,black]
    edgeB=Counter(); edgeK=Counter()
    for t in T:
        key=(min(t['mask'],t['flip']),max(t['mask'],t['flip']))
        if key in seen: continue
        seen.add(key)
        w=bym[t['flip']]; col='B' if t['gs'] else 'K'; a,b=t['mn'],w['mn']
        if a==b:
            selfloop[cat[a]][0 if col=='B' else 1]+=1
        else:
            ca,cb=sorted([cat[a],cat[b]])
            if col=='B': adjB[(ca,cb)]+=1; edgeB[a]+=1; edgeB[b]+=1
            else: adjK[(ca,cb)]+=1; edgeK[a]+=1; edgeK[b]+=1
    print(f"  BLUE edges by category-pair: {dict(adjB)}")
    print(f"  BLACK edges by category-pair: {dict(adjK)}")
    print(f"  self-loops by category [blue,black]: {dict(selfloop)}")
    # (3) tau parity by SC/NS
    par=defaultdict(Counter)
    for v in nodes: par['SC' if SC[v] else 'NS'][tau[v]%2]+=1
    print(f"  tau parity: SC nodes {dict(par['SC'])} (0=even,1=odd), NS nodes {dict(par['NS'])}  => SC odd & NS even? {par['SC'][0]==0 and par['NS'][1]==0}")
    # (4) blue-edge-degree odd on SC(pureBLUE/mixed) & even(0) on NS; black-edge-degree even everywhere
    beOdSC=all(edgeB[v]%2==1 for v in nodes if SC[v]); keEven=all(edgeK[v]%2==0 for v in nodes)
    print(f"  blue-edge-deg ODD on all SC nodes? {beOdSC}; black-edge-deg EVEN on ALL nodes? {keEven}")
