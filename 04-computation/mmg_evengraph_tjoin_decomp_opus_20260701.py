"""BLACK = even graph (Eulerian), BLUE = T-join (T=SC nodes). Verify + structure (components, cycle-rank,
bipartite), self-loop distribution, |SC| parity, recursion table."""
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
            A[xi][yi]=1 if bits[i]==0 else 0; A[yi][xi]=1-A[xi][yi]
        return A
    def canon(A):
        best=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]; A=bits_adj(bits); c=canon(A)
        cop=canon([[A[j][i] for j in range(n)] for i in range(n)])
        T.append(dict(mask=mask,c=c,cop=cop,gs=all(bits[i]==bits[TRANS[i]] for i in range(m) if TRANS[i]!=i),flip=mask^((1<<m)-1)))
    cls=sorted(set(t['c'] for t in T)); cid={c:i for i,c in enumerate(cls)}
    for t in T: t['ci']=cid[t['c']]; t['copi']=cid[t['cop']]
    tpt={t['ci']:t['copi'] for t in T}
    for t in T: t['mn']=min(t['ci'],tpt[t['ci']])
    return n,m,T,tpt
def components(nodes, edges):
    par={v:v for v in nodes}
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for a,b in edges: par[f(a)]=f(b)
    return len(set(f(v) for v in nodes))
def bipartite(nodes, edges):
    adj=defaultdict(list)
    for a,b in edges: adj[a].append(b); adj[b].append(a)
    col={}
    for s in nodes:
        if s in col: continue
        col[s]=0; st=[s]
        while st:
            u=st.pop()
            for w in adj[u]:
                if w not in col: col[w]=col[u]^1; st.append(w)
                elif col[w]==col[u]: return False
    return True
print(f"{'n':>2} {'B':>3} {'M':>3} {'K':>3} {'SC':>3} {'|SC|even':>8} {'blueL':>6} {'blackL':>7} {'blkEdges':>8} {'blkEvenGraph':>12} {'blk#comp':>8} {'blkCycRank':>10} {'blkBipart':>9} {'blueTjoin(T=SC)':>15}")
for n in [4,5,6]:
    n,m,T,tpt=build(n)
    nodes=sorted(set(t['mn'] for t in T))
    tau=Counter(t['mn'] for t in T); gs=Counter(t['mn'] for t in T if t['gs'])
    cat={v:('B' if gs[v]==tau[v] else 'K' if gs[v]==0 else 'M') for v in nodes}
    SC={v:(v==tpt[v]) for v in nodes}
    bym={t['mask']:t for t in T}; seen=set()
    blackedges=[]; blueedges=[]; selfB=Counter(); selfK=Counter()
    bdeg=Counter(); kdeg=Counter()
    for t in T:
        key=(min(t['mask'],t['flip']),max(t['mask'],t['flip']))
        if key in seen: continue
        seen.add(key); w=bym[t['flip']]; a,b=t['mn'],w['mn']; col='B' if t['gs'] else 'K'
        if a==b:
            (selfB if col=='B' else selfK)[a]+=1
        else:
            if col=='B': blueedges.append((a,b)); bdeg[a]+=1; bdeg[b]+=1
            else: blackedges.append((a,b)); kdeg[a]+=1; kdeg[b]+=1
    B=sum(1 for v in nodes if cat[v]=='B'); M=sum(1 for v in nodes if cat[v]=='M'); K=sum(1 for v in nodes if cat[v]=='K')
    SCn=B+M
    kn=sorted(set([v for e in blackedges for v in e]))
    blk_even=all(kdeg[v]%2==0 for v in nodes)
    comp=components(kn,blackedges) if kn else 0
    cycrank=len(blackedges)-len(kn)+comp if kn else 0
    bip=bipartite(kn,blackedges) if kn else True
    tjoin=all(bdeg[v]%2==1 for v in nodes if SC[v]) and all(bdeg[v]%2==0 for v in nodes if not SC[v])
    print(f"{n:>2} {B:>3} {M:>3} {K:>3} {SCn:>3} {str(SCn%2==0):>8} {2**((m+ (n-1)//2)//2 -1):>6} {2**(m-1)-2**((m+(n-1)//2)//2-1):>7} {len(blackedges):>8} {str(blk_even):>12} {comp:>8} {cycrank:>10} {str(bip):>9} {str(tjoin):>15}")
    print(f"    self-loops: blue on {dict(selfB)} (cats {[cat[v] for v in selfB]}); black on {dict(selfK)} (cats {[cat[v] for v in selfK]})")
