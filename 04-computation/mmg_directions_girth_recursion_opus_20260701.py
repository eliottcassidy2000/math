"""Work the other directions: (A) black even-graph odd-cycle girth; (B) recursion dim D(n)=(m+f)/2;
(C) realization-degeneracy: is the black even-graph the UNIQUE even graph with its degree seq on M-union-K?"""
from itertools import permutations, combinations
from collections import defaultdict, Counter
def build(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES)
    tileidx={t:i for i,t in enumerate(TILES)}; TRANS=[tileidx[(n-y+1,n-x+1)] for (x,y) in TILES]
    perms=list(permutations(range(n)))
    def bits_adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=n-xL; yi=n-yL; A[xi][yi]=1 if bits[i]==0 else 0; A[yi][xi]=1-A[xi][yi]
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
    return n,m,T,tpt,len(TILES),sum(1 for i in range(m) if TRANS[i]==i)
def girth_odd(adj, nodes):  # BFS shortest cycle + shortest odd cycle (simple graph)
    import collections
    def bfs(src, parity_target=None):
        best=99
        dist={src:0}; par={src:-1}; q=collections.deque([src])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w]=dist[u]+1; par[w]=u; q.append(w)
                elif par[u]!=w:
                    cyc=dist[u]+dist[w]+1
                    if parity_target is None or cyc%2==parity_target: best=min(best,cyc)
        return best
    g=min((bfs(s) for s in nodes),default=99); godd=min((bfs(s,1) for s in nodes),default=99)
    return g,godd
print("n | D=(m+f)/2 | black even-graph: nodes edges  girth oddGirth  |  blue Tjoin nodes edges")
for n in [4,5,6]:
    n,m,T,tpt,mm,f=build(n)
    nodes=sorted(set(t['mn'] for t in T))
    tau=Counter(t['mn'] for t in T); gs=Counter(t['mn'] for t in T if t['gs'])
    cat={v:('B' if gs[v]==tau[v] else 'K' if gs[v]==0 else 'M') for v in nodes}
    bym={t['mask']:t for t in T}; seen=set()
    Kadj=defaultdict(set); Badj=defaultdict(set); ke=0; be=0
    for t in T:
        key=(min(t['mask'],t['flip']),max(t['mask'],t['flip']))
        if key in seen: continue
        seen.add(key); w=bym[t['flip']]; a,b=t['mn'],w['mn']
        if a==b: continue
        if t['gs']: Badj[a].add(b); Badj[b].add(a); be+=1
        else: Kadj[a].add(b); Kadj[b].add(a); ke+=1
    knodes=[v for v in nodes if Kadj[v]]; bnodes=[v for v in nodes if Badj[v]]
    kg,kog=girth_odd(Kadj,knodes) if knodes else (0,0)
    D=(mm+f)//2
    print(f"{n} |   {D}     | {len(knodes):>4} {ke:>5}   {kg:>4} {kog:>7}    | {len(bnodes):>4} {be:>4}")
# (B) recursion dims and even-n SC=A000568(n-1)
A={3:2,4:4,5:12,6:56,7:456}
print("\nrecursion: D=(m+f)/2 = 2,4,6,9 (n=4..7); SC=2,8,12,88; even-n SC(n)=A000568(n-1)? ", 
      f"n4:2=={A[3]}({2==A[3]}) n6:12=={A[5]}({12==A[5]})")
print("half-tiling DIM D vs (n-1)-staircase C(n-2,2):", {n:( (C:=[0,0,0,1,3,6,10,15][n-1]) , 'match' if (([0,0,0,1,3,6,10,15][n-1])==( ([n2 for n2 in [n]][0]) )) else '') for n in [4,5,6,7]})
