"""
ADVERSARIAL check:
 (F) strong-min(m) = min H over STRONGLY-CONNECTED tournaments on m vertices.
     Worker claims 3,5,9,15 for m=3..6 (exhaustive) and 25 for m=7 (jumps 15->25,
     skipping 21), with strong value sets excluding 7 and 21.
 (G) Multiplicativity: H(T) = prod over strong components of H(component).
     Test on a few reducible tournaments.
 (H) Strong H-spectra at m<=6 exclude 7 and 21.
"""
import itertools

def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    m=len(pairs)
    for bits in range(1<<m):
        A=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def ham_paths(A):
    n=len(A); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            cur=dp[mask][v]
            if cur==0 or not (mask>>v)&1: continue
            for w in range(n):
                if (mask>>w)&1: continue
                if A[v][w]==1: dp[mask|(1<<w)][w]+=cur
    return sum(dp[full][v] for v in range(n))

def is_strong(A):
    """strongly connected: from vertex 0 reach all, and all reach 0."""
    n=len(A)
    # forward reach
    def reach(adj_func):
        seen={0}; stack=[0]
        while stack:
            u=stack.pop()
            for v in range(n):
                if v not in seen and adj_func(u,v):
                    seen.add(v); stack.append(v)
        return len(seen)==n
    fwd=reach(lambda u,v: A[u][v]==1)
    bwd=reach(lambda u,v: A[v][u]==1)
    return fwd and bwd

def strong_components(A):
    """Tarjan-ish via Kosaraju. Return list of components (each a set), in topological
    order along the condensation (which is itself transitive for tournaments)."""
    n=len(A)
    visited=[False]*n; order=[]
    def dfs1(u):
        visited[u]=True
        for v in range(n):
            if A[u][v]==1 and not visited[v]:
                dfs1(v)
        order.append(u)
    for u in range(n):
        if not visited[u]: dfs1(u)
    comp=[-1]*n
    def dfs2(u,c):
        comp[u]=c
        for v in range(n):
            if A[v][u]==1 and comp[v]==-1:
                dfs2(v,c)
    c=0
    for u in reversed(order):
        if comp[u]==-1:
            dfs2(u,c); c+=1
    comps=[set() for _ in range(c)]
    for u in range(n): comps[comp[u]].add(u)
    return comps

def sub_tournament(A, vs):
    vs=sorted(vs); idx={v:i for i,v in enumerate(vs)}
    B=[[0]*len(vs) for _ in vs]
    for v in vs:
        for w in vs:
            if v!=w: B[idx[v]][idx[w]]=A[v][w]
    return B

if __name__=="__main__":
    print("=== (F)/(H) strong-min and strong H-spectra ===")
    for n in range(3,7):
        strong_spectrum=set()
        smin=None
        for A in all_tournaments(n):
            if is_strong(A):
                H=ham_paths(A)
                strong_spectrum.add(H)
                if smin is None or H<smin: smin=H
        print(f"m={n}: strong-min={smin}, strong-spectrum={sorted(strong_spectrum)}, "
              f"7 in? {7 in strong_spectrum}, 21 in? {21 in strong_spectrum}")

    print()
    print("=== (G) Multiplicativity H = prod over strong components ===")
    # test on all n=5 tournaments: verify H == prod of H over strong components
    bad=0; tested=0
    for n in (4,5):
        for A in all_tournaments(n):
            tested+=1
            H=ham_paths(A)
            comps=strong_components(A)
            prod=1
            for cset in comps:
                B=sub_tournament(A,cset)
                prod*=ham_paths(B)
            if prod!=H:
                bad+=1
                if bad<=5:
                    print(f"  MULT FAIL n={n}: H={H} prod={prod} comps={[len(c) for c in comps]}")
    print(f"  multiplicativity: tested={tested}, failures={bad}")
