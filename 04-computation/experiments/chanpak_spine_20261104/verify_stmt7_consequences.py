# Verify (1) the equality family for e_a <= (h-1)/2, (2) the sole-exit law H(X+v) = 2h - e_a + sigma, sigma >= 0,
# and that the naive equality H(X+v) = 2h - e_a fails in general.
import itertools, random
def ham_paths(n, adj):
    # adj[i][j] True iff i -> j ; count Hamiltonian paths by DP over subsets, return (h, ends[], starts[])
    full=(1<<n)-1
    # dp[mask][last] = number of paths covering mask ending at last
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        row=dp[mask]
        for last in range(n):
            c=row[last]
            if not c: continue
            for w in range(n):
                if not (mask>>w)&1 and adj[last][w]:
                    dp[mask|(1<<w)][w]+=c
    ends=dp[full][:]
    h=sum(ends)
    # starts: run on the converse
    adjc=[[adj[j][i] for j in range(n)] for i in range(n)]
    dp2=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp2[1<<v][v]=1
    for mask in range(1<<n):
        row=dp2[mask]
        for last in range(n):
            c=row[last]
            if not c: continue
            for w in range(n):
                if not (mask>>w)&1 and adjc[last][w]:
                    dp2[mask|(1<<w)][w]+=c
    starts=dp2[full][:]
    return h, ends, starts
print("== equality family: TT_m with arc (1,m) reversed (vertices 1..m -> indices 0..m-1) ==")
for m in range(3,9):
    adj=[[i<j for j in range(m)] for i in range(m)]
    adj[0][m-1]=False; adj[m-1][0]=True
    h,ends,starts=ham_paths(m,adj)
    eq=[v+1 for v in range(m) if 2*ends[v]==h-1 and any(adj[v])]
    print(f"m={m}: h={h} (2^(m-2)+1={2**(m-2)+1}), e_v={ends}, non-sink vertices with e=(h-1)/2: {eq}, e_(m-1)={ends[m-2]} =2^(m-3)={2**(m-3)}")
print("== sole-exit law on all tournaments n<=6 (X has n vertices, v -> a only) ==")
random.seed(1)
viol_lower=0; eq_naive=0; tot=0; max_sigma=0
for n in range(3,7):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[False]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: adj[i][j]=True
            else: adj[j][i]=True
        h,ends,starts=ham_paths(n,adj)
        for a in range(n):
            if not any(adj[a]): continue   # a is a sink -> skip (statement is about non-sinks)
            # extend by v = n: v -> a, all others -> v
            adj2=[row+[False] for row in adj]+[[False]*(n+1)]
            adj2[n][a]=True
            for x in range(n):
                if x!=a: adj2[x][n]=True
            h2,_,_=ham_paths(n+1,adj2)
            sigma=h2-(2*h-ends[a])
            tot+=1
            if sigma<0: viol_lower+=1
            if sigma==0: eq_naive+=1
            max_sigma=max(max_sigma,sigma)
            assert 2*ends[a]<=h-1, (n,bits,a,h,ends)
print(f"cases={tot}, sigma<0 violations={viol_lower}, naive equality (sigma=0) holds in {eq_naive} cases, max sigma={max_sigma}; e_a<=(h-1)/2 held in all")
