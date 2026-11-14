# Verify the root-chain transfer (H,E,lam,h) against brute force Hamiltonian path counts.
import itertools, random
def ham_stats(n, adj, r):
    # adj[i][j]=1 iff i->j. Return H, E(r), lam(r), h=H(T-r)
    full=(1<<n)-1
    # dp[mask][v] = #paths covering mask ending at v
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            for w in range(n):
                if not (mask>>w)&1 and adj[v][w]: dp[mask|(1<<w)][w]+=c
    def paths(mask): return 1 if mask==0 else sum(dp[mask])
    # paths starting at r covering mask: use reversed dp
    st=[[0]*n for _ in range(1<<n)]
    for v in range(n): st[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            c=st[mask][v]
            if not c: continue
            for w in range(n):
                if not (mask>>w)&1 and adj[w][v]: st[mask|(1<<w)][w]+=c
    H=paths(full); E=dp[full][r]
    lam=0
    for m2 in range(1<<n):
        if (m2>>r)&1:
            lam+=st[m2][r]*paths(full^m2)
    h=paths(full^(1<<r))
    return H,E,lam,h
def build(moves):
    # moves: list of (beats_root(0/1), move_root(0/1)); vertex 0 initial root
    n=len(moves)+1; adj=[[0]*n for _ in range(n)]
    r=0
    for k,(b,mv) in enumerate(moves,1):
        for j in range(k): adj[j][k]=1
        if b: adj[r][k]=0; adj[k][r]=1
        if mv: r=k
    return n,adj,r
def step(s,b,mv):
    H,E,lam,h=s
    if b==0 and mv==0: return (H,0,2*lam,h)
    if b==0 and mv==1: return (H,H,H,H)
    if b==1 and mv==0: return (H-E+lam,h,2*lam-h,h)
    return (H-E+lam,H-E,H+lam,H)
random.seed(1); bad=0
for trial in range(400):
    L=random.randint(1,8); moves=[(random.randint(0,1),random.randint(0,1)) for _ in range(L)]
    n,adj,r=build(moves); s=(1,1,1,1)
    for b,mv in moves: s=step(s,b,mv)
    if s!=ham_stats(n,adj,r): bad+=1; print("MISMATCH",moves,s,ham_stats(n,adj,r))
print("checked 400, mismatches",bad)
