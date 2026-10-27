import sys, itertools
n=int(sys.argv[1])
def ham(adj):
    full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1,full+1):
        row=dp[m]
        for v in range(n):
            c=row[v]
            if not c: continue
            for w in range(n):
                if adj[v][w] and not (m>>w)&1: dp[m|1<<w][w]+=c
    return sum(dp[full])
def auts(adj):
    s=[sum(r) for r in adj]
    res=[]; perm=[-1]*n; used=[False]*n
    def bt(i):
        if i==n: res.append(perm[:]); return
        for x in range(n):
            if used[x]: continue
            if all(adj[i][j]==adj[x][perm[j]] for j in range(i)):
                perm[i]=x; used[x]=True; bt(i+1); used[x]=False
        perm[i]=-1
    bt(0); return res
cnt=0
for line in sys.stdin:
    line=line.strip()
    if not line or line[0] not in '01': continue
    A=[[0]*n for _ in range(n)]; k=0
    for i in range(n):
        for j in range(i+1,n):
            if line[k]=='1': A[i][j]=1
            else: A[j][i]=1
            k+=1
    G=auts(A)
    orb={0}
    for g in G: orb.add(g[0])
    vt = len(orb)==n
    H=ham(A)
    cnt+=1
    print(f"regular class {cnt}: |Aut|={len(G)} vertex-transitive={vt} H={H} tc={H//len(G)}")
