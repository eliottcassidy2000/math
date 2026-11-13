import collections, sys
sys.setrecursionlimit(10000)
def hamcycle(out):
    n=len(out); res=[]
    def dfs(path,used):
        if res: return
        if len(path)==n:
            if out[path[-1]]&1: res.append(list(path))
            return
        v=path[-1]
        for w in range(n):
            if not used>>w&1 and out[v]>>w&1:
                path.append(w); dfs(path,used|1<<w); path.pop()
    dfs([0],1); return res[0] if res else None
def canon(out):
    n=len(out); cyc=hamcycle(out); best=None
    for s in range(n):
        lab=[cyc[(s+i)%n] for i in range(n)]
        w=tuple(out[lab[i]]>>lab[j]&1 for i in range(n) for j in range(i+1,n))
        if best is None or w<best: best=w
    return best
def conv(out):
    n=len(out); return [sum(1<<v for v in range(n) if out[v]>>w&1) for w in range(n)]
for n in (10,11,12):
    rows=[list(map(int,l.split()[1:])) for l in open('c1_%d.txt'%n)]
    ms=[]
    for l in open('c1arc_%d.txt'%n):
        left,right=l.split('|'); a=left.split(); ms.append((int(a[1]),int(a[2]),tuple(map(int,right.split()))))
    can=[canon(r) for r in rows]; idx={c:i for i,c in enumerate(can)}
    cv=[idx[canon(conv(r))] for r in rows]
    nsc=sum(1 for i in range(len(rows)) if cv[i]==i); merged=(len(rows)+nsc)//2
    G=collections.defaultdict(list)
    for i,m in enumerate(ms): G[m].append(i)
    bad=[g for g in G.values() if not (len(g)==1 and cv[g[0]]==g[0]) and not (len(g)==2 and cv[g[0]]==g[1])]
    print('order',n,': C=1 classes',len(rows),' self-converse',nsc,' merged',merged,' distinct multisets',len(G),' groups not explained by converse:',len(bad))
