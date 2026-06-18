from math import gcd
from itertools import combinations, permutations
MOD=14
UNITS=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def depth(r): r%=MOD; return min(r,MOD-r)
def m4(res,tb):
    m=len(res); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            wi=wj=0
            for a in UNITS:
                di=depth(res[i]*a); dj=depth(res[j]*a)
                if di>dj: wi+=1
                elif dj>di: wj+=1
            if wi>wj: adj[i][j]=1
            elif wj>wi: adj[j][i]=1
            else:
                if tb[i]<tb[j]: adj[i][j]=1
                else: adj[j][i]=1
    return adj
def canon(adj,m):
    best=None
    for p in permutations(range(m)):
        b=tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or b<best: best=b
    return best
def H(adj,m):
    c=0
    for p in permutations(range(m)):
        if all(adj[p[k]][p[k+1]] for k in range(m-1)): c+=1
    return c
def c3(adj,m):
    c=0
    for a,b,cc in combinations(range(m),3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c+=1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c+=1
    return c
def sc(adj,m): return tuple(sorted(sum(adj[i][j] for j in range(m) if j!=i) for i in range(m)))
m=5
FORB={(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))}
sigs=set(); hits=[]
# distinct nonzero residues, all orderings
for res in combinations(range(1,MOD),m):
    for perm in permutations(res):
        adj=m4(list(perm),list(range(m)))
        sig=(H(adj,m),c3(adj,m),sc(adj,m))
        sigs.add(sig)
        if sig in FORB: hits.append((perm,sig))
print("distinct-nonzero-residue SDR signatures realized:",len(sigs))
print("forbidden hits:", sorted(set(s for _,s in hits)) if hits else "NONE")
print("all realized sigs:")
for s in sorted(sigs): print("  ",s)
