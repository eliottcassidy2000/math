from math import gcd
from itertools import combinations, permutations, product, combinations_with_replacement
MOD=14; UNITS=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def depth(r): r%=MOD; return min(r,MOD-r)
WIN=[[0]*MOD for _ in range(MOD)]
for x in range(MOD):
    for y in range(MOD):
        wi=wj=0
        for a in UNITS:
            di=depth(x*a); dj=depth(y*a)
            if di>dj: wi+=1
            elif dj>di: wj+=1
        WIN[x][y]=1 if wi>wj else (-1 if wj>wi else 0)
m=5
def canon(adj):
    best=None
    for p in permutations(range(m)):
        b=tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or b<best: best=b
    return best
def H(adj):
    c=0
    for p in permutations(range(m)):
        if all(adj[p[k]][p[k+1]] for k in range(m-1)): c+=1
    return c
def c3(adj):
    c=0
    for a,b,cc in combinations(range(m),3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c+=1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c+=1
    return c
def sco(adj): return tuple(sorted(sum(adj[i][j] for j in range(m) if j!=i) for i in range(m)))
def sig(adj): return (H(adj),c3(adj),sco(adj))

# all 12 free classes, group by signature, show canon keys for the (9,3,...) signature
free={}; pairs=list(combinations(range(m),2))
for bits in product([0,1],repeat=len(pairs)):
    adj=[[0]*m for _ in range(m)]
    for (i,j),bb in zip(pairs,bits):
        if bb: adj[i][j]=1
        else: adj[j][i]=1
    k=canon(adj)
    if k not in free: free[k]=sig(adj)
tgt=(9,3,(1,1,2,3,3))
twins=[k for k,s in free.items() if s==tgt]
print(f"signature {tgt} has {len(twins)} iso classes (canon keys):")
for k in twins: print("   key:",k)

# Now: build the EXACT M4 ceiling as a set of canon KEYS (not just sigs), allowing
# repetition + all tie-break TOTAL orders. Determine which canon keys are reached.
def build(res,order):
    adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            w=WIN[res[i]][res[j]]
            if w==1: adj[i][j]=1
            elif w==-1: adj[j][i]=1
            else: adj[i][j]=1 if order[i]<order[j] else 0; adj[j][i]=1-adj[i][j]
    return adj
reached=set()
PERMS=list(permutations(range(m)))
for res in combinations_with_replacement(range(MOD),m):
    tie=any(WIN[res[i]][res[j]]==0 for i in range(m) for j in range(i+1,m))
    orders=PERMS if tie else [tuple(range(m))]
    for order in orders:
        reached.add(canon(build(res,order)))
print(f"\nEXACT M4 ceiling reaches {len(reached)} of 12 iso classes (canon keys).")
# which forbidden keys are reached?
FORB_sig=[(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))]
for fs in FORB_sig:
    keys=[k for k,s in free.items() if s==fs]
    for k in keys:
        print(f"  forbidden-sig {fs} key reached? {k in reached}")
# overall: which free classes NOT reached
notreached=[ (free[k],k) for k in free if k not in reached]
print(f"\nfree classes NOT reached by exact M4 ceiling: {len(notreached)}")
for s,k in sorted(notreached): print("   ",s)
print(f"free classes reached: {len([k for k in free if k in reached])}/12")
