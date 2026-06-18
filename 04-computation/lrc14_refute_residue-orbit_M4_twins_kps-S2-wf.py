from math import gcd
from itertools import combinations, permutations, product
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
def aut(adj):  # |Aut| size to distinguish the twins
    cnt=0
    for p in permutations(range(m)):
        if all(adj[p[i]][p[j]]==adj[i][j] for i in range(m) for j in range(m) if i!=j): cnt+=1
    return cnt
# enumerate the 12 free classes
free={}; pairs=list(combinations(range(m),2))
for bits in product([0,1],repeat=len(pairs)):
    adj=[[0]*m for _ in range(m)]
    for (i,j),bb in zip(pairs,bits):
        if bb: adj[i][j]=1
        else: adj[j][i]=1
    k=canon(adj)
    if k not in free: free[k]=(sig(adj),aut(adj))
tgt=(9,3,(1,1,2,3,3))
twins=[(k,info) for k,info in free.items() if info[0]==tgt]
print(f"signature {tgt}: {len(twins)} iso classes")
for k,info in twins:
    print(f"   key={k}  |Aut|={info[1]}")
# witness canon key from genuine speeds (1,3,5,9,19):
wkey=(0,0,0,1,1,0,0,0,1,1,0,0,1,1,1,0,0,1,1,1)
# rebuild adj from key, get its canon (should already be canon) + aut + 3cyc-structure
def from_key(k):
    adj=[[0]*m for _ in range(m)]; idx=0
    for i in range(m):
        for j in range(m):
            if i!=j: adj[i][j]=k[idx]; idx+=1
    return adj
wadj=from_key(wkey); wk=canon(wadj)
print(f"\nwitness {{1,3,5,9,19}} canon key = {wk}")
print(f"  matches which twin? {[ (k==wk) for k,_ in twins]}")
print(f"  witness |Aut| = {aut(wadj)}")
# Is the witness's class among the 8 reachable by perfect-SDR M4? Recompute SDR-reachable KEYS.
MOD=14; UNITS=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def depth(r): r%=MOD; return min(r,MOD-r)
def m4res(res):
    adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            wi=wj=0
            for a in UNITS:
                di=depth(res[i]*a); dj=depth(res[j]*a)
                if di>dj: wi+=1
                elif dj>di: wj+=1
            if wi>wj: adj[i][j]=1
            elif wj>wi: adj[j][i]=1
            else: adj[i][j]=1 if res[i]<res[j] else 0; adj[j][i]=1-adj[i][j]  # tie-break by residue value (proxy)
    return adj
sdr_keys=set()
for res in combinations(range(1,MOD),m):
    for perm in permutations(res):
        sdr_keys.add(canon(m4res(list(perm))))
print(f"\nperfect-SDR (distinct nonzero residue) reachable KEYS: {len(sdr_keys)}")
for k,_ in twins:
    print(f"  twin {k}: reachable by perfect SDR? {k in sdr_keys}")
print(f"  witness class reachable by perfect SDR? {wk in sdr_keys}")
