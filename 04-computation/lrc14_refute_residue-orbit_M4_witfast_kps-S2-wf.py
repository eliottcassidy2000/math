from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product
MOD=14; UNITS=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def depth(r): r%=MOD; return min(r,MOD-r)
def m4(S):
    m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            wi=wj=0
            for a in UNITS:
                di=depth(S[i]*a); dj=depth(S[j]*a)
                if di>dj: wi+=1
                elif dj>di: wj+=1
            if wi>wj: adj[i][j]=1
            elif wj>wi: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
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
def sco(adj,m): return tuple(sorted(sum(adj[i][j] for j in range(m) if j!=i) for i in range(m)))
def sig(adj,m): return (H(adj,m),c3(adj,m),sco(adj,m))
def valid(adj,m): return all(adj[i][j]+adj[j][i]==1 for i in range(m) for j in range(i+1,m))
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def gm(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mg(S):
    b=F(0); at=None
    for t in cand(S):
        v=gm(S,t)
        if v>b: b=v; at=t
    return b,at
m=5
# uniqueness of forbidden signatures among the 12 free classes
free={}; pairs=list(combinations(range(m),2))
for bits in product([0,1],repeat=len(pairs)):
    adj=[[0]*m for _ in range(m)]
    for (i,j),bb in zip(pairs,bits):
        if bb: adj[i][j]=1
        else: adj[j][i]=1
    free[canon(adj,m)]=sig(adj,m)
by_sig={}
for k,s in free.items(): by_sig.setdefault(s,[]).append(k)
FORB=[(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))]
print("free classes total:",len(free))
for f in FORB: print(f"  forbidden {f}: #iso classes w/ this signature = {len(by_sig.get(f,[]))}")
# genuine-speed witness for (9,3,...)
for S in ([1,15,29,3,5],[3,5,1,15,29],[1,3,5,15,29]):
    adj=m4(S); print(f"\nS={S} residues={[v%MOD for v in S]} valid={valid(adj,m)} sig={sig(adj,m)}")
    g,at=Mg(S); print(f"  gap={g}={float(g):.6f}  >=1/14? {g>=F(1,14)}  tau={at}")
