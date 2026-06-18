from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
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
# also find a NON-trivially-lonely witness (gap not from all-odd 1/2) if any
def primitive(s):
    g=0
    for v in s: g=gcd(g,v)
    return g==1
print("PRIMARY WITNESS for forbidden signature (9,3,(1,1,2,3,3)):")
for S in [[1,3,5,9,19],[1,15,29,3,5]]:
    adj=m4(S)
    print(f"  S={S} res(mod14)={[v%14 for v in S]} sig=({H(adj)},{c3(adj)},{sco(adj)}) gap={Mg(S)[0]} tau={Mg(S)[1]}")
# Search for a forbidden-class witness whose loneliness is NOT just the all-odd 1/2 peak
# (i.e. has an even speed, so 1/2 doesn't trivially work) within 1..40.
print("\nSearching for forbidden-class witness that is NOT all-odd (loneliness off the 1/2 peak)...")
found=None
for s in combinations(range(1,41),m):
    if not primitive(s): continue
    if all(v%2==1 for v in s): continue   # require at least one even speed
    adj=m4(list(s))
    if (H(adj),c3(adj),sco(adj))==(9,3,(1,1,2,3,3)):
        g,at=Mg(list(s))
        if g>=F(1,14):
            found=(s,g,at); break
if found:
    print(f"  NON-all-odd lonely witness: S={found[0]} gap={found[1]}={float(found[1]):.5f} tau={found[2]}")
else:
    print("  none with an even speed in 1..40")
