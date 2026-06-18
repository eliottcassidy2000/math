"""
Targeted GENUINE-SPEED realization of the forbidden M4 classes, with LRC-loneliness check.
Mechanism: residue collisions mod 14 create depth-ties broken by raw speed. We construct
speed sets with controlled residues (allowing several speeds to share a residue mod 14),
across a broad window, and test the EXACT M4 map + exact LRC gap.
"""
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
def canon(adj,m):
    best=None
    for p in permutations(range(m)):
        b=tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or b<best: best=b
    return best
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
FORB=[(9,3,(1,1,2,3,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))]
FORBset=set(FORB)
def primitive(s):
    g=0
    for v in s: g=gcd(g,v)
    return g==1
# Broad but smart: speeds 1..40 only (matches/exceeds claim author's window), full sweep but
# fast (C(40,5)=658008).
VMAX=40
found={}; found_lonely={}; n=0
for s in combinations(range(1,VMAX+1),m):
    if not primitive(s): continue
    n+=1
    adj=m4(list(s))
    if not valid(adj,m): continue
    sg=sig(adj,m)
    if sg in FORBset:
        if sg not in found: found[sg]=(s,canon(adj,m))
        if sg not in found_lonely:
            gp,at=Mg(list(s))
            if gp>=F(1,14): found_lonely[sg]=(s,gp,at)
print(f"GENUINE SPEEDS 1..{VMAX}: used {n} primitive 5-sets")
print("forbidden classes realized by genuine speeds:")
for sg in sorted(found):
    s,key=found[sg]
    print(f"   {sg}: first speeds {s}")
print("forbidden classes realized by LRC-LONELY (M>=1/14) genuine speeds:")
if found_lonely:
    for sg in sorted(found_lonely):
        s,gp,at=found_lonely[sg]
        print(f"   {sg}: speeds {s}, gap={gp}={float(gp):.5f}, tau={at}")
else:
    print("   none")
# Also report the canon key of the (9,3,...) realization
if (9,3,(1,1,2,3,3)) in found:
    s,key=found[(9,3,(1,1,2,3,3))]
    print(f"\n(9,3,...) canon key from genuine speeds {s}:\n   {key}")
