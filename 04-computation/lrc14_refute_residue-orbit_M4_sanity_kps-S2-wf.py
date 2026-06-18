# Sanity: confirm M4 is a pure function of residues mod 14 (+ tie-break), and that the
# witness's loneliness is real. Also: does a NON-all-odd lonely witness exist at all
# (broader window 1..70, residues forced to (1,3,5,9,5)-type collision with an even speed)?
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
def sig(adj):
    H=sum(1 for p in permutations(range(m)) if all(adj[p[k]][p[k+1]] for k in range(m-1)))
    c3=0
    for a,b,cc in combinations(range(m),3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c3+=1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c3+=1
    sc=tuple(sorted(sum(adj[i][j] for j in range(m) if j!=i) for i in range(m)))
    return (H,c3,sc)
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
def prim(s):
    g=0
    for v in s: g=gcd(g,v)
    return g==1
# (a) residue-function check: same residues+order -> same tournament regardless of representative
import random
random.seed(1)
ok=True
for _ in range(2000):
    base=[random.randint(1,13) for _ in range(5)]
    S1=base[:]
    S2=[base[i]+14*random.randint(0,5) for i in range(5)]
    # tie-break uses raw speed -> must keep same RELATIVE order to be equal; enforce by sorting key
    # Instead: M4 equal iff residues equal AND speed order equal. Check residue dependence by fixing order:
    # use speeds = residue + 14*rank so order matches residue order
    order=sorted(range(5), key=lambda i: base[i])
    rank={order[i]:i for i in range(5)}
    Sa=[base[i] for i in range(5)]
    Sb=[base[i]+14*(rank[i]+1) for i in range(5)]   # preserves residue, and order (adds bigger offset to higher rank)
    if m4(Sa)!=m4(Sb): ok=False; break
print("M4 depends only on (residues mod 14, speed order):", ok)
# (b) broaden: any forbidden-sig (9,3,...) lonely witness with an EVEN speed, 1..70?
found=None
cnt=0
for s in combinations(range(1,71),5):
    if any(v%2==0 for v in s) and prim(s):
        cnt+=1
        if cnt>400000: break
        if sig(m4(list(s)))==(9,3,(1,1,2,3,3)):
            g,at=Mg(list(s))
            if g>=F(1,14):
                found=(s,g,at); break
print("non-all-odd lonely witness for (9,3,...) within scan:", found, "(scanned",cnt,"even-containing sets)")
