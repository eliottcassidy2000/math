"""e10_analysis.py -- mac-mini-2026-10-27-potentials
Structure of the unique order-10 strong class without a black descent (THM-4811 (C1)): scores, minimum cuts, the nine
descent words (all blue), automorphisms and anti-automorphisms, black moves by resulting lambda, and the check that
every lambda = 2 neighbour still has a black descent avoiding E10."""
import itertools
code="111001111101111111111111111111111111000111101"
n=10
A=[[0]*n for _ in range(n)]
k=0
for i in range(n):
    for j in range(i+1,n):
        if code[k]=='1': A[i][j]=1
        else: A[j][i]=1
        k+=1
s=[sum(r) for r in A]
print("scores",s)
for i in range(n): print(i, "->", [j for j in range(n) if A[i][j]])
full=(1<<n)-1
best=99; cuts=[]
for X in range(1,full):
    f=sum(A[x][y] for x in range(n) if X>>x&1 for y in range(n) if not X>>y&1)
    if f<best: best=f; cuts=[]
    if f==best: cuts.append(X)
print("lambda",best,"min cuts",[[v for v in range(n) if X>>v&1] for X in cuts])
# all ham paths
paths=[]
def dfs(p,used):
    if len(p)==n: paths.append(tuple(p)); return
    for w in range(n):
        if not used>>w&1 and A[p[-1]][w]: dfs(p+[w],used|1<<w)
for v in range(n): dfs([v],1<<v)
print("H",len(paths))
def black(Q):
    r={Q[i]:Q[n-1-i] for i in range(n)}
    return any(A[u][v] and not A[r[v]][r[u]] for u in range(n) for v in range(n))
cand=[Q for Q in paths if any((X>>Q[0]&1) and not (X>>Q[-1]&1) for X in cuts)]
print("descent candidates",len(cand),"black among them",sum(black(Q) for Q in cand))
for Q in cand: print("  cand",Q,"black" if black(Q) else "blue", "closed" if A[Q[-1]][Q[0]] else "open", "scores ends",s[Q[0]],s[Q[-1]])
print("black paths total",sum(black(Q) for Q in paths))
# automorphisms / anti-automorphisms
auts=0; antis=0
for p in itertools.permutations(range(n)):
    if any(s[p[i]]!=s[i] and s[p[i]]!=n-1-s[i] for i in range(n)): continue
    if all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j): auts+=1
    if all(A[i][j]==A[p[j]][p[i]] for i in range(n) for j in range(n) if i!=j): antis+=1
print("|Aut|",auts,"#anti-automorphisms",antis)

def kmove(A,Q):
    B=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u!=v: B[u][v]=1-A[u][v]
    for u,v in zip(Q,Q[1:]): B[u][v]=1; B[v][u]=0
    return B
def lam_cuts(B):
    best=99; cs=[]
    for X in range(1,full):
        f=sum(B[x][y] for x in range(n) if X>>x&1 for y in range(n) if not X>>y&1)
        if f<best: best=f; cs=[]
        if f==best: cs.append(X)
    return best,cs
def has_black_descent(B):
    l,cs=lam_cuts(B)
    if l==0: return None
    ps=[]
    def dfs2(p,used):
        if len(p)==n: ps.append(tuple(p)); return
        for w in range(n):
            if not used>>w&1 and B[p[-1]][w]: dfs2(p+[w],used|1<<w)
    for v in range(n): dfs2([v],1<<v)
    def blk(Q):
        r={Q[i]:Q[n-1-i] for i in range(n)}
        return any(B[u][v] and not B[r[v]][r[u]] for u in range(n) for v in range(n))
    return any(blk(Q) and any((X>>Q[0]&1) and not (X>>Q[-1]&1) for X in cs) for Q in ps)
from collections import Counter
stats=Counter(); witness=None
for Q in paths:
    if not black(Q): continue
    B=kmove(A,Q)
    l,_=lam_cuts(B)
    if l==1:
        hb=has_black_descent(B)
        stats[(l,hb)]+=1
        if hb and witness is None: witness=Q
    else:
        stats[(l,None)]+=1
print("black moves from E10 by (lambda', black descent exists):",dict(stats))
print("witness first move",witness)

# lambda=2 neighbours of E10: do they have a black descent to a class other than E10?
E10sorted=sorted(s)
def allpaths(B):
    ps=[]
    def dfs2(p,used):
        if len(p)==n: ps.append(tuple(p)); return
        for w in range(n):
            if not used>>w&1 and B[p[-1]][w]: dfs2(p+[w],used|1<<w)
    for v in range(n): dfs2([v],1<<v)
    return ps
def blk_in(B,Q):
    r={Q[i]:Q[n-1-i] for i in range(n)}
    return any(B[u][v] and not B[r[v]][r[u]] for u in range(n) for v in range(n))
def canon(B):
    # cheap complete invariant for this check: sorted list of (score, sorted out-neighbour scores, sorted in-neighbour scores) plus H
    sc=[sum(r) for r in B]
    inv=sorted((sc[v],tuple(sorted(sc[w] for w in range(n) if B[v][w])),tuple(sorted(sc[w] for w in range(n) if B[w][v]))) for v in range(n))
    return tuple(inv)
E10inv=canon(A)
bad2=0; checked=0
seenT=set()
for Q in paths:
    B=kmove(A,Q)
    l,cs=lam_cuts(B)
    if l!=2: continue
    key=tuple(tuple(r) for r in B)
    if key in seenT: continue
    seenT.add(key)
    checked+=1
    ok=False
    for P in allpaths(B):
        if not any((X>>P[0]&1) and not (X>>P[-1]&1) for X in cs): continue
        if not blk_in(B,P): continue
        C=kmove(B,P)
        if canon(C)!=E10inv: ok=True; break
    if not ok: bad2+=1
print(f"lambda=2 neighbours of E10 (any colour move): {checked}; with no black descent avoiding the E10 invariant: {bad2}")
