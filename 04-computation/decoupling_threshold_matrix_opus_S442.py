#!/usr/bin/env python3
"""
THE UNIVERSAL DECOUPLING THRESHOLD (opus-S442, invariants/monoids/orbits lens).
A "decoupling" (invariant X not determined by invariant Y) is an ORBIT-REFINEMENT: the Y-level-set
(a coarse orbit) splits into several X-values. THM-1865 (H not score-determined) and THM-1930
(var(lambda^2) not c3-determined) BOTH first break at n=5. Is n=5 universal?

For each ordered pair (X | Y): threshold = smallest n where some Y-value carries >=2 X-values
(i.e. X is NOT a function of Y). Compute over ALL tournaments n=3..6.
Invariants: H (Ham paths), c3 (3-cycles), score (sorted out-degrees), charS (char-poly coeff tuple),
var (=var of lambda^2 -> use tr(S^4) exact int proxy), scc (#strong components), fas (min feedback arc set).
"""
import itertools, numpy as np

def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            adj[i][j]=1 if (bits>>k)&1 else 0; adj[j][i]=1-adj[i][j]
        yield adj

def H(adj,n):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if c:
                for u in range(n):
                    if not (mask>>u)&1 and adj[v][u]: dp[mask|(1<<u)][u]+=c
    return sum(dp[size-1])

def c3(adj,n):
    t=0
    for i,j,k in itertools.combinations(range(n),3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]): t+=1
    return t

def score(adj,n): return tuple(sorted(sum(adj[i]) for i in range(n)))
def tr4(adj,n):
    S=np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])
    S2=S@S; return int(round(np.trace(S2@S2)))
def charS(adj,n):
    S=np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])
    return tuple(int(round(c)) for c in np.poly(S).real)
def scc(adj,n):
    # Tarjan-lite via reachability closure
    reach=[[adj[i][j] or i==j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]: reach[i][j]=True
    comp={}
    cid=0
    seen=[False]*n
    for i in range(n):
        if seen[i]: continue
        grp=[j for j in range(n) if reach[i][j] and reach[j][i]]
        for j in grp: seen[j]=True
        cid+=1
    return cid

INV={'H':H,'c3':c3,'score':score,'tr4':tr4,'charS':charS,'scc':scc}
pairs=[('H','score'),('H','c3'),('tr4','c3'),('tr4','score'),('charS','score'),
       ('H','charS'),('H','tr4'),('scc','score'),('c3','score')]

print("DECOUPLING THRESHOLD: smallest n where X is NOT a function of Y (orbit-refinement onset)")
print("="*74)
data={}
for n in range(3,7):
    vals={name:[f(adj,n) for adj in edges_iter(n)] for name,f in INV.items()}
    data[n]=vals
print(f"{'X | Y':<16} {'n=3':>5} {'n=4':>5} {'n=5':>5} {'n=6':>5}   threshold")
for X,Y in pairs:
    row=[]; thr=None
    for n in range(3,7):
        byY={}
        for xv,yv in zip(data[n][X],data[n][Y]):
            byY.setdefault(yv,set()).add(xv)
        decoupled = any(len(s)>1 for s in byY.values())
        row.append('SPLIT' if decoupled else 'det')
        if decoupled and thr is None: thr=n
    print(f"{X+' | '+Y:<16} {row[0]:>5} {row[1]:>5} {row[2]:>5} {row[3]:>5}   {thr if thr else '>6'}")
print("\n READING: threshold = the n at which the Y-orbit first refines into multiple X-values.")
print(" If spectral/path-from-combinatorial pairs cluster at n=5 => the 'quaternion wall' conjecture.")
