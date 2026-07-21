#!/usr/bin/env python3
"""
tournament_zoo_full_deathstar_S79.py

"Keep adding to the zoo": a COMPREHENSIVE reusable tournament-invariant battery
(~20 invariants), computed over all tournaments n<=6, saved as a dataset, and
auto-mined for comparabilities / tight bounds / forbidden values. Supersedes the
S78 spectral engine by unifying the classical (kind-pasteur THM-1845, klein
THM-1850) and spectral (death-star THM-1858) zoos into ONE battery + adding new
invariants (cycle-trace spectrum, dichromatic, fas, arborescence sum, ham cycles,
domination, scc). Output feeds the invariant-atlas gap analysis.
"""
from math import comb
from itertools import combinations, permutations
import json
import numpy as np

def all_t(n):
    P = list(combinations(range(n), 2))
    for b in range(1 << len(P)):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(P):
            if b>>k&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def ham_paths_and_cycles(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1<<n):
        for v in range(n):
            c=dp[m][v]
            if c:
                for w in range(n):
                    if not(m>>w&1) and out[v]>>w&1: dp[m|1<<w][w]+=c
    H=sum(dp[full][v] for v in range(n))
    # ham cycles: paths 0..v then v->0
    cyc=sum(dp[full][v] for v in range(n) if A[v][0])  # fix start=0 (each cycle counted once)
    return H, cyc

def largest_transitive(A,n):
    best=0
    for m in range(1<<n):
        S=[v for v in range(n) if m>>v&1]
        if S and sorted(sum(A[v][w] for w in S) for v in S)==list(range(len(S))): best=max(best,len(S))
    return best

def domination(A,n):
    # min S: every vertex in S or beaten by some s in S
    for k in range(1,n+1):
        for S in combinations(range(n),k):
            dom=set(S)
            for s in S:
                for w in range(n):
                    if A[s][w]: dom.add(w)
            if len(dom)==n: return k
    return n

def kings(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    k=0
    for v in range(n):
        r=out[v]|(1<<v)
        for w in range(n):
            if out[v]>>w&1: r|=out[w]
        if r==full: k+=1
    return k

def scc_count(A,n):
    # count strongly connected components via reachability closure
    reach=[[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n): reach[i][i]=1
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if reach[i][k] and reach[k][j]: reach[i][j]=1
    seen=set(); c=0
    for v in range(n):
        if v in seen: continue
        comp={w for w in range(n) if reach[v][w] and reach[w][v]}
        seen|=comp; c+=1
    return c

def fas(A,n):
    # min feedback arc set = min back-arcs over all vertex orderings
    best=comb(n,2)
    for perm in permutations(range(n)):
        pos={v:i for i,v in enumerate(perm)}
        back=sum(1 for i in range(n) for j in range(n) if A[i][j] and pos[i]>pos[j])
        best=min(best,back)
        if best==0: break
    return best

def dichromatic(A,n):
    # min colors so each color class induces a transitive (acyclic) subtournament
    def acyclic(S):
        return sorted(sum(A[v][w] for w in S) for v in S)==list(range(len(S)))
    for k in range(1,n+1):
        # try to k-color; brute over assignments
        stack=[(0,[[] for _ in range(k)])]
        # simple backtracking
        def bt(v,parts):
            if v==n: return True
            for p in parts:
                p.append(v)
                if acyclic(p) and bt(v+1,parts): return True
                p.pop()
            return False
        if bt(0,[[] for _ in range(k)]): return k
    return n

def sum_arborescences(A,n):
    # sum over roots r of #spanning out-arborescences rooted at r (Matrix-Tree)
    total=0
    for r in range(n):
        # Laplacian L = Din - A^T restricted; #out-arbs rooted r = det(L with row/col r removed)
        # out-arborescence rooted r: edges point away from r; use L = D_in - A (in-degree)
        indeg=[sum(A[i][j] for i in range(n)) for j in range(n)]
        L=[[ (indeg[j] if i==j else 0) - A[i][j] for j in range(n)] for i in range(n)]
        M=[[L[i][j] for j in range(n) if j!=r] for i in range(n) if i!=r]
        total+=round(abs(np.linalg.det(np.array(M,float)))) if n>1 else 1
    return total

def spectral(A,n):
    M=np.array(A,float); ev=np.linalg.eigvals(M)
    mods=sorted((abs(z) for z in ev),reverse=True)
    ndev=len({(round(z.real,4),round(z.imag,4)) for z in ev})
    disc=round(abs(np.linalg.det(np.eye(n)+M-M.T)))//2**(n-1)
    p=[round(np.trace(np.linalg.matrix_power(M,k)).real) for k in (3,4,5)]
    return mods[0], (mods[1] if n>1 else 0.0), ndev, disc, p

def invariants(A,n):
    s=[sum(A[i]) for i in range(n)]
    H,cyc=ham_paths_and_cycles(A,n)
    rho,rho2,ndev,disc,p=spectral(A,n)
    return dict(n=n, c3=comb(n,3)-sum(comb(x,2) for x in s),
        p3=p[0],p4=p[1],p5=p[2], H=H, hamcyc=cyc, beta=largest_transitive(A,n),
        gamma=domination(A,n), kings=kings(A,n), scc=scc_count(A,n),
        fas=fas(A,n), dichr=dichromatic(A,n), arb=sum_arborescences(A,n),
        smax=max(s), smin=min(s), srange=max(s)-min(s), ndistscore=len(set(s)),
        ndev=ndev, disc=disc, rho=round(rho,4), rho2=round(rho2,4))

print("building FULL zoo n=3..6 (this takes a few min for n=6 fas/dichr) ...",flush=True)
DB={}
for n in (3,4,5,6):
    rows=[invariants(A,n) for A in all_t(n)]
    # dedup by full invariant vector for reporting
    seen=set(); uniq=[]
    for r in rows:
        key=tuple(sorted(r.items()))
        if key not in seen: seen.add(key); uniq.append(r)
    DB[n]=rows
    print(f"  n={n}: {len(rows)} tournaments, {len(uniq)} distinct invariant-profiles",flush=True)

# save dataset
with open("05-knowledge/results/tournament_zoo_full_deathstar_S79_data.json","w") as f:
    json.dump({str(n):DB[n] for n in DB}, f)
print("saved dataset -> 05-knowledge/results/tournament_zoo_full_deathstar_S79_data.json")

ALL=[r for n in DB for r in DB[n]]
INT=["c3","p3","p4","p5","H","hamcyc","beta","gamma","kings","scc","fas","dichr",
     "arb","smax","smin","srange","ndistscore","ndev","disc"]

print("\n===== FORBIDDEN VALUES (integer invariants that skip values, n<=6) =====")
for K in INT:
    vals=sorted({r[K] for r in ALL}); rng=list(range(min(vals),max(vals)+1))
    miss=[v for v in rng if v not in vals]
    if miss: print(f"  {K:11s}: {vals}   FORBIDDEN {miss}")

print("\n===== COMPARABILITY (X<=Y for ALL; T=tight) =====")
for X in INT:
    for Y in INT:
        if X!=Y and all(r[X]<=r[Y] for r in ALL):
            t="T" if any(r[X]==r[Y] for r in ALL) else " "
            print(f"  [{t}] {X} <= {Y}")

print("\n===== TIGHT ADDITIVE X <= Y + c (1<=c<=2, equality achieved) =====")
for X in INT:
    for Y in INT:
        if X==Y: continue
        c=max(r[X]-r[Y] for r in ALL)
        if 1<=c<=2 and any(r[X]==r[Y]+c for r in ALL):
            print(f"  {X} <= {Y} + {c}")
print("\ndone.")
