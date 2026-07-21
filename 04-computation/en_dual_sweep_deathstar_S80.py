#!/usr/bin/env python3
"""
en_dual_sweep_deathstar_S80.py  (v2: n<=6 exact, flushed, WOWII-103 with real b)

THE E_n DUAL SWEEP -- even-graph metagraph E_n = Q_m/S_n (m=C(n-1,2)), the DUAL of
the tournament metagraph G_n. Cycle-space model: even graph = XOR of fundamental
triangles {0-i,i-j,0-j} over cotree edges (i,j), 1<=i<j<=n-1. Sanity: V=A002854.
Runs the invariant zoo + WOWII-103 (alpha <= floor(b - ln ecc_avg)) on E_n, to
compare against klein's alpha(G_n)=2,5,18 at n=4,5,6.
"""
import sys
from math import comb, log, floor
from itertools import combinations, product, permutations
from collections import defaultdict

def path_cycle(i,j,n):
    # fundamental cycle of cotree edge (i,j), j>i+1, w.r.t. the base PATH 0-1-...-(n-1):
    # the path i..j plus the chord (i,j). Matches the tournament base-path tiles (repo E_n).
    es=set()
    for k in range(i,j): es ^= {tuple(sorted((k,k+1)))}
    es ^= {(i,j)}
    return es

def even_graph(subset, n):
    ed={}
    for (i,j) in subset:
        for e in path_cycle(i,j,n):
            ed[e]=ed.get(e,0)^1
    return frozenset(e for e,v in ed.items() if v)

def canonical(edges,n):
    # canonical form = min adjacency bitstring over relabelings to canonical POSITIONS
    # (vertices grouped by degree -> position blocks; minimize over within-block orderings)
    deg=[0]*n
    for (i,j) in edges: deg[i]+=1; deg[j]+=1
    groups=defaultdict(list)
    for v in range(n): groups[deg[v]].append(v)
    classes=[groups[d] for d in sorted(groups)]
    best=None
    for parts in product(*[list(permutations(c)) for c in classes]):
        order=[v for p in parts for v in p]        # order[position] = old vertex
        pos={order[p]:p for p in range(n)}         # old vertex -> canonical position
        bits=0
        for (i,j) in edges:
            a,c=sorted((pos[i],pos[j])); bits|=1<<(a*n+c)
        if best is None or bits<best: best=bits
    return best

def build_En(n):
    cot=[(i,j) for i in range(n) for j in range(i+2,n)]; m=len(cot)   # path-tree cotree edges
    reps={}; lab=[0]*(1<<m)
    for bits in range(1<<m):
        c=canonical(even_graph([cot[k] for k in range(m) if bits>>k&1],n),n)
        if c not in reps: reps[c]=len(reps)
        lab[bits]=reps[c]
    V=len(reps); adj=[set() for _ in range(V)]
    for bits in range(1<<m):
        u=lab[bits]
        for k in range(m):
            v=lab[bits^(1<<k)]
            if v!=u: adj[u].add(v); adj[v].add(u)
    return V,adj

def bk_max(adj,V):
    best=[0]
    def bk(R,P,X):
        if not P and not X: best[0]=max(best[0],len(R)); return
        piv=max(P|X,key=lambda x:len(adj[x]&P)) if P|X else None
        for v in list(P-(adj[piv] if piv is not None else set())):
            bk(R|{v},P&adj[v],X&adj[v]); P=P-{v}; X=X|{v}
    bk(set(),set(range(V)),set()); return best[0]
def alpha_(adj,V): return bk_max([set(range(V))-adj[v]-{v} for v in range(V)],V)
def omega_(adj,V): return bk_max(adj,V)

def largest_bipartite(adj,V):
    best=0
    for mask in range(1<<V):
        S=[v for v in range(V) if mask>>v&1]
        if len(S)<=best: continue
        col={}; ok=True; Ss=set(S)
        for st in S:
            if st in col: continue
            col[st]=0; q=[st]
            while q and ok:
                x=q.pop()
                for w in adj[x]&Ss:
                    if w not in col: col[w]=col[x]^1; q.append(w)
                    elif col[w]==col[x]: ok=False; break
            if not ok: break
        if ok: best=len(S)
    return best

def ecc(adj,V):
    es=[]
    for s in range(V):
        d={s:0}; q=[s]
        while q:
            x=q.pop(0)
            for w in adj[x]:
                if w not in d: d[w]=d[x]+1; q.append(w)
        if len(d)<V: return None,None
        es.append(max(d.values()))
    return sum(es)/V,max(es)

def greedy_chi(adj,V):
    col={}
    for v in sorted(range(V),key=lambda x:-len(adj[x])):
        u={col[w] for w in adj[v] if w in col}; c=0
        while c in u: c+=1
        col[v]=c
    return max(col.values())+1 if col else 0

A002854={3:2,4:3,5:7,6:16}
ALPHA_Gn={4:2,5:5,6:18}   # klein's tournament metagraph, for comparison
print("E_n DUAL SWEEP (even-graph metagraph; dual of tournament G_n)\n"+"-"*66,flush=True)
for n in (3,4,5,6):
    V,adj=build_En(n)
    E=sum(len(a) for a in adj)//2; dens=E/comb(V,2) if V>1 else 0
    al=alpha_(adj,V); om=omega_(adj,V); ch=greedy_chi(adj,V)
    ea,dm=ecc(adj,V)
    b=largest_bipartite(adj,V) if V<=16 else None
    line=f"n={n}: V={V} (A002854={A002854[n]}{'OK' if V==A002854[n] else 'MISMATCH'}) "\
         f"E={E} dens={dens:.2f} alpha={al} omega={om} chi<={ch} "\
         f"ecc_avg={ea:.3f} diam={dm}"
    print(line,flush=True)
    if b is not None and ea is not None:
        rhs=floor(b-log(ea)); print(f"      WOWII-103: b(largest-induced-bipartite)={b}, "
              f"floor(b-ln ecc_avg)={rhs}; alpha={al} <= {rhs}? {al<=rhs}",flush=True)
    if n in ALPHA_Gn:
        print(f"      DUAL COMPARE: alpha(E_{n})={al}  vs  alpha(G_{n})={ALPHA_Gn[n]} (klein)",flush=True)
print("\nSanity: V(E_n)=2,3,7,16 must equal A002854. chi(E_n) grows faster than chi(G_n)=n-1.",flush=True)
