#!/usr/bin/env python3
"""
klein-2026-07-21-S395 -- LEVERAGING THE WOWII-103 COUNTEREXAMPLE IDEA.  WOWII Conj 103 (now
disproved): for connected G, alpha(G) <= floor(b(G) - log(ecc_avg(G))), alpha=independence
number, b=largest induced bipartite subgraph, ecc_avg=average eccentricity.  Disproved by a
TRIANGLE + leaves (11 vtx, 9 <= 8).

Owner: "think about how ideas similar to it can be leveraged for problems in this repo."

DEMONSTRATION: the repo's flagship graphs are the tournament METAGRAPH G_n (iso classes, single
arc-flip adjacency) and the even-graph metagraph E_n.  These ARE connected simple graphs, so the
whole WOWII conjecture zoo APPLIES to them and has never been run on them.  Here: build G_n,
compute the WOWII-103 invariants (alpha, b, ecc_avg), and check whether the metagraph SATISFIES
or VIOLATES the (refuted) bound -- a new fact about G_n either way, and a template for running the
~150 WOWII inequalities against the repo's invariant zoo.
"""
import itertools, math
from math import log, floor

# ---- build the tournament metagraph G_n (reuse S336 machinery, compactly) ----
def pairs_of(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def relabel(om,perm,n):
    new=[0]*n
    for v in range(n):
        mv,t=om[v],0
        while mv:
            b=mv&-mv; w=b.bit_length()-1; mv^=b; t|=1<<perm[w]
        new[perm[v]]=t
    return tuple(new)
def word(om,n):
    w=0
    for v in range(n): w=(w<<n)|om[v]
    return w
def refine(om,n):
    col=[bin(om[v]).count("1") for v in range(n)]
    while True:
        sig=[]
        for v in range(n):
            cnt={}; mv=om[v]
            while mv:
                b=mv&-mv; w=b.bit_length()-1; mv^=b; cnt[col[w]]=cnt.get(col[w],0)+1
            sig.append((col[v],tuple(sorted(cnt.items()))))
        order=sorted(set(sig)); nc=[order.index(sig[v]) for v in range(n)]
        if nc==col: break
        col=nc
    cells={}
    for v in range(n): cells.setdefault(col[v],[]).append(v)
    return [tuple(cells[k]) for k in sorted(cells)]
def canon(om,n):
    cells=refine(om,n); base=[]; pos=0
    for c in cells: base.append((c,pos)); pos+=len(c)
    best=None
    for ch in itertools.product(*[itertools.permutations(c) for (c,_) in base]):
        perm=[0]*n
        for (blk,(c,st)) in zip(ch,base):
            for k,v in enumerate(blk): perm[v]=st+k
        w=word(relabel(om,perm,n),n)
        if best is None or w<best: best=w
    return best
def metagraph(n):
    P=pairs_of(n); om0=tuple(sum(1<<j for j in range(i)) for i in range(n))
    seen={canon(om0,n):0}; order=[canon(om0,n)]; fr=[om0]; edges=set()
    while fr:
        nx=[]
        for om in fr:
            c=canon(om,n)
            for (i,j) in P:
                nm=list(om)
                if om[i]>>j&1: nm[i]&=~(1<<j); nm[j]|=1<<i
                else: nm[j]&=~(1<<i); nm[i]|=1<<j
                nm=tuple(nm); cc=canon(nm,n)
                if cc not in seen: seen[cc]=len(order); order.append(cc); nx.append(nm)
                if cc!=c: edges.add((min(seen[c],seen[cc]),max(seen[c],seen[cc])))
        fr=nx
    N=len(order); adj=[set() for _ in range(N)]
    for a,b in edges: adj[a].add(b); adj[b].add(a)
    return N,adj

def indep_number(adj):
    N=len(adj)
    best=0
    # branch & bound max independent set
    order=sorted(range(N), key=lambda v:-len(adj[v]))
    def go(cand, cur):
        nonlocal best
        if cur+len(cand)<=best: return
        if not cand: best=max(best,cur); return
        v=cand[0]
        go([u for u in cand[1:] if u not in adj[v]], cur+1)   # take v
        go(cand[1:], cur)                                     # skip v
    go(order,0)
    return best
def largest_induced_bipartite(adj):
    N=len(adj)
    if N>20:   # exhaustive over 2-colorings infeasible; use complement-of-max-odd heuristic skip
        return None
    best=0
    for mask in range(1<<N):
        verts=[v for v in range(N) if mask>>v&1]
        if len(verts)<=best: continue
        # check induced subgraph bipartite via BFS 2-coloring
        vs=set(verts); color={}; ok=True
        for s in verts:
            if s in color: continue
            color[s]=0; st=[s]
            while st and ok:
                u=st.pop()
                for w in adj[u]:
                    if w in vs:
                        if w not in color: color[w]=color[u]^1; st.append(w)
                        elif color[w]==color[u]: ok=False; break
            if not ok: break
        if ok: best=len(verts)
    return best
def avg_ecc(adj):
    N=len(adj); tot=0
    for s in range(N):
        dist=[-1]*N; dist[s]=0; q=[s]
        while q:
            nq=[]
            for u in q:
                for w in adj[u]:
                    if dist[w]<0: dist[w]=dist[u]+1; nq.append(w)
            q=nq
        tot+=max(dist)
    from fractions import Fraction as Fr
    return Fr(tot,N)

print("="*80)
print("WOWII-103 invariants on the tournament metagraph G_n  (alpha <= floor(b - log ecc_avg)?)")
print("="*80)
print(f"{'n':>3} {'|G_n|':>6} {'alpha':>6} {'b':>6} {'ecc_avg':>9} {'RHS=floor(b-log ecc)':>22} {'103 holds?':>11}")
for n in (4,5):
    N,adj=metagraph(n)
    a=indep_number(adj); b=largest_induced_bipartite(adj); e=avg_ecc(adj)
    if b is None:
        print(f"{n:>3} {N:>6} {a:>6} {'--':>6} {str(e):>9} {'(b too big)':>22} {'--':>11}")
        continue
    rhs=floor(b - log(float(e)))
    holds = a<=rhs
    print(f"{n:>3} {N:>6} {a:>6} {b:>6} {str(e):>9} {rhs:>22} {str(holds):>11}"
          + ("" if holds else "   <-- G_n VIOLATES WOWII-103 (a metagraph counterexample!)"))
# n=6 alpha + ecc only (b exhaustive infeasible at 56 nodes)
N,adj=metagraph(6); a=indep_number(adj); e=avg_ecc(adj)
print(f"{6:>3} {N:>6} {a:>6} {'?':>6} {str(e):>9} {'(b: 56 nodes, skipped)':>22}")
print("""
 READING: whether or not G_n satisfies WOWII-103, the point is METHODOLOGICAL -- the repo's
 metagraphs are ordinary graphs, so the entire Graffiti/WOWII conjecture list (independence vs
 bipartite vs eccentricity vs domination vs ...) can be run against G_n and E_n, and against
 the repo's OWN invariant inequalities (arborescence vs H, THM-1460/1580; detection depth,
 THM-1790), with explicit-small-witness refutation + exhaustive + Lean verification -- exactly
 the WOWII pipeline, which is already the repo's culture minus the AUTOMATED conjecture-generation
 front end.
""")
