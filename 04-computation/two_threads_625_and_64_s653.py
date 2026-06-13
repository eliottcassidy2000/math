#!/usr/bin/env python3
"""
S653 — TWO PARALLEL THREADS on the arc's 2-ADIC SEAM, worked back and forth.

THREAD A (cont. Erdos 625): zeta = the sigma-symmetric chromatic number (sigma=complement, order 2).
  New sigma-bound: zeta(G) <= min(chi(G), chi(G^c)), so chi - zeta >= chi(G) - chi(G^c) (the chromatic
  sigma-asymmetry). [formalized cochromColorable_of_(compl_)colorable]
THREAD B (Erdos 64 = Erdos-Gyarfas): min deg >=3 forces a cycle of length 2^k (a POWER OF 2). The "2^k"
  IS the doubling tower / 2-adic seam (S640). Via Sidon sequences + the summand graph + the cauldron:
  cycle length = additive-relation length; a 2^k cycle = a 2^k-term +/- relation; a SIDON set forbids
  the 4-term (=2^2) relation -> its Cayley graph has NO 4-cycle -> the smallest 2^k cycle is >= 8.
  Avoiding ALL 2^k cycles = a Sidon condition at EVERY doubling level = 'doubling-tower-Sidon'.

CROSS-CONNECTION: both threads are about the prime 2 -- 625 = the order-2 involution sigma; 64 = the
powers of 2 (doubling). They are the two faces of the arc's 2-ADIC SEAM.
No external libs.
"""
from itertools import combinations, product

# ============================ THREAD A: Erdos 625 ============================
def chromatic(adj,n):
    for k in range(1,n+1):
        col=[-1]*n
        def bt(v):
            if v==n: return True
            for c in range(k):
                if all(col[u]!=c for u in range(n) if adj[v][u]):
                    col[v]=c
                    if bt(v+1): return True
                    col[v]=-1
            return False
        if bt(0): return k
    return n
def cochromatic(adj,n):
    def ok(S):
        return all(not adj[i][j] for i,j in combinations(S,2)) or all(adj[i][j] for i,j in combinations(S,2))
    for k in range(1,n+1):
        mem=[[] for _ in range(k)]
        def bt(v):
            if v==n: return True
            for c in range(k):
                mem[c].append(v)
                if ok(mem[c]) and bt(v+1): return True
                mem[c].pop()
            return False
        if bt(0): return k
    return n
def complement(adj,n): return [[i!=j and not adj[i][j] for j in range(n)] for i in range(n)]

import random
random.seed(653)
print("="*72)
print("THREAD A (625): zeta <= min(chi, chi_c);  chi - zeta >= chi(G) - chi(G^c)")
print("="*72)
print("  n | chi(G) | chi(G^c) | zeta | zeta<=min? | chi-zeta | chi(G)-chi(G^c) | bound ok?")
for n in range(4,9):
    edges=[(i,j) for i in range(n) for j in range(i+1,n) if random.random()<0.5]
    A=[[False]*n for _ in range(n)]
    for i,j in edges: A[i][j]=A[j][i]=True
    Ac=complement(A,n)
    cG,cGc,z=chromatic(A,n),chromatic(Ac,n),cochromatic(A,n)
    le=z<=min(cG,cGc); bd=(cG-z)>=(cG-cGc)
    print(f"  {n} | {cG:6d} | {cGc:8d} | {z:4d} | {str(le):>9} | {cG-z:8d} | {cG-cGc:15d} | {bd}")
print("  -> zeta<=min(chi,chi_c) ALWAYS (formalized); so chi-zeta >= chi(G)-chi(G^c) = the chromatic")
print("     sigma-asymmetry. 625's gap is driven by how non-self-complementary chi is.")

# ============================ THREAD B: Erdos 64 ============================
def is_sidon(S):
    # STRICT Sidon (distinct pairs a<b have distinct sums) -- the convention the summand graph (a!=b)
    # detects. (Full Sidon would also forbid a+b=2c; we use the distinct-pair version to match.)
    sums=set()
    for a,b in combinations(sorted(set(S)),2):
        if a+b in sums: return False
        sums.add(a+b)
    return True

def cayley_cycles_lengths(n, gens, maxlen):
    # Cay(Z/n, gens) with gens symmetric closure; return set of cycle lengths present (<=maxlen)
    S=set()
    for g in gens: S.add(g%n); S.add((-g)%n)
    S.discard(0)
    adj=[[False]*n for _ in range(n)]
    for v in range(n):
        for s in S: adj[v][(v+s)%n]=True
    # find girth-ish: BFS shortest cycle through each vertex; collect lengths up to maxlen
    lengths=set()
    # simple: for each closed walk... use DFS for cycles up to maxlen (small n)
    import sys
    found=set()
    def dfs(start,v,visited,length):
        if length>maxlen: return
        for w in range(n):
            if adj[v][w]:
                if w==start and length>=3:
                    found.add(length+1 if False else length)  # we count edges
                elif w not in visited and w>start:  # avoid recount, canonical
                    visited.add(w); dfs(start,w,visited,length+1); visited.discard(w)
    # simpler robust: BFS for shortest even/2^k cycle existence via length-l closed walks is complex;
    # use: cycle of length L exists iff there's a simple cycle; detect via brute over small subsets is hard.
    # Instead: report girth and whether a 2^k cycle (4,8,16) exists, via shortest-cycle search per pair.
    return adj,S

def has_cycle_length(adj,n,L):
    # does the graph have a simple cycle of length exactly L? DFS (small n,L).
    if L>n: return False
    def dfs(start,v,visited,depth):
        if depth==L:
            return adj[v][start]
        for w in range(n):
            if adj[v][w] and w not in visited and w!=start:
                visited.add(w)
                if dfs(start,w,visited,depth+1): visited.discard(w); return True
                visited.discard(w)
        return False
    for start in range(n):
        if dfs(start,start,{start},1): return True
    return False

print("\n" + "="*72)
print("THREAD B (64=Erdos-Gyarfas): the cycle=additive-relation dictionary, done HONESTLY")
print("="*72)
print("  CAUTION/CORRECTION: a 4-cycle of length L = an L-term +/- relation among generators -- BUT an")
print("  abelian Cayley graph Cay(Z/n,{a,b,...}) ALWAYS has the PARALLELOGRAM 4-cycle (0,a,a+b,b), so it")
print("  is the WRONG model for 'no short cycle'. (Verified below.) The Sidon<->C4-free statement lives")
print("  in the BIPARTITE SUMMAND/incidence graph, not the Cayley graph.")
def cayley_has_c4(n,gens):
    Sg=set()
    for g in gens: Sg.add(g%n); Sg.add((-g)%n)
    Sg.discard(0)
    adj=[[False]*n for _ in range(n)]
    for v in range(n):
        for s in Sg: adj[v][(v+s)%n]=True
    return has_cycle_length(adj,n,4)
print("  abelian-Cayley parallelogram check: " +
      ", ".join(f"Z/{n}{gens}:C4={cayley_has_c4(n,gens)}" for n,gens in [(13,[1,5]),(31,[1,5,11]),(41,[1,8,14])]))
print("  -> always True (the parallelogram). So Cayley is NOT the place Sidon kills C4.")

print("\n  THE CORRECT OBJECT = the SUMMAND GRAPH (repo): S is SIDON  <=>  no two distinct pairs in S have")
print("  the same sum  <=>  the summand graph on S has NO 4-cycle. (a+b=c+d is the 4-cycle.)")
print("  set S                  | Sidon? | #repeated-sums (=#summand-C4s) | match?")
def summand_c4_count(S):
    from collections import Counter
    c=Counter()
    for a,b in combinations(sorted(set(S)),2): c[a+b]+=1
    return sum(v*(v-1)//2 for v in c.values())  # # pairs-of-pairs with equal sum = # 4-cycles
for S in [[1,2,4,8],[1,2,3,5],[1,2,5,11,22],[1,3,7,12],[0,1,3,7],[1,2,3,4,5]]:
    sid=is_sidon(S); c4=summand_c4_count(S)
    print(f"  {str(S):22s} | {str(sid):>5}  | {c4:30d} | {(c4==0)==sid}")
print("  -> Sidon <=> summand-graph C4-free, EXACTLY. The summand graph IS the Sidon/C4 detector.")
print("     Sidon sets give the EXTREMAL C4-free graphs (Erdos-Renyi/Brown) -- the densest graphs with")
print("     NO 4-cycle. THESE are the HARD instances of Erdos-Gyarfas: a C4-free min-deg-3 graph still")
print("     must (conjecturally) have an 8- or 16-cycle. Markstrom's near-counterexamples (only a 16-cycle)")
print("     live exactly here. To DISPROVE 64 you must avoid C4 AND C8 AND C16 AND ... simultaneously =")
print("     a Sidon-type condition at EVERY doubling level, at min degree >=3 -- conjecturally impossible.")

print("\n" + "="*72)
print("THE SUMMAND GRAPH & CAULDRON (repo) = the additive-relation detectors")
print("="*72)
print("  SUMMAND GRAPH (repo): a->n,b->n iff a+b=n. A 4-CYCLE there (two pairs a+b=c+d=n with {a,b}!={c,d})")
print("  IS a Sidon violation. So 'Sidon' = 'the summand graph has no parallel pairs' = no 4-cycle.")
print("  CAULDRON GAME (repo, S618): color [1,N] avoiding a monochromatic sum a+b (a Schur 'boil'). A boil")
print("  = a 3-term additive relation; Schur's theorem = you cannot avoid boils forever. Erdos-Gyarfas is")
print("  the cycle/2^k analogue: you cannot avoid power-of-2 additive relations forever at min degree 3.")

print("\n" + "="*72)
print("CROSS-CONNECTION: both threads live on the 2-ADIC SEAM")
print("="*72)
print("  625: sigma = COMPLEMENT, the ORDER-2 involution; zeta = the sigma-fixed chromatic number;")
print("       chi-zeta = the sigma-asymmetry. (the 2 = the involution.)")
print("  64 : the target cycle length is a POWER OF 2 = the DOUBLING tower (S640: ord_p(2), the 2-adic")
print("       seam); a 2^k cycle = a doubling-closed additive relation. (the 2 = the doubling.)")
print("  -> 625 is the order-2 (involution) face and 64 is the 2^k (doubling) face of the SAME prime 2.")
print("     The arc's 2-adic seam (half-turn sigma + doubling tower) is exactly where both Erdos problems")
print("     sit. Back-and-forth payoff: the cube-root/odd machinery (S637-651) is orthogonal to both;")
print("     these two are the project's PURE 2-ADIC problems.")
