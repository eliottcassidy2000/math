#!/usr/bin/env python3
"""
S652 — Erdos Problem 625: chi(G) - zeta(G) for G(n,1/2), through the repo's sigma = complement lens.

ERDOS 625 (Erdos-Gimbel ~1991, $100/$1000): for G ~ G(n,1/2), does chi(G) - zeta(G) -> infinity whp?
  chi = chromatic number (partition into INDEPENDENT sets).
  zeta = COCHROMATIC number (partition into independent sets OR cliques = the two 'pure phases').
Status: Heckel & Steiner (2024) YES for ~95% of n; conjectured chi - zeta ~ n/(log n)^3.

REPO LENS (sigma = complement, the arc's spine): sigma swaps clique<->independent, so
  zeta(G) = zeta(G^c)  (cochromatic is COMPLEMENT-INVARIANT = sigma-symmetric)  [formalized]
  chi(G) != chi(G^c) in general (chromatic is NOT sigma-symmetric).
So chi - zeta = the SIGMA-ASYMMETRY of coloring. On G(n,1/2) (sigma-self-complementary in
distribution), zeta is the sigma-symmetric core; chi-zeta is the gap from breaking sigma.

This script: (A) compute chi and zeta exactly for small graphs (incl. random G(n,1/2)); verify
zeta(G)=zeta(G^c) and that chi is NOT complement-invariant; (B) the chi-zeta gap & bounds; (C) the
Ising/two-ground-states and chi.alpha>=n (S634) reframes.
No external libs (exact small-graph chromatic/cochromatic by search).
"""
import random, itertools

def adj_from_edges(n, edges):
    A=[[False]*n for _ in range(n)]
    for (i,j) in edges: A[i][j]=A[j][i]=True
    return A
def complement(A,n):
    return [[ (i!=j) and (not A[i][j]) for j in range(n)] for i in range(n)]
def is_indep(A,S):
    return all(not A[i][j] for i in S for j in S if i<j)
def is_clique(A,S):
    return all(A[i][j] for i in S for j in S if i<j)
def chromatic(A,n):
    # min colors, each class independent. Try k=1,2,...
    verts=list(range(n))
    for k in range(1,n+1):
        # backtracking proper coloring
        color=[-1]*n
        def bt(v):
            if v==n: return True
            for c in range(k):
                if all(color[u]!=c for u in range(n) if A[v][u]):
                    color[v]=c
                    if bt(v+1): return True
                    color[v]=-1
            return False
        if bt(0): return k
    return n
def cochromatic(A,n):
    # min colors, each class clique OR independent. backtracking, but the 'valid class' test differs.
    for k in range(1,n+1):
        color=[-1]*n
        # a partial assignment is OK if each color class so far is (clique) or (independent) -- but a
        # class can switch between consistent-as-clique / consistent-as-indep until it has >=2 edges
        # decided. Track for each color: members. Validity: members form indep OR clique.
        def class_ok(members):
            return is_indep(A,members) or is_clique(A,members)
        members=[[] for _ in range(k)]
        def bt(v):
            if v==n: return True
            for c in range(k):
                members[c].append(v)
                if class_ok(members[c]):
                    if bt(v+1): return True
                members[c].pop()
            return False
        if bt(0): return k
    return n

print("="*70)
print("(A) sigma-symmetry: zeta(G)=zeta(G^c), but chi(G)!=chi(G^c) in general")
print("="*70)
random.seed(625)
print("  n | chi(G) | chi(G^c) | chi sym? | zeta(G) | zeta(G^c) | zeta sym? | chi-zeta")
for n in range(3,9):
    # one random G(n,1/2)
    edges=[(i,j) for i in range(n) for j in range(i+1,n) if random.random()<0.5]
    A=adj_from_edges(n,edges); Ac=complement(A,n)
    cG,cGc=chromatic(A,n),chromatic(Ac,n)
    zG,zGc=cochromatic(A,n),cochromatic(Ac,n)
    print(f"  {n} | {cG:6d} | {cGc:8d} | {str(cG==cGc):>7} | {zG:7d} | {zGc:9d} | {str(zG==zGc):>8} | {cG-zG}")
print("  -> zeta is ALWAYS complement-invariant (sigma-fixed); chi is NOT. [formalized cochrom_compl]")
print("     chi-zeta = the sigma-asymmetry of coloring.")

print("\n" + "="*70)
print("(B) the chi - zeta gap on G(n,1/2), averaged")
print("="*70)
print("  n | avg chi | avg zeta | avg(chi-zeta) | (Erdos 625: does chi-zeta -> oo ?)")
for n in range(4,10):
    T=40 if n<=8 else 12
    cs=zs=0
    for _ in range(T):
        edges=[(i,j) for i in range(n) for j in range(i+1,n) if random.random()<0.5]
        A=adj_from_edges(n,edges)
        cs+=chromatic(A,n); zs+=cochromatic(A,n)
    print(f"  {n} | {cs/T:7.2f} | {zs/T:8.2f} | {(cs-zs)/T:13.2f} |")
print("  -> small n: gap is tiny/noisy (it grows only slowly; conjectured ~ n/(log n)^3).")
print("     The gap is a LOWER-ORDER effect: chi ~ zeta ~ n/(2 log2 n); the difference is subtle.")

print("\n" + "="*70)
print("(C) repo reframes of zeta and chi-zeta")
print("="*70)
print("  ISING / two ground states: a clique = all-edges (ferromagnetic 'up'), an independent set =")
print("  no-edges ('down'). zeta = min #pieces each in a PURE PHASE (one of the two Ising ground")
print("  states). chi = min #pieces each 'down' only. chi-zeta = the cost of forbidding the 'up' phase.")
print("  Under sigma=complement, up<->down, so zeta is sigma-symmetric; chi is not. (Coloring-partition")
print("  unification, S626-633: zeta = partition into the two pure phases of the conflict graph.)")
print("")
print("  chi.alpha >= n (S634): chi >= n/alpha. zeta >= n/max(alpha,omega) (each pure piece <= max(a,w)).")
print("  For G(n,1/2): alpha ~ omega ~ 2 log2 n (sigma-symmetric!), chi ~ n/(2 log2 n). So zeta and chi")
print("  share the SAME leading term; Erdos 625 is about the SUBLEADING gap = where sigma-symmetry")
print("  (zeta) and sigma-asymmetry (chi) first diverge.")
print("")
print("  G(n,1/2) = the finite RADO graph (S638 did the random TOURNAMENT). zeta = its sigma-symmetric")
print("  coloring invariant. The cube in n/(log n)^3 (Heckel's conjecture) echoes the arc's 3 (cube root).")
