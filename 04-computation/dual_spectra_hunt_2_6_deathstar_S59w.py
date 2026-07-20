#!/usr/bin/env python3
"""
death-star-2026-07-20-S59w (HYP-8220) -- hunt for {2,6} as a FORBIDDEN pair in a
dual/even-sector spectrum (since H=OCF forbids {7,21} in the odd sector).
Compute several tournament invariant spectra n=3..6 and report forbidden values.
"""
from itertools import combinations, permutations

def all_tournaments(n):
    edges=list(combinations(range(n),2)); m=len(edges)
    for bits in range(2**m):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(edges):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def three_cycles(adj,n):
    c=0
    for i,j,k in combinations(range(n),3):
        # count directed 3-cycles among {i,j,k}
        for a,b,cc in [(i,j,k),(i,k,j)]:
            if adj[a][b] and adj[b][cc] and adj[cc][a]: c+=1
    return c
def score_seq(adj,n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))
def feedback_min(adj,n):
    # min feedback arc set size = arcs against the best linear order
    best=10**9
    for perm in permutations(range(n)):
        pos={v:i for i,v in enumerate(perm)}
        back=sum(1 for i in range(n) for j in range(n) if adj[i][j] and pos[i]>pos[j])
        best=min(best,back)
    return best

print("=== invariant spectra n=3..6, forbidden values (hunting {2,6}) ===")
for n in range(3,7):
    c3=set(); ss=set(); fb=set()
    for adj in all_tournaments(n):
        c3.add(three_cycles(adj,n)); ss.add(score_seq(adj,n)); fb.add(feedback_min(adj,n))
    c3=sorted(c3); fb=sorted(fb)
    forb_c3=[v for v in range(min(c3),max(c3)+1) if v not in c3]
    forb_fb=[v for v in range(min(fb),max(fb)+1) if v not in fb]
    print(f"  n={n}: #3-cycles achievable {c3}, forbidden interior {forb_c3}")
    print(f"        #score-sequences = {len(ss)}; min-feedback-arc achievable {fb}, forbidden {forb_fb}")

print("\n=== the number of NON-ISO tournaments and 3-cycle max ===")
# max #3-cycles in a regular tournament on n vertices = the 'most cyclic'
# n=3:1, n=5:5, n=7:? ; and 2,6 as small counts
from math import comb
for n in [3,4,5,6,7]:
    # max 3-cycles = C(n,3) - min; for regular tournament ~ C(n,3)*1/4 *...
    print(f"  n={n}: C(n,3)={comb(n,3)} (total triples); C(n,2)={comb(n,2)} edges")
print("  NOTE C(4,2)=6, C(7,2)=21 -> {2,6}? vs {7,21}: 6 and 21 are EDGES of K_4, K_7")
print("  and 2, 7: could be the CORE cyclic-group orders (C_2 vs C_7) of S_3, F_21")

print("\n=== {2,6} & {7,21} as {|reflection or core|, |group|} — the two exceptional")
print("    metacyclic symmetry groups of the repo's two special sizes (n=3, n=7) ===")
print("  n=3: the 3-cycle / JC fiber -> S_3 (order 6), reflection Z_2 (order 2): {2,6}")
print("  n=7: Paley T_7 / doubling freeze -> F_21 = C_7:C_3 (order 21), core C_7 (order 7): {7,21}")
print("  6 = 3! = |S_3|, 21 = |F_21|; both G = C_p:C_q with pq: S_3=3*2, F_21=7*3")
