#!/usr/bin/env python3
"""
death-star-2026-07-20-S59w (HYP-8220) -- locate {2,6} vs {7,21} empirically.
Compute tournament spectra n=3..6 exhaustively: # Hamiltonian PATHS (H, the OCF),
# Hamiltonian CYCLES, and test the {k,3k} / group-order / Omega_p hypotheses.
"""
from itertools import combinations, permutations
from math import comb

def all_tournaments(n):
    """yield adjacency (adj[i][j]=1 iff i->j) over all 2^C(n,2) labeled tournaments."""
    edges=list(combinations(range(n),2))
    m=len(edges)
    for bits in range(2**m):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(edges):
            if (bits>>k)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj

def ham_paths(adj,n):
    # Held-Karp count of Hamiltonian paths (directed)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1<<n):
        for last in range(n):
            v=dp[mask][last]
            if not v: continue
            for nxt in range(n):
                if not (mask>>nxt)&1 and adj[last][nxt]:
                    dp[mask|(1<<nxt)][nxt]+=v
    full=(1<<n)-1
    return sum(dp[full][i] for i in range(n))

def ham_cycles(adj,n):
    # count directed Hamiltonian cycles (fix start=0 to avoid rotation overcount)
    cnt=0
    for perm in permutations(range(1,n)):
        path=(0,)+perm
        if all(adj[path[i]][path[(i+1)%n]] for i in range(n)):
            cnt+=1
    return cnt

print("=== H-PATH spectrum (OCF) and H-CYCLE spectrum, n=3..6 ===")
for n in range(3,7):
    Hp=set(); Hc=set()
    for adj in all_tournaments(n):
        Hp.add(ham_paths(adj,n))
        Hc.add(ham_cycles(adj,n))
    Hp=sorted(Hp); Hc=sorted(Hc)
    # forbidden odd H-path values below max
    mx=max(Hp)
    forb_path=[v for v in range(1,mx+1,2) if v not in Hp]
    forb_cyc=[v for v in range(0,max(Hc)+1) if v not in Hc]
    print(f"  n={n}: H-path achievable {Hp}")
    print(f"        H-path forbidden odd < {mx}: {forb_path}")
    print(f"        H-cycle achievable {Hc}; forbidden < {max(Hc)}: {forb_cyc}")

print("\n=== the {k,3k} exceptional-pair family ===")
pairs={"{2,6}":(2,6),"{7,21}":(7,21),"{12,24}":(12,24)}
for name,(a,b) in pairs.items():
    print(f"  {name}: ratio {b//a if b%a==0 else b/a}, b={a}*{b//a}, a+b={a+b}, a*b={a*b}, "
          f"a=T? {a in [k*(k+1)//2 for k in range(9)]}, b=C(?,2)? {b in [comb(k,2) for k in range(9)]}")
# 6 = C(4,2) edges of K_4; 21 = C(7,2) edges of K_7; is 2,7 the vertex counts minus?
print("  6 = C(4,2)=|E(K_4)|, 21 = C(7,2)=|E(K_7)|; 2 and 7 as vertex-ish: 4 and 7 (K_4,K_7)")
print("  observer 2n+1: 2*3+1=7 (base 2 -> base 7 via 2*3+1); 6*3+... ")

print("\n=== group-order hypothesis: S_3 vs F_21 ===")
print("  S_3 = C_3 rtimes C_2, |S_3|=6, subgroups of order {1,2,3,6}; C_2 (reflection)=2")
print("  F_21 = C_7 rtimes C_3, |F_21|=21, core C_7=7")
print("  {2,6}: {|C_2 reflection|, |S_3|}? ; {7,21}: {|C_7 core|, |F_21|}")
print("  parallel: both G=C_p rtimes C_q with |G|=pq: (p,q)=(3,2)->6, (7,3)->21")

print("\n=== Omega_p / even-graph dims ===")
print("  Paley T_7 Omega_p: 7,21,42,63,63,42,21 -> first two = {7,21}")
print("  V(E_n)=2,3,7,16,54 (n=3..7): E_3=2, E_5=7 -> the '2' and '7' bases live at n=3,5")
