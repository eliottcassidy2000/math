#!/usr/bin/env python3
"""Verify the apex/atom H-spectrum hint: 49=7^2 and 75 as STRONG-CORE atoms at n=7,8
(irreducible since H=7 forbidden), and the APEX TILE (0,n-1) toggle role.
kind-mendel-2026-06-22-S6."""
import random
random.seed(0)
def H_count(adj, n):
    "number of directed Hamiltonian paths; adj[i] bitmask of out-neighbors"
    # f[mask] dict not needed; dp[mask][v]. Use arrays.
    size=1<<n
    # dp[v] over masks: iterate masks ascending
    import array
    dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        row=dp[mask]
        # extend
        for v in range(n):
            c=row[v]
            if not c: continue
            # from v go to w not in mask with v->w
            outs=adj[v] & ~mask
            w=0; o=outs
            while o:
                w=(o & -o).bit_length()-1
                dp[mask|(1<<w)][w]+=c
                o&=o-1
    full=size-1
    return sum(dp[full])
def is_strong(adj,n):
    "strongly connected?"
    def reach(start, A):
        seen=1<<start; stack=[start]
        while stack:
            u=stack.pop(); nb=A[u]&~seen
            while nb:
                w=(nb&-nb).bit_length()-1; seen|=1<<w; stack.append(w); nb&=nb-1
        return seen==(1<<n)-1
    radj=[0]*n
    for i in range(n):
        for j in range(n):
            if adj[i]>>j&1: radj[j]|=1<<i
    return reach(0,adj) and reach(0,radj)
def rand_tour(n):
    adj=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if random.random()<0.5: adj[i]|=1<<j
            else: adj[j]|=1<<i
    return adj

for n in [7,8]:
    seen={}; strong_vals=set()
    Nsamp = 40000 if n==7 else 120000
    for _ in range(Nsamp):
        adj=rand_tour(n); h=H_count(adj,n)
        seen[h]=seen.get(h,0)+1
        if is_strong(adj,n): strong_vals.add(h)
    print(f"n={n}: {len(seen)} distinct H sampled (max {max(seen)}); 49 in spectrum: {49 in seen} (strong:{49 in strong_vals}); 75: {75 in seen} (strong:{75 in strong_vals})")
    print(f"   strong-core (atom) H-values divisible by 7, <=200: {sorted(v for v in strong_vals if v%7==0 and v<=200)}")
    print(f"   is 49 a product of smaller H? 7 forbidden so 49=7*7 impossible => atom. (strong-realized:{49 in strong_vals})")

# APEX TILE toggle: base path n-1->...->0; apex tile = arc(0,n-1). Flip it, watch H & strong.
print("\n=== APEX TILE (0,n-1) toggle: transitive vs apex-up (master cycle) ===")
n=8
# transitive: i->j iff i>j (path n-1->...->0). adj[i]=bits j<i
adjT=[sum(1<<j for j in range(i)) for i in range(n)]
print(f"transitive n=8: H={H_count(adjT,n)} strong={is_strong(adjT,n)} (apex (0,{n-1}) is 0<-{n-1}, i.e. {n-1}->0 already in transitive)")
# flip apex: make 0->n-1 (close the master cycle 0->n-1->n-2->...->1->0? ) 
adjA=[r for r in adjT]
# remove (n-1)->0, add 0->(n-1)
adjA[n-1]&=~(1<<0); adjA[0]|=1<<(n-1)
print(f"apex-flipped: H={H_count(adjA,n)} strong={is_strong(adjA,n)}")
