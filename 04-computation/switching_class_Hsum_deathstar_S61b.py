#!/usr/bin/env python3
"""
death-star-2026-07-20-S61b (HYP-8300) -- the HIGH-LEVERAGE fiberwise odd/even bridge:
Sum over each SWITCHING CLASS of the odd-valued H = n! exactly. Switching classes = tilings
= (at odd n) even graphs / two-graphs. So the odd invariant H, summed over an even-graph
fiber, yields |S_n| -- the odd values and the even-graph structure meet on S_n, fiberwise
(refining the global Sum_T H = n!*2^C(n-1,2)). Verified exactly n=3..6.
"""
from itertools import permutations, combinations
from math import comb, factorial

def ham_paths(adj, n):
    return sum(1 for p in permutations(range(n)) if all(adj[p[k]][p[k+1]] for k in range(n-1)))

def switch_vertex(bits, v, n, idx):
    """reverse all arcs incident to v = flip bits of all edges touching v."""
    for u in range(n):
        if u!=v:
            e=(min(u,v),max(u,v)); bits ^= (1<<idx[e])
    return bits

def analyze(n):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(edges); idx={e:k for k,e in enumerate(edges)}
    def adj_of(bits):
        a=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(edges):
            if (bits>>k)&1: a[i][j]=1
            else: a[j][i]=1
        return a
    seen=set(); sums=[]; sizes=[]
    for start in range(1<<m):
        if start in seen: continue
        # BFS the switching class
        cls=set([start]); fr=[start]
        while fr:
            b=fr.pop()
            for v in range(n):
                nb=switch_vertex(b,v,n,idx)
                if nb not in cls: cls.add(nb); fr.append(nb)
        seen|=cls
        s=sum(ham_paths(adj_of(b),n) for b in cls)
        sums.append(s); sizes.append(len(cls))
    return sums, sizes, m

print("=== Sum of H over each SWITCHING CLASS (= tiling = even graph at odd n) ===\n")
for n in range(3,7):
    sums, sizes, m = analyze(n)
    nf = factorial(n)
    all_equal = all(s==nf for s in sums)
    print(f"n={n}: {len(sums)} switching classes (expect 2^C(n-1,2)={2**comb(n-1,2)}), each size {set(sizes)} (expect 2^(n-1)={2**(n-1)})")
    print(f"  every class-sum of H == n! = {nf}:  {all_equal}   (distinct sums observed: {sorted(set(sums))})")
    print()
print("THEOREM (verified n=3..6): Sum_{T in switching class} H(T) = n! for EVERY class.")
print("  => the odd-valued H, averaged over a switching class of size 2^(n-1), has mean n!/2^(n-1);")
print("     summed it is exactly |S_n|. Odd values + even (switching/even-graph) structure meet on S_n,")
print("     fiberwise. Proof: a fixed directed Ham path (n! of them) lies in exactly one tournament")
print("     per switching class (the n-1 path-arcs pin the switch coset), so each class contains each")
print("     directed Ham path exactly once => class-sum = #directed Ham paths = n!.")
