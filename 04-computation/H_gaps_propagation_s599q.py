"""Extend the H-gap study + propagation. (1) n=7 H-spectrum via heavy structured+random sampling
(test if 49=7^2 is a gap, refine the 7·odd pattern). (2) the unit-distance Harborth-value
spectrum h(n)=floor(3n-sqrt(12n-3)) and ITS gaps (integers never a max-UD value) — the parallel
gap phenomenon. (3) honest hunt for a ~1.014 exponent. opus-2026-06-03-S599q."""
from itertools import combinations
from math import sqrt, floor
import random
def Hcount(n,adj):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        row=dp[mask]
        for v in range(n):
            c=row[v]
            if not c: continue
            av=adj[v]
            for w in range(n):
                if not(mask>>w&1) and av>>w&1: dp[mask|1<<w][w]+=c
    return sum(dp[size-1][v] for v in range(n))
def rot_tour(n,Aset):  # rotational: i->j iff (j-i) mod n in Aset
    adj=[0]*n
    for i in range(n):
        for j in range(n):
            if i!=j and ((j-i)%n) in Aset: adj[i]|=1<<j
    return adj
def rand_tour(n,rng):
    adj=[0]*n
    for i,j in combinations(range(n),2):
        if rng.random()<.5: adj[i]|=1<<j
        else: adj[j]|=1<<i
    return adj
def main():
    rng=random.Random(11); n=7; vals=set()
    # heavy random
    for _ in range(120000): vals.add(Hcount(n,rand_tour(n,rng)))
    # all near-transitive up to 3 flips
    base=[0]*n
    for i in range(n):
        for j in range(i+1,n): base[i]|=1<<j
    pairs=list(combinations(range(n),2))
    for k in (1,2,3):
        for combo in combinations(pairs,k):
            adj=[x for x in base]
            for (i,j) in combo:
                if adj[i]>>j&1: adj[i]&=~(1<<j); adj[j]|=1<<i
                else: adj[j]&=~(1<<i); adj[i]|=1<<j
            vals.add(Hcount(n,adj))
    mx=max(vals)
    gaps=[h for h in range(1,mx+1,2) if h not in vals]
    print(f"n=7 (sampled+structured): max seen={mx}; #achievable seen={len(vals)}")
    print(f" apparent GAPS (sampled, may overstate): {[g for g in gaps if g<=70]}")
    print(f" 7 seen? {7 in vals}; 21? {21 in vals}; 35? {35 in vals}; 49? {49 in vals}; 63? {63 in vals}")
    print(f" => 7·odd pattern (7,21,35,49,63) status: {[(7*k, (7*k in vals)) for k in (1,3,5,7,9)]}")
    print("\nUNIT-DISTANCE Harborth-value spectrum h(n)=floor(3n-sqrt(12n-3)), its GAPS (skipped maxes):")
    hvals=set(); seq=[]
    for nn in range(2,400):
        h=floor(3*nn-sqrt(12*nn-3)); hvals.add(h); seq.append(h)
    hgaps=[v for v in range(1, max(hvals)+1) if v not in hvals]
    print(f" first Harborth values: {seq[:14]}")
    print(f" first skipped integers (impossible max-UD on lattice): {hgaps[:20]}")
    print("\nHONEST exponent hunt:")
    Pn={3:3,4:5,5:15,6:45}
    import math
    for n2 in (4,5,6):
        p=Pn[n2]
        print(f"  n={n2}: P={p}; P^(1/n^2)={p**(1/n2**2):.4f}; P^(1/C(n,2))={p**(1/(n2*(n2-1)//2)):.4f}; ln P/n={math.log(p)/n2:.4f}")
    print("  (looking for ~1.014; report what matches, honestly)")
if __name__=='__main__': main()
