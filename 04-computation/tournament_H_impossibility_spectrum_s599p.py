"""The H-IMPOSSIBILITY SPECTRUM: which Hamiltonian-path counts are achievable by tournaments?
Rédei => all odd. But some odd values are GAPS (impossible): 7, 21, ... Compute the spectrum
exhaustively n=3..6; test whether 7,21 EVER appear (incl near-transitive tournaments at n=7,8
where small H values live); growth of P(n)=max H and look for a ~1.014 exponent honestly.
opus-2026-06-03-S599p."""
from itertools import combinations, product
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
                if not (mask>>w&1) and av>>w&1:
                    dp[mask|1<<w][w]+=c
    full=size-1
    return sum(dp[full][v] for v in range(n))
def tour(n,bits):
    adj=[0]*n
    for (i,j),b in zip(combinations(range(n),2),bits):
        if b: adj[i]|=1<<j
        else: adj[j]|=1<<i
    return adj
def transitive(n):
    adj=[0]*n
    for i in range(n):
        for j in range(i+1,n): adj[i]|=1<<j
    return adj
def flip_edges(n,adj,pairs):
    adj=adj[:]
    for (i,j) in pairs:
        if adj[i]>>j&1: adj[i]&=~(1<<j); adj[j]|=1<<i
        else: adj[j]&=~(1<<i); adj[i]|=1<<j
    return adj
def main():
    print("FULL SPECTRUM (exhaustive) — achievable Hamiltonian-path counts; gaps = impossible odds")
    Pn={}
    allvals=set()
    for n in range(3,7):
        vals=set()
        for bits in product((0,1),repeat=n*(n-1)//2):
            vals.add(Hcount(n,tour(n,bits)))
        Pn[n]=max(vals); allvals|=vals
        mx=max(vals)
        gaps=[h for h in range(1,mx+1,2) if h not in vals]
        print(f" n={n}: max P(n)={mx}; #achievable={len(vals)}; GAPS (impossible odd ≤max)= {gaps}")
    print(f"\n union of achievable (n≤6): contains 7? {7 in allvals}  contains 21? {21 in allvals}")
    print("\nTEST: do 7 or 21 appear among NEAR-TRANSITIVE tournaments at n=7,8 (where small H live)?")
    found7=found21=False; seen=set()
    rng=random.Random(7)
    for n in (7,8):
        base=transitive(n); pairs=list(combinations(range(n),2))
        # all <=3-edge flips + random larger
        for k in (1,2,3):
            for combo in combinations(pairs,k):
                h=Hcount(n,flip_edges(n,base,list(combo))); seen.add(h)
                if h==7: found7=True
                if h==21: found21=True
        for _ in range(4000):
            k=rng.randint(4,8); combo=rng.sample(pairs,k)
            h=Hcount(n,flip_edges(n,base,combo)); seen.add(h)
            if h==7: found7=True
            if h==21: found21=True
    print(f" near-transitive n=7,8 scan: 7 found={found7}; 21 found={found21}; small odds seen={sorted(x for x in seen if x<=25)}")
    print("\nGROWTH: P(n)=max Ham paths, and exponents")
    import math
    ns=sorted(Pn); 
    for i in range(1,len(ns)):
        n=ns[i]; p=Pn[n]; pm=Pn[ns[i-1]]
        avg=math.factorial(n)/2**(n-1)
        print(f" n={n}: P={p}; P/P_prev={p/pm:.3f}; P^(1/n)={p**(1/n):.4f}; P/avg={p/avg:.3f}; (P/avg)^(1/n)={(p/avg)**(1/n):.4f}")
if __name__=='__main__': main()
