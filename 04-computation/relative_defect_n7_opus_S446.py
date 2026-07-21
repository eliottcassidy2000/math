#!/usr/bin/env python3
"""
The RELATIVE H-defect at n=7 (opus-S446) -- does the formula/#P edge stay real, or shrink?
defect_k(n)/Hbar(n): if it does NOT shrink, H lives on the far side of the edge (not a measure-zero
correction). Exact n=5,6 gave 0.53, 0.62 at k=3. Estimate n=7 by SAMPLING (exact enumeration of
2^21 is slow); partition by c3 (=the degree-3 census, poly), track max-min H per c3 and mean H.
"""
import itertools, random
def Hcount(adj,n):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if c:
                av=adj[v]
                for u in range(n):
                    if not (mask>>u)&1 and av[u]: dp[mask|(1<<u)][u]+=c
    return sum(dp[size-1])
def c3(adj,n):
    t=0
    for i,j,k in itertools.combinations(range(n),3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]): t+=1
    return t
def rand_tournament(n,rng):
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<.5: adj[i][j]=1
            else: adj[j][i]=1
    return adj

random.seed(2027)
print("relative H-defect (defect within same-c3 class)/mean-H, k=3 (degree-3 census)")
print("="*68)
# exact small n for calibration
def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            adj[i][j]=1 if (bits>>k)&1 else 0; adj[j][i]=1-adj[i][j]
        yield adj
for n in (5,6):
    byc3={}; tot=0; ssum=0
    for adj in edges_iter(n):
        h=Hcount(adj,n); c=c3(adj,n); byc3.setdefault(c,[]).append(h); ssum+=h; tot+=1
    defect=max(max(v)-min(v) for v in byc3.values()); mean=ssum/tot
    print(f"  n={n} (exact): max defect={defect}, mean H={mean:.2f}, RELATIVE={defect/mean:.3f}")
# n=7 sampled
for n,trials in [(7,400000)]:
    byc3={}; ssum=0
    rng=random.Random(2027)
    for _ in range(trials):
        adj=rand_tournament(n,rng); h=Hcount(adj,n); c=c3(adj,n)
        d=byc3.setdefault(c,[1<<60,-1])       # [min,max]
        if h<d[0]: d[0]=h
        if h>d[1]: d[1]=h
        ssum+=h
    defect=max(d[1]-d[0] for d in byc3.values() if d[1]>=0)
    mean=ssum/trials
    print(f"  n={n} (sampled {trials}): max defect>={defect}, mean H={mean:.2f}, RELATIVE>={defect/mean:.3f}")
print("\n READING: if the relative defect does NOT shrink (0.53, 0.62, ...), H lives on the far side")
print(" of the formula/#P edge -- the missing part has positive weight, not a measure-zero correction.")
